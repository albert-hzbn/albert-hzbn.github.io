// Phase-Field Lab: application controller. Owns the WebGPU device, the
// simulation, the renderers and all panel controls.
import { MODELS, resolveModel, defaults } from './models/index.js';
import { rng } from './models/init.js';
import { Simulation } from './gpu/simulation.js';
import { Renderer, planePolygon } from './render/renderer.js';
import { gradientCSS, CMAP_IDS } from './ui/colormaps.js';
import { exportCanvas, exportVTK, exportSliceCSV } from './ui/export.js';

const $ = (id) => document.getElementById(id);
const MODES = { surface: 0, volume: 1, iso: 2 };

const app = {
    device: null, sim: null, renderer: null, ctx3d: null, ctx2d: null,
    model: null, resolved: null, params: {}, init: {}, n: 96, seed: 1,
    view: null,
    display: { lo: 0, hi: 1, cmap: 0, mode: 'surface', iso: 0.5, opacity: 0.6, gamma: 1.5, quality: 256, showBox: true, autoRotate: false },
    plane: { n: [0, 0, 1], d: 0 },
    showPlane: true, clip: false, flip: false, smooth: true,
    running: false, speed: 'auto', autoSteps: 8, dirty: true, gpuBusy: false, needsReset: false,
    rate: { steps: 0, since: performance.now(), value: 0 },
};

// ---------------------------------------------------------------------------
// Small UI helpers
// ---------------------------------------------------------------------------

let toastTimer;
function toast(msg, isError = false) {
    const t = $('toast');
    t.textContent = msg;
    t.classList.toggle('error', isError);
    t.classList.add('show');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => t.classList.remove('show'), isError ? 6000 : 2600);
}

const stripTags = (html) => { const d = document.createElement('div'); d.innerHTML = html; return d.textContent; };
const decimals = (step) => Math.max(0, Math.min(6, -Math.floor(Math.log10(step) + 1e-9)));
const fmt = (v, step) => (Math.abs(v) < 1e-12 ? '0' : Number(v).toFixed(decimals(step)));

function segmented(el, onChange) {
    el.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
        el.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', String(x === b)));
        onChange(b.dataset.v, b);
    }));
}

function setSeg(el, v) {
    el.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', String(x.dataset.v === v)));
}

function busy(on, text = 'Preparing…') {
    $('busy-text').textContent = text;
    $('busy').classList.toggle('on', on);
}

// ---------------------------------------------------------------------------
// WebGPU set-up
// ---------------------------------------------------------------------------

async function initGPU() {
    if (!navigator.gpu) throw new Error('This browser does not expose WebGPU.');
    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) throw new Error('No WebGPU adapter was found (is hardware acceleration disabled?).');
    const lim = adapter.limits;
    const device = await adapter.requestDevice({
        requiredLimits: {
            maxStorageBufferBindingSize: Math.min(lim.maxStorageBufferBindingSize, 2 ** 31 - 4),
            maxBufferSize: Math.min(lim.maxBufferSize, 2 ** 31),
            maxStorageBuffersPerShaderStage: Math.min(lim.maxStorageBuffersPerShaderStage, 10),
        },
    });
    device.lost.then((info) => { if (info.reason !== 'destroyed') showNoGPU(`The GPU device was lost (${info.message}). Reload the page to continue.`); });
    device.addEventListener('uncapturederror', (e) => { console.error(e.error); toast('GPU error: ' + e.error.message.split('\n')[0], true); });
    app.device = device;
    const info = adapter.info || {};
    $('st-gpu').textContent = `GPU: ${[info.vendor, info.architecture || info.description].filter(Boolean).join(' ') || 'WebGPU'} · buffer limit ${(device.limits.maxStorageBufferBindingSize / 2 ** 20).toFixed(0)} MB`;
    const format = navigator.gpu.getPreferredCanvasFormat();
    app.ctx3d = $('c3d').getContext('webgpu');
    app.ctx2d = $('c2d').getContext('webgpu');
    for (const ctx of [app.ctx3d, app.ctx2d]) ctx.configure({ device, format, alphaMode: 'opaque' });
    app.sim = new Simulation(device);
    app.renderer = new Renderer(device, format);
}

function showNoGPU(msg) {
    $('nogpu').hidden = false;
    if (msg) $('nogpu-msg').textContent = msg;
    busy(false);
}

// ---------------------------------------------------------------------------
// Model panel
// ---------------------------------------------------------------------------

function buildModelTiles() {
    const box = $('model-tiles');
    box.innerHTML = MODELS.map((m) => `<button type="button" data-model="${m.id}" aria-pressed="false" title="${m.name}"><svg viewBox="0 0 24 24">${m.icon}</svg><span>${m.short}</span></button>`).join('');
    box.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => selectModel(b.dataset.model)));
}

function paramControl(p, value, onInput) {
    const row = document.createElement('div');
    row.className = 'param';
    row.innerHTML = `<label title="${p.hint ? stripTags(p.hint) : ''}">${p.label}</label>
        <input type="range" min="${p.min}" max="${p.max}" step="${p.step}" value="${value}" aria-label="${stripTags(p.label)}">
        <input type="number" min="${p.min}" max="${p.max}" step="${p.step}" value="${fmt(value, p.step)}" aria-label="${stripTags(p.label)} value">`;
    const [range, number] = row.querySelectorAll('input');
    range.addEventListener('input', () => { number.value = fmt(+range.value, p.step); onInput(+range.value); });
    number.addEventListener('change', () => { const v = parseFloat(number.value); if (Number.isFinite(v)) { range.value = v; onInput(v); } });
    row.set = (v) => { range.value = v; number.value = fmt(v, p.step); };
    return row;
}

// KaTeX is loaded on demand; until it arrives (or if the CDN is unreachable)
// the LaTeX source is shown instead.
let katexReady = null;
const loadKatex = () => (katexReady ??= import('https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.mjs')
    .then((mod) => mod.default).catch(() => null));

function renderEquations(box, m) {
    const note = m.equationNote ? `<p class="eq-note">${m.equationNote}</p>` : '';
    box.innerHTML = m.equations.map(() => '<div class="eq-line"></div>').join('') + note;
    const lines = box.querySelectorAll('.eq-line');
    loadKatex().then((katex) => {
        if (app.model !== m) return;
        m.equations.forEach((tex, k) => {
            if (katex) katex.render(tex, lines[k], { displayMode: true, throwOnError: false });
            else { lines[k].textContent = tex; lines[k].classList.add('raw'); }
        });
    });
}
loadKatex();

const controls = { params: {}, init: {} };

function renderModelPanel() {
    const m = app.model;
    $('model-icon').innerHTML = m.icon;
    $('model-name').textContent = m.name;
    $('model-short').textContent = `${m.params.length} parameters · ${m.views.length} views`;
    $('model-desc').innerHTML = m.description;
    renderEquations($('model-eq'), m);
    document.querySelectorAll('#model-tiles button').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.model === m.id)));

    // Presets
    const pre = $('presets');
    pre.innerHTML = (m.presets || []).map((p, k) => `<button type="button" data-k="${k}" aria-pressed="false">${p.label}</button>`).join('');
    pre.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => applyPreset(m.presets[+b.dataset.k], b)));

    // Initial-condition fields
    const initBox = $('init-fields');
    initBox.innerHTML = '';
    controls.init = {};
    for (const p of m.init) {
        const lab = document.createElement('label');
        lab.className = 'field';
        if (p.type === 'select') {
            lab.innerHTML = `<span>${p.label}</span><select>${p.options.map((o, k) => `<option value="${k}">${o}</option>`).join('')}</select>`;
            const sel = lab.querySelector('select');
            sel.value = app.init[p.key];
            sel.addEventListener('change', () => { app.init[p.key] = +sel.value; markNeedsReset(); });
            lab.set = (v) => { sel.value = v; };
        } else {
            lab.innerHTML = `<span title="${p.hint ? stripTags(p.hint) : ''}">${p.label}</span><input type="number" min="${p.min}" max="${p.max}" step="${p.step}" value="${app.init[p.key]}">`;
            const inp = lab.querySelector('input');
            inp.addEventListener('change', () => { const v = parseFloat(inp.value); if (Number.isFinite(v)) { app.init[p.key] = v; markNeedsReset(); } });
            lab.set = (v) => { inp.value = v; };
        }
        controls.init[p.key] = lab;
        initBox.appendChild(lab);
    }

    // Parameters, grouped
    const groups = $('param-groups');
    groups.innerHTML = '';
    controls.params = {};
    const byGroup = {};
    for (const p of m.params) (byGroup[p.group || 'Parameters'] = byGroup[p.group || 'Parameters'] || []).push(p);
    for (const [g, list] of Object.entries(byGroup)) {
        const t = document.createElement('div');
        t.className = 'group-title';
        t.textContent = g;
        groups.appendChild(t);
        const box = document.createElement('div');
        box.className = 'param-group';
        for (const p of list) {
            const row = paramControl(p, app.params[p.key], (v) => { app.params[p.key] = v; app.sim.setParams(app.params); app.dirty = true; });
            controls.params[p.key] = row;
            box.appendChild(row);
        }
        groups.appendChild(box);
    }
    $('grid').value = String(app.n);
    $('seed').value = app.seed;
    updateMemoryNote();
}

function applyPreset(preset, button) {
    $('presets').querySelectorAll('button').forEach((b) => b.setAttribute('aria-pressed', String(b === button)));
    for (const [k, v] of Object.entries(preset.params || {})) { app.params[k] = v; controls.params[k]?.set(v); }
    for (const [k, v] of Object.entries(preset.init || {})) { app.init[k] = v; controls.init[k]?.set(v); }
    app.sim.setParams(app.params);
    reset();
    toast(`Preset: ${preset.label}`);
}

function markNeedsReset() {
    app.needsReset = true;
    $('btn-apply').classList.add('pulse');
    updateMemoryNote();
}

function updateMemoryNote() {
    const n = +$('grid').value;
    const resolved = resolveModel(app.model, app.init);
    const bytes = Simulation.memoryEstimate(resolved, n);
    const limit = app.device ? app.device.limits.maxStorageBufferBindingSize : Infinity;
    const largest = Math.max(...resolved.fields.map((f) => f.comps)) * n ** 3 * 4;
    const tooBig = largest > limit;
    $('mem-note').innerHTML = `GPU memory ≈ <b>${(bytes / 2 ** 20).toFixed(0)} MB</b> for ${n}³ = ${(n ** 3).toLocaleString()} cells${tooBig ? ' · <span class="warn">too large for this GPU</span>' : ''}`;
    $('mem-note').classList.toggle('warn', tooBig);
}

function selectModel(id, { keepGrid = false } = {}) {
    const m = MODELS.find((x) => x.id === id);
    if (!m) return;
    app.model = m;
    app.params = defaults(m.params);
    app.init = defaults(m.init);
    if (!keepGrid) app.n = m.grid;
    app.view = null;
    renderModelPanel();
    reset();
    try { localStorage.setItem('pf-model', id); } catch { /* ignore */ }
}

// Build buffers and pipelines for the current settings and upload the initial state.
function reset() {
    const m = app.model;
    app.n = +$('grid').value || app.n;
    app.seed = parseInt($('seed').value, 10) || 1;
    app.valid = false;
    busy(true, `Initialising ${m.short} at ${app.n}³…`);
    setTimeout(async () => {
        try {
            const resolved = resolveModel(m, app.init);
            app.resolved = resolved;
            app.sim.load(resolved, app.n);
            const errs = await app.sim.compileErrors();
            if (errs.length) throw new Error('Shader error: ' + errs[0]);
            app.sim.setParams(app.params);
            app.sim.upload(resolved.initialState(app.n, app.init, rng(app.seed), app.params), app.seed);
            app.renderer.setTexture(app.sim.texture);
            const keep = app.view && resolved.views.find((v) => v.id === app.view.id);
            setView(keep ? keep.id : resolved.views[0].id, !keep);
            app.needsReset = false;
            $('btn-apply').classList.remove('pulse');
            app.dirty = true;
            updateStatus();
            app.valid = true;
        } catch (e) {
            // Never keep submitting a broken pipeline.
            app.valid = false;
            setRunning(false);
            console.error(e);
            toast(e.message, true);
        } finally {
            busy(false);
        }
    }, 30);
}

// ---------------------------------------------------------------------------
// Display panel
// ---------------------------------------------------------------------------

function fillViews() {
    $('view').innerHTML = app.resolved.views.map((v) => `<option value="${v.id}">${stripTags(v.label)}</option>`).join('');
}

function setView(id, applyDefaults = true) {
    fillViews();
    const v = app.resolved.views.find((x) => x.id === id) || app.resolved.views[0];
    app.view = v;
    $('view').value = v.id;
    if (applyDefaults) {
        const d = app.display;
        [d.lo, d.hi] = v.range || [0, 1];
        if (v.cmap) d.cmap = CMAP_IDS[v.cmap];
        const r = v.render || {};
        d.mode = r.mode || 'surface';
        d.iso = r.iso ?? 0.5;
        d.opacity = r.opacity ?? 0.6;
        d.gamma = r.gamma ?? 1.5;
        syncDisplayControls();
    }
    $('cmap').disabled = !v.cmap;
    $('range-auto').disabled = !v.cmap || v.autoRange === false;
    $('hud-view').innerHTML = v.label;
    updateColorbar();
    app.dirty = true;
}

function syncDisplayControls() {
    const d = app.display;
    $('cmap').value = String(d.cmap);
    $('range-lo').value = +d.lo.toPrecision(4);
    $('range-hi').value = +d.hi.toPrecision(4);
    setSeg($('render-mode'), d.mode);
    for (const [id, v, digits] of [['iso', d.iso, 2], ['opacity', d.opacity, 2], ['gamma', d.gamma, 1], ['quality', d.quality, 0]]) {
        $(id).value = v;
        $(id + '-out').value = Number(v).toFixed(digits);
    }
    $('iso-row').hidden = d.mode !== 'iso';
    $('opacity-row').hidden = d.mode !== 'volume';
    $('gamma-row').hidden = d.mode !== 'volume';
}

function updateColorbar() {
    const v = app.view;
    $('cbar').hidden = !v || !v.cmap;
    if (!v || !v.cmap) return;
    $('cbar-grad').style.background = gradientCSS(app.display.cmap);
    $('cbar-lo').textContent = +app.display.lo.toPrecision(3);
    $('cbar-hi').textContent = +app.display.hi.toPrecision(3);
}

function updateColorbarLabels() {
    $('cbar-lo').textContent = +app.display.lo.toPrecision(3);
    $('cbar-hi').textContent = +app.display.hi.toPrecision(3);
    $('range-hi').value = +app.display.hi.toPrecision(4);
}

async function autoRange() {
    const v = app.view;
    const [name, comp] = v.source || [app.resolved.fields[0].name, 0];
    const data = await app.sim.readField(name);
    const n3 = app.n ** 3;
    let lo = Infinity, hi = -Infinity;
    for (let i = comp * n3; i < (comp + 1) * n3; i++) {
        const x = data[i];
        if (!Number.isFinite(x) || (v.ignoreBelow !== undefined && x < v.ignoreBelow)) continue;
        if (x < lo) lo = x;
        if (x > hi) hi = x;
    }
    if (!(hi > lo)) { toast('The field is uniform; the range was not changed.'); return; }
    app.display.lo = lo;
    app.display.hi = hi;
    syncDisplayControls();
    updateColorbar();
    app.dirty = true;
    toast(`Colour range ${+lo.toPrecision(3)} … ${+hi.toPrecision(3)}`);
}

function initDisplayPanel() {
    const d = app.display;
    $('view').addEventListener('change', () => setView($('view').value));
    $('cmap').addEventListener('change', () => { d.cmap = +$('cmap').value; updateColorbar(); app.dirty = true; });
    $('range-lo').addEventListener('change', () => { d.lo = parseFloat($('range-lo').value); updateColorbar(); app.dirty = true; });
    $('range-hi').addEventListener('change', () => { d.hi = parseFloat($('range-hi').value); updateColorbar(); app.dirty = true; });
    $('range-auto').addEventListener('click', () => autoRange().catch((e) => toast(e.message, true)));
    segmented($('render-mode'), (v) => { d.mode = v; syncDisplayControls(); });
    for (const [id, key, digits] of [['iso', 'iso', 2], ['opacity', 'opacity', 2], ['gamma', 'gamma', 1], ['quality', 'quality', 0]]) {
        $(id).addEventListener('input', () => { d[key] = +$(id).value; $(id + '-out').value = d[key].toFixed(digits); });
    }
    $('show-box').addEventListener('change', (e) => { d.showBox = e.target.checked; });
    $('auto-rotate').addEventListener('change', (e) => { d.autoRotate = e.target.checked; });
    const dirs = { x: [1, 0.0001, 0.0001], y: [0.0001, 1, 0.0001], z: [0.0001, 0.001, 1], iso: [1, 1, 0.8] };
    document.querySelectorAll('[data-look]').forEach((b) => b.addEventListener('click', () => app.renderer.camera.setDirection(dirs[b.dataset.look])));
}

// ---------------------------------------------------------------------------
// Slice panel
// ---------------------------------------------------------------------------

function normalFromAngles(theta, phi) {
    const t = theta * Math.PI / 180, p = phi * Math.PI / 180;
    return [Math.sin(t) * Math.cos(p), Math.sin(t) * Math.sin(p), Math.cos(t)];
}

function setPlaneNormal(n, { fromAngles = false } = {}) {
    const l = Math.hypot(...n);
    if (!(l > 1e-9)) return toast('The normal cannot be zero.', true);
    app.plane.n = n.map((v) => v / l);
    if (!fromAngles) {
        const [x, y, z] = app.plane.n;
        $('pl-theta').value = (Math.acos(Math.max(-1, Math.min(1, z))) * 180 / Math.PI).toFixed(1);
        $('pl-phi').value = (Math.atan2(y, x) * 180 / Math.PI).toFixed(1);
    }
    $('pl-theta-out').value = `${(+$('pl-theta').value).toFixed(1)}°`;
    $('pl-phi-out').value = `${(+$('pl-phi').value).toFixed(1)}°`;
    document.querySelectorAll('#plane-presets button').forEach((b) => {
        const v = b.dataset.n.split(',').map(Number), lv = Math.hypot(...v);
        b.setAttribute('aria-pressed', String(v.every((c, k) => Math.abs(c / lv - app.plane.n[k]) < 1e-3)));
    });
    updatePlaneNote();
}

function setPlaneOffset(d) {
    app.plane.d = Math.max(-0.87, Math.min(0.87, d));
    $('pl-d').value = app.plane.d;
    $('pl-d-out').value = app.plane.d.toFixed(3);
    updatePlaneNote();
}

function updatePlaneNote() {
    const n = app.plane.n.map((v) => (Math.abs(v) < 5e-4 ? 0 : v));
    const poly = planePolygon(app.plane.n, app.plane.d);
    let area = 0;
    for (let k = 1; k + 1 < poly.length; k++) {
        const a = poly[k].map((v, i) => v - poly[0][i]), b = poly[k + 1].map((v, i) => v - poly[0][i]);
        area += Math.hypot(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]) / 2;
    }
    const txt = `n = (${n.map((v) => v.toFixed(3)).join(', ')}), offset ${(app.plane.d * app.n).toFixed(1)} cells`;
    $('plane-note').textContent = poly.length ? `${txt} · section ${poly.length}-gon, ${(area * app.n * app.n).toFixed(0)} cells²` : `${txt} · outside the box`;
    $('slice-label').textContent = `n = (${n.map((v) => v.toFixed(2)).join(', ')}) · d = ${(app.plane.d * app.n).toFixed(1)}`;
}

function initSlicePanel() {
    document.querySelectorAll('#plane-presets button').forEach((b) => b.addEventListener('click', () => {
        const v = b.dataset.n.split(',').map(Number);
        [$('pl-h').value, $('pl-k').value, $('pl-l').value] = v;
        setPlaneNormal(v);
    }));
    $('pl-apply').addEventListener('click', () => setPlaneNormal([+$('pl-h').value, +$('pl-k').value, +$('pl-l').value]));
    const fromAngles = () => setPlaneNormal(normalFromAngles(+$('pl-theta').value, +$('pl-phi').value), { fromAngles: true });
    $('pl-theta').addEventListener('input', fromAngles);
    $('pl-phi').addEventListener('input', fromAngles);
    $('pl-d').addEventListener('input', () => setPlaneOffset(+$('pl-d').value));
    $('show-plane').addEventListener('change', (e) => { app.showPlane = e.target.checked; });
    $('clip').addEventListener('change', (e) => { app.clip = e.target.checked; });
    $('clip-flip').addEventListener('change', (e) => { app.flip = e.target.checked; });
    $('smooth').addEventListener('change', (e) => { app.smooth = e.target.checked; });
    $('look-normal').addEventListener('click', () => app.renderer.camera.setDirection(app.plane.n.map((v, k) => v + (k === 1 ? 1e-4 : 0))));
    setPlaneNormal([0, 0, 1]);
    setPlaneOffset(0);
}

// ---------------------------------------------------------------------------
// Canvas interaction, triad and scale bar
// ---------------------------------------------------------------------------

function initInteraction() {
    const c3 = $('c3d'), c2 = $('c2d');
    let drag = null;
    c3.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY }; c3.setPointerCapture(e.pointerId); });
    c3.addEventListener('pointermove', (e) => {
        if (!drag) return;
        app.renderer.camera.rotate(e.clientX - drag.x, e.clientY - drag.y);
        drag = { x: e.clientX, y: e.clientY };
    });
    c3.addEventListener('pointerup', () => { drag = null; });
    c3.addEventListener('wheel', (e) => { e.preventDefault(); app.renderer.camera.zoom(Math.exp(e.deltaY * 0.0012)); }, { passive: false });
    c3.addEventListener('dblclick', () => app.renderer.camera.setDirection([1, 1, 0.8]));

    // Slice view: wheel moves the plane, drag tilts it.
    let sdrag = null;
    c2.addEventListener('wheel', (e) => { e.preventDefault(); setPlaneOffset(app.plane.d - e.deltaY * 0.0006); }, { passive: false });
    c2.addEventListener('pointerdown', (e) => { sdrag = { x: e.clientX, y: e.clientY, t: +$('pl-theta').value, p: +$('pl-phi').value }; c2.setPointerCapture(e.pointerId); });
    c2.addEventListener('pointermove', (e) => {
        if (!sdrag) return;
        const t = Math.max(0, Math.min(180, sdrag.t + (e.clientY - sdrag.y) * 0.25));
        const p = ((sdrag.p + (e.clientX - sdrag.x) * 0.25 + 540) % 360) - 180;
        $('pl-theta').value = t;
        $('pl-phi').value = p;
        setPlaneNormal(normalFromAngles(t, p), { fromAngles: true });
    });
    c2.addEventListener('pointerup', () => { sdrag = null; });
}

function drawTriad() {
    const cam = app.renderer.camera;
    const { view } = cam.matrices(1);
    const axes = [[1, 0, 0, '#e2574c', 'x'], [0, 1, 0, '#3aa35b', 'y'], [0, 0, 1, '#3e7fd6', 'z']].map(([x, y, z, col, lab]) => {
        const sx = view[0] * x + view[4] * y + view[8] * z;
        const sy = view[1] * x + view[5] * y + view[9] * z;
        const sz = view[2] * x + view[6] * y + view[10] * z;
        return { sx, sy, sz, col, lab };
    }).sort((a, b) => a.sz - b.sz);
    $('triad').innerHTML = axes.map((a) => {
        const x = a.sx * 30, y = -a.sy * 30, op = a.sz < -0.3 ? 0.45 : 1;
        return `<g opacity="${op}"><line x1="0" y1="0" x2="${x}" y2="${y}" stroke="${a.col}" stroke-width="3" stroke-linecap="round"/><circle cx="${x}" cy="${y}" r="9" fill="${a.col}"/><text x="${x}" y="${y + 0.5}">${a.lab}</text></g>`;
    }).join('');
}

function updateScalebar(info) {
    const bar = $('scalebar');
    if (!info) { bar.hidden = true; return; }
    const cssW = $('c2d').clientWidth;
    const cellsAcross = info.width * app.n;
    const target = cellsAcross * 0.22;
    const pow = 10 ** Math.floor(Math.log10(target));
    const nice = [1, 2, 5, 10].map((m) => m * pow).reduce((best, v) => (Math.abs(v - target) < Math.abs(best - target) ? v : best));
    bar.hidden = false;
    bar.querySelector('i').style.width = `${(nice / cellsAcross) * cssW}px`;
    const dx = app.params.dx;
    bar.querySelector('span').textContent = dx && dx !== 1 ? `${nice} cells = ${+(nice * dx).toPrecision(3)}` : `${nice} cells`;
}

// ---------------------------------------------------------------------------
// Frame loop
// ---------------------------------------------------------------------------

function resizeCanvas(c) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(1, Math.round(c.clientWidth * dpr)), h = Math.max(1, Math.round(c.clientHeight * dpr));
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
}

function themeColors() {
    const dark = document.documentElement.dataset.theme === 'dark';
    return {
        bg: dark ? [0.105, 0.118, 0.133] : [0.965, 0.970, 0.978],
        line: dark ? [0.72, 0.76, 0.8, 0.55] : [0.18, 0.2, 0.24, 0.55],
        accent: dark ? [0.35, 0.65, 0.93] : [0.18, 0.49, 0.88],
    };
}

function render3D() {
    const c = $('c3d');
    if (!c.clientWidth) return;
    resizeCanvas(c);
    const th = themeColors();
    const d = app.display;
    app.renderer.render3D(app.ctx3d, c, {
        mode: MODES[d.mode], iso: d.iso, opacity: d.opacity, gamma: d.gamma, steps: d.quality,
        clip: app.clip, flip: app.flip, showPlane: app.showPlane, plane: app.plane, showBox: d.showBox,
        bg: th.bg, accent: th.accent, lineColor: th.line,
    });
}

function renderSlice() {
    const c = $('c2d');
    if (!c.clientWidth) return null;
    resizeCanvas(c);
    return app.renderer.renderSlice(app.ctx2d, c, { plane: app.plane, bg: themeColors().bg, smooth: app.smooth });
}

function stepsPerFrame() {
    return app.speed === 'auto' ? Math.max(1, Math.round(app.autoSteps)) : +app.speed;
}

function frame() {
    requestAnimationFrame(frame);
    if (!app.valid || !app.sim?.texture || app.gpuBusy || $('busy').classList.contains('on')) return;
    const t0 = performance.now();
    let k = 0;
    if (app.running) {
        k = stepsPerFrame();
        app.sim.advance(k);
        app.dirty = true;
    }
    if (app.dirty && app.view) {
        // Views of "time of transformation" stretch their range with the run.
        if (app.view.followSteps && app.sim.step > 0 && app.display.lo === 0) {
            app.display.hi = Math.max(10, app.sim.step);
            if (app.running) updateColorbarLabels();
        }
        app.sim.display(app.view, app.display);
        app.dirty = false;
    }
    if (app.display.autoRotate) app.renderer.camera.yaw += 0.004;
    render3D();
    updateScalebar(renderSlice());
    drawTriad();
    app.gpuBusy = true;
    app.device.queue.onSubmittedWorkDone().then(() => {
        app.gpuBusy = false;
        const ms = performance.now() - t0;
        if (app.running && app.speed === 'auto' && k) {
            // Aim for ~12 ms of GPU work per frame.
            app.autoSteps = Math.max(1, Math.min(256, app.autoSteps * Math.max(0.6, Math.min(1.4, 12 / Math.max(ms, 1)))));
        }
        app.rate.steps += k;
        const now = performance.now();
        if (now - app.rate.since > 800) {
            app.rate.value = app.rate.steps * 1000 / (now - app.rate.since);
            app.rate.steps = 0;
            app.rate.since = now;
        }
        updateStatus();
    });
}

function updateStatus() {
    const s = app.sim;
    const t = s.time;
    const tStr = t < 1e-2 || t > 1e5 ? t.toExponential(2) : t.toPrecision(4);
    $('st-model').textContent = app.model.name;
    $('st-step').textContent = s.step.toLocaleString();
    $('st-time').textContent = tStr;
    $('st-rate').textContent = app.running ? `${Math.round(app.rate.value).toLocaleString()} steps/s` : 'paused';
    $('st-grid').textContent = `${app.n}³ grid`;
    $('hud-model').textContent = app.model.name;
    $('hud-stats').innerHTML = `t = <b>${tStr}</b><br>${s.step.toLocaleString()} steps${app.running ? `<br>${Math.round(app.rate.value).toLocaleString()} steps/s` : ''}`;
}

// ---------------------------------------------------------------------------
// Chrome: transport, layout, tabs, menus, export, keyboard, theme
// ---------------------------------------------------------------------------

function setRunning(on) {
    app.running = on;
    $('btn-run').setAttribute('aria-pressed', String(on));
    $('btn-run').querySelector('span').textContent = on ? 'Pause' : 'Run';
    if (on && app.needsReset) toast('Some settings apply only after “Reset with these settings”.');
    app.rate = { steps: 0, since: performance.now(), value: app.rate.value };
    updateStatus();
}

function initChrome() {
    $('btn-run').addEventListener('click', () => setRunning(!app.running));
    $('btn-step').addEventListener('click', () => { if (!app.valid) return; app.sim.advance(stepsPerFrame()); app.dirty = true; updateStatus(); });
    $('btn-reset').addEventListener('click', () => reset());
    $('btn-apply').addEventListener('click', () => reset());
    $('speed').addEventListener('change', (e) => { app.speed = e.target.value; });
    $('grid').addEventListener('change', markNeedsReset);
    $('seed').addEventListener('change', markNeedsReset);
    segmented($('layout'), (v) => { $('viewport').dataset.layout = v; });
    document.querySelectorAll('.tabs [data-tab]').forEach((b) => b.addEventListener('click', () => {
        document.querySelectorAll('.tabs [data-tab]').forEach((x) => x.setAttribute('aria-selected', String(x === b)));
        document.querySelectorAll('[data-body]').forEach((x) => { x.hidden = x.dataset.body !== b.dataset.tab; });
    }));
    $('panel-toggle').addEventListener('click', () => document.body.classList.toggle('panel-open'));
    document.querySelectorAll('.menu').forEach((m) => m.querySelector('[data-menu]').addEventListener('click', (e) => { e.stopPropagation(); m.classList.toggle('open'); }));
    document.addEventListener('click', () => document.querySelectorAll('.menu.open').forEach((o) => o.classList.remove('open')));
    document.querySelectorAll('[data-export]').forEach((b) => b.addEventListener('click', () => doExport(b.dataset.export).catch((e) => toast(e.message, true))));
    $('btn-theme').addEventListener('click', () => {
        const t = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark';
        document.documentElement.dataset.theme = t;
        try { localStorage.setItem('studio-theme', t); } catch { /* ignore */ }
    });
    window.addEventListener('keydown', (e) => {
        if (/INPUT|SELECT|TEXTAREA/.test(document.activeElement?.tagName) || e.ctrlKey || e.metaKey || e.altKey) return;
        if (e.key === ' ') { e.preventDefault(); setRunning(!app.running); }
        else if (e.key.toLowerCase() === 's') $('btn-step').click();
        else if (e.key.toLowerCase() === 'r') reset();
        else if (e.key === '1' || e.key === '2' || e.key === '3') $('layout').querySelectorAll('button')[+e.key - 1].click();
    });
}

async function doExport(kind) {
    const name = `${app.model.id}_${app.n}_step${app.sim.step}`;
    if (kind === 'png3d') return exportCanvas($('c3d'), render3D, `${name}_3d.png`);
    if (kind === 'pngslice') return exportCanvas($('c2d'), renderSlice, `${name}_slice.png`);
    if (kind === 'vtk') {
        busy(true, 'Reading fields from the GPU…');
        try { await exportVTK(app.sim, name, app.params.dx || 1); } finally { busy(false); }
        return toast('Saved VTK volume (open it in ParaView)');
    }
    if (kind === 'slicecsv') {
        const [field, comp] = app.view.source || [app.resolved.fields[0].name, 0];
        await exportSliceCSV(app.sim, field, comp, app.plane, name);
        return toast(`Saved slice values of ${field}${app.resolved.fields.find((f) => f.name === field)?.comps > 1 ? `[${comp}]` : ''}`);
    }
}

// ---------------------------------------------------------------------------
// Start
// ---------------------------------------------------------------------------

async function start() {
    busy(true, 'Starting WebGPU…');
    try {
        await initGPU();
    } catch (e) {
        showNoGPU(`${e.message} This simulator needs WebGPU: use a recent Chrome, Edge or Opera (113+), Safari 18+, or Firefox 141+ on Windows, with hardware acceleration enabled.`);
        return;
    }
    buildModelTiles();
    initDisplayPanel();
    initSlicePanel();
    initInteraction();
    initChrome();
    let first = MODELS[0].id;
    try { const saved = localStorage.getItem('pf-model'); if (MODELS.some((m) => m.id === saved)) first = saved; } catch { /* ignore */ }
    const q = new URLSearchParams(location.search).get('model');
    if (MODELS.some((m) => m.id === q)) first = q;
    selectModel(first);
    requestAnimationFrame(frame);
}

start();

window.pf = app;
