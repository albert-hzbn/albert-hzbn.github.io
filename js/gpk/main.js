// GPK Web: application controller. The simulation runs in worker.js; this
// module owns the panel, the 2D rendering, the statistics chart and exports.
import { MODELS, defaults } from './models.js';
import { COLORMAPS, gradientCSS } from '../phasefield/ui/colormaps.js';

const $ = (id) => document.getElementById(id);
const HALF_PI = Math.PI / 2;

const VIEWS = {
    grains: { label: 'Grains (labels)', note: 'Each label has its own colour. With α blending on, interface cells mix the resident and competing colours by α.' },
    phase: { label: 'Phases', note: 'Colour by phase label.' },
    orient: { label: 'Orientation g', note: 'Hue follows the in-plane orientation angle g in [0, 90°).' },
    conc: { label: 'Composition c', cmap: 0, note: 'Solute content of each cell.' },
    kappa: { label: 'Curvature κ (µm⁻¹)', cmap: 3, diverging: true, iface: true, note: 'Curvature of the resident label from the Gaussian-projected Hessian; positive for a convex grain.' },
    vel: { label: 'Velocity vₙ (µm s⁻¹)', cmap: 3, diverging: true, iface: true, note: 'Projected normal velocity; negative means the resident label is retreating.' },
    alpha: { label: 'Sub-cell position α', cmap: 0, iface: true, note: 'Fraction of the interface cell still held by its resident label.' },
};

const PHASE_COLORS = {
    solidification: [[0.86, 0.89, 0.93], [0.16, 0.36, 0.72]],
    precipitation: [[0.80, 0.81, 0.83], [0.20, 0.45, 0.82], [0.89, 0.45, 0.18]],
    wulff_shape: [[0.86, 0.87, 0.89], [0.10, 0.55, 0.55]],
};
const PHASE_NAMES = { solidification: ['liquid', 'solid'], precipitation: ['α matrix', 'β', 'γ'], wulff_shape: ['matrix', 'grain'] };
const SERIES_COLORS = ['#2f6fdb', '#e0662f', '#1f9d63'];

const app = {
    worker: null, model: null, params: {}, init: {}, n: 128, dx: 1, seed: 1, tEnd: 1,
    view: 'grains', cmap: 0, lo: 0, hi: 1, autoRange: true,
    boundary: 'projected', smooth: false, edges: false, normals: false, logChart: false,
    running: false, frame: null, history: [], needsReset: false,
    rate: { steps: 0, last: 0, since: performance.now(), value: 0 },
    hover: null, dirty: true,
};

// ---------------------------------------------------------------------------
// Helpers
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

function fmtTime(t) {
    if (!(t > 0)) return '0 s';
    if (t >= 1) return `${+t.toPrecision(4)} s`;
    if (t >= 1e-3) return `${+(t * 1e3).toPrecision(4)} ms`;
    return `${+(t * 1e6).toPrecision(4)} µs`;
}

function busy(on, text = 'Preparing…') {
    $('busy-text').textContent = text;
    $('busy').classList.toggle('on', on);
}

function download(name, data, type = 'application/octet-stream') {
    const blob = data instanceof Blob ? data : new Blob([data], { type });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}

function hsl(h, s, l) {
    const k = (n) => (n + h * 12) % 12, a = s * Math.min(l, 1 - l);
    const f = (n) => l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
    return [f(0), f(8), f(4)];
}

// Distinct colours for grain labels (golden-ratio hue walk).
const labelLUT = { n: 0, rgb: new Float32Array(0) };
function labelColors(maxLabel) {
    if (labelLUT.n > maxLabel) return labelLUT.rgb;
    const n = Math.max(maxLabel + 1, 2 * labelLUT.n, 64);
    const rgb = new Float32Array(3 * n);
    for (let l = 0; l < n; l++) {
        const c = hsl((l * 0.6180339887 + 0.13) % 1, 0.42 + 0.16 * ((l * 7) % 3) / 2, 0.58 + 0.14 * ((l * 5) % 4) / 3);
        rgb.set(c, 3 * l);
    }
    labelLUT.n = n; labelLUT.rgb = rgb;
    return rgb;
}

// ---------------------------------------------------------------------------
// Worker
// ---------------------------------------------------------------------------

function startWorker() {
    app.worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    app.worker.onmessage = (e) => {
        const m = e.data;
        if (m.type === 'frame') {
            if (m.history) app.history = m.history;
            else if (m.points) app.history.push(...m.points);
            if (app.history.length > 4000) app.history = app.history.filter((_, k) => k % 2 === 0);
            const prev = app.frame;
            app.frame = m;
            if (prev && m.step > prev.step) app.rate.steps += m.step - prev.step;
            app.dirty = true;
            busy(false);
        } else if (m.type === 'stopped') {
            setRunning(false, false);
            toast(m.reason === 'end' ? `Reached t = ${fmtTime(app.tEnd)}` : 'The process is complete.');
        } else if (m.type === 'note') {
            if (m.text) toast(m.text);
        } else if (m.type === 'error') {
            setRunning(false, false);
            busy(false);
            toast(m.message, true);
        }
    };
    app.worker.onerror = (e) => { busy(false); toast(`Worker error: ${e.message}`, true); };
}

const send = (msg) => app.worker.postMessage(msg);

// ---------------------------------------------------------------------------
// Model panel
// ---------------------------------------------------------------------------

function buildModelTiles() {
    const box = $('model-tiles');
    box.innerHTML = MODELS.map((m) => `<button type="button" data-model="${m.id}" aria-pressed="false" title="${m.name}"><svg viewBox="0 0 24 24">${m.icon}</svg><span>${m.short}</span></button>`).join('');
    box.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => selectModel(b.dataset.model)));
}

function paramControl(p, value, onInput) {
    if (p.type === 'check') {
        const lab = document.createElement('label');
        lab.className = 'check';
        lab.innerHTML = `<input type="checkbox" ${value ? 'checked' : ''}> <span>${p.label}</span>`;
        const cb = lab.querySelector('input');
        cb.addEventListener('change', () => onInput(cb.checked));
        lab.set = (v) => { cb.checked = !!v; };
        return lab;
    }
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

let katexReady = null;
const loadKatex = () => (katexReady ??= import('https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.mjs')
    .then((mod) => mod.default).catch(() => null));

function renderEquations(box, m) {
    box.innerHTML = m.equations.map(() => '<div class="eq-line"></div>').join('');
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
    $('model-short').textContent = `${m.params.length} parameters · 2D`;
    $('model-desc').innerHTML = m.description;
    renderEquations($('model-eq'), m);
    document.querySelectorAll('#model-tiles button').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.model === m.id)));

    const pre = $('presets');
    pre.innerHTML = (m.presets || []).map((p, k) => `<button type="button" data-k="${k}" aria-pressed="false">${p.label}</button>`).join('');
    pre.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => applyPreset(m.presets[+b.dataset.k], b)));

    const initBox = $('init-fields');
    initBox.innerHTML = '';
    controls.init = {};
    for (const p of m.init) {
        const lab = document.createElement('label');
        lab.className = 'field';
        lab.innerHTML = `<span>${p.label}</span><input type="number" min="${p.min}" max="${p.max}" step="${p.step}" value="${app.init[p.key]}">`;
        const inp = lab.querySelector('input');
        inp.addEventListener('change', () => {
            const v = parseFloat(inp.value);
            if (Number.isFinite(v)) { app.init[p.key] = Math.min(p.max, Math.max(p.min, v)); inp.value = app.init[p.key]; markNeedsReset(); }
        });
        lab.set = (v) => { inp.value = v; };
        controls.init[p.key] = lab;
        initBox.appendChild(lab);
    }

    const groups = $('param-groups');
    groups.innerHTML = '';
    controls.params = {};
    const byGroup = {};
    for (const p of m.params) (byGroup[p.group || 'Parameters'] ||= []).push(p);
    for (const [g, list] of Object.entries(byGroup)) {
        const t = document.createElement('div');
        t.className = 'group-title';
        t.textContent = g;
        groups.appendChild(t);
        const box = document.createElement('div');
        box.className = 'param-group';
        for (const p of list) {
            const row = paramControl(p, app.params[p.key], (v) => { app.params[p.key] = v; send({ type: 'params', params: app.params }); });
            controls.params[p.key] = row;
            box.appendChild(row);
        }
        groups.appendChild(box);
    }
    syncGridControls();
    $('click-hint').textContent = m.seed ? 'Click to seed' : '';
    $('click-hint').hidden = !m.seed;
}

function syncGridControls() {
    $('grid').value = String(app.n);
    $('dx').value = app.dx;
    $('seed').value = app.seed;
    $('tend').value = +app.tEnd.toPrecision(4);
    updateGridNote();
}

function updateGridNote() {
    const n = +$('grid').value, dx = +$('dx').value;
    $('grid-note').innerHTML = `Domain <b>${+(n * dx).toPrecision(4)} µm</b> wide, ${(n * n).toLocaleString()} cells, periodic in x and y. Gaussian width ℓ<sub>G</sub> = 2Δx = ${+(2 * dx).toPrecision(3)} µm.`;
}

function applyPreset(preset, button) {
    $('presets').querySelectorAll('button').forEach((b) => b.setAttribute('aria-pressed', String(b === button)));
    for (const [k, v] of Object.entries(preset.params || {})) { app.params[k] = v; controls.params[k]?.set(v); }
    for (const [k, v] of Object.entries(preset.init || {})) { app.init[k] = v; controls.init[k]?.set(v); }
    if (preset.n) app.n = preset.n; else app.n = app.model.n;
    syncGridControls();
    reset();
    toast(`Preset: ${preset.label}`);
}

function markNeedsReset() {
    app.needsReset = true;
    $('btn-apply').classList.add('pulse');
    updateGridNote();
}

function selectModel(id) {
    const m = MODELS.find((x) => x.id === id);
    if (!m) return;
    app.model = m;
    app.params = defaults(m.params);
    app.init = defaults(m.init);
    app.n = m.n; app.dx = m.dx; app.tEnd = m.tEnd;
    renderModelPanel();
    fillViews();
    setView(m.views[0]);
    reset();
    try { localStorage.setItem('gpk-model', id); } catch { /* ignore */ }
}

function reset() {
    setRunning(false, false);
    app.n = +$('grid').value || app.n;
    app.dx = Math.max(0.001, parseFloat($('dx').value) || app.dx);
    app.seed = parseInt($('seed').value, 10) || 1;
    app.tEnd = Math.max(1e-12, parseFloat($('tend').value) || app.tEnd);
    app.needsReset = false;
    $('btn-apply').classList.remove('pulse');
    app.frame = null;
    app.history = [];
    busy(true, `Initialising ${app.model.short} on ${app.n}²…`);
    send({ type: 'init', model: app.model.id, params: app.params, init: app.init, n: app.n, dx: app.dx, seed: app.seed, tEnd: app.tEnd });
    send({ type: 'speed', steps: +$('speed').value });
}

// ---------------------------------------------------------------------------
// Display panel
// ---------------------------------------------------------------------------

function fillViews() {
    $('view').innerHTML = app.model.views.map((v) => `<option value="${v}">${VIEWS[v].label}</option>`).join('');
}

function setView(id) {
    app.view = id;
    $('view').value = id;
    const v = VIEWS[id];
    if (v.cmap !== undefined) app.cmap = v.cmap;
    $('cmap').value = String(app.cmap);
    const hasMap = v.cmap !== undefined;
    for (const el of ['cmap-row', 'cbar', 'auto-row', 'range-row']) $(el).hidden = !hasMap;
    $('view-note').textContent = v.note;
    $('hud-view').textContent = v.label;
    app.autoRange = true;
    $('auto-range').checked = true;
    updateColorbar();
    app.dirty = true;
}

function updateColorbar() {
    $('cbar-grad').style.background = gradientCSS(app.cmap);
    $('cbar-lo').textContent = +app.lo.toPrecision(3);
    $('cbar-hi').textContent = +app.hi.toPrecision(3);
    if (document.activeElement !== $('range-lo')) $('range-lo').value = +app.lo.toPrecision(4);
    if (document.activeElement !== $('range-hi')) $('range-hi').value = +app.hi.toPrecision(4);
}

function initDisplayPanel() {
    $('view').addEventListener('change', (e) => setView(e.target.value));
    $('cmap').addEventListener('change', (e) => { app.cmap = +e.target.value; updateColorbar(); app.dirty = true; });
    $('auto-range').addEventListener('change', (e) => { app.autoRange = e.target.checked; app.dirty = true; });
    for (const id of ['range-lo', 'range-hi']) $(id).addEventListener('change', () => {
        const lo = parseFloat($('range-lo').value), hi = parseFloat($('range-hi').value);
        if (Number.isFinite(lo) && Number.isFinite(hi) && hi > lo) {
            app.lo = lo; app.hi = hi; app.autoRange = false; $('auto-range').checked = false; updateColorbar(); app.dirty = true;
        }
    });
    for (const k of ['smooth', 'edges', 'normals']) $(k).addEventListener('change', (e) => { app[k] = e.target.checked; app.dirty = true; });
    $('boundary-mode').addEventListener('change', (e) => { app.boundary = e.target.value; app.dirty = true; });
    $('log-chart').addEventListener('change', (e) => { app.logChart = e.target.checked; app.dirty = true; });
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

const off = document.createElement('canvas');
const offCtx = off.getContext('2d');
let image = null;

function scalarField(f, view) {
    const N = f.n * f.n;
    if (view === 'conc') return { data: f.c, scale: 1, all: true };
    if (view === 'kappa') return { data: f.kappa, scale: 1e-6 };
    if (view === 'vel') return { data: f.vn, scale: 1e6 };
    if (view === 'alpha') return { data: f.alpha, scale: 1 };
    return { data: new Float32Array(N), scale: 1 };
}

// Colour of each cell's resident label (own) and of its competitor (comp).
function cellColors(f) {
    const m = app.model, n = f.n, N = n * n;
    const own = new Float32Array(3 * N), comp = new Float32Array(3 * N);
    const phaseView = app.view === 'phase' || (VIEWS[app.view].iface && PHASE_COLORS[m.id]);
    const deformed = [0.62, 0.63, 0.66];
    const zeroIsMatrix = m.id === 'recrystallisation';
    const orient = app.view === 'orient';
    const lut = labelColors(f.maxLabel);
    const pc = PHASE_COLORS[m.id];
    const byLabel = (l) => (phaseView && pc ? pc[Math.min(l, pc.length - 1)] : [lut[3 * l], lut[3 * l + 1], lut[3 * l + 2]]);
    const byCell = (j) => {
        if (zeroIsMatrix && f.phi[j] === 0) return deformed;
        return orient ? hsl((f.g[j] / HALF_PI) % 1, 0.55, 0.58) : byLabel(f.phi[j]);
    };
    for (let i = 0; i < N; i++) {
        const c1 = byCell(i);
        own.set(c1, 3 * i);
        let c2 = c1;
        if (f.mask[i] && f.comp[i] !== f.phi[i]) {
            if (orient) {
                const x = i % n, y = (i / n) | 0;
                for (const j of [y * n + (x + 1) % n, y * n + (x + n - 1) % n, ((y + 1) % n) * n + x, ((y + n - 1) % n) * n + x])
                    if (f.phi[j] === f.comp[i]) { c2 = byCell(j); break; }
            } else {
                c2 = zeroIsMatrix && f.comp[i] === 0 ? deformed : byLabel(f.comp[i]);
            }
        }
        comp.set(c2, 3 * i);
    }
    return { own, comp };
}

// Cell-resolution colours with interface cells blended by alpha.
function blendedColors(f, cols) {
    const N = f.n * f.n, out = new Float32Array(3 * N);
    for (let i = 0; i < N; i++) {
        const a = f.mask[i] && f.comp[i] !== f.phi[i] ? Math.min(1, Math.max(0, f.alpha[i])) : 1;
        for (let k = 0; k < 3; k++) out[3 * i + k] = a * cols.own[3 * i + k] + (1 - a) * cols.comp[3 * i + k];
    }
    return out;
}

// Supersampled Gaussian-projected label map. Each pixel takes the label
// with the largest Gaussian-weighted vote of the 5x5 nearby cells; an
// interface cell votes alpha for its label and 1 - alpha for its competitor.
const PROJ_SIGMA = 0.9;
function projectedImage(f, cols, U) {
    const n = f.n, m = n * U, out = new Float32Array(3 * m * m);
    const R = 2, S = 2 * R + 1;
    const wt = new Float32Array(U * U * S * S);
    for (let sy = 0; sy < U; sy++) for (let sx = 0; sx < U; sx++) {
        const u = (sx + 0.5) / U - 0.5, v = (sy + 0.5) / U - 0.5;
        for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++)
            wt[((sy * U + sx) * S + dy + R) * S + dx + R] = Math.exp(-((u - dx) ** 2 + (v - dy) ** 2) / (2 * PROJ_SIGMA ** 2));
    }
    const labs = new Int32Array(64), ws = new Float32Array(64), src = new Int32Array(64), srcComp = new Uint8Array(64);
    for (let cy = 0; cy < n; cy++) for (let cx = 0; cx < n; cx++) {
        for (let sy = 0; sy < U; sy++) for (let sx = 0; sx < U; sx++) {
            let nl = 0;
            const base = (sy * U + sx) * S * S;
            for (let dy = -R; dy <= R; dy++) {
                const row = ((cy + dy + n) % n) * n;
                for (let dx = -R; dx <= R; dx++) {
                    const j = row + (cx + dx + n) % n, w = wt[base + (dy + R) * S + dx + R];
                    const split = f.mask[j] && f.comp[j] !== f.phi[j];
                    const a = split ? Math.min(1, Math.max(0, f.alpha[j])) : 1;
                    for (let part = 0; part < (split ? 2 : 1); part++) {
                        const l = part ? f.comp[j] : f.phi[j], wv = part ? (1 - a) * w : a * w;
                        let k = 0;
                        while (k < nl && labs[k] !== l) k++;
                        if (k === nl) { labs[nl] = l; ws[nl] = 0; src[nl] = j; srcComp[nl] = part; nl++; }
                        ws[k] += wv;
                    }
                }
            }
            let b1 = 0, b2 = -1;
            for (let k = 1; k < nl; k++) {
                if (ws[k] > ws[b1]) { b2 = b1; b1 = k; } else if (b2 < 0 || ws[k] > ws[b2]) b2 = k;
            }
            const o = 3 * ((cy * U + sy) * m + cx * U + sx);
            const c1 = srcComp[b1] ? cols.comp : cols.own, j1 = 3 * src[b1];
            if (b2 < 0) { out[o] = c1[j1]; out[o + 1] = c1[j1 + 1]; out[o + 2] = c1[j1 + 2]; continue; }
            const t = ws[b1] / (ws[b1] + ws[b2]);
            const s = Math.min(1, (t - 0.5) / 0.12) * 0.5 + 0.5;   // anti-alias over about one pixel
            const c2 = srcComp[b2] ? cols.comp : cols.own, j2 = 3 * src[b2];
            for (let k = 0; k < 3; k++) out[o + k] = s * c1[j1 + k] + (1 - s) * c2[j2 + k];
        }
    }
    return out;
}

function renderImage() {
    const f = app.frame;
    if (!f) return;
    const n = f.n, N = n * n;
    const v = VIEWS[app.view];
    let rgb, U = 1;
    if (v.cmap === undefined) {
        const cols = cellColors(f);
        if (app.boundary === 'projected') {
            U = n <= 96 ? 6 : n <= 128 ? 4 : n <= 192 ? 3 : 2;
            rgb = projectedImage(f, cols, U);
        } else {
            rgb = app.boundary === 'alpha' ? blendedColors(f, cols) : cols.own;
        }
    } else {
        const { data, scale, all } = scalarField(f, app.view);
        const use = (i) => all || f.mask[i];
        if (app.autoRange) {
            let lo = Infinity, hi = -Infinity;
            for (let i = 0; i < N; i++) if (use(i)) { const x = data[i] * scale; if (x < lo) lo = x; if (x > hi) hi = x; }
            if (!(hi >= lo)) { lo = 0; hi = 1; }
            if (v.diverging) {
                // Symmetric range from the 98th percentile of |x|, so a few
                // junction cells at the curvature limit do not wash out the rest.
                const mags = [];
                for (let i = 0; i < N; i++) if (use(i)) mags.push(Math.abs(data[i] * scale));
                mags.sort((p, q) => p - q);
                const a = Math.max(mags[Math.floor(0.98 * (mags.length - 1))] || 0, 1e-12);
                lo = -a; hi = a;
            }
            if (hi - lo < 1e-12) { hi = lo + 1e-6; }
            app.lo = lo; app.hi = hi;
            updateColorbar();
        }
        const cm = COLORMAPS[app.cmap];
        const bg = all ? null : blendedColors(f, cellColors(f));
        const dark = document.documentElement.dataset.theme === 'dark';
        const w = 0.35, b = dark ? 0.08 : 0.97;
        rgb = new Float32Array(3 * N);
        const span = app.hi - app.lo;
        for (let i = 0; i < N; i++) {
            if (use(i)) {
                const c = cm(Math.min(1, Math.max(0, (data[i] * scale - app.lo) / span)));
                rgb[3 * i] = c[0]; rgb[3 * i + 1] = c[1]; rgb[3 * i + 2] = c[2];
            } else {
                // Dim the microstructure behind interface-only fields.
                for (let k = 0; k < 3; k++) rgb[3 * i + k] = w * bg[3 * i + k] + (1 - w) * b;
            }
        }
    }
    const m = n * U;
    if (off.width !== m) { off.width = off.height = m; image = offCtx.createImageData(m, m); }
    const px = image.data;
    // Row 0 is the bottom of the domain (y up).
    for (let y = 0; y < m; y++) for (let x = 0; x < m; x++) {
        const i = y * m + x, o = 4 * ((m - 1 - y) * m + x);
        px[o] = Math.max(0, Math.min(255, rgb[3 * i] * 255));
        px[o + 1] = Math.max(0, Math.min(255, rgb[3 * i + 1] * 255));
        px[o + 2] = Math.max(0, Math.min(255, rgb[3 * i + 2] * 255));
        px[o + 3] = 255;
    }
    offCtx.putImageData(image, 0, 0);
    app.projected = U > 1;
}

function mapGeometry() {
    const c = $('map'), f = app.frame;
    const W = c.width, H = c.height, n = f ? f.n : app.n;
    const pad = 16 * devicePixelRatio;
    const size = Math.max(10, Math.min(W, H) - 2 * pad);
    return { x0: (W - size) / 2, y0: (H - size) / 2, size, cell: size / n, n };
}

function resizeCanvas(c) {
    const r = c.getBoundingClientRect();
    const w = Math.max(1, Math.round(r.width * devicePixelRatio)), h = Math.max(1, Math.round(r.height * devicePixelRatio));
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; return true; }
    return false;
}

function drawMap() {
    const c = $('map'), ctx = c.getContext('2d'), f = app.frame;
    const css = getComputedStyle(document.documentElement);
    ctx.fillStyle = css.getPropertyValue('--app-bg').trim() || '#eee';
    ctx.fillRect(0, 0, c.width, c.height);
    if (!f) return;
    const { x0, y0, size, cell, n } = mapGeometry();
    ctx.imageSmoothingEnabled = app.smooth || app.projected;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(off, x0, y0, size, size);

    if (app.edges && cell >= 2) {
        ctx.beginPath();
        for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
            const i = y * n + x, l = f.phi[i];
            const sx = x0 + x * cell, sy = y0 + (n - 1 - y) * cell;
            if (x < n - 1 && f.phi[i + 1] !== l) { ctx.moveTo(sx + cell, sy); ctx.lineTo(sx + cell, sy + cell); }
            if (y < n - 1 && f.phi[i + n] !== l) { ctx.moveTo(sx, sy); ctx.lineTo(sx + cell, sy); }
        }
        ctx.strokeStyle = 'rgba(0,0,0,0.55)';
        ctx.lineWidth = Math.max(1, devicePixelRatio * 0.8);
        ctx.stroke();
    }
    if (app.normals && cell >= 5) {
        // Outward normal of the resident label from the Gaussian projection.
        ctx.beginPath();
        for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
            const i = y * n + x;
            if (!f.mask[i] || (!f.nx[i] && !f.ny[i])) continue;
            const cx = x0 + (x + 0.5) * cell, cy = y0 + (n - 0.5 - y) * cell;
            ctx.moveTo(cx, cy);
            ctx.lineTo(cx + f.nx[i] * cell * 1.2, cy - f.ny[i] * cell * 1.2);
        }
        ctx.strokeStyle = getComputedStyle(document.documentElement).getPropertyValue('--ink').trim() || '#222';
        ctx.globalAlpha = 0.8;
        ctx.lineWidth = Math.max(1, devicePixelRatio);
        ctx.stroke();
        ctx.globalAlpha = 1;
        ctx.beginPath();
        ctx.stroke();
    }
    ctx.strokeStyle = getComputedStyle(document.documentElement).getPropertyValue('--line').trim() || '#ccc';
    ctx.lineWidth = devicePixelRatio;
    ctx.strokeRect(x0 - 0.5, y0 - 0.5, size + 1, size + 1);
    updateScalebar(cell / devicePixelRatio);
}

function updateScalebar(cellCss) {
    const umPerPx = app.dx / cellCss;
    const target = 90 * umPerPx;
    const p = 10 ** Math.floor(Math.log10(target));
    const nice = [1, 2, 5, 10].map((k) => k * p).filter((v) => v <= target).pop() || p;
    $('scalebar').querySelector('i').style.width = `${nice / umPerPx}px`;
    $('scalebar').querySelector('span').textContent = nice >= 1000 ? `${nice / 1000} mm` : `${+nice.toPrecision(3)} µm`;
}

// ---------------------------------------------------------------------------
// Chart
// ---------------------------------------------------------------------------

function drawChart() {
    const c = $('chart'), ctx = c.getContext('2d'), m = app.model;
    const css = getComputedStyle(document.documentElement);
    const ink = css.getPropertyValue('--ink-2').trim() || '#555', faint = css.getPropertyValue('--ink-3').trim() || '#999';
    const line = css.getPropertyValue('--line').trim() || '#ddd';
    const W = c.width, H = c.height, d = devicePixelRatio;
    ctx.clearRect(0, 0, W, H);
    const h = app.history;
    $('chart-title').textContent = m.chart.y;
    $('chart-legend').innerHTML = m.chart.series.length > 1
        ? m.chart.series.map((s, k) => `<span><i style="background:${s.color || SERIES_COLORS[k]}"></i>${s.label}</span>`).join('') : '';
    if (h.length < 2) {
        ctx.fillStyle = faint;
        ctx.font = `${12 * d}px Inter, sans-serif`;
        ctx.textAlign = 'center';
        ctx.fillText('Press Run to record statistics', W / 2, H / 2);
        $('chart-note').textContent = '';
        return;
    }
    const L = 52 * d, R = 12 * d, T = 10 * d, B = 30 * d;
    const log = app.logChart;
    const xs = h.map((p) => p.t);
    let xmin = log ? Math.max(xs.find((t) => t > 0) || 1e-12, 1e-30) : 0, xmax = Math.max(...xs, 1e-30);
    if (log && !(xmax > xmin)) xmax = xmin * 10;
    let ymin = Infinity, ymax = -Infinity;
    for (const s of m.chart.series) for (const p of h) { const y = p[s.key]; if (y < ymin) ymin = y; if (y > ymax) ymax = y; }
    if (m.chart.yMax !== undefined) { ymin = 0; ymax = m.chart.yMax; }
    else { ymin = Math.min(0, ymin); const pad = (ymax - ymin) * 0.06 || 1; ymax += pad; }
    const X = (t) => L + (W - L - R) * (log ? (Math.log10(Math.max(t, xmin)) - Math.log10(xmin)) / (Math.log10(xmax) - Math.log10(xmin) || 1) : t / xmax);
    const Y = (y) => H - B - (H - T - B) * (y - ymin) / (ymax - ymin || 1);

    // Grid and ticks
    ctx.font = `${10.5 * d}px Inter, sans-serif`;
    ctx.strokeStyle = line; ctx.lineWidth = d; ctx.fillStyle = faint;
    const yt = niceTicks(ymin, ymax, 5);
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const v of yt) {
        const y = Y(v);
        ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(W - R, y); ctx.stroke();
        ctx.fillText(+v.toPrecision(3), L - 6 * d, y);
    }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    const scaleT = xmax >= 1 ? [1, 's'] : xmax >= 1e-3 ? [1e3, 'ms'] : [1e6, 'µs'];
    const xt = log ? logTicks(xmin, xmax) : niceTicks(0, xmax, 5);
    for (const v of xt) {
        const x = X(v);
        if (x < L - 1 || x > W - R + 1) continue;
        ctx.beginPath(); ctx.moveTo(x, H - B); ctx.lineTo(x, H - B + 4 * d); ctx.stroke();
        ctx.fillText(+(v * scaleT[0]).toPrecision(3), x, H - B + 6 * d);
    }
    ctx.fillStyle = ink;
    ctx.fillText(`time (${scaleT[1]})`, (L + W - R) / 2, H - 13 * d);
    ctx.strokeStyle = faint;
    ctx.beginPath(); ctx.moveTo(L, T); ctx.lineTo(L, H - B); ctx.lineTo(W - R, H - B); ctx.stroke();

    m.chart.series.forEach((s, k) => {
        ctx.beginPath();
        h.forEach((p, j) => { const x = X(p.t), y = Y(p[s.key]); j ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
        ctx.strokeStyle = s.color || SERIES_COLORS[k]; ctx.lineWidth = 2 * d; ctx.lineJoin = 'round';
        ctx.stroke();
    });

    // Parabolic-law slope for grain growth: d<R>^2/dt from the second half.
    let note = '';
    if (m.id === 'grain_growth' && h.length > 8) {
        const half = h.slice(Math.floor(h.length / 2));
        const n = half.length, mx = half.reduce((a, p) => a + p.t, 0) / n, my = half.reduce((a, p) => a + p.R2, 0) / n;
        let sxy = 0, sxx = 0;
        for (const p of half) { sxy += (p.t - mx) * (p.R2 - my); sxx += (p.t - mx) ** 2; }
        if (sxx > 0) note = `d⟨R⟩²/dt ≈ ${+(sxy / sxx).toPrecision(3)} µm² s⁻¹`;
    }
    if (m.id === 'recrystallisation' && h.length > 8) {
        // JMAK exponent from the slope of ln(-ln(1-X)) against ln t, 5% < X < 95%.
        const pts = h.filter((p) => p.t > 0 && p.X > 0.05 && p.X < 0.95).map((p) => [Math.log(p.t), Math.log(-Math.log(1 - p.X))]);
        if (pts.length > 4) {
            const n = pts.length, mx = pts.reduce((a, p) => a + p[0], 0) / n, my = pts.reduce((a, p) => a + p[1], 0) / n;
            let sxy = 0, sxx = 0;
            for (const [x, y] of pts) { sxy += (x - mx) * (y - my); sxx += (x - mx) ** 2; }
            if (sxx > 0) note = `Avrami exponent n ≈ ${(sxy / sxx).toFixed(2)}`;
        }
    }
    $('chart-note').textContent = note;
}

function niceTicks(lo, hi, count) {
    const span = hi - lo;
    if (!(span > 0)) return [lo];
    const step0 = span / count, p = 10 ** Math.floor(Math.log10(step0));
    const step = [1, 2, 2.5, 5, 10].map((k) => k * p).find((s) => s >= step0) || p * 10;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return out;
}

function logTicks(lo, hi) {
    const out = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) out.push(10 ** e);
    return out;
}

// ---------------------------------------------------------------------------
// HUD, status and hover
// ---------------------------------------------------------------------------

function updateStatus() {
    const f = app.frame, m = app.model;
    $('st-model').textContent = m.name;
    $('hud-model').textContent = m.name;
    $('st-grid').textContent = `${app.n}² grid · Δx = ${app.dx} µm`;
    if (!f) return;
    $('st-step').textContent = f.step.toLocaleString();
    $('st-time').textContent = fmtTime(f.time);
    $('st-rate').textContent = app.running ? `${Math.round(app.rate.value).toLocaleString()} steps/s` : 'paused';
    const lines = f.stats ? m.hud(f.stats) : [];
    $('hud-stats').innerHTML = [`t = <b>${fmtTime(f.time)}</b>`, ...lines].join('<br>');
}

function updateHover() {
    const el = $('hover'), f = app.frame, h = app.hover;
    if (!f || !h) { el.hidden = true; return; }
    const i = h.y * f.n + h.x, m = app.model;
    const names = PHASE_NAMES[m.id];
    const label = names ? `${names[f.phi[i]] ?? f.phi[i]} (${f.phi[i]})` : (m.id === 'recrystallisation' && f.phi[i] === 0 ? 'deformed (0)' : `grain ${f.phi[i]}`);
    const rows = [`cell (${h.x}, ${h.y})`, label];
    if (f.mask[i]) {
        const comp = names ? (names[f.comp[i]] ?? f.comp[i]) : f.comp[i];
        rows.push(`α = ${f.alpha[i].toFixed(3)} → ${comp}`);
        rows.push(`κ = ${(f.kappa[i] * 1e-6).toFixed(4)} µm⁻¹`);
        rows.push(`vₙ = ${(f.vn[i] * 1e6).toPrecision(3)} µm/s`);
    }
    if (m.id === 'solidification' || m.id === 'precipitation') rows.push(`c = ${f.c[i].toFixed(4)}`);
    if (m.id === 'grain_growth' || m.id === 'recrystallisation') rows.push(`g = ${(f.g[i] * 180 / Math.PI).toFixed(1)}°`);
    el.innerHTML = rows.join('<br>');
    el.hidden = false;
}

function cellAt(e) {
    const c = $('map'), r = c.getBoundingClientRect(), f = app.frame;
    if (!f) return null;
    const { x0, y0, cell, n } = mapGeometry();
    const px = (e.clientX - r.left) * devicePixelRatio, py = (e.clientY - r.top) * devicePixelRatio;
    const x = Math.floor((px - x0) / cell), yv = Math.floor((py - y0) / cell);
    if (x < 0 || x >= n || yv < 0 || yv >= n) return null;
    return { x, y: n - 1 - yv };
}

function initInteraction() {
    const c = $('map');
    c.addEventListener('pointermove', (e) => { app.hover = cellAt(e); updateHover(); });
    c.addEventListener('pointerleave', () => { app.hover = null; updateHover(); });
    c.addEventListener('click', (e) => {
        const p = cellAt(e);
        if (!p) return;
        if (!app.model.seed) { toast('This model has no seeding; change the settings and reset instead.'); return; }
        send({ type: 'seed', x: p.x, y: p.y });
    });
}

// ---------------------------------------------------------------------------
// Main loop
// ---------------------------------------------------------------------------

function frame() {
    requestAnimationFrame(frame);
    const sized = resizeCanvas($('map')) | resizeCanvas($('chart'));
    if (app.dirty || sized) {
        if (app.dirty) renderImage();
        drawMap();
        drawChart();
        updateStatus();
        updateHover();
        app.dirty = false;
    }
    const now = performance.now();
    if (now - app.rate.since > 800) {
        app.rate.value = app.rate.steps * 1000 / (now - app.rate.since);
        app.rate.steps = 0;
        app.rate.since = now;
        if (app.running) updateStatus();
    }
}

// ---------------------------------------------------------------------------
// Chrome: transport, tabs, menus, export, keyboard, theme
// ---------------------------------------------------------------------------

function setRunning(on, notify = true) {
    app.running = on;
    $('btn-run').setAttribute('aria-pressed', String(on));
    $('btn-run').querySelector('span').textContent = on ? 'Pause' : 'Run';
    if (notify) send({ type: 'run', on });
    if (on && app.needsReset) toast('Grid and initial settings apply only after “Reset with these settings”.');
    if (on && app.frame && app.frame.time >= app.tEnd * (1 - 1e-9)) {
        app.tEnd *= 2;
        $('tend').value = +app.tEnd.toPrecision(4);
        send({ type: 'tEnd', tEnd: app.tEnd });
        toast(`Extended the stop time to ${fmtTime(app.tEnd)}`);
    }
    app.rate = { steps: 0, since: performance.now(), value: app.rate.value };
    updateStatus();
}

function initChrome() {
    $('btn-run').addEventListener('click', () => setRunning(!app.running));
    $('btn-step').addEventListener('click', () => { if (app.running) setRunning(false); send({ type: 'step', steps: 10 }); });
    $('btn-reset').addEventListener('click', () => reset());
    $('btn-apply').addEventListener('click', () => reset());
    $('speed').addEventListener('change', (e) => send({ type: 'speed', steps: +e.target.value }));
    $('grid').addEventListener('change', markNeedsReset);
    $('dx').addEventListener('change', markNeedsReset);
    $('seed').addEventListener('change', markNeedsReset);
    $('tend').addEventListener('change', () => {
        const v = parseFloat($('tend').value);
        if (v > 0) { app.tEnd = v; send({ type: 'tEnd', tEnd: v }); }
    });
    document.querySelectorAll('.tabs [data-tab]').forEach((b) => b.addEventListener('click', () => {
        document.querySelectorAll('.tabs [data-tab]').forEach((x) => x.setAttribute('aria-selected', String(x === b)));
        document.querySelectorAll('[data-body]').forEach((x) => { x.hidden = x.dataset.body !== b.dataset.tab; });
    }));
    $('panel-toggle').addEventListener('click', () => document.body.classList.toggle('panel-open'));
    document.querySelectorAll('.menu').forEach((m) => m.querySelector('[data-menu]').addEventListener('click', (e) => { e.stopPropagation(); m.classList.toggle('open'); }));
    document.addEventListener('click', () => document.querySelectorAll('.menu.open').forEach((o) => o.classList.remove('open')));
    document.querySelectorAll('[data-export]').forEach((b) => b.addEventListener('click', () => doExport(b.dataset.export)));
    $('btn-theme').addEventListener('click', () => {
        const t = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark';
        document.documentElement.dataset.theme = t;
        try { localStorage.setItem('studio-theme', t); } catch { /* ignore */ }
        app.dirty = true;
    });
    window.addEventListener('keydown', (e) => {
        if (/INPUT|SELECT|TEXTAREA/.test(document.activeElement?.tagName) || e.ctrlKey || e.metaKey || e.altKey) return;
        if (e.key === ' ') { e.preventDefault(); setRunning(!app.running); }
        else if (e.key.toLowerCase() === 's') $('btn-step').click();
        else if (e.key.toLowerCase() === 'r') reset();
    });
}

function doExport(kind) {
    const f = app.frame, m = app.model;
    if (!f) return;
    const name = `gpk_${m.id}_${f.n}_step${f.step}`;
    if (kind === 'png') {
        const scale = Math.max(1, Math.floor(1024 / off.width));
        const c = document.createElement('canvas');
        c.width = c.height = off.width * scale;
        const ctx = c.getContext('2d');
        ctx.imageSmoothingEnabled = app.smooth || app.projected;
        ctx.drawImage(off, 0, 0, c.width, c.height);
        c.toBlob((b) => b && download(`${name}.png`, b, 'image/png'), 'image/png');
    } else if (kind === 'chartpng') {
        $('chart').toBlob((b) => b && download(`${name}_chart.png`, b, 'image/png'), 'image/png');
    } else if (kind === 'csv') {
        const rows = ['x_um,y_um,label,alpha,competitor,interface,composition,orientation_rad,curvature_per_um,velocity_um_per_s'];
        const n = f.n;
        for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
            const i = y * n + x;
            rows.push([((x + 0.5) * app.dx).toFixed(4), ((y + 0.5) * app.dx).toFixed(4), f.phi[i], f.alpha[i].toFixed(5), f.comp[i], f.mask[i],
                f.c[i].toPrecision(6), f.g[i].toFixed(5), (f.kappa[i] * 1e-6).toPrecision(5), (f.vn[i] * 1e6).toPrecision(5)].join(','));
        }
        download(`${name}.csv`, `# GPK Web, ${m.name}, t = ${f.time} s, step ${f.step}\n` + rows.join('\n'), 'text/csv');
    } else if (kind === 'history') {
        const keys = Object.keys(app.history[0] || { t: 0 });
        download(`gpk_${m.id}_statistics.csv`, [keys.join(','), ...app.history.map((p) => keys.map((k) => p[k]).join(','))].join('\n'), 'text/csv');
    } else if (kind === 'ini') {
        const P = {};
        for (const p of m.params) P[p.key] = p.type === 'check' ? !!app.params[p.key] : app.params[p.key] * (p.scale || 1);
        const deck = m.deck({ n: app.n, dx: app.dx * 1e-6, tEnd: app.tEnd, seed: app.seed }, P, app.init);
        download(`${m.id}.ini`, deck, 'text/plain');
        toast('Saved a GPK input deck');
    }
}

// ---------------------------------------------------------------------------
// Start
// ---------------------------------------------------------------------------

function start() {
    startWorker();
    buildModelTiles();
    initDisplayPanel();
    initInteraction();
    initChrome();
    let first = MODELS[0].id;
    try { const saved = localStorage.getItem('gpk-model'); if (MODELS.some((m) => m.id === saved)) first = saved; } catch { /* ignore */ }
    const q = new URLSearchParams(location.search).get('model');
    if (MODELS.some((m) => m.id === q)) first = q;
    selectModel(first);
    const view = new URLSearchParams(location.search).get('view');
    if (app.model.views.includes(view)) setView(view);
    if (new URLSearchParams(location.search).has('run')) setTimeout(() => setRunning(true), 300);
    requestAnimationFrame(frame);
}

start();

window.gpk = app;
