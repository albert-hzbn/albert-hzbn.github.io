// Structure Builder application: UI wiring, state, undo/redo, rendering.
import {
    Structure, supercell, transformCell, cellParameters, cellVolume, latticeFromParameters, vecmat, inv3,
    add, sub, norm, dot, cross, unit, removeDuplicates, NeighborGrid, reciprocal,
} from './core.js';
import { element, normalizeSymbol, BY_SYMBOL } from './elements.js';
import {
    loadSpaceGroups, getSettings, settingsForNumber, defaultSetting, crystalSystem, buildCrystal, PRESETS, presetParameters,
} from './symmetry.js';
import {
    buildSlab, addVacuum, listTiltBoundaries, buildTiltGB, SHAPES, buildNanoparticle, buildWulff,
    buildPolycrystal, buildSolidSolution, buildStackingFault,
} from './builders.js';
import { findBonds, coordinationNumbers, radialDistribution } from './analysis.js';
import { parseStructure, WRITERS } from './io.js';
import { Viewer } from './viewer.js';

const $ = (id) => document.getElementById(id);
const num = (id) => parseFloat($(id).value);
const int = (id) => parseInt($(id).value, 10);
const fmt = (v, d = 3) => (Math.abs(v) < 1e-9 ? 0 : v).toFixed(d);

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

const state = {
    structure: new Structure(),
    undo: [],
    redo: [],
    selection: new Set(),
    mode: 'select',
    pending: [],            // picked instances for the current measurement
    measurements: [],
    customColors: {},
    hidden: new Set(),
    bonds: [],
    cn: [],
    display: [],
    spaceGroupLabel: '',
    view: {
        style: 'ball', atomScale: 1, bondRadius: 0.14, tolerance: 1.15, bonds: true, labels: false,
        color: 'element', cell: true, boundary: true, reps: [1, 1, 1],
        plane: null, direction: null,
    },
};

const viewer = new Viewer($('viewport'));

// ---------------------------------------------------------------------------
// Utilities: toasts, busy overlay, downloads, local storage
// ---------------------------------------------------------------------------

let toastTimer;
function toast(msg, isError = false) {
    const t = $('toast');
    t.textContent = msg;
    t.classList.toggle('error', isError);
    t.classList.add('show');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => t.classList.remove('show'), isError ? 5200 : 2600);
}

function withBusy(label, fn) {
    $('busy-text').textContent = label;
    $('busy').classList.add('on');
    return new Promise((resolve) => {
        requestAnimationFrame(() => setTimeout(() => {
            try { resolve(fn()); } catch (e) { console.error(e); toast(e.message, true); resolve(null); }
            finally { $('busy').classList.remove('on'); }
        }, 20));
    });
}

function download(name, content, type = 'text/plain') {
    const blob = content instanceof Blob ? content : new Blob([content], { type });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
}

function saveLocal() {
    try {
        if (state.structure.count <= 20000) localStorage.setItem('studio-last', WRITERS.xyz.fn(state.structure));
    } catch (e) { /* storage unavailable */ }
}

// ---------------------------------------------------------------------------
// Structure changes and undo/redo
// ---------------------------------------------------------------------------

function setStructure(s, { message = '', record = true, fit = true, keepSelection = false, spaceGroup = '' } = {}) {
    if (record) {
        state.undo.push({ s: state.structure, sg: state.spaceGroupLabel });
        if (state.undo.length > 60) state.undo.shift();
        state.redo = [];
    }
    state.structure = s;
    state.spaceGroupLabel = spaceGroup;
    if (!keepSelection) state.selection.clear();
    state.measurements = [];
    state.pending = [];
    state.hidden = new Set([...state.hidden].filter((el) => s.symbols.includes(el)));
    refreshAll({ fit });
    if (message) toast(message);
    saveLocal();
}

function undo() {
    if (!state.undo.length) return;
    state.redo.push({ s: state.structure, sg: state.spaceGroupLabel });
    const prev = state.undo.pop();
    state.structure = prev.s;
    state.spaceGroupLabel = prev.sg;
    state.selection.clear();
    refreshAll({ fit: false });
}

function redo() {
    if (!state.redo.length) return;
    state.undo.push({ s: state.structure, sg: state.spaceGroupLabel });
    const next = state.redo.pop();
    state.structure = next.s;
    state.spaceGroupLabel = next.sg;
    state.selection.clear();
    refreshAll({ fit: false });
}

function requireCrystal() {
    if (!state.structure.periodic || !state.structure.count) throw new Error('This builder needs a periodic crystal. Build or open one first.');
    return state.structure;
}

function requireCubic() {
    const s = requireCrystal();
    const p = cellParameters(s.cell);
    const ok = Math.abs(p.a - p.b) < 1e-3 && Math.abs(p.a - p.c) < 1e-3 && [p.alpha, p.beta, p.gamma].every((x) => Math.abs(x - 90) < 1e-2)
        && Math.abs(s.cell[0][1]) + Math.abs(s.cell[0][2]) + Math.abs(s.cell[1][0]) + Math.abs(s.cell[1][2]) < 1e-6;
    if (!ok) throw new Error('This builder needs a conventional cubic cell (a = b = c, 90° angles). Build one from a cubic space group first.');
    return s;
}

// ---------------------------------------------------------------------------
// Colours
// ---------------------------------------------------------------------------

const TAG_COLORS = ['#9aa5b1', '#4e79a7', '#f28e2b', '#59a14f', '#e15759', '#76b7b2', '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#86bcb6', '#d37295', '#a0cbe8', '#ffbe7d', '#8cd17d'];
const CN_COLORS = ['#6b7280', '#a855f7', '#6366f1', '#3b82f6', '#0ea5e9', '#14b8a6', '#22c55e', '#84cc16', '#eab308', '#f59e0b', '#f97316', '#ef4444', '#b91c1c', '#7f1d1d'];
const RAMP = [[0.00, [68, 1, 84]], [0.25, [59, 82, 139]], [0.5, [33, 145, 140]], [0.75, [94, 201, 98]], [1.0, [253, 231, 37]]];

function ramp(t) {
    t = Math.max(0, Math.min(1, t));
    for (let i = 1; i < RAMP.length; i++) {
        if (t <= RAMP[i][0]) {
            const [t0, c0] = RAMP[i - 1], [t1, c1] = RAMP[i];
            const u = (t - t0) / (t1 - t0);
            const c = c0.map((v, k) => Math.round(v + (c1[k] - v) * u));
            return '#' + c.map((v) => v.toString(16).padStart(2, '0')).join('');
        }
    }
    return '#fde725';
}

const elementColor = (sym) => state.customColors[sym] || element(sym).color;

// ---------------------------------------------------------------------------
// Rendering pipeline
// ---------------------------------------------------------------------------

function computeAnalysis() {
    const s = state.structure;
    state.bonds = s.count && s.count <= 250000 ? findBonds(s, { tolerance: state.view.tolerance }) : [];
    state.cn = coordinationNumbers(s, state.bonds);
}

function atomRadius(sym) {
    const r = element(sym).radius;
    const v = state.view;
    if (v.style === 'space') return r * 1.05 * v.atomScale;
    if (v.style === 'stick') return Math.max(v.bondRadius * 1.02, 0.05);
    return (0.28 + r * 0.3) * v.atomScale;
}

function atomColor(i) {
    const s = state.structure, v = state.view;
    if (v.color === 'tag') return TAG_COLORS[(s.tags[i] || 0) % TAG_COLORS.length];
    if (v.color === 'cn') return CN_COLORS[Math.min(state.cn[i] || 0, CN_COLORS.length - 1)];
    if (v.color === 'height') {
        const { lo, hi } = state.heightRange;
        return ramp((s.positions[i][2] - lo) / Math.max(1e-6, hi - lo));
    }
    return elementColor(s.symbols[i]);
}

function buildDisplay() {
    const s = state.structure, v = state.view;
    const out = [];
    const reps = s.periodic ? v.reps : [1, 1, 1];
    const fr = s.periodic ? s.fractionalPositions() : null;
    const zs = s.positions.map((p) => p[2]);
    state.heightRange = { lo: Math.min(...zs, 0), hi: Math.max(...zs, 1) };
    const limit = 600000;
    for (let i = 0; i < s.count; i++) {
        if (state.hidden.has(s.symbols[i])) continue;
        const color = atomColor(i), radius = atomRadius(s.symbols[i]);
        if (!s.periodic) { out.push({ index: i, position: s.positions[i], color, radius, ghost: false }); continue; }
        const shifts = [0, 1, 2].map((k) => {
            const list = [];
            for (let n = 0; n < reps[k]; n++) list.push(n);
            if (v.boundary && fr[i][k] < 1e-3) list.push(reps[k]);
            return list;
        });
        for (const a of shifts[0]) for (const b of shifts[1]) for (const c of shifts[2]) {
            const ghost = a >= reps[0] || b >= reps[1] || c >= reps[2];
            out.push({ index: i, position: vecmat([fr[i][0] + a, fr[i][1] + b, fr[i][2] + c], s.cell), color, radius, ghost, image: a + b + c > 0 });
            if (out.length > limit) break;
        }
    }
    return out;
}

function displayBonds(display) {
    if (!state.view.bonds || state.view.style === 'space' || display.length > 200000) return [];
    // Bonds between displayed atoms, using the same covalent criterion.
    const tmp = new Structure({ symbols: display.map((d) => state.structure.symbols[d.index]), positions: display.map((d) => d.position) });
    const bonds = findBonds(tmp, { tolerance: state.view.tolerance });
    const r = state.view.style === 'stick' ? state.view.bondRadius * 1.0 : state.view.bondRadius;
    return bonds.map((b) => ({ a: display[b.i].position, b: display[b.j].position, colorA: display[b.i].color, colorB: display[b.j].color, radius: r }));
}

function render({ fit = false } = {}) {
    const s = state.structure;
    const display = buildDisplay();
    state.display = display;
    let cell = null;
    if (s.periodic && state.view.cell) {
        const R = state.view.reps;
        cell = s.cell.map((row, k) => row.map((x) => x * R[k]));
    }
    viewer.setContent({ atoms: display, bonds: displayBonds(display), cell, bondRadius: state.view.bondRadius });
    renderHighlights();
    renderOverlays();
    if (s.periodic) viewer.setGizmoAxes(s.cell, ['a', 'b', 'c'], ['#e2574c', '#3aa35b', '#3e7fd6']);
    else viewer.setGizmoAxes([[1, 0, 0], [0, 1, 0], [0, 0, 1]], ['x', 'y', 'z'], ['#e2574c', '#3aa35b', '#3e7fd6']);
    if (fit) viewer.fit(boundsOf(display), !state.firstFitDone ? false : true);
    state.firstFitDone = true;
}

function boundsOf(display) {
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (const d of display) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], d.position[k] - d.radius); hi[k] = Math.max(hi[k], d.position[k] + d.radius); }
    const s = state.structure;
    if (s.periodic) {
        const R = state.view.reps;
        for (let i = 0; i < 8; i++) {
            const c = vecmat([(i & 1) * R[0], ((i >> 1) & 1) * R[1], ((i >> 2) & 1) * R[2]], s.cell);
            for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], c[k]); hi[k] = Math.max(hi[k], c[k]); }
        }
    }
    if (!isFinite(lo[0])) return { lo: [-5, -5, -5], hi: [5, 5, 5] };
    return { lo, hi };
}

function renderHighlights() {
    const items = [];
    const sel = state.selection;
    if (sel.size) for (const d of state.display) if (sel.has(d.index)) items.push({ position: d.position, radius: d.radius, color: '#2f7de1' });
    for (const p of state.pending) items.push({ position: p.position, radius: p.radius, color: '#e0662f' });
    viewer.setHighlights(items);
}

function renderOverlays() {
    viewer.clearOverlays();
    const s = state.structure;
    // Measurements
    for (const m of state.measurements) {
        viewer.addLine(m.points, '#e0662f', true);
        const mid = m.points.length === 2 ? m.points[0].map((v, k) => (v + m.points[1][k]) / 2) : m.points[1];
        viewer.addLabel(mid, m.text);
    }
    // Element labels for small structures
    if (state.view.labels && state.display.length <= 400) {
        for (const d of state.display) viewer.addLabel(d.position, s.symbols[d.index], 'el');
    }
    if (!s.periodic) return;
    const R = state.view.reps;
    if (state.view.plane) {
        const { hkl, d } = state.view.plane;
        const poly = planePolygon(hkl, d, R);
        if (poly.length >= 3) viewer.addPolygon(poly.map((f) => vecmat(f, s.cell)), '#2f7de1', 0.25);
        else toast('That plane does not cut the displayed cell; change d.', true);
    }
    if (state.view.direction) {
        const v = vecmat(state.view.direction, s.cell);
        viewer.addArrow([0, 0, 0], v, norm(v), '#c2410c');
        viewer.addLabel(v, `[${state.view.direction.join(' ')}]`);
    }
}

// Intersection of the plane h·f1 + k·f2 + l·f3 = d with the box [0,R]³ (fractional).
function planePolygon([h, k, l], d, R) {
    const corners = [];
    for (let i = 0; i < 8; i++) corners.push([(i & 1) * R[0], ((i >> 1) & 1) * R[1], ((i >> 2) & 1) * R[2]]);
    const edges = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
    const f = (p) => h * p[0] + k * p[1] + l * p[2] - d;
    const pts = [];
    for (const [i, j] of edges) {
        const a = corners[i], b = corners[j], fa = f(a), fb = f(b);
        if (Math.abs(fa) < 1e-9) pts.push(a);
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0)) {
            const t = fa / (fa - fb);
            pts.push(a.map((v, n) => v + (b[n] - v) * t));
        }
    }
    const uniq = [];
    for (const p of pts) if (!uniq.some((q) => Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]) < 1e-6)) uniq.push(p);
    if (uniq.length < 3) return uniq;
    // Sort around the centroid in the plane.
    const c = uniq.reduce((acc, p) => acc.map((v, n) => v + p[n] / uniq.length), [0, 0, 0]);
    const n = [h, k, l];
    const u = unit(sub(uniq[0], c));
    const w = unit(cross(n, u));
    return uniq.sort((p, q) => Math.atan2(dot(sub(p, c), w), dot(sub(p, c), u)) - Math.atan2(dot(sub(q, c), w), dot(sub(q, c), u)));
}

function refreshAll({ fit = false } = {}) {
    computeAnalysis();
    render({ fit });
    refreshPanels();
}

// ---------------------------------------------------------------------------
// Panels: HUD, status bar, legend, info, selectors
// ---------------------------------------------------------------------------

function refreshPanels() {
    const s = state.structure;
    $('hud-title').textContent = s.title || s.formula() || 'Empty structure';
    const sub1 = [s.formula(), state.spaceGroupLabel, `${s.count.toLocaleString()} atoms`].filter(Boolean);
    $('hud-sub').textContent = sub1.join(' · ');
    $('st-atoms').textContent = s.count.toLocaleString();
    $('st-formula').textContent = s.formula();
    if (s.periodic) {
        const p = cellParameters(s.cell);
        $('st-cell').textContent = `a ${fmt(p.a)} · b ${fmt(p.b)} · c ${fmt(p.c)} Å · ${fmt(p.alpha, 1)}° ${fmt(p.beta, 1)}° ${fmt(p.gamma, 1)}°`;
    } else $('st-cell').textContent = 'no cell (cluster)';
    $('btn-undo').disabled = !state.undo.length;
    $('btn-redo').disabled = !state.redo.length;
    refreshSelectionUI();
    refreshLegend();
    refreshInfo();
    refreshSpeciesSelects();
    refreshCellInputs();
    refreshCoordination();
}

function refreshLegend() {
    const s = state.structure, box = $('legend');
    box.innerHTML = '';
    const counts = s.speciesCounts();
    if (state.view.color === 'element') {
        for (const [sym, n] of Object.entries(counts)) {
            const li = document.createElement('div');
            li.className = 'li';
            li.innerHTML = `<input type="color" value="${elementColor(sym)}" aria-label="Colour of ${sym}"><span><b>${sym}</b> <span class="note">${element(sym).name}</span></span><span class="count">${n.toLocaleString()}</span><button class="eye" aria-pressed="${!state.hidden.has(sym)}" title="Show / hide ${sym}">${eyeIcon}</button>`;
            li.querySelector('input').addEventListener('input', (e) => { state.customColors[sym] = e.target.value; render(); });
            li.querySelector('.eye').addEventListener('click', () => {
                if (state.hidden.has(sym)) state.hidden.delete(sym); else state.hidden.add(sym);
                render(); refreshLegend();
            });
            box.appendChild(li);
        }
    } else if (state.view.color === 'tag') {
        const tags = [...new Set(s.tags)].sort((a, b) => a - b);
        box.innerHTML = tags.slice(0, 24).map((t) => `<div class="li"><span style="width:14px;height:14px;border-radius:50%;background:${TAG_COLORS[t % TAG_COLORS.length]}"></span><span>${t === 0 ? 'Untagged' : 'Region ' + t}</span><span class="count">${s.tags.filter((x) => x === t).length}</span><span></span></div>`).join('')
            + (tags.length > 24 ? `<p class="note">… ${tags.length - 24} more</p>` : '');
    } else if (state.view.color === 'cn') {
        const hist = {};
        state.cn.forEach((c) => { hist[c] = (hist[c] || 0) + 1; });
        box.innerHTML = Object.keys(hist).map(Number).sort((a, b) => a - b).map((c) => `<div class="li"><span style="width:14px;height:14px;border-radius:50%;background:${CN_COLORS[Math.min(c, CN_COLORS.length - 1)]}"></span><span>CN ${c}</span><span class="count">${hist[c]}</span><span></span></div>`).join('');
    } else {
        const { lo, hi } = state.heightRange;
        box.innerHTML = `<div style="height:12px;border-radius:6px;background:linear-gradient(90deg,${[0, 0.25, 0.5, 0.75, 1].map(ramp).join(',')})"></div><div class="row c2" style="font-family:var(--mono);font-size:11.5px"><span>${fmt(lo, 2)} Å</span><span style="text-align:right">${fmt(hi, 2)} Å</span></div>`;
    }
}

const eyeIcon = '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="1.8"><path d="M2 12s3.6-7 10-7 10 7 10 7-3.6 7-10 7S2 12 2 12z"/><circle cx="12" cy="12" r="3"/></svg>';

function refreshInfo() {
    const s = state.structure;
    const rows = [['Formula', s.formula() || '—'], ['Atoms', s.count.toLocaleString()]];
    if (state.spaceGroupLabel) rows.push(['Space group', state.spaceGroupLabel]);
    if (s.periodic) {
        const p = cellParameters(s.cell);
        rows.push(['a, b, c (Å)', `${fmt(p.a)}, ${fmt(p.b)}, ${fmt(p.c)}`]);
        rows.push(['α, β, γ (°)', `${fmt(p.alpha, 2)}, ${fmt(p.beta, 2)}, ${fmt(p.gamma, 2)}`]);
        rows.push(['Volume (Å³)', fmt(cellVolume(s.cell), 3)]);
        rows.push(['Density (g/cm³)', fmt(s.density(), 4)]);
        rows.push(['Volume / atom (Å³)', s.count ? fmt(cellVolume(s.cell) / s.count, 3) : '—']);
    } else if (s.count) {
        const { lo, hi } = s.bounds();
        rows.push(['Extent (Å)', [0, 1, 2].map((k) => fmt(hi[k] - lo[k], 2)).join(' × ')]);
    }
    rows.push(['Mass (u)', fmt(s.mass(), 3)]);
    rows.push(['Bonds shown', state.bonds.length.toLocaleString()]);
    for (const [sym, n] of Object.entries(s.speciesCounts())) rows.push([`${sym} fraction`, `${(100 * n / s.count).toFixed(2)} %`]);
    $('info').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
}

function refreshCoordination() {
    const s = state.structure;
    const per = {};
    s.symbols.forEach((sym, i) => { (per[sym] = per[sym] || []).push(state.cn[i] || 0); });
    $('cn-table').innerHTML = Object.entries(per).map(([sym, list]) => {
        const mean = list.reduce((a, b) => a + b, 0) / list.length;
        const hist = {};
        list.forEach((c) => { hist[c] = (hist[c] || 0) + 1; });
        const top = Object.entries(hist).sort((a, b) => b[1] - a[1]).slice(0, 3).map(([c, n]) => `${c}:${n}`).join(' ');
        return `<dt>${sym} (mean ${mean.toFixed(2)})</dt><dd>${top}</dd>`;
    }).join('') || '<dt>—</dt><dd></dd>';
}

function refreshSpeciesSelects() {
    const species = Object.keys(state.structure.speciesCounts());
    const fill = (id, withAll = false) => {
        const el = $(id), cur = el.value;
        el.innerHTML = (withAll ? '<option value="">All</option>' : '') + species.map((x) => `<option>${x}</option>`).join('');
        if ([...el.options].some((o) => o.value === cur)) el.value = cur;
    };
    fill('sel-el'); fill('rn-from'); fill('ss-host'); fill('rdf-a', true); fill('rdf-b', true);
}

function refreshCellInputs() {
    const s = state.structure;
    const ids = ['ce-a', 'ce-b', 'ce-c', 'ce-al', 'ce-be', 'ce-ga'];
    if (!s.periodic) { ids.forEach((id) => { $(id).value = ''; }); return; }
    const p = cellParameters(s.cell);
    [p.a, p.b, p.c, p.alpha, p.beta, p.gamma].forEach((v, k) => { $(ids[k]).value = +v.toFixed(k < 3 ? 5 : 4); });
    $('title-in').value = s.title || '';
}

function refreshSelectionUI() {
    const s = state.structure, sel = [...state.selection];
    $('sel-count').textContent = sel.length ? `${sel.length} selected` : 'none';
    $('st-sel').textContent = sel.length ? `${sel.length} selected` : '';
    const box = $('atom-editor');
    if (sel.length === 1) {
        const i = sel[0];
        const p = s.positions[i];
        const f = s.periodic ? s.toFractional(p) : null;
        box.innerHTML = `<p class="note">Atom #${i + 1}: <b>${s.symbols[i]}</b> (${element(s.symbols[i]).name}), CN ${state.cn[i] ?? '—'}</p>
            <div class="row c3">
              <label class="field"><span>x (Å)</span><input type="number" step="0.01" value="${fmt(p[0], 5)}" data-c="0"></label>
              <label class="field"><span>y (Å)</span><input type="number" step="0.01" value="${fmt(p[1], 5)}" data-c="1"></label>
              <label class="field"><span>z (Å)</span><input type="number" step="0.01" value="${fmt(p[2], 5)}" data-c="2"></label>
            </div>
            ${f ? `<div class="row c3">
              <label class="field"><span>a (frac)</span><input type="number" step="0.01" value="${fmt(f[0], 5)}" data-f="0"></label>
              <label class="field"><span>b (frac)</span><input type="number" step="0.01" value="${fmt(f[1], 5)}" data-f="1"></label>
              <label class="field"><span>c (frac)</span><input type="number" step="0.01" value="${fmt(f[2], 5)}" data-f="2"></label>
            </div>` : ''}`;
        box.querySelectorAll('input[data-c]').forEach((inp) => inp.addEventListener('change', () => {
            const ns = s.clone();
            ns.positions[i][+inp.dataset.c] = parseFloat(inp.value);
            setStructure(ns, { keepSelection: true, fit: false, spaceGroup: '' });
        }));
        box.querySelectorAll('input[data-f]').forEach((inp) => inp.addEventListener('change', () => {
            const ns = s.clone();
            const fr = ns.toFractional(ns.positions[i]);
            fr[+inp.dataset.f] = parseFloat(inp.value);
            ns.positions[i] = ns.toCartesian(fr);
            setStructure(ns, { keepSelection: true, fit: false, spaceGroup: '' });
        }));
    } else if (sel.length > 1) {
        const counts = {};
        sel.forEach((i) => { counts[s.symbols[i]] = (counts[s.symbols[i]] || 0) + 1; });
        box.innerHTML = `<p class="note">${sel.length} atoms: ${Object.entries(counts).map(([k, v]) => `${v} ${k}`).join(', ')}</p>`;
    } else box.innerHTML = '<p class="note">Nothing selected.</p>';
}

// ---------------------------------------------------------------------------
// Build tab: bulk crystal
// ---------------------------------------------------------------------------

const sgState = { number: 225, setting: null };

function fillPresetSelect() {
    $('cr-preset').innerHTML = '<option value="">Custom…</option>' + PRESETS.map((p) => `<option value="${p.id}">${p.label}</option>`).join('');
}

function fillSpaceGroups(system) {
    const settings = getSettings();
    const nums = [...new Set(settings.filter((s) => s.system === system).map((s) => s.n))];
    $('cr-sg').innerHTML = nums.map((n) => {
        const st = defaultSetting(n);
        return `<option value="${n}">${n} · ${st.short || st.hm}</option>`;
    }).join('');
}

function fillSettings(n) {
    const list = settingsForNumber(n);
    const def = defaultSetting(n);
    $('cr-setting').innerHTML = list.map((s) => `<option value="${s.h}" ${s === def ? 'selected' : ''}>${s.hm}${s.choice ? ` (${settingLabel(s.choice)})` : ''}</option>`).join('');
    $('cr-setting').disabled = list.length < 2;
    sgState.number = n;
}

function settingLabel(c) {
    if (c === '1' || c === '2') return `origin choice ${c}`;
    if (c === 'H') return 'hexagonal axes';
    if (c === 'R') return 'rhombohedral axes';
    return c;
}

function currentSetting() {
    const h = int('cr-setting');
    return getSettings().find((s) => s.h === h);
}

// Enforce lattice constraints of the crystal system on the inputs.
function applyLatticeConstraints() {
    const st = currentSetting();
    if (!st) return;
    const sys = st.system;
    const rh = st.choice === 'R';
    const set = (id, v, dis) => { if (v !== null) $(id).value = v; $(id).disabled = dis; };
    const a = num('cr-a');
    switch (sys) {
        case 'cubic': set('cr-b', a, true); set('cr-c', a, true); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 90, true); break;
        case 'tetragonal': set('cr-b', a, true); set('cr-c', null, false); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 90, true); break;
        case 'hexagonal':
        case 'trigonal':
            if (rh) { set('cr-b', a, true); set('cr-c', a, true); set('cr-al', null, false); set('cr-be', num('cr-al'), true); set('cr-ga', num('cr-al'), true); }
            else { set('cr-b', a, true); set('cr-c', null, false); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 120, true); }
            break;
        case 'orthorhombic': ['cr-b', 'cr-c'].forEach((id) => set(id, null, false)); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 90, true); break;
        case 'monoclinic': ['cr-b', 'cr-c', 'cr-be'].forEach((id) => set(id, null, false)); set('cr-al', 90, true); set('cr-ga', 90, true); break;
        default: ['cr-b', 'cr-c', 'cr-al', 'cr-be', 'cr-ga'].forEach((id) => set(id, null, false));
    }
}

function siteRow(sym = 'Cu', x = 0, y = 0, z = 0) {
    const tr = document.createElement('div');
    tr.className = 'tr';
    const f = (v) => +(+v).toFixed(5);
    tr.innerHTML = `<input value="${sym}" aria-label="Element"><input type="number" step="0.01" value="${f(x)}" aria-label="x"><input type="number" step="0.01" value="${f(y)}" aria-label="y"><input type="number" step="0.01" value="${f(z)}" aria-label="z"><button class="x-btn" title="Remove site" aria-label="Remove site">×</button>`;
    tr.querySelector('button').addEventListener('click', () => tr.remove());
    return tr;
}

function setSites(sites) {
    const t = $('cr-sites');
    t.innerHTML = '<div class="tr head"><span>Element</span><span>x</span><span>y</span><span>z</span><span></span></div>';
    sites.forEach((s) => t.appendChild(siteRow(...s)));
}

function readSites() {
    return [...$('cr-sites').querySelectorAll('.tr:not(.head)')].map((tr) => {
        const [e, x, y, z] = tr.querySelectorAll('input');
        const frac = (v) => {
            const s = String(v.value).trim();
            if (s.includes('/')) { const [p, q] = s.split('/'); return parseFloat(p) / parseFloat(q); }
            return parseFloat(s);
        };
        return { symbol: normalizeSymbol(e.value), x: frac(x), y: frac(y), z: frac(z) };
    }).filter((s) => [s.x, s.y, s.z].every(Number.isFinite));
}

function loadPreset(id) {
    const p = PRESETS.find((x) => x.id === id);
    if (!p) return;
    const sys = crystalSystem(p.n);
    $('cr-system').value = sys;
    fillSpaceGroups(sys);
    $('cr-sg').value = p.n;
    fillSettings(p.n);
    const prm = presetParameters(p);
    $('cr-a').value = prm.a; $('cr-b').value = prm.b; $('cr-c').value = prm.c;
    $('cr-al').value = prm.alpha; $('cr-be').value = prm.beta; $('cr-ga').value = prm.gamma;
    applyLatticeConstraints();
    setSites(p.sites);
}

function buildFromPanel({ quiet = false } = {}) {
    const setting = currentSetting();
    applyLatticeConstraints();
    const sites = readSites();
    if (!sites.length) throw new Error('Add at least one site.');
    const prm = { a: num('cr-a'), b: num('cr-b'), c: num('cr-c'), alpha: num('cr-al'), beta: num('cr-be'), gamma: num('cr-ga') };
    if (![prm.a, prm.b, prm.c].every((v) => v > 0)) throw new Error('Lattice lengths must be positive.');
    const preset = PRESETS.find((p) => p.id === $('cr-preset').value);
    const { structure, multiplicities } = buildCrystal({ setting, ...prm, sites });
    structure.title = preset ? preset.label : `${structure.formula()} (${setting.short || setting.hm})`;
    const label = `${setting.hm} (No. ${setting.n})`;
    $('cr-note').textContent = `Generated ${structure.count} atoms. Site multiplicities: ${sites.map((s, i) => `${s.symbol} ${multiplicities[i]}`).join(', ')}.`;
    setStructure(structure, { message: quiet ? '' : `Built ${structure.formula()} in ${setting.hm}`, spaceGroup: label });
}

// ---------------------------------------------------------------------------
// Build tab: other builders
// ---------------------------------------------------------------------------

let gbChoice = null;
function refreshGBList() {
    const axis = $('gb-axis').value.split(',').map(Number);
    const list = listTiltBoundaries(axis, 7, int('gb-maxsig') || 51);
    const box = $('gb-list');
    box.innerHTML = list.map((g, k) => `<button data-k="${k}" aria-pressed="${k === 0}"><span>Σ${g.sigma}</span><span>(${g.plane.join(' ')})</span><span>${g.angle.toFixed(2)}°</span></button>`).join('')
        || '<p class="note" style="padding:8px">No boundaries for this axis.</p>';
    gbChoice = list[0] || null;
    box.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
        box.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', 'false'));
        b.setAttribute('aria-pressed', 'true');
        gbChoice = list[+b.dataset.k];
    }));
}

function facetRow(h = 1, k = 1, l = 1, e = 1) {
    const tr = document.createElement('div');
    tr.className = 'tr';
    tr.innerHTML = `<input type="number" value="${h}" aria-label="h"><input type="number" value="${k}" aria-label="k"><input type="number" value="${l}" aria-label="l"><input type="number" step="0.01" value="${e}" aria-label="Surface energy"><button class="x-btn" aria-label="Remove facet">×</button>`;
    tr.querySelector('button').addEventListener('click', () => tr.remove());
    return tr;
}

function soluteRow(sym = 'Ni', pct = 25) {
    const tr = document.createElement('div');
    tr.className = 'tr';
    tr.innerHTML = `<input value="${sym}" aria-label="Solute element"><input type="number" step="0.5" value="${pct}" aria-label="Percent of host sites"><button class="x-btn" aria-label="Remove solute">×</button>`;
    tr.querySelector('button').addEventListener('click', () => tr.remove());
    return tr;
}

function initBuilders() {
    fillPresetSelect();
    $('cr-preset').addEventListener('change', (e) => { if (e.target.value) loadPreset(e.target.value); });
    $('cr-system').addEventListener('change', (e) => {
        fillSpaceGroups(e.target.value);
        fillSettings(int('cr-sg'));
        $('cr-preset').value = '';
        if (e.target.value === 'hexagonal' || e.target.value === 'trigonal') { $('cr-ga').value = 120; }
        applyLatticeConstraints();
    });
    $('cr-sg').addEventListener('change', () => { fillSettings(int('cr-sg')); $('cr-preset').value = ''; applyLatticeConstraints(); });
    $('cr-setting').addEventListener('change', applyLatticeConstraints);
    ['cr-a', 'cr-al'].forEach((id) => $(id).addEventListener('input', applyLatticeConstraints));
    $('cr-add-site').addEventListener('click', () => $('cr-sites').appendChild(siteRow('O', 0, 0, 0)));
    $('cr-build').addEventListener('click', () => { try { buildFromPanel(); } catch (e) { toast(e.message, true); } });

    // Supercell / transformation
    $('sc-apply').addEventListener('click', () => withBusy('Building supercell…', () => {
        const s = requireCrystal();
        const n = [int('sc-a'), int('sc-b'), int('sc-c')];
        if (n.some((v) => !(v >= 1))) throw new Error('Repeats must be ≥ 1.');
        const total = s.count * n[0] * n[1] * n[2];
        if (total > 1500000) throw new Error(`${total.toLocaleString()} atoms is too many for the browser.`);
        setStructure(supercell(s, ...n), { message: `Supercell ${n.join('×')}: ${total.toLocaleString()} atoms` });
    }));
    $('tm-apply').addEventListener('click', () => withBusy('Transforming…', () => {
        const s = requireCrystal();
        const v = [...$('tm').querySelectorAll('input')].map((i) => {
            const t = i.value.trim();
            if (t.includes('/')) { const [p, q] = t.split('/'); return parseFloat(p) / parseFloat(q); }
            return parseFloat(t);
        });
        const P = [v.slice(0, 3), v.slice(3, 6), v.slice(6, 9)];
        setStructure(transformCell(s, P), { message: 'Cell transformed' });
    }));
    $('tm-prim').addEventListener('click', () => {
        let lattice;
        try { lattice = detectCentring(requireCrystal()); } catch (e) { return toast(e.message, true); }
        if (!lattice) return toast('The current cell is not F- or I-centred, so it is already primitive (or needs a custom matrix).', true);
        const M = lattice === 'I'
            ? ['-1/2', '1/2', '1/2', '1/2', '-1/2', '1/2', '1/2', '1/2', '-1/2']
            : ['0', '1/2', '1/2', '1/2', '0', '1/2', '1/2', '1/2', '0'];
        $('tm').querySelectorAll('input').forEach((inp, k) => { inp.value = M[k]; });
        toast(`${lattice}-centred cell detected: matrix filled in. Press Transform to apply.`);
    });

    // Slab
    $('sl-build').addEventListener('click', () => withBusy('Cutting slab…', () => {
        const s = requireCrystal();
        const hkl = [int('sl-h'), int('sl-k'), int('sl-l')];
        let slab = buildSlab(s, hkl, int('sl-layers'), num('sl-vac'));
        const rep = $('sl-rep').value.trim().split(/\s+/).map(Number);
        if (rep.length === 2 && rep.every((x) => x >= 1) && (rep[0] > 1 || rep[1] > 1)) slab = supercell(slab, rep[0], rep[1], 1);
        setStructure(slab, { message: `(${hkl.join(' ')}) slab: ${slab.count} atoms` });
    }));

    // Grain boundary
    $('gb-axis').addEventListener('change', refreshGBList);
    $('gb-maxsig').addEventListener('change', refreshGBList);
    refreshGBList();
    $('gb-build').addEventListener('click', () => withBusy('Building bicrystal…', () => {
        const s = requireCubic();
        if (!gbChoice) throw new Error('Choose a boundary from the list.');
        const axis = $('gb-axis').value.split(',').map(Number);
        const gb = buildTiltGB(s, axis, gbChoice.plane, {
            reps: [int('gb-n1'), int('gb-n2'), int('gb-n3')], shift: [num('gb-sy') || 0, num('gb-sz') || 0], overlap: num('gb-ov'),
        });
        state.view.color = 'tag'; $('v-color').value = 'tag';
        setStructure(gb, { message: `${gb.title}: ${gb.count} atoms, ${gb.removed} overlapping atoms removed` });
    }));

    // Nanoparticle
    $('np-shape').innerHTML = Object.entries(SHAPES).map(([k, v]) => `<option value="${k}">${v.label}</option>`).join('');
    segmented($('np-mode'), (v) => { $('np-shape-opts').hidden = v !== 'shape'; $('np-wulff-opts').hidden = v !== 'wulff'; });
    const facets = $('wf-facets');
    facets.innerHTML = '<div class="tr head"><span>h</span><span>k</span><span>l</span><span>γ (J/m²)</span><span></span></div>';
    facets.appendChild(facetRow(1, 1, 1, 1.0));
    facets.appendChild(facetRow(1, 0, 0, 1.15));
    facets.appendChild(facetRow(1, 1, 0, 1.25));
    $('wf-add').addEventListener('click', () => facets.appendChild(facetRow(2, 1, 1, 1.3)));
    $('np-build').addEventListener('click', () => withBusy('Cutting nanoparticle…', () => {
        const s = requireCrystal();
        const mode = $('np-mode').querySelector('[aria-pressed="true"]').dataset.v;
        let np;
        if (mode === 'wulff') {
            const fs = [...facets.querySelectorAll('.tr:not(.head)')].map((tr) => {
                const v = [...tr.querySelectorAll('input')].map((i) => parseFloat(i.value));
                return { hkl: v.slice(0, 3), energy: v[3] };
            }).filter((f) => f.energy > 0 && f.hkl.some((x) => x));
            np = buildWulff(s, fs, { radius: num('wf-r'), cubicSymmetry: $('wf-cubic').checked, centerOn: $('np-center').value, vacuum: num('np-vac') });
        } else {
            np = buildNanoparticle(s, { shape: $('np-shape').value, radius: num('np-r'), height: num('np-h'), radiusY: num('np-ry'), centerOn: $('np-center').value, vacuum: num('np-vac') });
        }
        setStructure(np, { message: `${np.title}: ${np.count.toLocaleString()} atoms` });
    }));

    // Polycrystal
    $('pc-build').addEventListener('click', () => withBusy('Building polycrystal…', () => {
        const s = requireCrystal();
        const pc = buildPolycrystal(s, { box: [num('pc-x'), num('pc-y'), num('pc-z')], grains: int('pc-n'), seed: int('pc-seed'), overlap: num('pc-ov') });
        state.view.color = 'tag'; $('v-color').value = 'tag';
        state.view.bonds = pc.count < 60000; $('v-bonds').checked = state.view.bonds;
        setStructure(pc, { message: `${pc.title}: ${pc.count.toLocaleString()} atoms` });
    }));

    // Stacking fault
    const sfOut = () => { $('sf-u-out').value = num('sf-u').toFixed(2); $('sf-presets').querySelectorAll('button').forEach((b) => b.setAttribute('aria-pressed', String(Math.abs(+b.dataset.u - num('sf-u')) < 1e-6))); };
    $('sf-u').addEventListener('input', sfOut);
    $('sf-presets').querySelectorAll('button').forEach((b) => b.addEventListener('click', () => { $('sf-u').value = b.dataset.u; sfOut(); }));
    $('sf-build').addEventListener('click', () => withBusy('Building stacking fault…', () => {
        const s = requireCubic();
        if (s.count !== 4 || !isFcc(s)) throw new Error('Needs a 4-atom FCC conventional cell (e.g. the FCC preset).');
        const sf = buildStackingFault(s, { reps: [int('sf-nx'), int('sf-ny'), int('sf-nz')], u: num('sf-u') });
        state.view.color = 'tag'; $('v-color').value = 'tag';
        setStructure(sf, { message: `${sf.title}: ${sf.count} atoms` });
    }));

    // Solid solution
    const sol = $('ss-solutes');
    sol.innerHTML = '<div class="tr head"><span>Solute</span><span>% of host sites</span><span></span></div>';
    sol.appendChild(soluteRow('Ni', 25));
    $('ss-add').addEventListener('click', () => sol.appendChild(soluteRow('Al', 10)));
    $('ss-build').addEventListener('click', () => withBusy('Substituting…', () => {
        const s = state.structure;
        if (!s.count) throw new Error('Build or open a structure first.');
        const solutes = [...sol.querySelectorAll('.tr:not(.head)')].map((tr) => {
            const [e, p] = tr.querySelectorAll('input');
            return { symbol: normalizeSymbol(e.value), fraction: parseFloat(p.value) / 100 };
        }).filter((x) => x.fraction > 0);
        const total = solutes.reduce((a, b) => a + b.fraction, 0);
        if (total > 1) throw new Error('Solute fractions add up to more than 100%.');
        const { structure, report } = buildSolidSolution(s, $('ss-host').value, solutes, int('ss-seed'));
        setStructure(structure, { message: `Substituted ${report}`, fit: false });
    }));
}

// Detect F or I centring by checking whether the centring translations map the
// structure onto itself.
function detectCentring(s) {
    const fr = s.fractionalPositions();
    const key = (sym, f) => sym + f.map((v) => (((Math.round(v * 1000) / 1000) % 1) + 1) % 1).map((v) => (v > 0.9995 ? 0 : v).toFixed(3)).join(',');
    const set = new Set(fr.map((f, i) => key(s.symbols[i], f)));
    const maps = (t) => fr.every((f, i) => set.has(key(s.symbols[i], [f[0] + t[0], f[1] + t[1], f[2] + t[2]])));
    if ([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]].every(maps)) return 'F';
    if (maps([0.5, 0.5, 0.5])) return 'I';
    return null;
}

function isFcc(s) {
    const f = s.fractionalPositions().map((v) => v.map((x) => Math.round(x * 2) / 2 % 1));
    const keys = new Set(f.map((v) => v.join(',')));
    return ['0,0,0', '0,0.5,0.5', '0.5,0,0.5', '0.5,0.5,0'].every((k) => keys.has(k));
}

function segmented(el, onChange) {
    el.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
        el.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', String(x === b)));
        onChange(b.dataset.v, b);
    }));
}

// ---------------------------------------------------------------------------
// Edit tab
// ---------------------------------------------------------------------------

function initEdit() {
    $('sel-all').addEventListener('click', () => selectAll());
    $('sel-none').addEventListener('click', () => setSelection([]));
    $('sel-invert').addEventListener('click', () => setSelection([...Array(state.structure.count).keys()].filter((i) => !state.selection.has(i))));
    $('sel-by-el').addEventListener('click', () => {
        const el = $('sel-el').value;
        setSelection(state.structure.symbols.map((x, i) => (x === el ? i : -1)).filter((i) => i >= 0));
    });
    $('sel-grow').addEventListener('click', () => {
        const s = state.structure, r = num('sel-r');
        if (!state.selection.size) return toast('Select at least one atom first.', true);
        const grid = new NeighborGrid(s, r);
        const next = new Set(state.selection);
        for (const i of state.selection) grid.forEachNeighbor(i, r, (j) => next.add(j));
        setSelection([...next]);
    });

    $('ed-set-el').addEventListener('click', () => {
        const el = normalizeSymbol($('ed-el').value);
        if (!BY_SYMBOL[el]) return toast(`Unknown element "${$('ed-el').value}".`, true);
        if (!state.selection.size) return toast('Select atoms first.', true);
        const ns = state.structure.clone();
        for (const i of state.selection) ns.symbols[i] = el;
        setStructure(ns, { keepSelection: true, fit: false, message: `${state.selection.size} atom(s) changed to ${el}` });
    });
    $('ed-move').addEventListener('click', () => {
        if (!state.selection.size) return toast('Select atoms first.', true);
        const d = [num('ed-dx') || 0, num('ed-dy') || 0, num('ed-dz') || 0];
        const ns = state.structure.clone();
        for (const i of state.selection) ns.positions[i] = add(ns.positions[i], d);
        setStructure(ns, { keepSelection: true, fit: false });
    });
    $('ed-dup').addEventListener('click', () => {
        if (!state.selection.size) return toast('Select atoms first.', true);
        const ns = state.structure.clone();
        const d = [num('ed-dx') || 0, num('ed-dy') || 0, num('ed-dz') || 1];
        const start = ns.count;
        for (const i of state.selection) ns.push(ns.symbols[i], add(ns.positions[i], d), ns.tags[i]);
        state.selection = new Set([...Array(ns.count - start).keys()].map((k) => start + k));
        setStructure(ns, { keepSelection: true, fit: false, message: `Duplicated ${ns.count - start} atom(s), shifted by (${d.join(', ')}) Å` });
    });
    $('ed-del').addEventListener('click', deleteSelection);

    $('add-frac').addEventListener('change', (e) => {
        const f = e.target.checked;
        $('add-lx').textContent = f ? 'a' : 'x'; $('add-ly').textContent = f ? 'b' : 'y'; $('add-lz').textContent = f ? 'c' : 'z';
    });
    $('add-btn').addEventListener('click', () => {
        const el = normalizeSymbol($('add-el').value);
        if (!BY_SYMBOL[el]) return toast(`Unknown element "${$('add-el').value}".`, true);
        const v = [num('add-x'), num('add-y'), num('add-z')];
        const ns = state.structure.clone();
        if ($('add-frac').checked && !ns.periodic) return toast('Fractional coordinates need a cell.', true);
        ns.push(el, $('add-frac').checked ? ns.toCartesian(v) : v);
        state.selection = new Set([ns.count - 1]);
        setStructure(ns, { keepSelection: true, fit: false, message: `Added ${el}` });
    });

    $('ce-apply').addEventListener('click', () => {
        const v = ['ce-a', 'ce-b', 'ce-c', 'ce-al', 'ce-be', 'ce-ga'].map(num);
        if (v.some((x) => !(x > 0))) return toast('Enter all six cell parameters.', true);
        const s = state.structure, ns = s.clone();
        const newCell = latticeFromParameters(...v);
        if (s.periodic && $('ce-scale').checked) {
            const fr = s.fractionalPositions();
            ns.cell = newCell;
            ns.setFractional(fr);
        } else ns.cell = newCell;
        setStructure(ns, { message: 'Cell updated', fit: true });
    });
    $('ce-vac-btn').addEventListener('click', () => { try { setStructure(addVacuum(state.structure, num('ce-vac')), { message: 'Vacuum added' }); } catch (e) { toast(e.message, true); } });
    $('ce-wrap').addEventListener('click', () => { if (!state.structure.periodic) return toast('No cell to wrap into.', true); setStructure(state.structure.clone().wrap(), { fit: false, message: 'Atoms wrapped into the cell' }); });
    $('ce-center').addEventListener('click', () => {
        const s = state.structure.clone();
        if (!s.count) return;
        const c = s.center();
        const target = s.periodic ? vecmat([0.5, 0.5, 0.5], s.cell) : [0, 0, 0];
        const d = sub(target, c);
        s.positions = s.positions.map((p) => add(p, d));
        setStructure(s, { fit: false, message: 'Centred' });
    });
    $('ce-remove').addEventListener('click', () => { const s = state.structure.clone(); s.cell = null; setStructure(s, { message: 'Cell removed: structure is now a cluster' }); });

    $('rn-apply').addEventListener('click', () => {
        const from = $('rn-from').value, to = normalizeSymbol($('rn-to').value);
        if (!BY_SYMBOL[to]) return toast(`Unknown element "${$('rn-to').value}".`, true);
        const s = state.structure.clone();
        s.symbols = s.symbols.map((x) => (x === from ? to : x));
        setStructure(s, { fit: false, message: `${from} → ${to}`, spaceGroup: state.spaceGroupLabel });
    });
    $('title-apply').addEventListener('click', () => {
        const s = state.structure.clone(); s.title = $('title-in').value;
        setStructure(s, { fit: false, spaceGroup: state.spaceGroupLabel });
    });
    $('ed-dedupe').addEventListener('click', () => {
        const s = state.structure.clone(); const n0 = s.count;
        removeDuplicates(s, 0.5);
        setStructure(s, { fit: false, message: `Removed ${n0 - s.count} overlapping atom(s)` });
    });
}

function setSelection(list) {
    state.selection = new Set(list);
    renderHighlights();
    refreshSelectionUI();
}

function selectAll() { setSelection([...Array(state.structure.count).keys()]); }

function deleteSelection() {
    if (!state.selection.size) return toast('Select atoms first.', true);
    const n = state.selection.size;
    const ns = state.structure.clone().removeIndices([...state.selection]);
    setStructure(ns, { fit: false, message: `Deleted ${n} atom(s)` });
}

// ---------------------------------------------------------------------------
// View tab
// ---------------------------------------------------------------------------

function initView() {
    segmented($('v-style'), (v) => { state.view.style = v; render(); });
    const bindRange = (id, key, digits = 2) => {
        const out = $(id + '-out');
        $(id).addEventListener('input', () => {
            state.view[key] = num(id);
            out.value = num(id).toFixed(digits);
            if (key === 'tolerance') computeAnalysis();
            render();
            if (key === 'tolerance') { refreshCoordination(); refreshInfo(); }
        });
    };
    bindRange('v-atom', 'atomScale');
    bindRange('v-bond', 'bondRadius');
    bindRange('v-tol', 'tolerance');
    $('v-bonds').addEventListener('change', (e) => { state.view.bonds = e.target.checked; render(); });
    $('v-labels').addEventListener('change', (e) => { state.view.labels = e.target.checked; renderOverlays(); });
    $('v-color').addEventListener('change', (e) => { state.view.color = e.target.value; render(); refreshLegend(); });
    $('v-cell').addEventListener('change', (e) => { state.view.cell = e.target.checked; render(); });
    $('v-boundary').addEventListener('change', (e) => { state.view.boundary = e.target.checked; render(); });
    ['v-ra', 'v-rb', 'v-rc'].forEach((id, k) => $(id).addEventListener('change', () => {
        const v = Math.max(1, Math.min(10, int(id) || 1));
        $(id).value = v;
        state.view.reps[k] = v;
        render({ fit: true });
    }));
    const updPlane = () => {
        state.view.plane = $('pl-on').checked ? { hkl: [int('pl-h'), int('pl-k'), int('pl-l')], d: num('pl-d') } : null;
        if (state.view.plane && state.view.plane.hkl.every((x) => !x)) { state.view.plane = null; toast('Miller indices cannot all be zero.', true); }
        renderOverlays();
    };
    ['pl-on', 'pl-h', 'pl-k', 'pl-l', 'pl-d'].forEach((id) => $(id).addEventListener('change', updPlane));
    const updDir = () => {
        state.view.direction = $('dir-on').checked ? [int('dir-u'), int('dir-v'), int('dir-w')] : null;
        if (state.view.direction && state.view.direction.every((x) => !x)) state.view.direction = null;
        renderOverlays();
    };
    ['dir-on', 'dir-u', 'dir-v', 'dir-w'].forEach((id) => $(id).addEventListener('change', updDir));
    segmented($('v-proj'), (v) => viewer.setProjection(v));
    $('v-rotate').addEventListener('change', (e) => { viewer.autoRotate = e.target.checked; });
    $('va-go').addEventListener('click', () => {
        const s = state.structure;
        const uvw = [int('va-u'), int('va-v'), int('va-w')];
        const dir = s.periodic ? vecmat(uvw, s.cell) : uvw;
        viewer.viewAlong(dir, s.periodic ? s.cell[1] : [0, 1, 0]);
    });
    document.querySelectorAll('[data-view]').forEach((b) => b.addEventListener('click', () => {
        const s = state.structure, v = b.dataset.view;
        if (v === 'fit') return viewer.fit(boundsOf(state.display));
        const idx = { a: 0, b: 1, c: 2 }[v];
        if (s.periodic) {
            // Look down the axis: along the reciprocal direction so the other two lie in the screen plane.
            const rec = reciprocal(s.cell);
            const up = idx === 2 ? s.cell[1] : s.cell[2];
            viewer.viewAlong(rec[idx], up);
        } else viewer.viewAlong([[1, 0, 0], [0, 1, 0], [0, 0, 1]][idx], idx === 2 ? [0, 1, 0] : [0, 0, 1]);
    }));
}

// ---------------------------------------------------------------------------
// Analyze tab: RDF
// ---------------------------------------------------------------------------

function initAnalyze() {
    $('rdf-go').addEventListener('click', () => withBusy('Computing g(r)…', () => {
        const s = state.structure;
        if (s.count < 2) throw new Error('Need at least two atoms.');
        const res = radialDistribution(s, { rmax: num('rdf-r'), bins: int('rdf-bins'), a: $('rdf-a').value || null, b: $('rdf-b').value || null });
        drawPlot($('rdf-plot'), res.r, res.g, { xlabel: 'r (Å)', ylabel: 'g(r)' });
        $('rdf-summary').innerHTML = [
            ['First peak', res.firstPeak ? `${fmt(res.firstPeak)} Å` : '—'],
            ['First minimum', res.firstMin ? `${fmt(res.firstMin)} Å` : '—'],
            ['First-shell CN', res.firstShellCN ? fmt(res.firstShellCN, 2) : '—'],
            ['Boundary', res.periodic ? 'periodic' : 'open (cluster)'],
        ].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    }));
    $('measure-clear').addEventListener('click', () => { state.measurements = []; state.pending = []; refreshMeasurements(); renderOverlays(); renderHighlights(); });
}

function drawPlot(canvas, xs, ys, { xlabel, ylabel }) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    const css = getComputedStyle(document.documentElement);
    const ink = css.getPropertyValue('--ink-2').trim(), line = css.getPropertyValue('--line').trim(), acc = css.getPropertyValue('--accent-2').trim();
    ctx.clearRect(0, 0, W, H);
    const m = { l: 58, r: 16, t: 16, b: 50 };
    const xmax = xs[xs.length - 1], ymax = Math.max(1.2, ...ys) * 1.08;
    const X = (x) => m.l + (x / xmax) * (W - m.l - m.r);
    const Y = (y) => H - m.b - (y / ymax) * (H - m.t - m.b);
    ctx.font = '22px Inter, sans-serif';
    ctx.fillStyle = ink; ctx.strokeStyle = line; ctx.lineWidth = 1.5;
    // grid + ticks
    for (let i = 0; i <= 4; i++) {
        const y = (ymax / 4) * i;
        ctx.beginPath(); ctx.moveTo(m.l, Y(y)); ctx.lineTo(W - m.r, Y(y)); ctx.stroke();
        ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(y.toFixed(y < 10 ? 1 : 0), m.l - 8, Y(y));
    }
    const step = xmax > 10 ? 2 : 1;
    for (let x = 0; x <= xmax + 1e-9; x += step) {
        ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(x.toString(), X(x), H - m.b + 8);
    }
    ctx.textAlign = 'center'; ctx.fillText(xlabel, (m.l + W - m.r) / 2, H - 24);
    ctx.save(); ctx.translate(16, (m.t + H - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(ylabel, 0, 0); ctx.restore();
    // g = 1 reference
    ctx.setLineDash([8, 6]); ctx.beginPath(); ctx.moveTo(m.l, Y(1)); ctx.lineTo(W - m.r, Y(1)); ctx.stroke(); ctx.setLineDash([]);
    ctx.strokeStyle = acc; ctx.lineWidth = 3; ctx.beginPath();
    xs.forEach((x, i) => { const px = X(x), py = Y(ys[i]); if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py); });
    ctx.stroke();
}

function refreshMeasurements() {
    $('measure-list').innerHTML = state.measurements.map((m, k) => `<div><span>${k + 1}. ${m.label}</span><b>${m.text}</b></div>`).join('') || '<p class="note">No measurements yet.</p>';
}

// ---------------------------------------------------------------------------
// Viewport interaction: hover, click-select, box-select, measuring
// ---------------------------------------------------------------------------

function initInteraction() {
    const vp = $('viewport'), canvas = viewer.renderer.domElement;
    const card = $('hover-card'), rect = $('select-rect');
    let down = null, boxing = false, hoverRaf = 0, lastMove = null;

    canvas.addEventListener('pointerdown', (e) => {
        down = { x: e.clientX, y: e.clientY, shift: e.shiftKey };
        if (e.shiftKey && state.mode === 'select' && e.button === 0) {
            boxing = true;
            viewer.controls.enabled = false;
            rect.style.display = 'block';
            updateRect(e);
        }
    });
    const updateRect = (e) => {
        const r = vp.getBoundingClientRect();
        const x0 = Math.min(down.x, e.clientX) - r.left, y0 = Math.min(down.y, e.clientY) - r.top;
        Object.assign(rect.style, { left: x0 + 'px', top: y0 + 'px', width: Math.abs(e.clientX - down.x) + 'px', height: Math.abs(e.clientY - down.y) + 'px' });
    };
    window.addEventListener('pointermove', (e) => {
        if (boxing) { updateRect(e); return; }
        if (e.target !== canvas) { card.style.display = 'none'; return; }
        lastMove = e;
        if (!hoverRaf) hoverRaf = requestAnimationFrame(() => { hoverRaf = 0; hover(lastMove); });
    });
    window.addEventListener('pointerup', (e) => {
        if (boxing) {
            boxing = false;
            viewer.controls.enabled = true;
            rect.style.display = 'none';
            const ids = viewer.atomsInRect(down.x, down.y, e.clientX, e.clientY);
            const next = e.ctrlKey || e.metaKey ? new Set([...state.selection, ...ids]) : new Set(ids);
            setSelection([...next]);
            down = null;
            return;
        }
        if (!down || e.target !== canvas) { down = null; return; }
        const moved = Math.hypot(e.clientX - down.x, e.clientY - down.y);
        down = null;
        if (moved > 4 || e.button !== 0) return;
        click(e);
    });
    canvas.addEventListener('pointerleave', () => { card.style.display = 'none'; });

    function hover(e) {
        const hit = viewer.pick(e.clientX, e.clientY);
        if (!hit) { card.style.display = 'none'; canvas.style.cursor = ''; return; }
        canvas.style.cursor = 'pointer';
        const s = state.structure, i = hit.index;
        const el = element(s.symbols[i]);
        const f = s.periodic ? s.toFractional(s.positions[i]) : null;
        card.innerHTML = `<b>${el.symbol}</b> <span class="note">${el.name} · #${i + 1}${hit.image ? ' · image' : ''}</span><br>
            <span class="mono">xyz ${hit.position.map((v) => fmt(v)).join('  ')}</span><br>
            ${f ? `<span class="mono">frac ${f.map((v) => fmt(v, 4)).join('  ')}</span><br>` : ''}
            <span class="mono">CN ${state.cn[i] ?? '—'}${s.tags[i] ? ` · region ${s.tags[i]}` : ''}</span>`;
        const r = vp.getBoundingClientRect();
        let x = e.clientX - r.left + 16, y = e.clientY - r.top + 16;
        card.style.display = 'block';
        if (x + card.offsetWidth > r.width - 8) x = e.clientX - r.left - card.offsetWidth - 16;
        if (y + card.offsetHeight > r.height - 8) y = e.clientY - r.top - card.offsetHeight - 16;
        card.style.left = x + 'px'; card.style.top = y + 'px';
    }

    function click(e) {
        const hit = viewer.pick(e.clientX, e.clientY);
        if (state.mode === 'select') {
            if (!hit) { if (!(e.ctrlKey || e.metaKey)) setSelection([]); return; }
            const next = new Set(e.ctrlKey || e.metaKey ? state.selection : []);
            if ((e.ctrlKey || e.metaKey) && next.has(hit.index)) next.delete(hit.index); else next.add(hit.index);
            setSelection([...next]);
            if (next.size === 1) activateTab('edit', false);
            return;
        }
        if (!hit) return;
        const need = { distance: 2, angle: 3, dihedral: 4 }[state.mode];
        state.pending.push(hit);
        renderHighlights();
        if (state.pending.length === need) {
            const pts = state.pending.map((p) => p.position);
            const s = state.structure;
            const names = state.pending.map((p) => `${s.symbols[p.index]}${p.index + 1}`).join('–');
            let text;
            if (need === 2) text = `${fmt(norm(sub(pts[1], pts[0])))} Å`;
            else if (need === 3) {
                const u = sub(pts[0], pts[1]), v = sub(pts[2], pts[1]);
                text = `${fmt(Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (norm(u) * norm(v))))) * 180 / Math.PI, 2)}°`;
            } else {
                const b1 = sub(pts[1], pts[0]), b2 = sub(pts[2], pts[1]), b3 = sub(pts[3], pts[2]);
                const n1 = cross(b1, b2), n2 = cross(b2, b3);
                const m1 = cross(n1, unit(b2));
                text = `${fmt(Math.atan2(dot(m1, n2), dot(n1, n2)) * 180 / Math.PI, 2)}°`;
            }
            state.measurements.push({ points: pts, text, label: names });
            state.pending = [];
            renderOverlays(); renderHighlights(); refreshMeasurements();
            toast(`${names}: ${text}`);
        }
    }
}

function setMode(mode) {
    state.mode = mode;
    state.pending = [];
    document.querySelectorAll('#modebar [data-mode]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.mode === mode)));
    const hints = {
        select: 'Click to select · Ctrl-click to add · Shift-drag for box selection · Delete removes the selection',
        distance: 'Click two atoms to measure their distance',
        angle: 'Click three atoms: the angle is measured at the second',
        dihedral: 'Click four atoms to measure the dihedral angle',
    };
    $('st-hint').textContent = hints[mode];
    renderHighlights();
}

// ---------------------------------------------------------------------------
// Tabs, menus, files, keyboard, theme
// ---------------------------------------------------------------------------

function activateTab(name, openPanel = true) {
    document.querySelectorAll('.tabs [data-tab]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.tab === name)));
    document.querySelectorAll('[data-body]').forEach((b) => { b.hidden = b.dataset.body !== name; });
    if (openPanel && window.innerWidth <= 860) document.body.classList.add('panel-open');
}

const EXAMPLES = [
    { label: 'Cu (FCC)', run: () => presetExample('fcc') },
    { label: 'Fe (BCC)', run: () => presetExample('bcc') },
    { label: 'Mg (HCP)', run: () => presetExample('hcp') },
    { label: 'Si (diamond)', run: () => presetExample('diamond') },
    { label: 'NaCl (rock salt)', run: () => presetExample('nacl') },
    { label: 'SrTiO₃ (perovskite)', run: () => presetExample('perovskite') },
    { label: 'Al₂O₃ (corundum)', run: () => presetExample('corundum') },
    { label: 'MgAl₂O₄ (spinel)', run: () => presetExample('spinel') },
    { sep: true },
    { label: 'Cu Σ5 (310)[001] grain boundary', run: () => { presetExample('fcc', true); $('gb-axis').value = '0,0,1'; refreshGBList(); gbChoice = listTiltBoundaries([0, 0, 1], 7, 51).find((g) => g.sigma === 5 && g.angle < 40); $('gb-n1').value = 4; $('gb-n3').value = 3; $('gb-build').click(); } },
    { label: 'Pt (111) slab', run: () => { presetExample('fcc', true, { Cu: 'Pt', a: 3.924 }); $('sl-h').value = 1; $('sl-k').value = 1; $('sl-l').value = 1; $('sl-layers').value = 4; $('sl-rep').value = '3 3'; $('sl-build').click(); } },
    { label: 'Au Wulff nanoparticle', run: () => { presetExample('fcc', true, { Cu: 'Au', a: 4.078 }); $('np-mode').querySelector('[data-v="wulff"]').click(); $('np-build').click(); } },
    { label: 'Cu polycrystal (8 grains)', run: () => { presetExample('fcc', true); $('pc-x').value = $('pc-y').value = $('pc-z').value = 45; $('pc-build').click(); } },
    { label: 'Cu intrinsic stacking fault', run: () => { presetExample('fcc', true); $('sf-u').value = 1; $('sf-build').click(); } },
    { label: 'Cu₃Au-type random alloy', run: () => { presetExample('fcc', true); setStructure(supercell(state.structure, 4, 4, 4), { message: '' }); $('ss-solutes').querySelector('.tr:not(.head) input').value = 'Au'; $('ss-build').click(); } },
];

function presetExample(id, quiet = false, subst = null) {
    $('cr-preset').value = id;
    loadPreset(id);
    if (subst) {
        const row = $('cr-sites').querySelector('.tr:not(.head) input');
        for (const [from, to] of Object.entries(subst)) if (from !== 'a' && row.value === from) row.value = to;
        if (subst.a) { $('cr-a').value = subst.a; applyLatticeConstraints(); }
        $('cr-preset').value = '';
    }
    buildFromPanel({ quiet });
    if (subst) state.structure.title = `${state.structure.formula().replace(/\d+/g, '')} (FCC)`;
}

function initChrome() {
    document.querySelectorAll('.tabs [data-tab]').forEach((b) => b.addEventListener('click', () => activateTab(b.dataset.tab)));
    $('panel-toggle').addEventListener('click', () => document.body.classList.toggle('panel-open'));
    // Menus
    document.querySelectorAll('.menu').forEach((m) => {
        m.querySelector('[data-menu]').addEventListener('click', (e) => {
            e.stopPropagation();
            document.querySelectorAll('.menu.open').forEach((o) => { if (o !== m) o.classList.remove('open'); });
            m.classList.toggle('open');
        });
    });
    document.addEventListener('click', () => document.querySelectorAll('.menu.open').forEach((o) => o.classList.remove('open')));
    const ex = $('examples-list');
    for (const item of EXAMPLES) {
        if (item.sep) { ex.appendChild(document.createElement('hr')); continue; }
        const b = document.createElement('button');
        b.textContent = item.label;
        b.addEventListener('click', () => { try { item.run(); } catch (e) { toast(e.message, true); } });
        ex.appendChild(b);
    }
    document.querySelectorAll('[data-export]').forEach((b) => b.addEventListener('click', () => exportAs(b.dataset.export)));

    // Files
    $('btn-open').addEventListener('click', () => $('file-input').click());
    $('file-input').addEventListener('change', (e) => { if (e.target.files[0]) openFile(e.target.files[0]); e.target.value = ''; });
    const vp = $('viewport');
    let dragDepth = 0;
    vp.addEventListener('dragenter', (e) => { e.preventDefault(); dragDepth++; vp.classList.add('dragging'); });
    vp.addEventListener('dragover', (e) => e.preventDefault());
    vp.addEventListener('dragleave', () => { if (--dragDepth <= 0) { dragDepth = 0; vp.classList.remove('dragging'); } });
    vp.addEventListener('drop', (e) => {
        e.preventDefault(); dragDepth = 0; vp.classList.remove('dragging');
        if (e.dataTransfer.files[0]) openFile(e.dataTransfer.files[0]);
    });

    $('btn-undo').addEventListener('click', undo);
    $('btn-redo').addEventListener('click', redo);
    $('btn-theme').addEventListener('click', () => setTheme(document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark'));
    document.querySelectorAll('#modebar [data-mode]').forEach((b) => b.addEventListener('click', () => setMode(b.dataset.mode)));

    window.addEventListener('keydown', (e) => {
        const typing = /INPUT|SELECT|TEXTAREA/.test(document.activeElement?.tagName);
        const mod = e.ctrlKey || e.metaKey;
        if (mod && e.key.toLowerCase() === 'z' && !typing) { e.preventDefault(); e.shiftKey ? redo() : undo(); return; }
        if (mod && e.key.toLowerCase() === 'y' && !typing) { e.preventDefault(); redo(); return; }
        if (mod && e.key.toLowerCase() === 'o') { e.preventDefault(); $('file-input').click(); return; }
        if (typing) return;
        if (mod && e.key.toLowerCase() === 'a') { e.preventDefault(); selectAll(); return; }
        if (e.key === 'Delete' || e.key === 'Backspace') { if (state.selection.size) { e.preventDefault(); deleteSelection(); } return; }
        if (e.key === 'Escape') { setSelection([]); setMode('select'); document.body.classList.remove('panel-open'); return; }
        if (mod || e.altKey) return;
        const k = e.key.toLowerCase();
        if (k === 'f') viewer.fit(boundsOf(state.display));
        else if (k === 's') setMode('select');
        else if (k === 'd') setMode('distance');
        else if (k === 'a') setMode('angle');
        else if (k === '1' || k === '2' || k === '3') document.querySelector(`[data-view="${'abc'[+k - 1]}"]`).click();
    });
}

function setTheme(t) {
    document.documentElement.dataset.theme = t;
    try { localStorage.setItem('studio-theme', t); } catch (e) { /* ignore */ }
    viewer.setTheme(t === 'dark');
}

async function openFile(file) {
    if (file.size > 80 * 1024 * 1024) return toast('File is larger than 80 MB.', true);
    const text = await file.text();
    await withBusy(`Reading ${file.name}…`, () => {
        const { structure, format } = parseStructure(file.name, text);
        if (!structure.count) throw new Error('No atoms found in the file.');
        state.view.color = 'element'; $('v-color').value = 'element';
        state.view.bonds = structure.count < 80000; $('v-bonds').checked = state.view.bonds;
        setStructure(structure, { message: `Opened ${file.name} (${format}, ${structure.count.toLocaleString()} atoms${structure.frames > 1 ? `, last of ${structure.frames} frames` : ''})` });
    });
}

function exportAs(kind) {
    const s = state.structure;
    const base = (s.formula() || 'structure').replace(/[^\w]+/g, '');
    try {
        if (kind === 'png' || kind === 'png-t') {
            const url = viewer.screenshot(kind === 'png' ? 2 : 4, kind === 'png-t');
            const a = document.createElement('a'); a.href = url; a.download = `${base}.png`; a.click();
            return;
        }
        const w = WRITERS[kind];
        const name = kind === 'poscar' ? `${base}.vasp` : `${base}.${w.ext}`;
        download(name, w.fn(s));
        toast(`Saved ${name}`);
    } catch (e) { toast(e.message, true); }
}

// ---------------------------------------------------------------------------
// Start-up
// ---------------------------------------------------------------------------

async function start() {
    viewer.setTheme(document.documentElement.dataset.theme === 'dark');
    $('hud-title').textContent = 'Loading space groups…';
    await loadSpaceGroups();
    initChrome();
    initBuilders();
    initEdit();
    initView();
    initAnalyze();
    initInteraction();
    refreshMeasurements();
    setMode('select');
    // Restore the last session if there is one, else start with FCC copper.
    let restored = false;
    try {
        const last = localStorage.getItem('studio-last');
        if (last) {
            const { structure } = parseStructure('last.xyz', last);
            if (structure.count) { setStructure(structure, { record: false, message: 'Restored your last structure' }); restored = true; }
        }
    } catch (e) { /* ignore */ }
    loadPreset('fcc');
    $('cr-preset').value = 'fcc';
    if (!restored) {
        buildFromPanel({ quiet: true });
        state.undo = [];
        refreshPanels();
    }
}

start().catch((e) => { console.error(e); toast('Could not start: ' + e.message, true); });

// Expose for debugging in the console.
window.studio = { state, viewer };
