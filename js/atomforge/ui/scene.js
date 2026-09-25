// Turns the current structure and view settings into viewer content:
// display atoms (with repeats and face images), bonds, colours, overlays.
import { Structure, vecmat, norm, element } from '../core/index.js';
import { findBonds, coordinationNumbers } from '../analysis/bonds.js';
import { planePolygon } from '../render/geometry.js';
import { state, viewer } from './context.js';
import { toast } from './dom.js';

export const TAG_COLORS = ['#9aa5b1', '#4e79a7', '#f28e2b', '#59a14f', '#e15759', '#76b7b2', '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#86bcb6', '#d37295', '#a0cbe8', '#ffbe7d', '#8cd17d'];
export const CN_COLORS = ['#6b7280', '#a855f7', '#6366f1', '#3b82f6', '#0ea5e9', '#14b8a6', '#22c55e', '#84cc16', '#eab308', '#f59e0b', '#f97316', '#ef4444', '#b91c1c', '#7f1d1d'];
const RAMP = [[0.00, [68, 1, 84]], [0.25, [59, 82, 139]], [0.5, [33, 145, 140]], [0.75, [94, 201, 98]], [1.0, [253, 231, 37]]];

export function ramp(t) {
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


export const elementColor = (sym) => state.customColors[sym] || element(sym).color;

export function computeAnalysis() {
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


export function render({ fit = false } = {}) {
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


export function boundsOf(display) {
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


export function renderHighlights() {
    const items = [];
    const sel = state.selection;
    if (sel.size) for (const d of state.display) if (sel.has(d.index)) items.push({ position: d.position, radius: d.radius, color: '#2f7de1' });
    for (const p of state.pending) items.push({ position: p.position, radius: p.radius, color: '#e0662f' });
    viewer.setHighlights(items);
}


export function renderOverlays() {
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

