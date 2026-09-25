// Structure changes with undo/redo, plus guards used by the builders.
// Every change emits a 'structure' event; the panels re-render on it.
import { cellParameters } from '../core/index.js';
import { WRITERS } from '../io/index.js';
import { state, bus } from './context.js';
import { toast } from './dom.js';

export function saveLocal() {
    try {
        if (state.structure.count <= 20000) localStorage.setItem('studio-last', WRITERS.xyz.fn(state.structure));
    } catch { /* storage unavailable */ }
}


export function setStructure(s, { message = '', record = true, fit = true, keepSelection = false, spaceGroup = '' } = {}) {
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
    bus.emit('structure', { fit });
    if (message) toast(message);
    saveLocal();
}


export function undo() {
    if (!state.undo.length) return;
    state.redo.push({ s: state.structure, sg: state.spaceGroupLabel });
    const prev = state.undo.pop();
    state.structure = prev.s;
    state.spaceGroupLabel = prev.sg;
    state.selection.clear();
    bus.emit('structure', { fit: false });
}


export function redo() {
    if (!state.redo.length) return;
    state.undo.push({ s: state.structure, sg: state.spaceGroupLabel });
    const next = state.redo.pop();
    state.structure = next.s;
    state.spaceGroupLabel = next.sg;
    state.selection.clear();
    bus.emit('structure', { fit: false });
}


// The current structure, which must be periodic (supercells, transforms).
export function requireCrystal() {
    if (!state.structure.periodic || !state.structure.count) throw new Error('This needs a periodic structure. Build or open one first.');
    return state.structure;
}

// ---------- Source crystal for the derived builders ----------

export function setSource(s, label = '') {
    if (!s || !s.periodic || !s.count) return;
    state.source = s;
    state.sourceLabel = label || s.title || s.formula();
    bus.emit('source');
}

export function requireSource() {
    if (!state.source) throw new Error('Build or open a periodic crystal first; it becomes the source for this builder.');
    return state.source;
}

export function isCubicCell(s) {
    const p = cellParameters(s.cell);
    return Math.abs(p.a - p.b) < 1e-3 && Math.abs(p.a - p.c) < 1e-3 && [p.alpha, p.beta, p.gamma].every((x) => Math.abs(x - 90) < 1e-2)
        && Math.abs(s.cell[0][1]) + Math.abs(s.cell[0][2]) + Math.abs(s.cell[1][0]) + Math.abs(s.cell[1][2]) < 1e-6;
}

export function requireCubicSource() {
    const s = requireSource();
    if (!isCubicCell(s)) throw new Error('Needs a conventional cubic source crystal (a = b = c, 90°). Build one from a cubic space group.');
    return s;
}
