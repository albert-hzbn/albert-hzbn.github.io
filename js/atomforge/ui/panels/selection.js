// Selection state and the selected-atom editor.
import { element } from '../../core/index.js';
import { state } from '../context.js';
import { $, fmt, toast } from '../dom.js';
import { setStructure } from '../history.js';
import { renderHighlights } from '../scene.js';

export function refreshSelectionUI() {
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


export function setSelection(list) {
    state.selection = new Set(list);
    renderHighlights();
    refreshSelectionUI();
}


export function selectAll() { setSelection([...Array(state.structure.count).keys()]); }

export function deleteSelection() {
    if (!state.selection.size) return toast('Select atoms first.', true);
    const n = state.selection.size;
    const ns = state.structure.clone().removeIndices([...state.selection]);
    setStructure(ns, { fit: false, message: `Deleted ${n} atom(s)` });
}

