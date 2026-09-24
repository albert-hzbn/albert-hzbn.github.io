// Edit panel: selection tools, atom editing, add atom, cell and species.
import { normalizeSymbol, BY_SYMBOL, NeighborGrid, add, sub, vecmat, latticeFromParameters, removeDuplicates } from '../../core/index.js';
import { addVacuum } from '../../builders/index.js';
import { state } from '../context.js';
import { $, num, toast } from '../dom.js';
import { setStructure } from '../history.js';
import { setSelection, selectAll, deleteSelection } from './selection.js';

export function initEdit() {
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

