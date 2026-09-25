// Supercell and general cell-transformation card (acts on the current structure).
import { supercell, transformCell } from '../../core/index.js';
import { $, int, toast } from '../dom.js';
import { setStructure, requireCrystal } from '../history.js';
import { registerCard, PREVIEW_LIMIT } from './card.js';

const readRepeats = () => {
    const n = [int('sc-a'), int('sc-b'), int('sc-c')];
    if (n.some((v) => !(v >= 1))) throw new Error('Repeats must be whole numbers ≥ 1.');
    return n;
};

// Entries may be fractions such as "1/2".
const readMatrix = () => {
    const v = [...$('tm').querySelectorAll('input')].map((i) => {
        const t = i.value.trim();
        if (t.includes('/')) { const [p, q] = t.split('/'); return parseFloat(p) / parseFloat(q); }
        return parseFloat(t);
    });
    if (v.some((x) => !Number.isFinite(x))) throw new Error('Fill in all nine matrix entries.');
    return [v.slice(0, 3), v.slice(3, 6), v.slice(6, 9)];
};

const isIdentity = (P) => P.every((row, i) => row.every((x, j) => Math.abs(x - (i === j ? 1 : 0)) < 1e-9));

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

function transformed() {
    const s = requireCrystal();
    const P = readMatrix();
    const base = isIdentity(P) ? s : transformCell(s, P);
    const n = readRepeats();
    return { base, n, total: base.count * n[0] * n[1] * n[2] };
}

export function initTransformCard() {
    $('tm-prim').addEventListener('click', () => {
        let lattice;
        try { lattice = detectCentring(requireCrystal()); } catch (e) { return toast(e.message, true); }
        if (!lattice) return toast('The current cell is not F- or I-centred, so it is already primitive (or needs a custom matrix).', true);
        const M = lattice === 'I'
            ? ['-1/2', '1/2', '1/2', '1/2', '-1/2', '1/2', '1/2', '1/2', '-1/2']
            : ['0', '1/2', '1/2', '1/2', '0', '1/2', '1/2', '1/2', '0'];
        $('tm').querySelectorAll('input').forEach((inp, k) => { inp.value = M[k]; });
        $('tm').dispatchEvent(new Event('change', { bubbles: true }));
        toast(`${lattice}-centred cell: primitive matrix filled in.`);
    });
    $('tm-reset').addEventListener('click', () => {
        $('tm').querySelectorAll('input').forEach((inp, k) => { inp.value = k % 4 === 0 ? '1' : '0'; });
        $('tm').dispatchEvent(new Event('change', { bubbles: true }));
    });
    registerCard('transform', {
        busy: 'Building supercell…',
        preview: () => {
            const { base, n, total } = transformed();
            // Show as many repeats as fit the preview budget.
            const shown = n.slice();
            while (base.count * shown[0] * shown[1] * shown[2] > PREVIEW_LIMIT && shown.some((v) => v > 1)) {
                const k = shown.indexOf(Math.max(...shown));
                shown[k]--;
            }
            const note = shown.join() !== n.join() ? `<br><span>Preview shows ${shown.join('×')}</span>` : '';
            return { structure: supercell(base, ...shown), caption: `<b>${total.toLocaleString()} atoms</b> · ${n.join(' × ')} of a ${base.count}-atom cell${note}` };
        },
        build: () => {
            const { base, n, total } = transformed();
            if (total > 1500000) throw new Error(`${total.toLocaleString()} atoms is too many for the browser.`);
            setStructure(supercell(base, ...n), { message: `${total.toLocaleString()} atoms (${n.join('×')})` });
        },
    });
}
