// PDB reader (CRYST1 cell and ATOM / HETATM records).
import { Structure, latticeFromParameters } from '../core/index.js';
import { lines } from './text.js';

export function parsePDB(text) {
    const L = lines(text);
    let cell = null;
    const s = new Structure();
    for (const l of L) {
        if (l.startsWith('CRYST1')) {
            const v = [l.slice(6, 15), l.slice(15, 24), l.slice(24, 33), l.slice(33, 40), l.slice(40, 47), l.slice(47, 54)].map(Number);
            if (v[0] > 1.5) cell = latticeFromParameters(...v);
        } else if (l.startsWith('ATOM') || l.startsWith('HETATM')) {
            const sym = (l.slice(76, 78).trim() || l.slice(12, 16).trim());
            s.push(sym, [+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)]);
        }
    }
    s.cell = cell;
    return s;
}
