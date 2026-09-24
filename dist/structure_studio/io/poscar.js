// VASP POSCAR / CONTCAR (VASP 4 and 5, selective dynamics, scale factors).
import { Structure, vecmat, inv3, dot, cross, ELEMENTS } from '../core/index.js';
import { lines, toks } from './text.js';

export function parsePOSCAR(text) {
    const L = lines(text).map((l) => l.trim());
    const title = L[0];
    const sf = toks(L[1]).map(Number);
    let cell = [toks(L[2]), toks(L[3]), toks(L[4])].map((r) => r.slice(0, 3).map(Number));
    if (sf.length === 3) cell = cell.map((r) => r.map((v, k) => v * sf[k]));
    else if (sf[0] < 0) {
        const vol = Math.abs(dot(cell[0], cross(cell[1], cell[2])));
        const f = Math.cbrt(-sf[0] / vol);
        cell = cell.map((r) => r.map((v) => v * f));
    } else cell = cell.map((r) => r.map((v) => v * sf[0]));

    let li = 5;
    let species = toks(L[li]);
    let counts;
    if (species.every((t) => /^\d+$/.test(t))) {
        counts = species.map(Number);
        // VASP 4: species may be in the title line
        const guess = toks(title).filter((t) => /^[A-Z][a-z]?$/.test(t));
        species = guess.length === counts.length ? guess : counts.map((_, i) => ELEMENTS[i].symbol);
    } else {
        species = species.map((t) => t.split(/[_/]/)[0]);
        li++;
        counts = toks(L[li]).map(Number);
    }
    li++;
    if (/^s/i.test(L[li])) li++; // selective dynamics
    const cartesian = /^[ck]/i.test(L[li]);
    li++;
    const s = new Structure({ cell, title });
    let k = 0;
    species.forEach((sp, si) => {
        for (let n = 0; n < counts[si]; n++, k++) {
            const v = toks(L[li + k]).slice(0, 3).map(Number);
            const p = cartesian ? v.map((x) => x * (sf.length === 1 && sf[0] > 0 ? sf[0] : 1)) : vecmat(v, cell);
            s.push(sp, p);
        }
    });
    return s;
}

export function writePOSCAR(s, { direct = true } = {}) {
    if (!s.periodic) throw new Error('POSCAR needs a periodic cell. Add a box first (Modify → Cell).');
    const order = [...new Set(s.symbols)];
    const idx = order.flatMap((sp) => s.symbols.map((x, i) => (x === sp ? i : -1)).filter((i) => i >= 0));
    const inv = inv3(s.cell);
    const f = (v) => v.toFixed(10).padStart(16);
    let out = `${s.title || s.formula()}\n1.0\n`;
    out += s.cell.map((r) => r.map(f).join(' ')).join('\n') + '\n';
    out += order.map((x) => x.padStart(5)).join('') + '\n';
    out += order.map((sp) => String(s.symbols.filter((x) => x === sp).length).padStart(5)).join('') + '\n';
    out += direct ? 'Direct\n' : 'Cartesian\n';
    out += idx.map((i) => (direct ? vecmat(s.positions[i], inv) : s.positions[i]).map(f).join(' ')).join('\n') + '\n';
    return out;
}
