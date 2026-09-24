// Supercells, general cell transformations and duplicate removal.
import { add, scale, vecmat, inv3, matmul, det3 } from './math.js';
import { Structure } from './structure.js';
import { NeighborGrid } from './neighbors.js';

// Integer supercell along a, b, c.
export function supercell(s, na, nb, nc) {
    if (!s.periodic) throw new Error('Supercells need a periodic cell.');
    const out = new Structure({ cell: [scale(s.cell[0], na), scale(s.cell[1], nb), scale(s.cell[2], nc)], title: s.title });
    for (let i = 0; i < na; i++) for (let j = 0; j < nb; j++) for (let k = 0; k < nc; k++) {
        const t = vecmat([i, j, k], s.cell);
        for (let n = 0; n < s.count; n++) out.push(s.symbols[n], add(s.positions[n], t), s.tags[n]);
    }
    return out;
}

// General transformation: new cell rows = P · old cell rows (P integer, det ≠ 0).
export function transformCell(s, P) {
    if (!s.periodic) throw new Error('Transformations need a periodic cell.');
    const d = det3(P);
    if (Math.abs(d) < 1e-6) throw new Error('Transformation matrix must have a non-zero determinant.');
    const newCell = matmul(P, s.cell);
    const inv = inv3(newCell);
    // Range of old-cell translations needed to cover the new cell.
    const lo = [0, 0, 0], hi = [0, 0, 0];
    for (let i = 0; i < 8; i++) {
        const corner = [i & 1, (i >> 1) & 1, (i >> 2) & 1];
        const f = vecmat(corner, P);
        for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], Math.floor(f[k])); hi[k] = Math.max(hi[k], Math.ceil(f[k])); }
    }
    const out = new Structure({ cell: newCell, title: s.title });
    const expected = Math.round(Math.abs(d) * s.count);
    const eps = 1e-6;
    const fr = s.fractionalPositions();
    for (let i = lo[0] - 1; i <= hi[0]; i++) for (let j = lo[1] - 1; j <= hi[1]; j++) for (let k = lo[2] - 1; k <= hi[2]; k++) {
        for (let n = 0; n < s.count; n++) {
            const p = vecmat(add(fr[n], [i, j, k]), s.cell);
            const f = vecmat(p, inv);
            if (f.every((v) => v >= -eps && v < 1 - eps)) out.push(s.symbols[n], p, s.tags[n]);
        }
    }
    if (out.count !== expected) console.warn(`transformCell: expected ${expected} atoms, got ${out.count}`);
    return out.wrap();
}

// Remove atoms closer than tol (Å) to an earlier atom; periodic-aware.
export function removeDuplicates(s, tol = 0.1) {
    const drop = new Set();
    const grid = new NeighborGrid(s, tol);
    for (let i = 0; i < s.count; i++) {
        if (drop.has(i)) continue;
        grid.forEachNeighbor(i, tol, (j) => { if (j > i) drop.add(j); });
    }
    return s.removeIndices([...drop]);
}
