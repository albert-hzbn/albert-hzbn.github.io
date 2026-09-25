// Helpers shared by the builders: titles, lattice translations, oriented
// orthogonal boxes, vacuum padding and overlap removal.
import {
    NeighborGrid, add, scale, dot, cross, unit, vecmat, det3, transformCell, supercell,
} from '../core/index.js';

// "Cu", "NaCl", "SrTiO": species of a crystal without unit-cell counts, for titles.
export const species = (s) => [...new Set(s.symbols)].join('');

// Is shift vector v (Cartesian) a translation symmetry of periodic structure s?
export function isTranslation(s, v, tol = 0.05) {
    const probe = s.clone();
    probe.positions = s.positions.map((p) => add(p, v));
    // For every shifted atom there must be an atom of the same species within tol.
    const merged = s.clone();
    const n0 = s.count;
    for (let i = 0; i < probe.count; i++) merged.push(probe.symbols[i], probe.positions[i]);
    const g = new NeighborGrid(merged, 0.5);
    for (let i = n0; i < merged.count; i++) {
        let ok = false;
        g.forEachNeighbor(i, tol, (j, d) => { if (j < n0 && merged.symbols[j] === merged.symbols[i] && d <= tol) ok = true; });
        if (!ok) return false;
    }
    return true;
}

// rows: three mutually orthogonal lattice directions [uvw] (in units of the
// crystal's cell vectors). The shortest period along each direction is found
// automatically, and the result is rotated so the rows lie along x, y, z.
export function orientedBox(crystal, rows, reps = [1, 1, 1]) {
    const cart = rows.map((r) => vecmat(r, crystal.cell));
    for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++) {
        if (Math.abs(dot(unit(cart[i]), unit(cart[j]))) > 1e-6) throw new Error('Box directions must be mutually orthogonal.');
    }
    // Shortest lattice period along each direction.
    const P = rows.map((r, k) => {
        for (const div of [6, 4, 3, 2]) {
            const v = scale(cart[k], 1 / div);
            if (isTranslation(crystal, v)) return r.map((x) => x / div);
        }
        return r.slice();
    });
    let box = transformCell(crystal, P);
    box = supercell(box, reps[0], reps[1], reps[2]);
    // Rotate into the x, y, z frame.
    const R = box.cell.map(unit);                  // rows: new axes in old frame
    const rot = (p) => [dot(R[0], p), dot(R[1], p), dot(R[2], p)];
    if (det3(R) < 0) R[2] = scale(R[2], -1);
    box.positions = box.positions.map(rot);
    box.cell = box.cell.map(rot).map((v) => v.map((x) => (Math.abs(x) < 1e-9 ? 0 : x)));
    return box.wrap();
}

// Add vacuum along c (periodic) or build a padded box (non-periodic).
export function addVacuum(s, vacuum, { center = true } = {}) {
    s = s.clone();
    if (!s.periodic) {
        const { lo, hi } = s.bounds();
        const L = [0, 1, 2].map((k) => hi[k] - lo[k] + vacuum);
        s.cell = [[L[0], 0, 0], [0, L[1], 0], [0, 0, L[2]]];
        const shift = [0, 1, 2].map((k) => -lo[k] + vacuum / 2);
        s.positions = s.positions.map((p) => add(p, shift));
        return s;
    }
    const n = unit(cross(s.cell[0], s.cell[1]));
    const heights = s.positions.map((p) => dot(p, n));
    const zmin = Math.min(...heights), zmax = Math.max(...heights);
    const thickness = zmax - zmin;
    const cn = dot(s.cell[2], n);
    const newLen = thickness + vacuum;
    s.cell[2] = scale(s.cell[2], newLen / cn);
    if (center) {
        const shift = scale(n, vacuum / 2 - zmin);
        s.positions = s.positions.map((p) => add(p, shift));
    }
    return s;
}

// Remove atoms of grain `tag` that sit closer than dmin to any other atom.
export function removeCloseInGrain(s, dmin, tag) {
    const grid = new NeighborGrid(s, dmin);
    const drop = new Set();
    for (let i = 0; i < s.count; i++) {
        if (s.tags[i] !== tag || drop.has(i)) continue;
        grid.forEachNeighbor(i, dmin, (j) => { if (!drop.has(j) && j !== i) drop.add(i); });
    }
    s.removeIndices([...drop]);
    return drop.size;
}

export function nearestNeighbourDistance(crystal) {
    const grid = new NeighborGrid(crystal, 6);
    let best = Infinity;
    for (let i = 0; i < crystal.count; i++) grid.forEachNeighbor(i, 6, (j, d) => { if (d > 0.1) best = Math.min(best, d); });
    return isFinite(best) ? best : 2.5;
}
