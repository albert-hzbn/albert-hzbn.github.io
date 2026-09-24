// Symmetric tilt CSL grain boundaries in cubic crystals (from AtomForge's
// CSLComputation, restricted to symmetric tilt boundaries).
import { Structure, dot, cross, unit, gcd, removeDuplicates } from '../core/index.js';
import { species, orientedBox, removeCloseInGrain, nearestNeighbourDistance } from './common.js';

export function sigmaOf([h, k, l]) {
    let s = h * h + k * k + l * l;
    while (s % 2 === 0 && s > 0) s /= 2;
    return s;
}

// Tilt angles are measured from the nearest low-index mirror plane that contains
// the axis ({100} for [001], (001) for [110], {110} for [111]), which reproduces
// the usual literature values (e.g. Σ5(310)[001] = 36.87°, Σ3(112)[110] = 70.53°).
function mirrorReferences(t, maxIndex) {
    let best = Infinity, refs = [];
    for (let h = -maxIndex; h <= maxIndex; h++) for (let k = -maxIndex; k <= maxIndex; k++) for (let l = -maxIndex; l <= maxIndex; l++) {
        if ((!h && !k && !l) || h * t[0] + k * t[1] + l * t[2] !== 0) continue;
        const n2 = h * h + k * k + l * l;
        if (n2 < best) { best = n2; refs = []; }
        if (n2 === best) refs.push([h, k, l]);
    }
    return refs;
}

export function listTiltBoundaries(axis, maxIndex = 7, maxSigma = 99) {
    const t = axis;
    const refs = mirrorReferences(t, 3);
    const out = new Map();
    for (let h = -maxIndex; h <= maxIndex; h++) for (let k = -maxIndex; k <= maxIndex; k++) for (let l = -maxIndex; l <= maxIndex; l++) {
        if (!h && !k && !l) continue;
        if (h * t[0] + k * t[1] + l * t[2] !== 0) continue;
        if (gcd(gcd(h, k), l) !== 1) continue;
        const sig = sigmaOf([h, k, l]);
        if (sig > maxSigma || sig === 1) continue;
        const angle = Math.min(...refs.map((r) => 2 * Math.acos(Math.min(1, Math.abs(dot(unit([h, k, l]), unit(r))))) * 180 / Math.PI));
        const key = `${sig}|${angle.toFixed(2)}`;
        const cand = { sigma: sig, plane: [h, k, l], angle };
        const prev = out.get(key);
        const score = (p) => p.plane.filter((x) => x >= 0).length * 10 - p.plane.reduce((a, b) => a + Math.abs(b), 0);
        if (!prev || score(cand) > score(prev)) out.set(key, cand);
    }
    return [...out.values()].sort((a, b) => a.sigma - b.sigma || a.angle - b.angle);
}

// crystal: conventional cubic cell with axes along x, y, z.
// overlap: atoms of grain 2 closer than overlap × nearest-neighbour distance to
// another atom are removed.
export function buildTiltGB(crystal, axis, plane, { reps = [3, 1, 1], shift = [0, 0], overlap = 0.7 } = {}) {
    if (dot(axis, plane) !== 0) throw new Error('The boundary plane must contain the tilt axis.');
    const n = plane, t = axis;
    const y = cross(t, n);
    const gy = gcd(gcd(y[0], y[1]), y[2]);
    const rows = [n, y.map((v) => v / gy), t];
    const grain = orientedBox(crystal, rows, reps);
    const Lx = grain.cell[0][0];
    const out = new Structure({
        cell: [[2 * Lx, 0, 0], grain.cell[1].slice(), grain.cell[2].slice()],
        title: `${species(crystal)} Σ${sigmaOf(plane)} (${plane.join(' ')})[${axis.join(' ')}] tilt boundary`,
    });
    const Ly = grain.cell[1][1], Lz = grain.cell[2][2];
    for (let i = 0; i < grain.count; i++) {
        const p = grain.positions[i];
        out.push(grain.symbols[i], p.slice(), 1);
        // Mirror image across x = Lx gives the second grain, with optional rigid-body shift.
        out.push(grain.symbols[i], [2 * Lx - p[0], p[1] + shift[0] * Ly, p[2] + shift[1] * Lz], 2);
    }
    out.wrap();
    removeDuplicates(out, 0.05);
    out.removed = removeCloseInGrain(out, overlap * nearestNeighbourDistance(crystal), 2);
    return out;
}

