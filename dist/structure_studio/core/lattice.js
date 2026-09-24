// Lattice parameters, volumes and reciprocal vectors.
import { dot, norm, det3, inv3, transpose } from './math.js';

// Same convention as AtomForge's buildLatticeFromParameters: a along x, b in xy.
export function latticeFromParameters(a, b, c, alpha, beta, gamma) {
    const d = Math.PI / 180;
    const ca = Math.cos(alpha * d), cb = Math.cos(beta * d), cg = Math.cos(gamma * d), sg = Math.sin(gamma * d);
    const cx = c * cb;
    const cy = Math.abs(sg) > 1e-10 ? c * (ca - cb * cg) / sg : 0;
    const cz2 = c * c - cx * cx - cy * cy;
    return [[a, 0, 0], [b * cg, b * sg, 0], [cx, cy, cz2 > 0 ? Math.sqrt(cz2) : 0]];
}

export function cellParameters(cell) {
    const [a, b, c] = cell.map(norm);
    const ang = (u, v) => Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (norm(u) * norm(v))))) * 180 / Math.PI;
    return { a, b, c, alpha: ang(cell[1], cell[2]), beta: ang(cell[0], cell[2]), gamma: ang(cell[0], cell[1]) };
}

export const cellVolume = (cell) => Math.abs(det3(cell));

// Reciprocal vectors without the 2π factor (rows).
export function reciprocal(cell) {
    return transpose(inv3(cell));
}
