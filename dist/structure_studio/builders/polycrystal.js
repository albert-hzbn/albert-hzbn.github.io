// Voronoi polycrystals with random grain orientations (from AtomForge's
// PolyCrystalBuilder).
import { Structure, NeighborGrid, add, dot, norm, vecmat, reciprocal, randomRotation, mulberry32, cellVolume } from '../core/index.js';
import { species, nearestNeighbourDistance } from './common.js';

export function buildPolycrystal(crystal, { box = [60, 60, 60], grains = 8, seed = 1, overlap = 0.7, maxAtoms = 400000 } = {}) {
    if (!crystal.periodic) throw new Error('Polycrystals are built from a periodic crystal.');
    const rand = mulberry32(seed);
    const [Lx, Ly, Lz] = box;
    const volAtom = cellVolume(crystal.cell) / crystal.count;
    const estimate = Lx * Ly * Lz / volAtom;
    if (estimate > maxAtoms) throw new Error(`About ${Math.round(estimate).toLocaleString()} atoms: reduce the box (limit ${maxAtoms.toLocaleString()}).`);
    const seeds = Array.from({ length: grains }, () => [rand() * Lx, rand() * Ly, rand() * Lz]);
    const rots = seeds.map(() => randomRotation(rand));
    const micD2 = (a, b) => {
        let dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        dx -= Lx * Math.round(dx / Lx); dy -= Ly * Math.round(dy / Ly); dz -= Lz * Math.round(dz / Lz);
        return dx * dx + dy * dy + dz * dz;
    };
    const nearest = (p) => {
        let best = 0, bd = Infinity;
        for (let g = 0; g < grains; g++) { const d = micD2(p, seeds[g]); if (d < bd) { bd = d; best = g; } }
        return best;
    };
    // Extent of each Voronoi cell, sampled on a grid.
    const extent = new Array(grains).fill(0);
    const G = 24;
    for (let i = 0; i < G; i++) for (let j = 0; j < G; j++) for (let k = 0; k < G; k++) {
        const p = [(i + 0.5) * Lx / G, (j + 0.5) * Ly / G, (k + 0.5) * Lz / G];
        const g = nearest(p);
        extent[g] = Math.max(extent[g], Math.sqrt(micD2(p, seeds[g])));
    }
    const pad = Math.hypot(Lx, Ly, Lz) / G + 2;
    const s = new Structure({ cell: [[Lx, 0, 0], [0, Ly, 0], [0, 0, Lz]], title: `${grains}-grain ${species(crystal)} polycrystal` });
    const rec = reciprocal(crystal.cell);
    const spacing = rec.map((r) => 1 / norm(r));
    for (let g = 0; g < grains; g++) {
        const R = rots[g];
        const Rg = extent[g] + pad;
        const n = spacing.map((d) => Math.ceil(Rg / d) + 1);
        for (let i = -n[0]; i <= n[0]; i++) for (let j = -n[1]; j <= n[1]; j++) for (let k = -n[2]; k <= n[2]; k++) {
            const T = vecmat([i, j, k], crystal.cell);
            for (let a = 0; a < crystal.count; a++) {
                const local = add(crystal.positions[a], T);
                if (dot(local, local) > Rg * Rg) continue;
                const r = [dot(R[0], local), dot(R[1], local), dot(R[2], local)];
                const p = [
                    pyModF(seeds[g][0] + r[0], Lx), pyModF(seeds[g][1] + r[1], Ly), pyModF(seeds[g][2] + r[2], Lz),
                ];
                if (nearest(p) === g) s.push(crystal.symbols[a], p, g + 1);
            }
        }
    }
    // Resolve overlaps at grain boundaries.
    const nn = nearestNeighbourDistance(crystal);
    const grid = new NeighborGrid(s, overlap * nn);
    const drop = new Set();
    for (let i = 0; i < s.count; i++) {
        if (drop.has(i)) continue;
        grid.forEachNeighbor(i, overlap * nn, (j) => { if (j > i && !drop.has(j)) drop.add(j); });
    }
    s.removeIndices([...drop]);
    return s;
}

const pyModF = (a, L) => a - L * Math.floor(a / L);

