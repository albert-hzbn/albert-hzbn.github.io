// Covalent bonds (periodic aware) and coordination numbers.
import { NeighborGrid, add, element } from '../core/index.js';

// Bonds between atoms whose distance is below tolerance × (r_i + r_j).
// Returns pairs within the displayed structure (bonds across periodic
// boundaries are returned with their image vector so they can be drawn as stubs).
export function findBonds(s, { tolerance = 1.15, maxBonds = 400000 } = {}) {
    if (!s.count) return [];
    const radii = s.symbols.map((x) => element(x).radius);
    const rmax = Math.max(...radii);
    const cutoff = 2 * rmax * tolerance;
    const grid = new NeighborGrid(s, cutoff);
    const bonds = [];
    for (let i = 0; i < s.count && bonds.length < maxBonds; i++) {
        grid.forEachNeighbor(i, cutoff, (j, d, vec) => {
            if (j < i && !vec) return;
            if (d < 0.4) return;
            const lim = tolerance * (radii[i] + radii[j]);
            if (d > lim) return;
            if (vec) {
                // Periodic: keep each pair once; decide whether it crosses the cell.
                const direct = [s.positions[j][0] - s.positions[i][0], s.positions[j][1] - s.positions[i][1], s.positions[j][2] - s.positions[i][2]];
                const crosses = Math.abs(direct[0] - vec[0]) + Math.abs(direct[1] - vec[1]) + Math.abs(direct[2] - vec[2]) > 1e-4;
                if (!crosses && j < i) return;
                if (crosses) bonds.push({ i, j, d, end: add(s.positions[i], vec), crosses: true });
                else bonds.push({ i, j, d });
            } else bonds.push({ i, j, d });
        });
    }
    return bonds;
}

export function coordinationNumbers(s, bonds) {
    const cn = new Array(s.count).fill(0);
    for (const b of bonds) {
        cn[b.i]++;
        if (!b.crosses) cn[b.j]++;
    }
    return cn;
}
