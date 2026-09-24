// Structure analysis: bonds, coordination, radial distribution function.
// The RDF follows AtomForge's RadialDistributionAnalysis (pair-resolved,
// periodic, normalised to g(r)).
import { NeighborGrid, cellVolume, add } from './core.js';
import { element } from './elements.js';

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

// g(r) for all pairs or a chosen pair of species.
export function radialDistribution(s, { rmax = 8, bins = 200, a = null, b = null } = {}) {
    const n = s.count;
    const dr = rmax / bins;
    const hist = new Float64Array(bins);
    const grid = new NeighborGrid(s, rmax);
    const isA = (i) => !a || s.symbols[i] === a;
    const isB = (j) => !b || s.symbols[j] === b;
    let nA = 0, nB = 0;
    for (let i = 0; i < n; i++) { if (isA(i)) nA++; if (isB(i)) nB++; }
    for (let i = 0; i < n; i++) {
        if (!isA(i)) continue;
        grid.forEachNeighbor(i, rmax, (j, d) => {
            if (!isB(j) || d < 1e-6) return;
            const k = Math.floor(d / dr);
            if (k < bins) hist[k]++;
        });
    }
    const vol = s.periodic ? cellVolume(s.cell) : boundingVolume(s);
    const rhoB = nB / vol;
    const r = [], g = [], coord = [];
    let cum = 0;
    for (let k = 0; k < bins; k++) {
        const r0 = k * dr, r1 = r0 + dr;
        const shell = (4 / 3) * Math.PI * (r1 ** 3 - r0 ** 3);
        const val = nA ? hist[k] / (nA * rhoB * shell) : 0;
        cum += nA ? hist[k] / nA : 0;
        r.push(r0 + dr / 2); g.push(val); coord.push(cum);
    }
    // First peak and following minimum, as in AtomForge's RDF summary.
    let peak = 0;
    for (let k = 1; k < bins; k++) if (g[k] > g[peak]) peak = k;
    let first = -1;
    for (let k = 1; k < bins - 1; k++) if (g[k] > 1 && g[k] >= g[k - 1] && g[k] >= g[k + 1]) { first = k; break; }
    let minAfter = -1;
    if (first >= 0) for (let k = first + 1; k < bins - 1; k++) if (g[k] <= g[k - 1] && g[k] <= g[k + 1]) { minAfter = k; break; }
    return {
        r, g, coord,
        firstPeak: first >= 0 ? r[first] : null,
        firstMin: minAfter >= 0 ? r[minAfter] : null,
        firstShellCN: minAfter >= 0 ? coord[minAfter] : null,
        periodic: s.periodic,
    };
}

function boundingVolume(s) {
    const { lo, hi } = s.bounds();
    return Math.max(1, (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]));
}
