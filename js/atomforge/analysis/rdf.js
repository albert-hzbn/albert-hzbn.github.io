// Radial distribution function g(r), following AtomForge's
// RadialDistributionAnalysis (pair-resolved, periodic, normalised).
import { NeighborGrid, cellVolume } from '../core/index.js';

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
