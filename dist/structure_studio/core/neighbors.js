// Periodic-aware cell-list neighbour search and minimum-image vectors.
import { vecmat, inv3, norm, sub } from './math.js';
import { reciprocal } from './lattice.js';

export class NeighborGrid {
    constructor(s, cutoff) {
        this.s = s;
        this.cutoff = Math.max(cutoff, 0.5);
        this.periodic = s.periodic;
        const n = s.count;
        if (this.periodic) {
            // Work in fractional space; bin widths chosen from plane spacings.
            this.inv = inv3(s.cell);
            const rec = reciprocal(s.cell);
            const spacing = rec.map((r) => 1 / norm(r));
            this.nb = spacing.map((d) => Math.max(1, Math.floor(d / this.cutoff)));
            this.frac = s.positions.map((p) => vecmat(p, this.inv).map((v) => v - Math.floor(v)));
            // How many image shells are needed per direction.
            this.reach = spacing.map((d, k) => Math.ceil(this.cutoff / (d / this.nb[k])));
        } else {
            const { lo, hi } = s.bounds();
            this.lo = lo;
            this.nb = [0, 1, 2].map((k) => Math.max(1, Math.ceil((hi[k] - lo[k] + 1e-6) / this.cutoff)));
        }
        this.bins = new Map();
        for (let i = 0; i < n; i++) {
            const key = this.keyOf(this.binOf(i));
            let list = this.bins.get(key);
            if (!list) this.bins.set(key, list = []);
            list.push(i);
        }
    }

    binOf(i) {
        if (this.periodic) return this.frac[i].map((f, k) => Math.min(this.nb[k] - 1, Math.floor(f * this.nb[k])));
        const p = this.s.positions[i];
        return [0, 1, 2].map((k) => Math.min(this.nb[k] - 1, Math.floor((p[k] - this.lo[k]) / this.cutoff)));
    }

    keyOf(b) { return (b[0] * 73856093) ^ (b[1] * 19349663) ^ (b[2] * 83492791); }

    // Calls fn(j, distance, shiftVector) for each neighbour j of atom i within r.
    forEachNeighbor(i, r, fn) {
        const s = this.s, pi = s.positions[i], r2 = r * r;
        const b = this.binOf(i);
        if (!this.periodic) {
            for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
                const list = this.bins.get(this.keyOf([b[0] + dx, b[1] + dy, b[2] + dz]));
                if (!list) continue;
                for (const j of list) {
                    if (j === i) continue;
                    const pj = s.positions[j];
                    const d2 = (pj[0] - pi[0]) ** 2 + (pj[1] - pi[1]) ** 2 + (pj[2] - pi[2]) ** 2;
                    if (d2 <= r2) fn(j, Math.sqrt(d2), null);
                }
            }
            return;
        }
        const fi = this.frac[i];
        const [ra, rb, rc] = this.reach;
        const seen = new Set();
        for (let dx = -ra; dx <= ra; dx++) for (let dy = -rb; dy <= rb; dy++) for (let dz = -rc; dz <= rc; dz++) {
            const cb = [b[0] + dx, b[1] + dy, b[2] + dz];
            const wrapped = cb.map((v, k) => ((v % this.nb[k]) + this.nb[k]) % this.nb[k]);
            const img = cb.map((v, k) => Math.floor(v / this.nb[k]));
            const tag = wrapped.join(',') + '|' + img.join(',');
            if (seen.has(tag)) continue;
            seen.add(tag);
            const list = this.bins.get(this.keyOf(wrapped));
            if (!list) continue;
            for (const j of list) {
                if (j === i && img[0] === 0 && img[1] === 0 && img[2] === 0) continue;
                const df = [this.frac[j][0] + img[0] - fi[0], this.frac[j][1] + img[1] - fi[1], this.frac[j][2] + img[2] - fi[2]];
                const d = vecmat(df, s.cell);
                const d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                if (d2 <= r2) fn(j, Math.sqrt(d2), d);
            }
        }
    }
}

// Minimum-image distance vector from atom i to atom j (brute force over ±1 images).
export function mic(s, i, j) {
    const d0 = sub(s.positions[j], s.positions[i]);
    if (!s.periodic) return d0;
    const f = vecmat(d0, inv3(s.cell)).map((v) => v - Math.round(v));
    let best = null, bestN = Infinity;
    for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++) {
        const d = vecmat([f[0] + a, f[1] + b, f[2] + c], s.cell);
        const n = norm(d);
        if (n < bestN) { bestN = n; best = d; }
    }
    return best;
}
