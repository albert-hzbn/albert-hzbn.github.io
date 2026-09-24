// Core structure model and lattice maths.
// Cells are stored as three row vectors [a, b, c] in Angstrom, so a Cartesian
// position is r = f · cell for fractional coordinates f.
import { element, normalizeSymbol } from './elements.js';

// ---------- 3-vector / 3x3 helpers ----------

export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const norm = (a) => Math.hypot(a[0], a[1], a[2]);
export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
export const unit = (a) => { const n = norm(a) || 1; return [a[0] / n, a[1] / n, a[2] / n]; };

export function det3(m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

export function inv3(m) {
    const d = det3(m);
    if (Math.abs(d) < 1e-12) throw new Error('Singular matrix');
    const r = [
        [m[1][1] * m[2][2] - m[1][2] * m[2][1], m[0][2] * m[2][1] - m[0][1] * m[2][2], m[0][1] * m[1][2] - m[0][2] * m[1][1]],
        [m[1][2] * m[2][0] - m[1][0] * m[2][2], m[0][0] * m[2][2] - m[0][2] * m[2][0], m[0][2] * m[1][0] - m[0][0] * m[1][2]],
        [m[1][0] * m[2][1] - m[1][1] * m[2][0], m[0][1] * m[2][0] - m[0][0] * m[2][1], m[0][0] * m[1][1] - m[0][1] * m[1][0]],
    ];
    return r.map((row) => row.map((v) => v / d));
}

export function matmul(a, b) {
    return a.map((row) => [0, 1, 2].map((j) => row[0] * b[0][j] + row[1] * b[1][j] + row[2] * b[2][j]));
}

export const transpose = (m) => [0, 1, 2].map((i) => [m[0][i], m[1][i], m[2][i]]);

// Row vector times matrix: v · M
export const vecmat = (v, m) => [
    v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
    v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
    v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
];

// Matrix times column vector: M · v
export const matvec = (m, v) => [dot(m[0], v), dot(m[1], v), dot(m[2], v)];

export function gcd(a, b) {
    a = Math.abs(a); b = Math.abs(b);
    while (b) [a, b] = [b, a % b];
    return a;
}

// Rotation matrix (acting on column vectors) about unit axis by angle (rad).
export function rotationMatrix(axis, angle) {
    const [x, y, z] = unit(axis), c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
    return [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ];
}

// Uniformly distributed random rotation (Shoemake's method).
export function randomRotation(rand = Math.random) {
    const u1 = rand(), u2 = rand() * 2 * Math.PI, u3 = rand() * 2 * Math.PI;
    const a = Math.sqrt(1 - u1), b = Math.sqrt(u1);
    const [w, x, y, z] = [a * Math.sin(u2), a * Math.cos(u2), b * Math.sin(u3), b * Math.cos(u3)];
    return [
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ];
}

// Small seeded PRNG so builders are reproducible.
export function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// ---------- Lattice parameters ----------

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

// ---------- Structure ----------

export class Structure {
    constructor({ cell = null, symbols = [], positions = [], tags = null, title = '' } = {}) {
        this.cell = cell;                       // [[a],[b],[c]] or null for molecules/clusters
        this.symbols = symbols.map(normalizeSymbol);
        this.positions = positions;             // Cartesian, array of [x,y,z]
        this.tags = tags || new Array(symbols.length).fill(0);
        this.title = title;
    }

    get count() { return this.symbols.length; }
    get periodic() { return !!this.cell && cellVolume(this.cell) > 1e-8; }

    clone() {
        return new Structure({
            cell: this.cell ? this.cell.map((r) => r.slice()) : null,
            symbols: this.symbols.slice(),
            positions: this.positions.map((p) => p.slice()),
            tags: this.tags.slice(),
            title: this.title,
        });
    }

    toFractional(p) { return vecmat(p, inv3(this.cell)); }
    toCartesian(f) { return vecmat(f, this.cell); }

    fractionalPositions() {
        const inv = inv3(this.cell);
        return this.positions.map((p) => vecmat(p, inv));
    }

    setFractional(fracs) {
        this.positions = fracs.map((f) => vecmat(f, this.cell));
    }

    push(symbol, position, tag = 0) {
        this.symbols.push(normalizeSymbol(symbol));
        this.positions.push(position);
        this.tags.push(tag);
    }

    wrap() {
        if (!this.periodic) return this;
        const fr = this.fractionalPositions().map((f) => f.map((v) => {
            let w = v - Math.floor(v);
            if (w > 1 - 1e-9) w = 0;
            return w;
        }));
        this.setFractional(fr);
        return this;
    }

    formula() {
        const counts = {};
        for (const s of this.symbols) counts[s] = (counts[s] || 0) + 1;
        // Object keys keep insertion order, i.e. order of first appearance (SrTiO3, not O3SrTi).
        return Object.keys(counts).map((k) => k + (counts[k] > 1 ? counts[k] : '')).join('');
    }

    speciesCounts() {
        const counts = {};
        for (const s of this.symbols) counts[s] = (counts[s] || 0) + 1;
        return counts;
    }

    mass() { return this.symbols.reduce((m, s) => m + element(s).mass, 0); }

    density() {
        if (!this.periodic) return null;
        return this.mass() * 1.66053906660 / cellVolume(this.cell); // g/cm^3
    }

    center() {
        const n = this.count || 1;
        const c = [0, 0, 0];
        for (const p of this.positions) { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; }
        return c.map((v) => v / n);
    }

    bounds() {
        const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
        for (const p of this.positions) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k]); hi[k] = Math.max(hi[k], p[k]); }
        if (this.periodic) {
            for (let i = 0; i < 8; i++) {
                const v = vecmat([i & 1, (i >> 1) & 1, (i >> 2) & 1], this.cell);
                for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], v[k]); hi[k] = Math.max(hi[k], v[k]); }
            }
        }
        if (!isFinite(lo[0])) return { lo: [0, 0, 0], hi: [0, 0, 0] };
        return { lo, hi };
    }

    removeIndices(indices) {
        const drop = new Set(indices);
        const keep = (_, i) => !drop.has(i);
        this.symbols = this.symbols.filter(keep);
        this.positions = this.positions.filter(keep);
        this.tags = this.tags.filter(keep);
        return this;
    }
}

// ---------- Cell transformations ----------

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

// ---------- Neighbour search (cell list, periodic aware) ----------

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
