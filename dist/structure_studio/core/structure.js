// The Structure model: symbols, Cartesian positions, per-atom tags and an
// optional cell (rows a, b, c in Angstrom). r = f · cell for fractional f.
import { vecmat, inv3 } from './math.js';
import { cellVolume } from './lattice.js';
import { element, normalizeSymbol } from './elements.js';

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
