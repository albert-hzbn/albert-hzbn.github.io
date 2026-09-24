// Amorphous structures by Random Sequential Addition (port of AtomForge's
// AmorphousBuilder). Atoms are placed one at a time at uniformly random
// positions in an orthogonal periodic box; a candidate is accepted only if it
// respects the minimum separation to every atom already placed. A cell list
// makes each check O(1).
import { Structure, element, normalizeSymbol, mulberry32 } from '../core/index.js';

const AMU_A3_TO_G_CM3 = 1.66054;

// Box edge (Å) of a cube whose volume gives the target density (g/cm³).
export function boxForDensity(composition, density) {
    const mass = composition.reduce((m, c) => m + c.count * element(c.symbol).mass, 0);
    return Math.cbrt(mass * AMU_A3_TO_G_CM3 / Math.max(1e-3, density));
}

/**
 * composition: [{ symbol, count }]
 * box: [a, b, c] in Å, or null to derive a cube from `density`
 * pairDistances: { 'Si-O': 1.4, ... } explicit minimum separations (Å)
 * tolerance: fraction of (r_cov,i + r_cov,j) used for other pairs
 * scaleFactor: multiplies the box before packing (> 1 eases placement)
 */
export function buildAmorphous({
    composition, box = null, density = 2.2, scaleFactor = 1, pairDistances = {},
    tolerance = 0.75, seed = 42, maxAttempts = 1000,
}) {
    const comp = composition
        .map((c) => ({ symbol: normalizeSymbol(c.symbol), count: Math.max(0, Math.round(c.count)) }))
        .filter((c) => c.count > 0);
    const total = comp.reduce((n, c) => n + c.count, 0);
    if (!total) throw new Error('Add at least one element with a non-zero count.');
    if (total > 200000) throw new Error('At most 200,000 atoms can be placed in the browser.');

    let [A, B, C] = box || Array(3).fill(boxForDensity(comp, density));
    [A, B, C] = [A * scaleFactor, B * scaleFactor, C * scaleFactor];
    if (![A, B, C].every((v) => Number.isFinite(v) && v >= 0.1)) throw new Error('Box edges must be at least 0.1 Å.');

    // Minimum separation for each species pair.
    const species = comp.map((c) => c.symbol);
    const key = (x, y) => (x < y ? `${x}-${y}` : `${y}-${x}`);
    const explicit = {};
    for (const [k, v] of Object.entries(pairDistances)) {
        const [x, y] = k.split('-').map(normalizeSymbol);
        explicit[key(x, y)] = v;
    }
    const minDist = {};
    let maxMin = 0.01;
    for (const x of species) for (const y of species) {
        const d = explicit[key(x, y)] ?? tolerance * (element(x).radius + element(y).radius);
        minDist[`${x}|${y}`] = d * d;
        maxMin = Math.max(maxMin, d);
    }

    // Shuffled list of atoms to place.
    const rand = mulberry32(seed || Date.now());
    const list = comp.flatMap((c) => Array(c.count).fill(c.symbol));
    for (let i = list.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [list[i], list[j]] = [list[j], list[i]]; }

    // Cell list with periodic wrapping; cell edge ≥ largest minimum distance.
    const nA = Math.max(1, Math.floor(A / maxMin)), nB = Math.max(1, Math.floor(B / maxMin)), nC = Math.max(1, Math.floor(C / maxMin));
    const cells = new Map();
    const cellKey = (i, j, k) => (i * nB + j) * nC + k;
    const wrap = (v, n) => ((v % n) + n) % n;
    const xs = [], ys = [], zs = [], syms = [];
    let skipped = 0;

    for (const sym of list) {
        let placed = false;
        for (let attempt = 0; attempt < maxAttempts && !placed; attempt++) {
            const x = rand() * A, y = rand() * B, z = rand() * C;
            const gi = Math.min(nA - 1, Math.floor(x / A * nA)), gj = Math.min(nB - 1, Math.floor(y / B * nB)), gk = Math.min(nC - 1, Math.floor(z / C * nC));
            let conflict = false;
            for (let di = -1; di <= 1 && !conflict; di++) for (let dj = -1; dj <= 1 && !conflict; dj++) for (let dk = -1; dk <= 1 && !conflict; dk++) {
                const bucket = cells.get(cellKey(wrap(gi + di, nA), wrap(gj + dj, nB), wrap(gk + dk, nC)));
                if (!bucket) continue;
                for (const idx of bucket) {
                    let dx = x - xs[idx], dy = y - ys[idx], dz = z - zs[idx];
                    dx -= Math.round(dx / A) * A; dy -= Math.round(dy / B) * B; dz -= Math.round(dz / C) * C;
                    if (dx * dx + dy * dy + dz * dz < minDist[`${sym}|${syms[idx]}`]) { conflict = true; break; }
                }
            }
            if (conflict) continue;
            const idx = xs.length;
            xs.push(x); ys.push(y); zs.push(z); syms.push(sym);
            const k = cellKey(gi, gj, gk);
            let bucket = cells.get(k);
            if (!bucket) cells.set(k, bucket = []);
            bucket.push(idx);
            placed = true;
        }
        if (!placed) skipped++;
    }

    const s = new Structure({ cell: [[A, 0, 0], [0, B, 0], [0, 0, C]], title: `Amorphous ${species.join('')}` });
    for (let i = 0; i < xs.length; i++) s.push(syms[i], [xs[i], ys[i], zs[i]]);
    const actualDensity = s.mass() * AMU_A3_TO_G_CM3 / (A * B * C);
    return {
        structure: s,
        requested: total,
        placed: xs.length,
        skipped,
        density: actualDensity,
        box: [A, B, C],
        message: skipped
            ? `Placed ${xs.length} of ${total} atoms; ${skipped} could not be placed. Increase the box, scale factor or attempts, or lower the tolerance.`
            : `Placed all ${total} atoms, density ${actualDensity.toFixed(3)} g/cm³.`,
    };
}
