// Structure builders, ported from AtomForge's algorithms
// (BulkCrystalBuilder, NanoCrystalBuilder, PolyCrystalBuilder, CSLComputation,
// SubstitutionalSolidSolutionBuilder, StackingFaultBuilder) and simplified for
// the browser.
import {
    Structure, NeighborGrid, add, sub, scale, dot, cross, norm, unit, vecmat, inv3, matmul, det3,
    reciprocal, transformCell, removeDuplicates, randomRotation, mulberry32, gcd, supercell, cellVolume,
} from './core.js';

const EPS = 1e-6;

// "Cu", "NaCl", "SrTiO": species of a crystal without unit-cell counts, for titles.
const species = (s) => [...new Set(s.symbols)].join('');

// Python-style integer division and modulo (floor semantics), needed for the
// slab basis algorithm.
const pyDiv = (a, b) => Math.floor(a / b);
const pyMod = (a, b) => ((a % b) + b) % b;

function extGcd(a, b) {
    if (b === 0) return [1, 0];
    if (pyMod(a, b) === 0) return [0, 1];
    const [x, y] = extGcd(b, pyMod(a, b));
    return [y, x - y * pyDiv(a, b)];
}

// Is shift vector v (Cartesian) a translation symmetry of periodic structure s?
function isTranslation(s, v, tol = 0.05) {
    const probe = s.clone();
    probe.positions = s.positions.map((p) => add(p, v));
    // For every shifted atom there must be an atom of the same species within tol.
    const merged = s.clone();
    const n0 = s.count;
    for (let i = 0; i < probe.count; i++) merged.push(probe.symbols[i], probe.positions[i]);
    const g = new NeighborGrid(merged, 0.5);
    for (let i = n0; i < merged.count; i++) {
        let ok = false;
        g.forEachNeighbor(i, tol, (j, d) => { if (j < n0 && merged.symbols[j] === merged.symbols[i] && d <= tol) ok = true; });
        if (!ok) return false;
    }
    return true;
}

// ---------- Oriented orthogonal box from a crystal ----------

// rows: three mutually orthogonal lattice directions [uvw] (in units of the
// crystal's cell vectors). The shortest period along each direction is found
// automatically, and the result is rotated so the rows lie along x, y, z.
export function orientedBox(crystal, rows, reps = [1, 1, 1]) {
    const cart = rows.map((r) => vecmat(r, crystal.cell));
    for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++) {
        if (Math.abs(dot(unit(cart[i]), unit(cart[j]))) > 1e-6) throw new Error('Box directions must be mutually orthogonal.');
    }
    // Shortest lattice period along each direction.
    const P = rows.map((r, k) => {
        for (const div of [6, 4, 3, 2]) {
            const v = scale(cart[k], 1 / div);
            if (isTranslation(crystal, v)) return r.map((x) => x / div);
        }
        return r.slice();
    });
    let box = transformCell(crystal, P);
    box = supercell(box, reps[0], reps[1], reps[2]);
    // Rotate into the x, y, z frame.
    const R = box.cell.map(unit);                  // rows: new axes in old frame
    const rot = (p) => [dot(R[0], p), dot(R[1], p), dot(R[2], p)];
    if (det3(R) < 0) R[2] = scale(R[2], -1);
    box.positions = box.positions.map(rot);
    box.cell = box.cell.map(rot).map((v) => v.map((x) => (Math.abs(x) < 1e-9 ? 0 : x)));
    return box.wrap();
}

// ---------- Surface slab (algorithm of ASE's ase.build.surface) ----------

export function buildSlab(crystal, [h, k, l], layers = 4, vacuum = 10, { center = true } = {}) {
    if (!crystal.periodic) throw new Error('Slabs need a periodic bulk structure.');
    if (!h && !k && !l) throw new Error('Miller indices cannot all be zero.');
    const g = gcd(gcd(h, k), l);
    [h, k, l] = [h / g, k / g, l / g];
    const [a1, a2, a3] = crystal.cell;
    let c1, c2, c3;
    const h0 = h === 0, k0 = k === 0, l0 = l === 0;
    if ((h0 && k0) || (h0 && l0) || (k0 && l0)) {
        if (!h0) [c1, c2, c3] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]];
        if (!k0) [c1, c2, c3] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];
        if (!l0) [c1, c2, c3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    } else {
        let [p, q] = extGcd(k, l);
        const t1 = sub(scale(a1, k), scale(a2, h)), t2 = sub(scale(a1, l), scale(a3, h)), t3 = sub(scale(a2, l), scale(a3, k));
        const k1 = dot(add(scale(t1, p), scale(t2, q)), t3);
        const k2 = dot(sub(scale(t1, l), scale(t2, k)), t3);
        if (Math.abs(k2) > 1e-10) {
            const i = -Math.round(k1 / k2);
            p += i * l; q -= i * k;
        }
        const [a, b] = extGcd(p * k + q * l, h);
        c1 = [p * k + q * l, -p * h, -q * h];
        const gg = Math.abs(gcd(l, k));
        c2 = [0, l / gg, -k / gg];
        c3 = [b, a * p, a * q];
    }
    const basis = [c1, c2, c3];
    // Express atoms in the new basis and fold into [0,1).
    const newCell = matmul(basis, crystal.cell);
    const inv = inv3(newCell);
    let s = new Structure({ cell: newCell, title: `${species(crystal)} (${h}${k}${l}) slab` });
    for (let i = 0; i < crystal.count; i++) {
        const f = vecmat(crystal.positions[i], inv).map((v) => v - Math.floor(v + 1e-8));
        s.push(crystal.symbols[i], vecmat(f, newCell), crystal.tags[i]);
    }
    s = supercell(s, 1, 1, layers);
    // Make c perpendicular to the surface.
    const [b1, b2, b3] = s.cell;
    const n = cross(b1, b2);
    const c = scale(n, dot(b3, n) / dot(n, n));
    s.cell = [b1, b2, c];
    // Standard orientation: a along x, surface normal along z.
    const fr = s.fractionalPositions();
    const A = norm(b1), bx = dot(b1, b2) / A, by = Math.sqrt(Math.max(0, dot(b2, b2) - bx * bx));
    s.cell = [[A, 0, 0], [bx, by, 0], [0, 0, norm(c)]];
    s.setFractional(fr.map((f) => [f[0] - Math.floor(f[0] + 1e-9), f[1] - Math.floor(f[1] + 1e-9), f[2]]));
    return addVacuum(s, vacuum, { center });
}

// Add vacuum along c (periodic) or build a padded box (non-periodic).
export function addVacuum(s, vacuum, { center = true } = {}) {
    s = s.clone();
    if (!s.periodic) {
        const { lo, hi } = s.bounds();
        const L = [0, 1, 2].map((k) => hi[k] - lo[k] + vacuum);
        s.cell = [[L[0], 0, 0], [0, L[1], 0], [0, 0, L[2]]];
        const shift = [0, 1, 2].map((k) => -lo[k] + vacuum / 2);
        s.positions = s.positions.map((p) => add(p, shift));
        return s;
    }
    const n = unit(cross(s.cell[0], s.cell[1]));
    const heights = s.positions.map((p) => dot(p, n));
    const zmin = Math.min(...heights), zmax = Math.max(...heights);
    const thickness = zmax - zmin;
    const cn = dot(s.cell[2], n);
    const newLen = thickness + vacuum;
    s.cell[2] = scale(s.cell[2], newLen / cn);
    if (center) {
        const shift = scale(n, vacuum / 2 - zmin);
        s.positions = s.positions.map((p) => add(p, shift));
    }
    return s;
}

// ---------- Symmetric tilt grain boundaries (cubic lattices) ----------

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

// Remove atoms of grain `tag` that sit closer than dmin to any other atom.
function removeCloseInGrain(s, dmin, tag) {
    const grid = new NeighborGrid(s, dmin);
    const drop = new Set();
    for (let i = 0; i < s.count; i++) {
        if (s.tags[i] !== tag || drop.has(i)) continue;
        grid.forEachNeighbor(i, dmin, (j) => { if (!drop.has(j) && j !== i) drop.add(i); });
    }
    s.removeIndices([...drop]);
    return drop.size;
}

// ---------- Nanoparticles ----------

function cubicFamily([h, k, l]) {
    const out = new Map();
    const perms = [[h, k, l], [h, l, k], [k, h, l], [k, l, h], [l, h, k], [l, k, h]];
    for (const p of perms) for (let sx = -1; sx <= 1; sx += 2) for (let sy = -1; sy <= 1; sy += 2) for (let sz = -1; sz <= 1; sz += 2) {
        const v = [p[0] * sx, p[1] * sy, p[2] * sz];
        out.set(v.join(','), v);
    }
    return [...out.values()];
}

export const SHAPES = {
    sphere: { label: 'Sphere', test: (r, R) => dot(r, r) <= R * R },
    cube: { label: 'Cube', test: (r, R) => Math.max(Math.abs(r[0]), Math.abs(r[1]), Math.abs(r[2])) <= R },
    octahedron: { label: 'Octahedron', test: (r, R) => Math.abs(r[0]) + Math.abs(r[1]) + Math.abs(r[2]) <= R },
    cuboctahedron: { label: 'Cuboctahedron', test: (r, R) => Math.max(Math.abs(r[0]), Math.abs(r[1]), Math.abs(r[2])) <= R && Math.abs(r[0]) + Math.abs(r[1]) + Math.abs(r[2]) <= 2 * R },
    truncoct: { label: 'Truncated octahedron', test: (r, R) => Math.max(Math.abs(r[0]), Math.abs(r[1]), Math.abs(r[2])) <= R && Math.abs(r[0]) + Math.abs(r[1]) + Math.abs(r[2]) <= 1.5 * R },
    cylinder: { label: 'Cylinder (along z)', test: (r, R, H) => r[0] * r[0] + r[1] * r[1] <= R * R && Math.abs(r[2]) <= H / 2 },
    ellipsoid: { label: 'Ellipsoid', test: (r, R, H, E) => (r[0] / R) ** 2 + (r[1] / (E || R)) ** 2 + (r[2] / (H / 2 || R)) ** 2 <= 1 },
};

function tileAround(crystal, center, radius, fn) {
    const rec = reciprocal(crystal.cell);
    const spacing = rec.map((r) => 1 / norm(r));
    const n = spacing.map((d) => Math.ceil(radius / d) + 1);
    const fc = vecmat(center, inv3(crystal.cell)).map(Math.floor);
    for (let i = fc[0] - n[0]; i <= fc[0] + n[0]; i++) for (let j = fc[1] - n[1]; j <= fc[1] + n[1]; j++) for (let k = fc[2] - n[2]; k <= fc[2] + n[2]; k++) {
        const T = vecmat([i, j, k], crystal.cell);
        for (let a = 0; a < crystal.count; a++) fn(crystal.symbols[a], add(crystal.positions[a], T));
    }
}

export function buildNanoparticle(crystal, { shape = 'sphere', radius = 12, height = 20, radiusY = 0, centerOn = 'atom', vacuum = 8 } = {}) {
    if (!crystal.periodic) throw new Error('Nanoparticles are cut from a periodic crystal.');
    const center = centerOn === 'atom' ? crystal.positions[0].slice() : vecmat([0.5, 0.5, 0.5], crystal.cell);
    const test = SHAPES[shape].test;
    const reach = Math.max(radius * 2, height, radiusY) + 2;
    const s = new Structure({ title: `${SHAPES[shape].label} ${species(crystal)} nanoparticle` });
    tileAround(crystal, center, reach, (sym, p) => {
        const r = sub(p, center);
        if (test(r, radius, height, radiusY)) s.push(sym, r);
    });
    removeDuplicates(s, 0.05);
    return vacuum > 0 ? addVacuum(s, vacuum * 2) : s;
}

// facets: [{ hkl: [h,k,l], energy }]; cubicSymmetry expands each family.
export function buildWulff(crystal, facets, { radius = 15, cubicSymmetry = true, centerOn = 'atom', vacuum = 8 } = {}) {
    if (!crystal.periodic) throw new Error('Nanoparticles are cut from a periodic crystal.');
    if (!facets.length) throw new Error('Add at least one facet family.');
    const rec = reciprocal(crystal.cell);
    const gmin = Math.min(...facets.map((f) => f.energy));
    const planes = [];
    for (const f of facets) {
        const members = cubicSymmetry ? cubicFamily(f.hkl) : [f.hkl, f.hkl.map((v) => -v)];
        for (const m of members) {
            const nrm = unit(vecmat(m, rec));
            planes.push({ n: nrm, d: radius * f.energy / gmin });
        }
    }
    const center = centerOn === 'atom' ? crystal.positions[0].slice() : vecmat([0.5, 0.5, 0.5], crystal.cell);
    const reach = Math.max(...planes.map((p) => p.d)) * 1.8 + 2;
    const s = new Structure({ title: `Wulff ${species(crystal)} nanoparticle` });
    tileAround(crystal, center, reach, (sym, p) => {
        const r = sub(p, center);
        if (dot(r, r) > reach * reach) return;
        for (const pl of planes) if (dot(r, pl.n) > pl.d + EPS) return;
        s.push(sym, r);
    });
    removeDuplicates(s, 0.05);
    return vacuum > 0 ? addVacuum(s, vacuum * 2) : s;
}

// ---------- Voronoi polycrystal ----------

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

export function nearestNeighbourDistance(crystal) {
    const grid = new NeighborGrid(crystal, 6);
    let best = Infinity;
    for (let i = 0; i < crystal.count; i++) grid.forEachNeighbor(i, 6, (j, d) => { if (d > 0.1) best = Math.min(best, d); });
    return isFinite(best) ? best : 2.5;
}

// ---------- Substitutional solid solution ----------

// solutes: [{ symbol, fraction }] replacing randomly chosen host atoms.
export function buildSolidSolution(s, host, solutes, seed = 1) {
    const out = s.clone();
    const rand = mulberry32(seed);
    const sites = out.symbols.map((x, i) => (x === host ? i : -1)).filter((i) => i >= 0);
    if (!sites.length) throw new Error(`No ${host} atoms in the structure.`);
    // Fisher–Yates shuffle
    for (let i = sites.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [sites[i], sites[j]] = [sites[j], sites[i]]; }
    let k = 0;
    const report = [];
    for (const sol of solutes) {
        const n = Math.round(sol.fraction * sites.length);
        for (let m = 0; m < n && k < sites.length; m++, k++) out.symbols[sites[k]] = sol.symbol;
        report.push(`${n} ${sol.symbol}`);
    }
    out.title = `${out.formula()} solid solution`;
    return { structure: out, report: report.join(', ') + ` on ${sites.length} ${host} sites` };
}

// ---------- FCC stacking fault (generalised stacking-fault path) ----------

// Displacement u is in units of the Shockley partial b_p = a/√6 along [11-2].
// The cell's c vector is tilted by the same displacement, so the cell holds a
// single fault plane and stays periodic (the usual GSFE set-up).
export function buildStackingFault(fccCrystal, { reps = [2, 1, 4], u = 1 } = {}) {
    const box = orientedBox(fccCrystal, [[1, -1, 0], [1, 1, -2], [1, 1, 1]], reps);
    const a = fccCrystal.cell[0][0];
    const bp = a / Math.sqrt(6);
    const zs = [...new Set(box.positions.map((p) => +p[2].toFixed(4)))].sort((x, y) => x - y);
    const mid = zs[Math.floor(zs.length / 2)] - 1e-3;
    const shift = [0, u * bp, 0];
    box.positions = box.positions.map((p) => (p[2] >= mid ? add(p, shift) : p));
    box.tags = box.positions.map((p) => (p[2] >= mid ? 2 : 1));
    box.cell[2] = add(box.cell[2], shift);
    box.title = `${species(fccCrystal)} stacking fault, u = ${u.toFixed(2)} bₚ`;
    return box.wrap();
}
