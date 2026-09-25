// Surface slabs for any (hkl), following the basis algorithm of ASE's
// ase.build.surface, with vacuum along the surface normal.
import { Structure, add, sub, scale, dot, cross, norm, vecmat, inv3, matmul, gcd, supercell } from '../core/index.js';
import { species, addVacuum } from './common.js';

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
