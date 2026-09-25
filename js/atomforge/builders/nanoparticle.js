// Nanoparticles: geometric shapes and Wulff constructions (from AtomForge's
// NanoCrystalBuilder).
import { Structure, add, sub, dot, norm, unit, vecmat, inv3, reciprocal, removeDuplicates } from '../core/index.js';
import { convexPolyhedron } from '../core/polyhedron.js';
import { species, addVacuum } from './common.js';

const EPS = 1e-6;

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

// Half-space planes of a Wulff construction: facet family f sits at distance
// radius · γ_f / γ_min from the centre. Each plane carries its family label.
export function wulffPlanes(crystal, facets, { radius = 15, cubicSymmetry = true } = {}) {
    if (!facets.length) throw new Error('Add at least one facet family.');
    const rec = reciprocal(crystal.cell);
    const gmin = Math.min(...facets.map((f) => f.energy));
    const planes = [];
    facets.forEach((f, family) => {
        const members = cubicSymmetry ? cubicFamily(f.hkl) : [f.hkl, f.hkl.map((v) => -v)];
        for (const m of members) planes.push({ n: unit(vecmat(m, rec)), d: radius * f.energy / gmin, family, hkl: m });
    });
    return planes;
}

// The Wulff polyhedron itself, with the surface-area fraction of each family.
export function wulffShape(crystal, facets, options = {}) {
    const planes = wulffPlanes(crystal, facets, options);
    const poly = convexPolyhedron(planes);
    const total = poly.faces.reduce((a, f) => a + f.area, 0) || 1;
    const fractions = facets.map((_, family) => poly.faces.filter((f) => f.plane.family === family).reduce((a, f) => a + f.area, 0) / total);
    return { ...poly, fractions, closed: poly.faces.length >= 4 };
}

// Outline of a geometric shape for display: half-spaces for polyhedral shapes,
// analytic descriptions for the curved ones.
export function shapeOutline(shape, R, H = 2 * R, E = R) {
    const cube = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
    const oct = [];
    for (const x of [-1, 1]) for (const y of [-1, 1]) for (const z of [-1, 1]) oct.push([x, y, z].map((v) => v / Math.sqrt(3)));
    const planes = (list, d, family) => list.map((n) => ({ n, d, family }));
    switch (shape) {
        case 'cube': return { kind: 'poly', ...convexPolyhedron(planes(cube, R, 0)) };
        case 'octahedron': return { kind: 'poly', ...convexPolyhedron(planes(oct, R / Math.sqrt(3), 1)) };
        case 'cuboctahedron': return { kind: 'poly', ...convexPolyhedron([...planes(cube, R, 0), ...planes(oct, 2 * R / Math.sqrt(3), 1)]) };
        case 'truncoct': return { kind: 'poly', ...convexPolyhedron([...planes(cube, R, 0), ...planes(oct, 1.5 * R / Math.sqrt(3), 1)]) };
        case 'cylinder': return { kind: 'cylinder', radius: R, height: H };
        case 'ellipsoid': return { kind: 'ellipsoid', radii: [R, E || R, H / 2 || R] };
        default: return { kind: 'sphere', radius: R };
    }
}

// facets: [{ hkl: [h,k,l], energy }]; cubicSymmetry expands each family.
export function buildWulff(crystal, facets, { radius = 15, cubicSymmetry = true, centerOn = 'atom', vacuum = 8 } = {}) {
    if (!crystal.periodic) throw new Error('Nanoparticles are cut from a periodic crystal.');
    const planes = wulffPlanes(crystal, facets, { radius, cubicSymmetry });
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
