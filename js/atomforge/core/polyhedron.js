// Convex polyhedron from half-spaces n·r ≤ d (e.g. a Wulff construction).
// Vertices are intersections of plane triples that satisfy every constraint;
// each plane that touches ≥ 3 vertices becomes a face.
import { dot, cross, sub, norm, unit, det3 } from './math.js';

const TOL = 1e-6;

/**
 * planes: [{ n: unit normal, d: distance from origin, ...extra }]
 * Returns { vertices, faces: [{ plane, points (ordered), area }] }.
 * Planes that end up with no face (lying outside the shape) are simply absent.
 */
export function convexPolyhedron(planes) {
    const verts = [];
    const addVertex = (v) => {
        for (const w of verts) if (Math.abs(w[0] - v[0]) + Math.abs(w[1] - v[1]) + Math.abs(w[2] - v[2]) < 1e-6 * (1 + norm(v))) return;
        verts.push(v);
    };
    const n = planes.length;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) for (let k = j + 1; k < n; k++) {
        const A = [planes[i].n, planes[j].n, planes[k].n];
        const D = det3(A);
        if (Math.abs(D) < 1e-9) continue;
        // Cramer's rule for A x = d.
        const b = [planes[i].d, planes[j].d, planes[k].d];
        const x = [0, 1, 2].map((c) => det3(A.map((row, r) => row.map((v, cc) => (cc === c ? b[r] : v)))) / D);
        if (planes.every((p) => dot(p.n, x) <= p.d + TOL * (1 + Math.abs(p.d)))) addVertex(x);
    }
    const faces = [];
    for (const plane of planes) {
        const on = verts.filter((v) => Math.abs(dot(plane.n, v) - plane.d) < 1e-5 * (1 + Math.abs(plane.d)));
        if (on.length < 3) continue;
        const c = on.reduce((acc, v) => [acc[0] + v[0] / on.length, acc[1] + v[1] / on.length, acc[2] + v[2] / on.length], [0, 0, 0]);
        const u = unit(sub(on[0], c));
        const w = cross(plane.n, u);
        const points = on
            .map((v) => ({ v, a: Math.atan2(dot(sub(v, c), w), dot(sub(v, c), u)) }))
            .sort((p, q) => p.a - q.a)
            .map((p) => p.v);
        let area = 0;
        for (let t = 1; t + 1 < points.length; t++) area += norm(cross(sub(points[t], points[0]), sub(points[t + 1], points[0]))) / 2;
        faces.push({ plane, points, area });
    }
    return { vertices: verts, faces };
}
