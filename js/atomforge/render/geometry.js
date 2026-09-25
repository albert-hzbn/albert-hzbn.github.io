// Geometry helpers for overlays.
import { sub, dot, cross, unit } from '../core/index.js';

// Intersection of the plane h·f1 + k·f2 + l·f3 = d with the box [0,R]³ (fractional).
export function planePolygon([h, k, l], d, R) {
    const corners = [];
    for (let i = 0; i < 8; i++) corners.push([(i & 1) * R[0], ((i >> 1) & 1) * R[1], ((i >> 2) & 1) * R[2]]);
    const edges = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
    const f = (p) => h * p[0] + k * p[1] + l * p[2] - d;
    const pts = [];
    for (const [i, j] of edges) {
        const a = corners[i], b = corners[j], fa = f(a), fb = f(b);
        if (Math.abs(fa) < 1e-9) pts.push(a);
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0)) {
            const t = fa / (fa - fb);
            pts.push(a.map((v, n) => v + (b[n] - v) * t));
        }
    }
    const uniq = [];
    for (const p of pts) if (!uniq.some((q) => Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]) < 1e-6)) uniq.push(p);
    if (uniq.length < 3) return uniq;
    // Sort around the centroid in the plane.
    const c = uniq.reduce((acc, p) => acc.map((v, n) => v + p[n] / uniq.length), [0, 0, 0]);
    const n = [h, k, l];
    const u = unit(sub(uniq[0], c));
    const w = unit(cross(n, u));
    return uniq.sort((p, q) => Math.atan2(dot(sub(p, c), w), dot(sub(p, c), u)) - Math.atan2(dot(sub(q, c), w), dot(sub(q, c), u)));
}

