// Exports: PNG images of the views, VTK volumes for ParaView, and CSV slices.
import { planeBasis, planePolygon } from '../render/renderer.js';

export function download(name, blobOrText, type = 'application/octet-stream') {
    const blob = blobOrText instanceof Blob ? blobOrText : new Blob([blobOrText], { type });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}

// Render, then capture in the same task so the WebGPU canvas still holds the frame.
export function exportCanvas(canvas, renderFn, name) {
    renderFn();
    return new Promise((resolve) => canvas.toBlob((blob) => { if (blob) download(name, blob, 'image/png'); resolve(!!blob); }, 'image/png'));
}

// Legacy VTK (binary, big-endian) STRUCTURED_POINTS with one scalar array per
// component of every state field.
export async function exportVTK(sim, name, spacing = 1) {
    const n = sim.n, n3 = n ** 3;
    const parts = [`# vtk DataFile Version 3.0\n${name} (Phase-Field Lab, step ${sim.step})\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS ${n} ${n} ${n}\nORIGIN 0 0 0\nSPACING ${spacing} ${spacing} ${spacing}\nPOINT_DATA ${n3}\n`];
    for (const f of sim.fields) {
        const data = await sim.readField(f.name);
        for (let k = 0; k < f.comps; k++) {
            parts.push(`SCALARS ${f.name}${f.comps > 1 ? '_' + k : ''} float 1\nLOOKUP_TABLE default\n`);
            const buf = new ArrayBuffer(n3 * 4);
            const view = new DataView(buf);
            for (let i = 0; i < n3; i++) view.setFloat32(i * 4, data[k * n3 + i], false);
            parts.push(buf, '\n');
        }
    }
    download(`${name}.vtk`, new Blob(parts));
}

// Sample a field component on the slice plane (trilinear, periodic) and save CSV.
export async function exportSliceCSV(sim, fieldName, comp, plane, name, res = 200) {
    const n = sim.n, n3 = n ** 3;
    const data = await sim.readField(fieldName);
    const { e1, e2 } = planeBasis(plane.n);
    const origin = plane.n.map((v) => v * plane.d);
    const poly = planePolygon(plane.n, plane.d);
    if (!poly.length) throw new Error('The plane does not cut the box.');
    const proj = (p, e) => (p[0] - origin[0]) * e[0] + (p[1] - origin[1]) * e[1] + (p[2] - origin[2]) * e[2];
    const us = poly.map((p) => proj(p, e1)), vs = poly.map((p) => proj(p, e2));
    const [u0, u1, v0, v1] = [Math.min(...us), Math.max(...us), Math.min(...vs), Math.max(...vs)];
    const at = (x, y, z) => data[comp * n3 + ((((z % n) + n) % n) * n + (((y % n) + n) % n)) * n + (((x % n) + n) % n)];
    const rows = ['u,v,x,y,z,value'];
    for (let b = 0; b < res; b++) for (let a = 0; a < res; a++) {
        const u = u0 + (u1 - u0) * (a + 0.5) / res, v = v0 + (v1 - v0) * (b + 0.5) / res;
        const q = [0, 1, 2].map((k) => origin[k] + u * e1[k] + v * e2[k]);
        if (q.some((c) => Math.abs(c) > 0.5)) continue;
        const g = q.map((c) => (c + 0.5) * n - 0.5);
        const i0 = g.map(Math.floor), t = g.map((c, k) => c - i0[k]);
        let val = 0;
        for (let dz = 0; dz < 2; dz++) for (let dy = 0; dy < 2; dy++) for (let dx = 0; dx < 2; dx++) {
            const w = (dx ? t[0] : 1 - t[0]) * (dy ? t[1] : 1 - t[1]) * (dz ? t[2] : 1 - t[2]);
            val += w * at(i0[0] + dx, i0[1] + dy, i0[2] + dz);
        }
        rows.push([u * n, v * n, g[0], g[1], g[2], val].map((x) => +x.toFixed(5)).join(','));
    }
    download(`${name}_slice.csv`, rows.join('\n'), 'text/csv');
}
