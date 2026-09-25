// WebGPU renderers for the microstructure: a 3D ray-marching view (box
// surfaces, volume or isosurface, with an optional cut along the slice plane)
// and a 2D view of an arbitrary slice through the volume.
import { OrbitCamera } from './camera.js';

const VOLUME_WGSL = /* wgsl */`
struct U {
    invViewProj: mat4x4<f32>,
    eye: vec4<f32>,
    plane: vec4<f32>,          // xyz normal, w offset (box coordinates, centred)
    opts: vec4<f32>,           // x mode (0 surface, 1 volume, 2 iso), y iso, z opacity, w gamma
    flags: vec4<f32>,          // x clip, y show plane, z steps, w flip clip side
    bg: vec4<f32>,
    accent: vec4<f32>,
};
@group(0) @binding(0) var<uniform> u: U;
@group(0) @binding(1) var vol: texture_3d<f32>;
@group(0) @binding(2) var samp: sampler;

struct VSOut { @builtin(position) pos: vec4<f32>, @location(0) ndc: vec2<f32> };

@vertex fn vs(@builtin(vertex_index) i: u32) -> VSOut {
    var p = array<vec2<f32>, 3>(vec2<f32>(-1.0, -1.0), vec2<f32>(3.0, -1.0), vec2<f32>(-1.0, 3.0));
    var o: VSOut;
    o.pos = vec4<f32>(p[i], 0.0, 1.0);
    o.ndc = p[i];
    return o;
}

fn sampleAt(q: vec3<f32>) -> vec4<f32> { return textureSampleLevel(vol, samp, q + vec3<f32>(0.5), 0.0); }

fn gradA(q: vec3<f32>, h: f32) -> vec3<f32> {
    return vec3<f32>(
        sampleAt(q + vec3<f32>(h, 0.0, 0.0)).a - sampleAt(q - vec3<f32>(h, 0.0, 0.0)).a,
        sampleAt(q + vec3<f32>(0.0, h, 0.0)).a - sampleAt(q - vec3<f32>(0.0, h, 0.0)).a,
        sampleAt(q + vec3<f32>(0.0, 0.0, h)).a - sampleAt(q - vec3<f32>(0.0, 0.0, h)).a);
}

fn shade(c: vec3<f32>, n: vec3<f32>, view: vec3<f32>) -> vec3<f32> {
    let l = normalize(vec3<f32>(0.45, 0.35, 0.82));
    let diff = abs(dot(n, l));
    let h = normalize(l + view);
    let spec = pow(max(abs(dot(n, h)), 0.0), 40.0);
    return c * (0.32 + 0.72 * diff) + vec3<f32>(0.18 * spec);
}

// Face normal of the unit box at a surface point.
fn boxNormal(q: vec3<f32>) -> vec3<f32> {
    let a = abs(q);
    if (a.x > a.y && a.x > a.z) { return vec3<f32>(sign(q.x), 0.0, 0.0); }
    if (a.y > a.z) { return vec3<f32>(0.0, sign(q.y), 0.0); }
    return vec3<f32>(0.0, 0.0, sign(q.z));
}

@fragment fn fs(in: VSOut) -> @location(0) vec4<f32> {
    let near = u.invViewProj * vec4<f32>(in.ndc, 0.0, 1.0);
    let far = u.invViewProj * vec4<f32>(in.ndc, 1.0, 1.0);
    let ro = near.xyz / near.w;
    let rd = normalize(far.xyz / far.w - ro);
    // Ray / box [-0.5, 0.5]^3
    let inv = 1.0 / rd;
    let t0 = (vec3<f32>(-0.5) - ro) * inv;
    let t1 = (vec3<f32>(0.5) - ro) * inv;
    let tmin = min(t0, t1);
    let tmax = max(t0, t1);
    var tn = max(max(tmin.x, tmin.y), tmin.z);
    var tf = min(min(tmax.x, tmax.y), tmax.z);
    let bg = u.bg.rgb;
    if (tf <= max(tn, 0.0)) { return vec4<f32>(bg, 1.0); }
    tn = max(tn, 0.0);
    var cutNormal = vec3<f32>(0.0);
    var entryIsCut = false;
    // Cut away the half-space in front of the slice plane.
    let pn = u.plane.xyz * select(1.0, -1.0, u.flags.w > 0.5);
    let pd = u.plane.w * select(1.0, -1.0, u.flags.w > 0.5);
    if (u.flags.x > 0.5) {
        let denom = dot(pn, rd);
        let s0 = dot(pn, ro + rd * tn) - pd;
        let tp = (pd - dot(pn, ro)) / denom;
        if (s0 > 0.0) {
            // Entry point is in the removed half: start at the plane.
            if (denom >= 0.0 || tp > tf) { return vec4<f32>(bg, 1.0); }
            tn = max(tn, tp);
            entryIsCut = true;
            cutNormal = -pn;
        } else if (denom > 0.0 && tp < tf) {
            tf = tp;
        }
    }
    let view = -rd;
    let mode = u32(u.opts.x + 0.5);
    let steps = u.flags.z;
    let dt = 1.7320508 / steps;
    var col = vec3<f32>(0.0);
    var alpha = 0.0;

    if (mode == 0u) {
        // Coloured box faces (and the cut face).
        let q = ro + rd * (tn + 1e-4);
        let s = sampleAt(q);
        let n = select(boxNormal(q), cutNormal, entryIsCut);
        col = shade(s.rgb, n, view);
        alpha = 1.0;
    } else if (mode == 1u) {
        // Emission-absorption volume rendering, front to back.
        var t = tn + dt * 0.5;
        loop {
            if (t > tf || alpha > 0.985) { break; }
            let s = sampleAt(ro + rd * t);
            let a = clamp(pow(s.a, u.opts.w) * u.opts.z * dt * 60.0, 0.0, 1.0);
            col += (1.0 - alpha) * a * s.rgb;
            alpha += (1.0 - alpha) * a;
            t += dt;
        }
        col = col + (1.0 - alpha) * bg;
        alpha = 1.0;
    } else {
        // Isosurface of the scalar channel, refined by bisection.
        let iso = u.opts.y;
        var t = tn;
        var prev = sampleAt(ro + rd * t).a;
        if (entryIsCut && prev >= iso) {
            let s = sampleAt(ro + rd * (t + 1e-4));
            col = shade(s.rgb, cutNormal, view);
            alpha = 1.0;
        } else {
            loop {
                t += dt;
                if (t > tf) { break; }
                let cur = sampleAt(ro + rd * t).a;
                if (cur >= iso) {
                    var a = t - dt;
                    var b = t;
                    for (var k = 0; k < 6; k++) {
                        let m = 0.5 * (a + b);
                        if (sampleAt(ro + rd * m).a >= iso) { b = m; } else { a = m; }
                    }
                    let q = ro + rd * b;
                    let g = gradA(q, 1.0 / 128.0);
                    let n = select(-normalize(g), boxNormal(q), length(g) < 1e-5);
                    col = shade(sampleAt(q).rgb, n, view);
                    alpha = 1.0;
                    break;
                }
                prev = cur;
            }
        }
        if (alpha < 0.5) { col = bg; alpha = 1.0; }
    }

    // Translucent slice plane where the ray crosses it inside the box.
    if (u.flags.y > 0.5 && u.flags.x < 0.5) {
        let denom = dot(u.plane.xyz, rd);
        if (abs(denom) > 1e-5) {
            let tp = (u.plane.w - dot(u.plane.xyz, ro)) / denom;
            let q = ro + rd * tp;
            if (tp > 0.0 && all(abs(q) <= vec3<f32>(0.5))) {
                col = mix(col, u.accent.rgb, 0.22);
            }
        }
    }
    return vec4<f32>(col, 1.0);
}
`;

const LINES_WGSL = /* wgsl */`
struct U { viewProj: mat4x4<f32> };
@group(0) @binding(0) var<uniform> u: U;
struct VIn { @location(0) pos: vec3<f32>, @location(1) color: vec4<f32> };
struct VOut { @builtin(position) pos: vec4<f32>, @location(0) color: vec4<f32> };
@vertex fn vs(v: VIn) -> VOut {
    var o: VOut;
    o.pos = u.viewProj * vec4<f32>(v.pos, 1.0);
    o.color = v.color;
    return o;
}
@fragment fn fs(v: VOut) -> @location(0) vec4<f32> { return v.color; }
`;

const SLICE_WGSL = /* wgsl */`
struct U {
    origin: vec4<f32>,
    e1: vec4<f32>,
    e2: vec4<f32>,
    extent: vec4<f32>,         // u0, u1, v0, v1 in plane coordinates
    bg: vec4<f32>,
    opts: vec4<f32>,           // x: grid lines (1 = show box outline)
};
@group(0) @binding(0) var<uniform> u: U;
@group(0) @binding(1) var vol: texture_3d<f32>;
@group(0) @binding(2) var samp: sampler;
struct VSOut { @builtin(position) pos: vec4<f32>, @location(0) uv: vec2<f32> };
@vertex fn vs(@builtin(vertex_index) i: u32) -> VSOut {
    var p = array<vec2<f32>, 3>(vec2<f32>(-1.0, -1.0), vec2<f32>(3.0, -1.0), vec2<f32>(-1.0, 3.0));
    var o: VSOut;
    o.pos = vec4<f32>(p[i], 0.0, 1.0);
    o.uv = p[i] * 0.5 + vec2<f32>(0.5);
    return o;
}
@fragment fn fs(in: VSOut) -> @location(0) vec4<f32> {
    let a = mix(u.extent.x, u.extent.y, in.uv.x);
    let b = mix(u.extent.z, u.extent.w, in.uv.y);
    let q = u.origin.xyz + a * u.e1.xyz + b * u.e2.xyz;
    if (any(abs(q) > vec3<f32>(0.5))) { return vec4<f32>(u.bg.rgb, 1.0); }
    return vec4<f32>(textureSampleLevel(vol, samp, q + vec3<f32>(0.5), 0.0).rgb, 1.0);
}
`;

// Plane through the box: unit normal n and offset d along it (box centred at 0).
export function planeBasis(n) {
    const a = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    // e1 ⟂ n, as "horizontal" as possible; e2 = n × e1.
    let e1 = [a[1] * n[2] - a[2] * n[1], a[2] * n[0] - a[0] * n[2], a[0] * n[1] - a[1] * n[0]];
    const l1 = Math.hypot(...e1); e1 = e1.map((v) => v / l1);
    const e2 = [n[1] * e1[2] - n[2] * e1[1], n[2] * e1[0] - n[0] * e1[2], n[0] * e1[1] - n[1] * e1[0]];
    return { e1, e2 };
}

// Polygon where the plane cuts the box (ordered), in 3D.
export function planePolygon(n, d) {
    const corners = [];
    for (let i = 0; i < 8; i++) corners.push([(i & 1) - 0.5, ((i >> 1) & 1) - 0.5, ((i >> 2) & 1) - 0.5]);
    const edges = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
    const f = (p) => n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - d;
    const pts = [];
    for (const [i, j] of edges) {
        const A = corners[i], B = corners[j], fa = f(A), fb = f(B);
        if ((fa <= 0 && fb > 0) || (fa > 0 && fb <= 0)) {
            const t = fa / (fa - fb);
            pts.push(A.map((v, k) => v + (B[k] - v) * t));
        }
    }
    if (pts.length < 3) return [];
    const c = pts.reduce((acc, p) => acc.map((v, k) => v + p[k] / pts.length), [0, 0, 0]);
    const { e1, e2 } = planeBasis(n);
    const ang = (p) => Math.atan2((p[0] - c[0]) * e2[0] + (p[1] - c[1]) * e2[1] + (p[2] - c[2]) * e2[2], (p[0] - c[0]) * e1[0] + (p[1] - c[1]) * e1[1] + (p[2] - c[2]) * e1[2]);
    return pts.sort((p, q) => ang(p) - ang(q));
}

export class Renderer {
    constructor(device, format) {
        this.device = device;
        this.format = format;
        this.camera = new OrbitCamera();
        this.sampler = device.createSampler({ magFilter: 'linear', minFilter: 'linear', addressModeU: 'clamp-to-edge', addressModeV: 'clamp-to-edge', addressModeW: 'clamp-to-edge' });
        this.nearest = device.createSampler({ magFilter: 'nearest', minFilter: 'nearest' });
        this.volUniform = device.createBuffer({ size: 176, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.lineUniform = device.createBuffer({ size: 64, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.sliceUniform = device.createBuffer({ size: 96, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.lineBuf = device.createBuffer({ size: 64 * 28 * 4, usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });

        const volModule = device.createShaderModule({ code: VOLUME_WGSL, label: 'volume' });
        this.volPipeline = device.createRenderPipeline({
            layout: 'auto',
            vertex: { module: volModule, entryPoint: 'vs' },
            fragment: { module: volModule, entryPoint: 'fs', targets: [{ format }] },
            primitive: { topology: 'triangle-list' },
        });
        const lineModule = device.createShaderModule({ code: LINES_WGSL, label: 'lines' });
        this.linePipeline = device.createRenderPipeline({
            layout: 'auto',
            vertex: {
                module: lineModule, entryPoint: 'vs',
                buffers: [{ arrayStride: 28, attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x3' }, { shaderLocation: 1, offset: 12, format: 'float32x4' }] }],
            },
            fragment: {
                module: lineModule, entryPoint: 'fs',
                targets: [{ format, blend: { color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } } }],
            },
            primitive: { topology: 'line-list' },
        });
        this.lineGroup = device.createBindGroup({ layout: this.linePipeline.getBindGroupLayout(0), entries: [{ binding: 0, resource: { buffer: this.lineUniform } }] });
        const sliceModule = device.createShaderModule({ code: SLICE_WGSL, label: 'slice' });
        this.slicePipeline = device.createRenderPipeline({
            layout: 'auto',
            vertex: { module: sliceModule, entryPoint: 'vs' },
            fragment: { module: sliceModule, entryPoint: 'fs', targets: [{ format }] },
            primitive: { topology: 'triangle-list' },
        });
    }

    setTexture(texture) {
        this.texture = texture;
        const view = texture.createView();
        this.volGroup = this.device.createBindGroup({
            layout: this.volPipeline.getBindGroupLayout(0),
            entries: [{ binding: 0, resource: { buffer: this.volUniform } }, { binding: 1, resource: view }, { binding: 2, resource: this.sampler }],
        });
        this.sliceGroups = [this.sampler, this.nearest].map((s) => this.device.createBindGroup({
            layout: this.slicePipeline.getBindGroupLayout(0),
            entries: [{ binding: 0, resource: { buffer: this.sliceUniform } }, { binding: 1, resource: view }, { binding: 2, resource: s }],
        }));
    }

    // 3D view. o: { mode, iso, opacity, gamma, clip, flip, showPlane, plane:{n,d}, bg, accent, lineColor, steps }
    render3D(context, canvas, o) {
        if (!this.volGroup) return;
        const aspect = canvas.width / canvas.height;
        const m = this.camera.matrices(aspect);
        const u = new Float32Array(44);
        u.set(m.invViewProj, 0);
        u.set([...m.eye, 1], 16);
        u.set([...o.plane.n, o.plane.d], 20);
        u.set([o.mode, o.iso, o.opacity, o.gamma], 24);
        u.set([o.clip ? 1 : 0, o.showPlane ? 1 : 0, o.steps || 256, o.flip ? 1 : 0], 28);
        u.set([...o.bg, 1], 32);
        u.set([...o.accent, 1], 36);
        this.device.queue.writeBuffer(this.volUniform, 0, u);
        this.device.queue.writeBuffer(this.lineUniform, 0, m.viewProj);

        // Box edges and the slice outline.
        const verts = [];
        const C = [];
        for (let i = 0; i < 8; i++) C.push([(i & 1) - 0.5, ((i >> 1) & 1) - 0.5, ((i >> 2) & 1) - 0.5]);
        const E = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
        const lc = o.lineColor;
        if (o.showBox) for (const [a, b] of E) verts.push(...C[a], ...lc, ...C[b], ...lc);
        if (o.showPlane) {
            const poly = planePolygon(o.plane.n, o.plane.d);
            for (let k = 0; k < poly.length; k++) verts.push(...poly[k], ...o.accent, 1, ...poly[(k + 1) % poly.length], ...o.accent, 1);
        }
        const lineData = new Float32Array(verts);
        this.device.queue.writeBuffer(this.lineBuf, 0, lineData);

        const enc = this.device.createCommandEncoder();
        const pass = enc.beginRenderPass({ colorAttachments: [{ view: context.getCurrentTexture().createView(), loadOp: 'clear', storeOp: 'store', clearValue: { r: o.bg[0], g: o.bg[1], b: o.bg[2], a: 1 } }] });
        pass.setPipeline(this.volPipeline);
        pass.setBindGroup(0, this.volGroup);
        pass.draw(3);
        if (lineData.length) {
            pass.setPipeline(this.linePipeline);
            pass.setBindGroup(0, this.lineGroup);
            pass.setVertexBuffer(0, this.lineBuf);
            pass.draw(lineData.length / 7);
        }
        pass.end();
        this.device.queue.submit([enc.finish()]);
    }

    // 2D slice view; returns the plane extent actually shown (for rulers).
    renderSlice(context, canvas, { plane, bg, smooth = true }) {
        if (!this.sliceGroups) return null;
        const { e1, e2 } = planeBasis(plane.n);
        const origin = plane.n.map((v) => v * plane.d);
        const poly = planePolygon(plane.n, plane.d);
        let u0 = -0.5, u1 = 0.5, v0 = -0.5, v1 = 0.5;
        if (poly.length) {
            const us = poly.map((p) => (p[0] - origin[0]) * e1[0] + (p[1] - origin[1]) * e1[1] + (p[2] - origin[2]) * e1[2]);
            const vs = poly.map((p) => (p[0] - origin[0]) * e2[0] + (p[1] - origin[1]) * e2[1] + (p[2] - origin[2]) * e2[2]);
            [u0, u1, v0, v1] = [Math.min(...us), Math.max(...us), Math.min(...vs), Math.max(...vs)];
        }
        // Keep the aspect ratio of the canvas; add a small margin.
        const pad = 0.04;
        let w = (u1 - u0) * (1 + pad), h = (v1 - v0) * (1 + pad);
        const cu = (u0 + u1) / 2, cv = (v0 + v1) / 2;
        const aspect = canvas.width / canvas.height;
        if (w / h > aspect) h = w / aspect; else w = h * aspect;
        const ext = [cu - w / 2, cu + w / 2, cv - h / 2, cv + h / 2];
        const data = new Float32Array(24);
        data.set([...origin, 0, ...e1, 0, ...e2, 0, ...ext, ...bg, 1, 0, 0, 0, 0]);
        this.device.queue.writeBuffer(this.sliceUniform, 0, data);
        const enc = this.device.createCommandEncoder();
        const pass = enc.beginRenderPass({ colorAttachments: [{ view: context.getCurrentTexture().createView(), loadOp: 'clear', storeOp: 'store', clearValue: { r: bg[0], g: bg[1], b: bg[2], a: 1 } }] });
        pass.setPipeline(this.slicePipeline);
        pass.setBindGroup(0, this.sliceGroups[smooth ? 0 : 1]);
        pass.draw(3);
        pass.end();
        this.device.queue.submit([enc.finish()]);
        return { extent: ext, width: w };
    }
}
