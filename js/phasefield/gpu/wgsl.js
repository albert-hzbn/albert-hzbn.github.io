// WGSL code generation for phase-field models.
//
// A model declares its state fields (double-buffered), auxiliary fields and
// parameters; this module turns that into shader modules with typed helpers:
//   <f>(k, i)          value of component k of field f at cell i
//   <f>_at(k, p)       value at integer position p (periodic)
//   <f>_lap(k, p)      7-point Laplacian (divided by dx²)
//   <f>_grad(k, p)     central-difference gradient (divided by 2dx)
//   <f>_set(k, i, v)   write the new value (state: next buffer; aux: in place)
// Parameters appear as `let` constants inside every kernel. Names starting
// with pf_ are reserved for generated code.

export const PARAM_SLOTS = 60;          // user parameters; slot 63 holds 1/dx²
export const WORKGROUP = 4;             // 4×4×4 threads

// ---------------------------------------------------------------------------
// Shared prelude: grid indexing, random numbers, colour maps
// ---------------------------------------------------------------------------

const PRELUDE = /* wgsl */`
struct Params { v: array<vec4<f32>, 16> };
struct Grid { n: u32, step: u32, seed: u32, pad: u32 };

fn n3() -> u32 { return G.n * G.n * G.n; }
fn inv_dx2() -> f32 { return P.v[15].w; }
fn inv_2dx() -> f32 { return 0.5 * sqrt(P.v[15].w); }

fn cell(p: vec3<i32>) -> u32 {
    let n = i32(G.n);
    let q = ((p % n) + n) % n;
    return u32((q.z * n + q.y) * n + q.x);
}

// PCG hash; rand(i, salt) is uniform in [0, 1) and changes every step.
fn pcg(v: u32) -> u32 {
    let s = v * 747796405u + 2891336453u;
    let w = ((s >> ((s >> 28u) + 4u)) ^ s) * 277803737u;
    return (w >> 22u) ^ w;
}
fn rand(i: u32, salt: u32) -> f32 {
    return f32(pcg(i ^ pcg(G.step * 7919u + salt * 104729u + G.seed))) / 4294967296.0;
}
`;

const COLORMAPS = /* wgsl */`
fn cm_viridis(t: f32) -> vec3<f32> {
    let x = clamp(t, 0.0, 1.0);
    let c0 = vec3<f32>(0.2777273272, 0.0054073445, 0.3340998053);
    let c1 = vec3<f32>(0.1050930431, 1.4046135299, 1.3845901626);
    let c2 = vec3<f32>(-0.3308618287, 0.2148475595, 0.0950951630);
    let c3 = vec3<f32>(-4.6342304990, -5.7991009734, -19.3324409563);
    let c4 = vec3<f32>(6.2282699363, 14.1799333668, 56.6905526007);
    let c5 = vec3<f32>(4.7763849977, -13.7451453777, -65.3530326334);
    let c6 = vec3<f32>(-5.4354558559, 4.6458526122, 26.3124352496);
    return clamp(c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * c6))))), vec3<f32>(0.0), vec3<f32>(1.0));
}
fn cm_magma(t: f32) -> vec3<f32> {
    let x = clamp(t, 0.0, 1.0);
    let c0 = vec3<f32>(-0.0021364851, -0.0007496551, -0.0053861279);
    let c1 = vec3<f32>(0.2516605407, 0.6775232437, 2.4940265993);
    let c2 = vec3<f32>(8.3537172792, -3.5777195150, 0.3144679030);
    let c3 = vec3<f32>(-27.6687330858, 14.2647307810, -13.6492131881);
    let c4 = vec3<f32>(52.1761398123, -27.9436060717, 12.9441694424);
    let c5 = vec3<f32>(-50.7685253647, 29.0465828213, 4.2341529938);
    let c6 = vec3<f32>(18.6557050659, -11.4897735200, -5.6019615087);
    return clamp(c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * c6))))), vec3<f32>(0.0), vec3<f32>(1.0));
}
fn cm_turbo(t: f32) -> vec3<f32> {
    let x = clamp(t, 0.0, 1.0);
    let v4 = vec4<f32>(1.0, x, x * x, x * x * x);
    let v2 = v4.zw * v4.z;
    return clamp(vec3<f32>(
        dot(v4, vec4<f32>(0.13572138, 4.61539260, -42.66032258, 132.13108234)) + dot(v2, vec2<f32>(-152.94239396, 59.28637943)),
        dot(v4, vec4<f32>(0.09140261, 2.19418839, 4.84296658, -14.18503333)) + dot(v2, vec2<f32>(4.27729857, 2.82956604)),
        dot(v4, vec4<f32>(0.10667330, 12.64194608, -60.58204836, 110.36276771)) + dot(v2, vec2<f32>(-89.90310912, 27.34824973))
    ), vec3<f32>(0.0), vec3<f32>(1.0));
}
fn cm_coolwarm(t: f32) -> vec3<f32> {
    let x = clamp(t, 0.0, 1.0);
    let cold = vec3<f32>(0.230, 0.299, 0.754);
    let mid = vec3<f32>(0.865, 0.865, 0.865);
    let warm = vec3<f32>(0.706, 0.016, 0.150);
    if (x < 0.5) { return mix(cold, mid, smoothstep(0.0, 1.0, x * 2.0)); }
    return mix(mid, warm, smoothstep(0.0, 1.0, x * 2.0 - 1.0));
}
fn cm_copper(t: f32) -> vec3<f32> {
    let x = clamp(t, 0.0, 1.0);
    return clamp(vec3<f32>(1.25 * x, 0.7812 * x, 0.4975 * x), vec3<f32>(0.0), vec3<f32>(1.0));
}
fn cm_gray(t: f32) -> vec3<f32> { return vec3<f32>(clamp(t, 0.0, 1.0)); }

// Colour map chosen in the Display panel: D.a.z = 0 viridis, 1 magma, 2 turbo,
// 3 coolwarm, 4 copper, 5 grey. Values are first normalised to [D.a.x, D.a.y].
fn scalar01(v: f32) -> f32 { return clamp((v - D.a.x) / max(D.a.y - D.a.x, 1e-6), 0.0, 1.0); }
fn cmap(t: f32) -> vec3<f32> {
    let id = u32(D.a.z + 0.5);
    switch id {
        case 1u: { return cm_magma(t); }
        case 2u: { return cm_turbo(t); }
        case 3u: { return cm_coolwarm(t); }
        case 4u: { return cm_copper(t); }
        case 5u: { return cm_gray(t); }
        default: { return cm_viridis(t); }
    }
}
fn hsv2rgb(h: f32, s: f32, v: f32) -> vec3<f32> {
    let k = vec3<f32>(1.0, 2.0 / 3.0, 1.0 / 3.0);
    let p = abs(fract(vec3<f32>(h) + k) * 6.0 - vec3<f32>(3.0));
    return v * mix(vec3<f32>(1.0), clamp(p - vec3<f32>(1.0), vec3<f32>(0.0), vec3<f32>(1.0)), s);
}
// Distinct, soft colours for integer labels (grains, variants, phases).
fn label_color(id: u32) -> vec3<f32> {
    let h = fract(f32(pcg(id * 2654435761u + 17u) % 1000u) / 1000.0 + f32(id) * 0.61803398875);
    return hsv2rgb(h, 0.55, 0.92);
}
`;

// ---------------------------------------------------------------------------
// Helpers per field
// ---------------------------------------------------------------------------

function fieldHelpers(name, readVar, writeVar) {
    let s = `
fn ${name}(k: u32, i: u32) -> f32 { return ${readVar}[k * n3() + i]; }
fn ${name}_at(k: u32, p: vec3<i32>) -> f32 { return ${name}(k, cell(p)); }
fn ${name}_lap(k: u32, p: vec3<i32>) -> f32 {
    let c = ${name}_at(k, p);
    let s = ${name}_at(k, p + vec3<i32>(1, 0, 0)) + ${name}_at(k, p - vec3<i32>(1, 0, 0))
          + ${name}_at(k, p + vec3<i32>(0, 1, 0)) + ${name}_at(k, p - vec3<i32>(0, 1, 0))
          + ${name}_at(k, p + vec3<i32>(0, 0, 1)) + ${name}_at(k, p - vec3<i32>(0, 0, 1));
    return (s - 6.0 * c) * inv_dx2();
}
fn ${name}_grad(k: u32, p: vec3<i32>) -> vec3<f32> {
    return vec3<f32>(
        ${name}_at(k, p + vec3<i32>(1, 0, 0)) - ${name}_at(k, p - vec3<i32>(1, 0, 0)),
        ${name}_at(k, p + vec3<i32>(0, 1, 0)) - ${name}_at(k, p - vec3<i32>(0, 1, 0)),
        ${name}_at(k, p + vec3<i32>(0, 0, 1)) - ${name}_at(k, p - vec3<i32>(0, 0, 1))) * inv_2dx();
}
`;
    if (writeVar) s += `fn ${name}_set(k: u32, i: u32, v: f32) { ${writeVar}[k * n3() + i] = v; }\n`;
    return s;
}

// `let` declarations for all parameters, inserted at the top of every kernel.
export function paramLets(model) {
    return model.params.map((p, j) => `    let ${p.key} = P.v[${j >> 2}][${j & 3}];`).join('\n');
}

export function paramIndex(model, key) {
    return model.params.findIndex((p) => p.key === key);
}

// Binding layout shared by the step module (see Simulation for the JS side):
// 0 P, 1 G (dynamic offset), then per state field [in, out], then aux fields.
export function stepModule(model) {
    let b = 2;
    let decl = `@group(0) @binding(0) var<uniform> P: Params;\n@group(0) @binding(1) var<uniform> G: Grid;\n`;
    let helpers = '';
    for (const f of model.fields) {
        decl += `@group(0) @binding(${b++}) var<storage, read> ${f.name}_in: array<f32>;\n`;
        decl += `@group(0) @binding(${b++}) var<storage, read_write> ${f.name}_out: array<f32>;\n`;
        helpers += fieldHelpers(f.name, `${f.name}_in`, `${f.name}_out`);
    }
    for (const a of model.aux || []) {
        decl += `@group(0) @binding(${b++}) var<storage, read_write> ${a.name}_buf: array<f32>;\n`;
        helpers += fieldHelpers(a.name, `${a.name}_buf`, `${a.name}_buf`);
    }
    // Colour maps reference D; declare a dummy in the step module.
    const kernels = model.passes.map((ps) => kernel(ps.entry, model, ps.wgsl)).join('\n');
    return PRELUDE + decl + helpers + (model.wgslCommon || '') + kernels;
}

// Display module: 0 P, 1 G, 2 D (display settings), fields (read), aux (read), texture.
export function displayModule(model, view) {
    let b = 3;
    let decl = `struct Display { a: vec4<f32>, b: vec4<f32> };
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<uniform> G: Grid;
@group(0) @binding(2) var<uniform> D: Display;\n`;
    let helpers = '';
    for (const f of model.fields) {
        decl += `@group(0) @binding(${b++}) var<storage, read> ${f.name}_in: array<f32>;\n`;
        helpers += fieldHelpers(f.name, `${f.name}_in`, null);
    }
    for (const a of model.aux || []) {
        decl += `@group(0) @binding(${b++}) var<storage, read> ${a.name}_buf: array<f32>;\n`;
        helpers += fieldHelpers(a.name, `${a.name}_buf`, null);
    }
    decl += `@group(0) @binding(${b}) var out_tex: texture_storage_3d<rgba8unorm, write>;\n`;
    const fn = `
fn view_value(p: vec3<i32>, i: u32) -> vec4<f32> {
${paramLets(model)}
${view.wgsl}
}
@compute @workgroup_size(${WORKGROUP}, ${WORKGROUP}, ${WORKGROUP})
fn display(@builtin(global_invocation_id) pf_id: vec3<u32>) {
    if (any(pf_id >= vec3<u32>(G.n))) { return; }
    let p = vec3<i32>(pf_id);
    textureStore(out_tex, pf_id, clamp(view_value(p, cell(p)), vec4<f32>(0.0), vec4<f32>(1.0)));
}
`;
    return PRELUDE + COLORMAPS + decl + helpers + (model.wgslCommon || '') + fn;
}

function kernel(entry, model, body) {
    return `
@compute @workgroup_size(${WORKGROUP}, ${WORKGROUP}, ${WORKGROUP})
fn ${entry}(@builtin(global_invocation_id) pf_id: vec3<u32>) {
    if (any(pf_id >= vec3<u32>(G.n))) { return; }
    let p = vec3<i32>(pf_id);
    let i = cell(p);
${paramLets(model)}
${body}
}
`;
}
