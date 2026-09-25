// Nucleation and growth: a metastable parent transforms to a product phase by
// random nucleation (constant rate) and interface-controlled growth, the
// situation described by Johnson–Mehl–Avrami–Kolmogorov (JMAK) kinetics.
// Each transformed cell records when it transformed and which nucleus it grew
// from, so the product can be shown as a Johnson–Mehl grain structure or
// coloured by transformation time.
export default {
    id: 'nucleation',
    name: 'Nucleation & growth (JMAK)',
    short: 'Nucleation',
    icon: '<circle cx="12" cy="12" r="8" opacity=".45"/><circle cx="12" cy="12" r="4"/><circle cx="12" cy="12" r="1.2" class="fill"/>',
    description: 'Supercritical nuclei appear at random in the untransformed parent and grow at constant velocity until they impinge. The transformed fraction follows the Avrami law X = 1 − exp(−k t⁴) for constant nucleation and growth rates in 3D.',
    equations: [
        String.raw`\frac{\partial \phi}{\partial t} = -L\left(\frac{\partial f}{\partial \phi} - \kappa\nabla^2\phi\right)`,
        String.raw`\begin{aligned} \frac{\partial f}{\partial \phi} = {}& 2W\phi(1-\phi)(1-2\phi) \\ & - 6\phi(1-\phi)\,\Delta G \end{aligned}`,
    ],
    equationNote: 'Nuclei are seeded at rate J per unit volume.',
    grid: 128,
    params: [
        { key: 'dG', label: 'Driving force ΔG', value: 0.25, min: 0.05, max: 0.6, step: 0.005, group: 'Thermodynamics', hint: 'Critical radius 2σ/ΔG must stay below the nucleus radius' },
        { key: 'W', label: 'Barrier W', value: 1, min: 0.3, max: 2, step: 0.05, group: 'Thermodynamics' },
        { key: 'kappa', label: 'Gradient energy κ', value: 1, min: 0.3, max: 2, step: 0.05, group: 'Thermodynamics' },
        { key: 'J', label: 'Nucleation rate J', value: 0.00001, min: 0, max: 0.0005, step: 0.00001, group: 'Nucleation', hint: 'Per 8³-cell block per step' },
        { key: 'R0', label: 'Nucleus radius', value: 3, min: 2, max: 6, step: 0.5, group: 'Nucleation' },
        { key: 'L', label: 'Mobility L', value: 1, min: 0.2, max: 2, step: 0.05, group: 'Kinetics' },
        { key: 'dt', label: 'Time step Δt', value: 0.1, min: 0.01, max: 0.15, step: 0.01, group: 'Numerics' },
    ],
    init: [],
    presets: [
        { label: 'Many small grains', params: { J: 0.0001, dG: 0.25 } },
        { label: 'Few large grains', params: { J: 0.000003, dG: 0.25 } },
        { label: 'Slow growth', params: { J: 0.00005, dG: 0.15 } },
    ],
    layout: () => ({ fields: [{ name: 'phi', comps: 1 }, { name: 'tt', comps: 1 }, { name: 'gid', comps: 1 }], aux: [] }),
    wgslCommon: `
const BLOCK: i32 = 8;
// Nucleation event in block b during this step: centre in xyz, and in w a
// grain label + 1 (0 when there is no event).
fn nucleus(b: vec3<i32>, J: f32) -> vec4<f32> {
    let nb = i32(G.n) / BLOCK;
    let q = ((b % nb) + nb) % nb;
    let id = u32((q.z * nb + q.y) * nb + q.x);
    if (rand(id, 11u) >= J) { return vec4<f32>(0.0); }
    let off = vec3<f32>(rand(id, 12u), rand(id, 13u), rand(id, 14u)) * f32(BLOCK);
    let label = f32(pcg(id * 7919u + G.step) & 0xFFFFFu) + 1.0;
    return vec4<f32>(vec3<f32>(b * BLOCK) + off, label);
}
`,
    passes: [
        {
            entry: 'update',
            wgsl: `
    var f = phi(0u, i);
    var id = gid(0u, i);
    // Seed nuclei in untransformed regions.
    if (f < 0.1 && J > 0.0) {
        let b = p / BLOCK;
        for (var dz = -1; dz <= 1; dz++) { for (var dy = -1; dy <= 1; dy++) { for (var dx2 = -1; dx2 <= 1; dx2++) {
            let nu = nucleus(b + vec3<i32>(dx2, dy, dz), J);
            if (nu.w > 0.5) {
                let d = vec3<f32>(p) + vec3<f32>(0.5) - nu.xyz;
                if (dot(d, d) < R0 * R0) { f = 1.0; id = nu.w - 1.0; }
            }
        }}}
    }
    let df = 2.0 * W * f * (1.0 - f) * (1.0 - 2.0 * f) - 6.0 * f * (1.0 - f) * dG - kappa * phi_lap(0u, p);
    let fn2 = clamp(f - dt * L * df, 0.0, 1.0);
    phi_set(0u, i, fn2);
    var t = tt(0u, i);
    if (t < 0.0 && fn2 > 0.5) { t = f32(G.step); }
    tt_set(0u, i, t);
    // A newly transformed cell joins the grain of its most transformed neighbour.
    if (id < 0.0 && fn2 > 0.5) {
        let offs = array<vec3<i32>, 6>(vec3<i32>(1, 0, 0), vec3<i32>(-1, 0, 0), vec3<i32>(0, 1, 0), vec3<i32>(0, -1, 0), vec3<i32>(0, 0, 1), vec3<i32>(0, 0, -1));
        var best = -1.0;
        for (var o = 0; o < 6; o++) {
            let q = p + offs[o];
            let g = gid_at(0u, q);
            let pf = phi_at(0u, q);
            if (g >= 0.0 && pf > best) { best = pf; id = g; }
        }
    }
    gid_set(0u, i, id);`,
        },
    ],
    views: [
        {
            id: 'grains', label: 'Grains (Johnson–Mehl structure)', range: [0, 1], render: { mode: 'surface', iso: 0.5 },
            wgsl: `    let f = phi(0u, i);
    let g = gid(0u, i);
    let c = select(vec3<f32>(0.16, 0.18, 0.22), label_color(u32(max(g, 0.0))), g >= 0.0);
    return vec4<f32>(c * (0.45 + 0.55 * f), f);`,
        },
        {
            id: 'time', label: 'Product, coloured by transformation time', source: ['tt', 0], ignoreBelow: 0, followSteps: true, range: [0, 3000], cmap: 'turbo', render: { mode: 'surface', iso: 0.5 },
            wgsl: `    let f = phi(0u, i);
    let t = tt(0u, i);
    let c = select(vec3<f32>(0.16, 0.18, 0.22), cmap(scalar01(t)), t >= 0.0);
    return vec4<f32>(c, f);`,
        },
        {
            id: 'phi', label: 'Product fraction φ', range: [0, 1], cmap: 'viridis', render: { mode: 'iso', iso: 0.5 },
            wgsl: `    let f = phi(0u, i);
    return vec4<f32>(cmap(f), f);`,
        },
    ],
    initialState: (n) => ({ phi: new Float32Array(n ** 3), tt: new Float32Array(n ** 3).fill(-1), gid: new Float32Array(n ** 3).fill(-1) }),
};
