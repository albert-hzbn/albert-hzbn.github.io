// Multi-phase-field model (Steinbach-type, double-obstacle potential) in the
// pairwise form: each i–j interface is driven by its own tension σ_ij, plus a
// higher-order triple-junction term that suppresses spurious third phases.
// (The full triple sum Σ_k(σ_jk − σ_ik)I_k is ill-conditioned for strongly
// unequal energies and was found to collapse into phase mixtures in 3D.)
// Every grain is its own phase field; grains are of two phase types, α and β,
// with separate interface energies, so dihedral angles and wetting follow from
// the energy balance at triple junctions.
import { voronoi } from './init.js';

export default {
    id: 'multiphase',
    name: 'Multi-phase field',
    short: 'Multi-phase',
    icon: '<path d="M12 3v7M12 10l-7 7M12 10l7 7"/><circle cx="12" cy="10" r="1.6" class="fill"/><path d="M3 3h18v18H3z" opacity=".5"/>',
    description: 'N phase fields φ<sub>i</sub> that sum to one, each a grain of phase α or β. Unequal interface energies set the dihedral angles at triple junctions; when σ<sub>αα</sub> > 2σ<sub>αβ</sub> the β phase wets the α grain boundaries. A driving force ΔG lets β grow into α.',
    equations: [
        String.raw`\begin{aligned} \frac{\partial \phi_i}{\partial t} = {}& \sum_{j\neq i}\frac{M}{N_a}\Big[\sigma_{ij}(I_i - I_j) \\ & + T_j - T_i \\ & + \frac{\pi}{\eta}\sqrt{\phi_i\phi_j}\,\Delta G_{ij}\Big] \end{aligned}`,
        String.raw`I_k = \nabla^2\phi_k + \frac{\pi^2}{\eta^2}\phi_k`,
        String.raw`T_i = \frac{\pi^2}{\eta^2}\,\sigma_3\!\!\sum_{j<k;\; j,k\neq i}\!\!\phi_j\phi_k`,
    ],
    equationNote: 'Pairwise double-obstacle form; N<sub>a</sub> is the number of phases present at a point.',
    grid: 64,
    params: [
        { key: 'saa', label: 'Energy σ<sub>αα</sub>', value: 1, min: 0.2, max: 3, step: 0.05, group: 'Interfaces' },
        { key: 'sab', label: 'Energy σ<sub>αβ</sub>', value: 0.6, min: 0.2, max: 3, step: 0.05, group: 'Interfaces' },
        { key: 'sbb', label: 'Energy σ<sub>ββ</sub>', value: 1, min: 0.2, max: 3, step: 0.05, group: 'Interfaces' },
        { key: 'width', label: 'Interface width η (cells)', value: 5, min: 3, max: 9, step: 0.5, group: 'Interfaces' },
        { key: 'st', label: 'Triple-junction penalty σ₃', value: 5, min: 0, max: 15, step: 0.5, group: 'Interfaces', hint: 'Higher-order term σ₃Σφᵢφⱼφₖ; suppresses spurious third phases' },
        { key: 'dG', label: 'Driving force ΔG (β over α)', value: 0, min: -0.3, max: 0.3, step: 0.01, group: 'Thermodynamics' },
        { key: 'M', label: 'Interface mobility M', value: 1, min: 0.1, max: 2, step: 0.05, group: 'Kinetics' },
        { key: 'dt', label: 'Time step Δt', value: 0.05, min: 0.005, max: 0.1, step: 0.005, group: 'Numerics' },
    ],
    init: [
        { key: 'NP', label: 'Grains (phase fields)', value: 12, min: 3, max: 16, step: 1 },
        { key: 'NB', label: 'of which β grains', value: 4, min: 0, max: 16, step: 1 },
    ],
    presets: [
        { label: 'Wetting β (σαβ = 0.4)', params: { saa: 1, sab: 0.4, sbb: 1, dG: 0 } },
        { label: 'Equal energies', params: { saa: 1, sab: 1, sbb: 1, dG: 0 } },
        { label: 'β grows (ΔG > 0)', params: { saa: 1, sab: 0.8, sbb: 1, dG: 0.15 } },
    ],
    layout: (v) => ({
        fields: [{ name: 'phi', comps: v.NP }],
        aux: [],
        common: `const NP: u32 = ${v.NP}u;\nconst NA: u32 = ${Math.max(0, v.NP - Math.min(v.NB, v.NP))}u;\n`,
    }),
    // σ(a, b) reads saa, sab, sbb, which are parameters 0, 1 and 2.
    wgslCommon: `
fn isBeta(k: u32) -> bool { return k >= NA; }
fn sigma(a: u32, b: u32) -> f32 {
    if (a == b) { return 0.0; }
    let ba = isBeta(a);
    let bb = isBeta(b);
    if (ba && bb) { return P.v[0].z; }
    if (!ba && !bb) { return P.v[0].x; }
    return P.v[0].y;
}
`,
    passes: [
        {
            entry: 'update',
            wgsl: `
    let PI = 3.14159265;
    let c2 = (PI / width) * (PI / width);
    var f: array<f32, 16>;
    var lap: array<f32, 16>;
    var act: array<bool, 16>;
    var nact = 0.0;
    let offs = array<vec3<i32>, 6>(vec3<i32>(1, 0, 0), vec3<i32>(-1, 0, 0), vec3<i32>(0, 1, 0), vec3<i32>(0, -1, 0), vec3<i32>(0, 0, 1), vec3<i32>(0, 0, -1));
    for (var k = 0u; k < NP; k++) {
        f[k] = phi(k, i);
        var a = f[k] > 1e-5;
        if (!a) { for (var o = 0; o < 6; o++) { if (phi_at(k, p + offs[o]) > 1e-5) { a = true; } } }
        act[k] = a;
        if (a) { nact += 1.0; lap[k] = phi_lap(k, p) + c2 * f[k]; } else { lap[k] = 0.0; }
    }
    // Higher-order triple-junction term: T_i = c2·σ₃·Σ_{j<k, j,k≠i} φ_j φ_k.
    var S = 0.0;
    var S2 = 0.0;
    for (var k = 0u; k < NP; k++) { S += f[k]; S2 += f[k] * f[k]; }
    let pairs = 0.5 * (S * S - S2);
    var next: array<f32, 16>;
    var total = 0.0;
    for (var ii = 0u; ii < NP; ii++) {
        var d = 0.0;
        if (act[ii] && nact > 1.5) {
            for (var jj = 0u; jj < NP; jj++) {
                if (jj == ii || !act[jj]) { continue; }
                // Pairwise form: each i–j interface is driven by its own tension σ_ij.
                var s = sigma(ii, jj) * (lap[ii] - lap[jj]);
                let Ti = c2 * st * (pairs - f[ii] * (S - f[ii]));
                let Tj = c2 * st * (pairs - f[jj] * (S - f[jj]));
                s += Tj - Ti;
                var g = 0.0;
                if (isBeta(ii) && !isBeta(jj)) { g = dG; }
                if (!isBeta(ii) && isBeta(jj)) { g = -dG; }
                d += (M / nact) * (s + (PI / width) * sqrt(max(f[ii] * f[jj], 0.0)) * g);
            }
        }
        next[ii] = clamp(f[ii] + dt * d, 0.0, 1.0);
        total += next[ii];
    }
    for (var k = 0u; k < NP; k++) { phi_set(k, i, next[k] / max(total, 1e-6)); }`,
        },
    ],
    views: [
        {
            id: 'phases', label: 'Grains and phases', range: [0, 1], render: { mode: 'surface', iso: 0.5 },
            wgsl: `    var best = 0u;
    var bv = -1.0;
    var s2 = 0.0;
    for (var k = 0u; k < NP; k++) { let v = phi(k, i); s2 += v * v; if (v > bv) { bv = v; best = k; } }
    let hue = select(0.55 + 0.1 * fract(f32(best) * 0.618), 0.03 + 0.08 * fract(f32(best) * 0.618), isBeta(best));
    let edge = clamp((s2 - 0.45) * 2.2, 0.0, 1.0);
    return vec4<f32>(hsv2rgb(hue, select(0.35, 0.7, isBeta(best)), 0.95) * (0.3 + 0.7 * edge), clamp(bv, 0.0, 1.0));`,
        },
        {
            id: 'beta', label: 'β phase fraction', autoRange: false, range: [0, 1], cmap: 'magma', render: { mode: 'iso', iso: 0.5 },
            wgsl: `    var b = 0.0;
    for (var k = NA; k < NP; k++) { b += phi(k, i); }
    return vec4<f32>(cmap(b), clamp(b, 0.0, 1.0));`,
        },
        {
            id: 'interfaces', label: 'Interfaces', autoRange: false, range: [0, 1], cmap: 'viridis', render: { mode: 'volume', opacity: 0.5, gamma: 1.4 },
            wgsl: `    var s2 = 0.0;
    for (var k = 0u; k < NP; k++) { let v = phi(k, i); s2 += v * v; }
    let t = clamp((1.0 - s2) * 2.0, 0.0, 1.0);
    return vec4<f32>(cmap(t), t);`,
        },
    ],
    initialState: (n, v, rand) => {
        const NP = v.NP;
        const { labels } = voronoi(n, NP, rand);
        const n3 = n ** 3;
        const phi = new Float32Array(NP * n3);
        for (let i = 0; i < n3; i++) phi[labels[i] * n3 + i] = 1;
        return { phi };
    },
};

