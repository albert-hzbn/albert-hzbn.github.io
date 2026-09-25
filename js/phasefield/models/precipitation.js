// Precipitation and coarsening: conserved composition c coupled to four
// non-conserved order parameters (antiphase variants of an ordered precipitate,
// as for γ′ in Ni-base superalloys). Variants repel each other, so touching
// precipitates of different variants do not merge.
import { noise, spheres, forSphere } from './init.js';

const V = 4;

export default {
    id: 'precipitation',
    name: 'Precipitation',
    short: 'Precipitation',
    icon: '<circle cx="7" cy="8" r="3"/><circle cx="16" cy="7" r="2"/><circle cx="15" cy="16" r="4"/><circle cx="6" cy="17" r="1.6"/>',
    description: 'Ordered precipitates nucleate in a supersaturated matrix, grow by drawing solute from it and then coarsen: large particles grow at the expense of small ones (Ostwald ripening).',
    equations: [
        String.raw`\begin{aligned} f = {}& (1-H)\,A(c-c_\alpha)^2 \\ & + H A(c-c_\beta)^2 \\ & + W\sum_p \eta_p^2(1-\eta_p)^2 \\ & + \omega W\sum_{p<q}\eta_p^2\eta_q^2 \end{aligned}`,
        String.raw`H = \min\Big(1, \sum_p h(\eta_p)\Big)`,
        String.raw`h(\eta) = \eta^3(6\eta^2 - 15\eta + 10)`,
        String.raw`\frac{\partial c}{\partial t} = \nabla\cdot M\nabla\left(\frac{\partial f}{\partial c} - \kappa_c\nabla^2 c\right)`,
        String.raw`\frac{\partial \eta_p}{\partial t} = -L\left(\frac{\partial f}{\partial \eta_p} - \kappa_\eta\nabla^2\eta_p\right)`,
    ],
    equationNote: 'Sums run over the ordered variants p.',
    grid: 96,
    params: [
        { key: 'A', label: 'Chemical energy A', value: 2, min: 0.5, max: 4, step: 0.1, group: 'Free energy' },
        { key: 'ca', label: 'Matrix composition c<sub>α</sub>', value: 0.05, min: 0, max: 0.3, step: 0.01, group: 'Free energy' },
        { key: 'cb', label: 'Precipitate composition c<sub>β</sub>', value: 0.95, min: 0.6, max: 1, step: 0.01, group: 'Free energy' },
        { key: 'Wb', label: 'Order barrier W', value: 0.5, min: 0.1, max: 2, step: 0.05, group: 'Free energy' },
        { key: 'omega', label: 'Variant repulsion ω', value: 2, min: 0, max: 5, step: 0.1, group: 'Free energy' },
        { key: 'kc', label: 'Gradient κ<sub>c</sub>', value: 0.5, min: 0.1, max: 1, step: 0.05, group: 'Free energy' },
        { key: 'ke', label: 'Gradient κ<sub>η</sub>', value: 0.5, min: 0.1, max: 1.5, step: 0.05, group: 'Free energy' },
        { key: 'M', label: 'Solute mobility M', value: 1, min: 0.1, max: 1.5, step: 0.05, group: 'Kinetics' },
        { key: 'L', label: 'Order mobility L', value: 1, min: 0.1, max: 3, step: 0.05, group: 'Kinetics' },
        { key: 'dt', label: 'Time step Δt', value: 0.01, min: 0.001, max: 0.015, step: 0.001, group: 'Numerics' },
    ],
    init: [
        { key: 'c0', label: 'Alloy composition c₀', value: 0.22, min: 0.08, max: 0.5, step: 0.01 },
        { key: 'nuclei', label: 'Nuclei', value: 60, min: 1, max: 400, step: 1 },
        { key: 'r0', label: 'Nucleus radius (cells)', value: 3.5, min: 2, max: 8, step: 0.5 },
    ],
    presets: [
        { label: 'Supersaturated', init: { c0: 0.22, nuclei: 60 } },
        { label: 'High fraction', init: { c0: 0.35, nuclei: 120 } },
        { label: 'Few, large', init: { c0: 0.2, nuclei: 12, r0: 6 } },
    ],
    layout: () => ({ fields: [{ name: 'conc', comps: 1 }, { name: 'eta', comps: V }], aux: [{ name: 'mu', comps: 1 }], common: `const V: u32 = ${V}u;\n` }),
    wgslCommon: `
fn hfun(e: f32) -> f32 { let x = clamp(e, 0.0, 1.0); return x * x * x * (6.0 * x * x - 15.0 * x + 10.0); }
fn dhfun(e: f32) -> f32 { let x = clamp(e, 0.0, 1.0); return 30.0 * x * x * (1.0 - x) * (1.0 - x); }
`,
    passes: [
        {
            entry: 'chemical_potential',
            wgsl: `
    var H = 0.0;
    for (var k = 0u; k < V; k++) { H += hfun(eta(k, i)); }
    H = min(H, 1.0);
    let cc = conc(0u, i);
    mu_set(0u, i, 2.0 * A * ((1.0 - H) * (cc - ca) + H * (cc - cb)) - kc * conc_lap(0u, p));`,
        },
        {
            entry: 'update',
            wgsl: `
    let cc = conc(0u, i);
    conc_set(0u, i, cc + dt * M * mu_lap(0u, p));
    let df = A * ((cc - cb) * (cc - cb) - (cc - ca) * (cc - ca));
    var s2 = 0.0;
    for (var k = 0u; k < V; k++) { let e = eta(k, i); s2 += e * e; }
    for (var k = 0u; k < V; k++) {
        let e = eta(k, i);
        let g = 2.0 * e * (1.0 - e) * (1.0 - 2.0 * e) + 2.0 * omega * e * (s2 - e * e);
        eta_set(k, i, e - dt * L * (dhfun(e) * df + Wb * g - ke * eta_lap(k, p)));
    }`,
        },
    ],
    views: [
        {
            id: 'variants', label: 'Precipitate variants', range: [0, 1], render: { mode: 'iso', iso: 0.5 },
            wgsl: `    var best = 0u;
    var bv = 0.0;
    for (var k = 0u; k < V; k++) { let e = eta(k, i); if (e > bv) { bv = e; best = k; } }
    let t = clamp(bv, 0.0, 1.0);
    return vec4<f32>(mix(vec3<f32>(0.62, 0.66, 0.72), label_color(best + 3u), smoothstep(0.2, 0.6, t)), t);`,
        },
        {
            id: 'c', label: 'Composition c', source: ['conc', 0], range: [0, 1], cmap: 'viridis', render: { mode: 'surface', iso: 0.5 },
            wgsl: `    let t = scalar01(conc(0u, i));
    return vec4<f32>(cmap(t), t);`,
        },
    ],
    initialState: (n, v, rand) => {
        const n3 = n ** 3;
        const conc = noise(n, v.c0, 0.01, rand);
        const eta = new Float32Array(V * n3);
        for (const c of spheres(n, v.nuclei, v.r0, rand, 1)) {
            const variant = Math.floor(rand() * V);
            forSphere(n, c, v.r0, (i) => { eta[variant * n3 + i] = 1; conc[i] = 0.95; });
        }
        return { conc, eta };
    },
};
