// Dendritic solidification of a pure undercooled melt (Kobayashi 1993),
// extended to 3D with cubic anisotropy of the interface energy.
import { forSphere } from './init.js';

export default {
    id: 'dendrite',
    name: 'Dendritic solidification',
    short: 'Dendrite',
    icon: '<path d="M12 2v20M2 12h20"/><path d="M12 6l-3-2M12 6l3-2M12 18l-3 2M12 18l3 2M6 12l-2-3M6 12l-2 3M18 12l2-3M18 12l2 3"/>',
    description: 'A solid seed grows into an undercooled melt. Latent heat released at the front must diffuse away, which destabilises the interface; cubic anisotropy of the interface energy selects six ⟨100⟩ primary arms, and noise triggers side branches.',
    equations: [
        String.raw`\begin{aligned} \tau\frac{\partial \phi}{\partial t} = {}& \nabla\cdot\Big[\varepsilon^2\nabla\phi \\ & + \varepsilon\,|\nabla\phi|^2\frac{\partial \varepsilon}{\partial(\nabla\phi)}\Big] \\ & + \phi(1-\phi)\big(\phi - \tfrac12 + m\big) \end{aligned}`,
        String.raw`m = \frac{\alpha}{\pi}\arctan\!\big[\gamma(T_{\mathrm{eq}} - T)\big]`,
        String.raw`\frac{\partial T}{\partial t} = \nabla^2 T + K\frac{\partial \phi}{\partial t}`,
        String.raw`\varepsilon(\mathbf{n}) = \bar\varepsilon\Big[1 + \delta\Big(4\sum_i n_i^4 - 3\Big)\Big]`,
    ],
    grid: 128,
    params: [
        { key: 'K', label: 'Latent heat K', value: 1.8, min: 0.8, max: 2.4, step: 0.05, group: 'Thermal', hint: 'Higher K: slower, thinner dendrites' },
        { key: 'Teq', label: 'Melting temperature T<sub>eq</sub>', value: 1, min: 0.5, max: 1.5, step: 0.01, group: 'Thermal' },
        { key: 'delta', label: 'Anisotropy δ', value: 0.04, min: 0, max: 0.065, step: 0.0025, group: 'Interface', hint: 'Above 1/15 the interface stiffness turns negative' },
        { key: 'epsb', label: 'Gradient ε̄', value: 0.01, min: 0.005, max: 0.015, step: 0.0005, group: 'Interface' },
        { key: 'tau', label: 'Relaxation time τ', value: 0.0003, min: 0.0002, max: 0.001, step: 0.00005, group: 'Kinetics' },
        { key: 'alphaK', label: 'Coupling α', value: 0.9, min: 0.5, max: 1, step: 0.01, group: 'Kinetics' },
        { key: 'gammaK', label: 'Coupling γ', value: 10, min: 5, max: 20, step: 0.5, group: 'Kinetics' },
        { key: 'noiseA', label: 'Noise amplitude', value: 0.01, min: 0, max: 0.1, step: 0.005, group: 'Kinetics' },
        { key: 'dx', label: 'Grid spacing Δx', value: 0.03, min: 0.02, max: 0.05, step: 0.001, group: 'Numerics' },
        { key: 'dt', label: 'Time step Δt', value: 0.0001, min: 0.00002, max: 0.00015, step: 0.00001, group: 'Numerics' },
    ],
    init: [
        { key: 'r0', label: 'Seed radius (cells)', value: 4, min: 2, max: 10, step: 0.5 },
        { key: 'T0', label: 'Initial melt temperature', value: 0, min: 0, max: 0.8, step: 0.05 },
    ],
    presets: [
        { label: 'Cubic dendrite', params: { K: 1.8, delta: 0.04 } },
        { label: 'Weak anisotropy', params: { K: 1.6, delta: 0.015 } },
        { label: 'Fast growth', params: { K: 1.4, delta: 0.04 } },
    ],
    // Face fluxes J = ε²∇φ + ε|∇φ|² ∂ε/∂(∇φ) on the +x, +y, +z faces of each cell
    // (the second term carries the anisotropic stiffness that selects ⟨100⟩ arms).
    layout: () => ({ fields: [{ name: 'phi', comps: 1 }, { name: 'temp', comps: 1 }], aux: [{ name: 'flux', comps: 3 }] }),
    wgslCommon: `
fn face_flux(p: vec3<i32>, axis: i32, epsb: f32, delta: f32) -> f32 {
    let ex = vec3<i32>(1, 0, 0);
    let ey = vec3<i32>(0, 1, 0);
    let ez = vec3<i32>(0, 0, 1);
    var ea = ex; var eb = ey; var ec = ez;
    if (axis == 1) { ea = ey; eb = ez; ec = ex; }
    if (axis == 2) { ea = ez; eb = ex; ec = ey; }
    let q = p + ea;
    let s = sqrt(inv_dx2());
    // Gradient at the face: normal component by a two-point difference,
    // tangential components averaged over the two cells.
    let gn = (phi_at(0u, q) - phi_at(0u, p)) * s;
    let gb = (phi_at(0u, p + eb) - phi_at(0u, p - eb) + phi_at(0u, q + eb) - phi_at(0u, q - eb)) * 0.25 * s;
    let gc = (phi_at(0u, p + ec) - phi_at(0u, p - ec) + phi_at(0u, q + ec) - phi_at(0u, q - ec)) * 0.25 * s;
    let g2 = gn * gn + gb * gb + gc * gc;
    if (g2 < 1e-10) { return epsb * epsb * gn; }
    let gl = sqrt(g2);
    let nn = gn / gl;
    let nb = gb / gl;
    let nc = gc / gl;
    let s4 = nn * nn * nn * nn + nb * nb * nb * nb + nc * nc * nc * nc;
    let e = epsb * (1.0 + delta * (4.0 * s4 - 3.0));
    // ∂ε/∂g_n = 16 ε̄ δ (n_n³ − n_n Σn⁴) / |g|
    let de = 16.0 * epsb * delta * (nn * nn * nn - nn * s4) / gl;
    return e * e * gn + e * g2 * de;
}
`,
    passes: [
        {
            entry: 'fluxes',
            wgsl: `
    flux_set(0u, i, face_flux(p, 0, epsb, delta));
    flux_set(1u, i, face_flux(p, 1, epsb, delta));
    flux_set(2u, i, face_flux(p, 2, epsb, delta));`,
        },
        {
            entry: 'update',
            wgsl: `
    let f0 = phi(0u, i);
    let s = sqrt(inv_dx2());
    let div = (flux(0u, i) - flux_at(0u, p - vec3<i32>(1, 0, 0))
             + flux(1u, i) - flux_at(1u, p - vec3<i32>(0, 1, 0))
             + flux(2u, i) - flux_at(2u, p - vec3<i32>(0, 0, 1))) * s;
    let T = temp(0u, i);
    let m = (alphaK / 3.14159265) * atan(gammaK * (Teq - T));
    let react = f0 * (1.0 - f0) * (f0 - 0.5 + m + noiseA * (rand(i, 1u) - 0.5));
    let f1 = clamp(f0 + dt / tau * (div + react), 0.0, 1.0);
    phi_set(0u, i, f1);
    temp_set(0u, i, T + dt * temp_lap(0u, p) + K * (f1 - f0));`,
        },
    ],
    views: [
        {
            id: 'solid', label: 'Solid (coloured by interface temperature)', source: ['temp', 0], range: [0.6, 1], cmap: 'turbo', render: { mode: 'iso', iso: 0.5 },
            wgsl: `    let t = scalar01(temp(0u, i));
    return vec4<f32>(cmap(t), phi(0u, i));`,
        },
        {
            id: 'temperature', label: 'Temperature field', source: ['temp', 0], range: [0, 1], cmap: 'magma', render: { mode: 'volume', opacity: 0.35, gamma: 1.2 },
            wgsl: `    let t = scalar01(temp(0u, i));
    return vec4<f32>(cmap(t), t);`,
        },
    ],
    initialState: (n, v) => {
        const n3 = n ** 3;
        const phi = new Float32Array(n3);
        const temp = new Float32Array(n3).fill(v.T0);
        forSphere(n, [n / 2, n / 2, n / 2], v.r0, (i) => { phi[i] = 1; temp[i] = Math.max(v.T0, 0.5); });
        return { phi, temp };
    },
};
