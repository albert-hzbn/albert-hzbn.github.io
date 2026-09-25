// Swift–Hohenberg pattern formation: a finite-wavelength instability that
// produces lamellae, or with a quadratic term, close-packed spheres.
import { noise } from './init.js';

export default {
    id: 'swift-hohenberg',
    name: 'Swift–Hohenberg patterns',
    short: 'Patterns',
    icon: '<circle cx="6" cy="6" r="2"/><circle cx="12" cy="6" r="2"/><circle cx="18" cy="6" r="2"/><circle cx="9" cy="12" r="2"/><circle cx="15" cy="12" r="2"/><circle cx="6" cy="18" r="2"/><circle cx="12" cy="18" r="2"/><circle cx="18" cy="18" r="2"/>',
    description: 'The generic model of pattern selection (convection rolls, Turing patterns). A band of wavelengths near 2π/q₀ grows; the cubic term saturates it. With g = 0 the pattern is lamellar; a quadratic term g favours spots, which order into a BCC-like lattice in 3D.',
    equations: [
        String.raw`\frac{\partial \psi}{\partial t} = r\psi - \left(q_0^2 + \nabla^2\right)^2\psi + g\psi^2 - \psi^3`,
    ],
    grid: 96,
    params: [
        { key: 'r', label: 'Control parameter r', value: 0.25, min: 0.01, max: 0.8, step: 0.01, group: 'Instability' },
        { key: 'g', label: 'Quadratic term g', value: 0, min: 0, max: 1.5, step: 0.05, group: 'Instability' },
        { key: 'q0', label: 'Wavenumber q₀', value: 1, min: 0.7, max: 1.3, step: 0.01, group: 'Instability' },
        { key: 'dx', label: 'Grid spacing Δx', value: 0.785, min: 0.6, max: 1, step: 0.005, group: 'Numerics', hint: 'π/4: eight cells per wavelength' },
        { key: 'dt', label: 'Time step Δt', value: 0.004, min: 0.0005, max: 0.006, step: 0.0005, group: 'Numerics' },
    ],
    init: [
        { key: 'amp', label: 'Noise amplitude', value: 0.1, min: 0.01, max: 0.5, step: 0.01 },
    ],
    presets: [
        { label: 'Lamellae', params: { r: 0.25, g: 0 } },
        { label: 'Spots', params: { r: 0.15, g: 1 } },
        { label: 'Weakly unstable', params: { r: 0.05, g: 0 } },
    ],
    layout: () => ({ fields: [{ name: 'psi', comps: 1 }], aux: [{ name: 'lp', comps: 1 }] }),
    passes: [
        { entry: 'laplacian', wgsl: `\n    lp_set(0u, i, psi_lap(0u, p));` },
        {
            entry: 'update',
            wgsl: `
    let s = psi(0u, i);
    let lin = q0 * q0 * q0 * q0 * s + 2.0 * q0 * q0 * lp(0u, i) + lp_lap(0u, p);
    psi_set(0u, i, s + dt * (r * s - lin + g * s * s - s * s * s));`,
        },
    ],
    views: [
        {
            id: 'psi', label: 'Field ψ', range: [-0.8, 0.8], cmap: 'viridis', render: { mode: 'iso', iso: 0.55 },
            wgsl: `    let t = scalar01(psi(0u, i));
    return vec4<f32>(cmap(t), t);`,
        },
    ],
    initialState: (n, v, rand) => ({ psi: noise(n, 0, v.amp, rand) }),
};
