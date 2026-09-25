// Block-copolymer microphase separation: Ohta–Kawasaki (Cahn–Hilliard with a
// long-range term that limits domain size to a finite period).
import { noise } from './init.js';

export default {
    id: 'copolymer',
    name: 'Block copolymer',
    short: 'Copolymer',
    icon: '<path d="M5 4v16M10 4v16M15 4v16M20 4v16"/>',
    description: 'Diblock copolymers cannot separate on a large scale because the two blocks are bonded, so they self-assemble into periodic lamellae, cylinders or spheres depending on the block fraction.',
    equations: [
        String.raw`\frac{\partial \psi}{\partial t} = M\nabla^2\left(\psi^3 - \psi - \kappa\nabla^2\psi\right) - \alpha(\psi - \bar\psi)`,
    ],
    grid: 96,
    params: [
        { key: 'mean', label: 'Block fraction ψ̄', value: 0, min: -0.5, max: 0.5, step: 0.01, group: 'Composition', hint: 'Also sets the initial mean' },
        { key: 'alpha', label: 'Long-range strength α', value: 0.04, min: 0.005, max: 0.15, step: 0.005, group: 'Free energy' },
        { key: 'kappa', label: 'Gradient energy κ', value: 1, min: 0.3, max: 1.4, step: 0.05, group: 'Free energy' },
        { key: 'M', label: 'Mobility M', value: 1, min: 0.1, max: 1.2, step: 0.05, group: 'Kinetics' },
        { key: 'dt', label: 'Time step Δt', value: 0.008, min: 0.001, max: 0.012, step: 0.001, group: 'Numerics' },
    ],
    init: [
        { key: 'amp', label: 'Noise amplitude', value: 0.1, min: 0.01, max: 0.4, step: 0.01 },
    ],
    presets: [
        { label: 'Lamellae', params: { mean: 0 } },
        { label: 'Cylinders', params: { mean: -0.25 } },
        { label: 'Spheres', params: { mean: -0.38 } },
    ],
    layout: () => ({ fields: [{ name: 'psi', comps: 1 }], aux: [{ name: 'mu', comps: 1 }] }),
    passes: [
        {
            entry: 'chemical_potential',
            wgsl: `
    let s = psi(0u, i);
    mu_set(0u, i, s * s * s - s - kappa * psi_lap(0u, p));`,
        },
        {
            entry: 'update',
            wgsl: `
    let s = psi(0u, i);
    psi_set(0u, i, s + dt * (M * mu_lap(0u, p) - alpha * (s - mean)));`,
        },
    ],
    views: [
        {
            id: 'psi', label: 'Order parameter ψ', range: [-1, 1], cmap: 'coolwarm', render: { mode: 'iso', iso: 0.5 },
            wgsl: `    let t = scalar01(psi(0u, i));
    return vec4<f32>(cmap(t), t);`,
        },
    ],
    // The mean is a model parameter here, so the initial state uses it too.
    initialState: (n, v, rand, params) => ({ psi: noise(n, params.mean, v.amp, rand) }),
};
