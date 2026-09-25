// Antiphase-domain coarsening: Allen–Cahn equation for a non-conserved order
// parameter (two ordered variants), optionally with thermal noise.
import { noise } from './init.js';

export default {
    id: 'allen-cahn',
    name: 'Antiphase domains (Allen–Cahn)',
    short: 'Allen–Cahn',
    icon: '<path d="M3 3h18v18H3z"/><path d="M3 14c4-1 5-6 9-6s6 5 9 4"/>',
    description: 'After a quench into an ordered state, two equivalent variants form domains separated by antiphase boundaries. The boundaries move by mean curvature, so the domain size grows as t<sup>1/2</sup> without any long-range diffusion.',
    equations: [
        String.raw`\frac{\partial \phi}{\partial t} = -L\left(\phi^3 - \phi - h - \kappa\nabla^2\phi\right) + \xi`,
    ],
    grid: 128,
    params: [
        { key: 'L', label: 'Mobility L', value: 1, min: 0.1, max: 2, step: 0.05, group: 'Kinetics' },
        { key: 'kappa', label: 'Gradient energy κ', value: 1, min: 0.2, max: 1.5, step: 0.05, group: 'Free energy' },
        { key: 'bias', label: 'Bias field h', value: 0, min: -0.3, max: 0.3, step: 0.01, group: 'Free energy', hint: 'Favours one variant' },
        { key: 'temp', label: 'Noise strength', value: 0, min: 0, max: 0.5, step: 0.01, group: 'Kinetics' },
        { key: 'dt', label: 'Time step Δt', value: 0.1, min: 0.01, max: 0.15, step: 0.01, group: 'Numerics' },
    ],
    init: [
        { key: 'amp', label: 'Initial noise', value: 0.1, min: 0.01, max: 1, step: 0.01 },
        { key: 'mean', label: 'Initial mean', value: 0, min: -0.5, max: 0.5, step: 0.01 },
    ],
    presets: [
        { label: 'Symmetric quench', params: { bias: 0, temp: 0 }, init: { mean: 0 } },
        { label: 'With noise', params: { bias: 0, temp: 0.2 } },
        { label: 'Biased', params: { bias: 0.08, temp: 0 } },
    ],
    layout: () => ({ fields: [{ name: 'phi', comps: 1 }], aux: [] }),
    passes: [
        {
            entry: 'update',
            wgsl: `
    let f = phi(0u, i);
    let xi = temp * sqrt(dt) * (rand(i, 3u) - 0.5) * 3.4641;
    phi_set(0u, i, f - dt * L * (f * f * f - f - bias - kappa * phi_lap(0u, p)) + xi);`,
        },
    ],
    views: [
        {
            id: 'phi', label: 'Order parameter φ', range: [-1, 1], cmap: 'coolwarm', render: { mode: 'iso', iso: 0.5 },
            wgsl: `    let t = scalar01(phi(0u, i));
    return vec4<f32>(cmap(t), t);`,
        },
        {
            id: 'apb', label: 'Antiphase boundaries', autoRange: false, range: [0, 1], cmap: 'magma', render: { mode: 'volume', opacity: 0.55, gamma: 1.5 },
            wgsl: `    let f = phi(0u, i);
    let t = clamp(1.0 - f * f, 0.0, 1.0);
    return vec4<f32>(cmap(t), t);`,
        },
    ],
    initialState: (n, v, rand) => ({ phi: noise(n, v.mean, v.amp, rand) }),
};
