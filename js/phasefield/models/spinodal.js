// Spinodal decomposition: Cahn–Hilliard equation for a conserved composition.
import { noise } from './init.js';

export default {
    id: 'spinodal',
    name: 'Spinodal decomposition',
    short: 'Spinodal',
    icon: '<path d="M4 12c2-6 5-6 8 0s6 6 8 0"/><path d="M4 18c2-3 5-3 8 0s6 3 8 0" opacity=".6"/><path d="M4 6c2-3 5-3 8 0s6 3 8 0" opacity=".6"/>',
    description: 'A quenched binary alloy separates into two phases. At the critical composition the phases form an interpenetrating network; off-critical compositions give droplets that coarsen by Ostwald ripening.',
    equations: [
        String.raw`\frac{\partial c}{\partial t} = \nabla\cdot M\nabla\mu`,
        String.raw`\mu = f'(c) - \kappa\nabla^2 c, \qquad f = W c^2(1-c)^2`,
    ],
    grid: 96,
    params: [
        { key: 'W', label: 'Barrier height W', value: 1, min: 0.2, max: 3, step: 0.05, group: 'Free energy' },
        { key: 'kappa', label: 'Gradient energy κ', value: 0.5, min: 0.1, max: 1.2, step: 0.05, group: 'Free energy' },
        { key: 'M', label: 'Mobility M', value: 1, min: 0.1, max: 1.5, step: 0.05, group: 'Kinetics' },
        { key: 'dt', label: 'Time step Δt', value: 0.01, min: 0.001, max: 0.02, step: 0.001, group: 'Numerics' },
    ],
    init: [
        { key: 'c0', label: 'Mean composition c₀', value: 0.5, min: 0.1, max: 0.9, step: 0.01 },
        { key: 'amp', label: 'Noise amplitude', value: 0.05, min: 0.001, max: 0.2, step: 0.001 },
    ],
    presets: [
        { label: '50 / 50 network', init: { c0: 0.5 } },
        { label: '30 % droplets', init: { c0: 0.3 } },
        { label: '15 % dilute', init: { c0: 0.15 } },
    ],
    layout: () => ({ fields: [{ name: 'conc', comps: 1 }], aux: [{ name: 'mu', comps: 1 }] }),
    passes: [
        {
            entry: 'chemical_potential',
            wgsl: `
    let cc = conc(0u, i);
    mu_set(0u, i, 2.0 * W * cc * (1.0 - cc) * (1.0 - 2.0 * cc) - kappa * conc_lap(0u, p));`,
        },
        {
            entry: 'update',
            wgsl: `
    conc_set(0u, i, conc(0u, i) + dt * M * mu_lap(0u, p));`,
        },
    ],
    views: [
        {
            id: 'c', label: 'Composition c', range: [0, 1], cmap: 'viridis', render: { mode: 'surface', iso: 0.5 },
            wgsl: `    let v = conc(0u, i);
    return vec4<f32>(cmap(scalar01(v)), scalar01(v));`,
        },
        {
            id: 'mu', label: 'Chemical potential μ', source: ['mu', 0], range: [-0.3, 0.3], cmap: 'coolwarm', render: { mode: 'surface' },
            wgsl: `    let v = mu(0u, i);
    return vec4<f32>(cmap(scalar01(v)), scalar01(v));`,
        },
    ],
    initialState: (n, v, rand) => ({ conc: noise(n, v.c0, v.amp, rand) }),
};
