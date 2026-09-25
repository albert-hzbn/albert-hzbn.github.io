// Normal grain growth: Fan–Chen multi-order-parameter model. Each grain is
// represented by one of Q order parameters; grains are coloured so touching
// grains never share a parameter.
import { voronoi, colorLabels } from './init.js';

export default {
    id: 'grain-growth',
    name: 'Grain growth',
    short: 'Grain growth',
    icon: '<path d="M3 3h18v18H3z"/><path d="M3 10l6-2 3-5M9 8l3 7 9-2M12 15l-3 6M12 15l-9-1"/>',
    description: 'Curvature-driven coarsening of a polycrystal. Grain boundaries migrate towards their centres of curvature, small grains shrink and vanish, and the mean grain size grows as √t.',
    equations: [
        String.raw`\frac{\partial \eta_i}{\partial t} = -L\left(\frac{\partial f}{\partial \eta_i} - \kappa\nabla^2\eta_i\right)`,
        String.raw`\frac{\partial f}{\partial \eta_i} = -\alpha\eta_i + \beta\eta_i^3 + 2\gamma\,\eta_i\sum_{j\neq i}\eta_j^2`,
    ],
    grid: 96,
    params: [
        { key: 'L', label: 'Mobility L', value: 1, min: 0.1, max: 2, step: 0.05, group: 'Kinetics' },
        { key: 'kappa', label: 'Gradient energy κ', value: 0.5, min: 0.1, max: 1.5, step: 0.05, group: 'Free energy' },
        { key: 'gamma', label: 'Interaction γ', value: 1.5, min: 0.6, max: 3, step: 0.05, group: 'Free energy' },
        { key: 'alpha', label: 'α', value: 1, min: 0.5, max: 2, step: 0.05, group: 'Free energy' },
        { key: 'beta', label: 'β', value: 1, min: 0.5, max: 2, step: 0.05, group: 'Free energy' },
        { key: 'dt', label: 'Time step Δt', value: 0.1, min: 0.01, max: 0.25, step: 0.01, group: 'Numerics' },
    ],
    init: [
        { key: 'grains', label: 'Initial grains', value: 300, min: 10, max: 2000, step: 10 },
        { key: 'Q', label: 'Order parameters Q', value: 16, min: 4, max: 24, step: 1, hint: 'Memory grows with Q' },
    ],
    presets: [
        { label: '300 grains', init: { grains: 300 } },
        { label: '1000 fine grains', init: { grains: 1000 } },
        { label: '60 coarse grains', init: { grains: 60 } },
    ],
    layout: (v) => ({
        fields: [{ name: 'eta', comps: v.Q }],
        aux: [],
        common: `const Q: u32 = ${v.Q}u;\n`,
    }),
    passes: [
        {
            entry: 'update',
            wgsl: `
    var sum2 = 0.0;
    for (var k = 0u; k < Q; k++) { let e = eta(k, i); sum2 += e * e; }
    for (var k = 0u; k < Q; k++) {
        let e = eta(k, i);
        let dfdn = -alpha * e + beta * e * e * e + 2.0 * gamma * e * (sum2 - e * e);
        eta_set(k, i, e - dt * L * (dfdn - kappa * eta_lap(k, p)));
    }`,
        },
    ],
    views: [
        {
            id: 'grains', label: 'Grains', range: [0, 1], render: { mode: 'surface', iso: 0.6 },
            wgsl: `    var best = 0u;
    var bv = -1.0;
    var sum2 = 0.0;
    for (var k = 0u; k < Q; k++) { let e = eta(k, i); sum2 += e * e; if (e > bv) { bv = e; best = k; } }
    let edge = clamp((sum2 - 0.5) * 2.0, 0.0, 1.0);
    return vec4<f32>(label_color(best) * (0.35 + 0.65 * edge), clamp(bv, 0.0, 1.0));`,
        },
        {
            id: 'boundaries', label: 'Grain boundaries', autoRange: false, range: [0, 1], cmap: 'magma', render: { mode: 'volume', opacity: 0.6, gamma: 1.6 },
            wgsl: `    var sum2 = 0.0;
    for (var k = 0u; k < Q; k++) { let e = eta(k, i); sum2 += e * e; }
    let b = clamp(1.0 - sum2, 0.0, 1.0) * 1.6;
    return vec4<f32>(cmap(b), clamp(b, 0.0, 1.0));`,
        },
    ],
    initialState: (n, v, rand) => {
        const { labels } = voronoi(n, v.grains, rand);
        const color = colorLabels(n, labels, v.grains, v.Q);
        const n3 = n ** 3;
        const eta = new Float32Array(v.Q * n3);
        for (let i = 0; i < n3; i++) eta[color[labels[i]] * n3 + i] = 1;
        return { eta };
    },
};
