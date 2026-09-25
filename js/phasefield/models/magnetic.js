// Ferromagnetic domains: Landau–Lifshitz–Gilbert dynamics of the unit
// magnetisation m with exchange, uniaxial anisotropy, an applied field and
// optional bulk Dzyaloshinskii–Moriya interaction (helices, skyrmion tubes).
// Long-range magnetostatics are not included.

export default {
    id: 'magnetic',
    name: 'Magnetic domains',
    short: 'Magnetic',
    icon: '<path d="M6 20V9a6 6 0 0 1 12 0v11"/><path d="M6 16h3M15 16h3"/><path d="M12 3v4"/>',
    description: 'The magnetisation relaxes by precession and damping (LLG). Uniaxial anisotropy gives up/down domains separated by Bloch walls; a bulk DMI twists m into helices, and a field along z turns them into skyrmion tubes. Demagnetising fields are neglected.',
    equations: [
        String.raw`\frac{\partial \mathbf{m}}{\partial t} = -\frac{\mathbf{m}\times\mathbf{H} + \alpha\,\mathbf{m}\times(\mathbf{m}\times\mathbf{H})}{1+\alpha^2}`,
        String.raw`\begin{aligned} \mathbf{H} = {}& A\nabla^2\mathbf{m} + 2K(\mathbf{m}\cdot\hat{\mathbf{z}})\,\hat{\mathbf{z}} \\ & + \mathbf{H}_{\mathrm{ext}} - 2D\,\nabla\times\mathbf{m} \end{aligned}`,
        String.raw`|\mathbf{m}| = 1`,
    ],
    grid: 96,
    params: [
        { key: 'Aex', label: 'Exchange A', value: 1, min: 0.2, max: 2, step: 0.05, group: 'Energy' },
        { key: 'Ku', label: 'Anisotropy K (easy z)', value: 0.2, min: -0.3, max: 0.6, step: 0.01, group: 'Energy', hint: 'Negative: easy plane' },
        { key: 'Dm', label: 'DMI strength D', value: 0, min: 0, max: 0.5, step: 0.01, group: 'Energy' },
        { key: 'Hx', label: 'Field H<sub>x</sub>', value: 0, min: -0.5, max: 0.5, step: 0.005, group: 'Applied field' },
        { key: 'Hy', label: 'Field H<sub>y</sub>', value: 0, min: -0.5, max: 0.5, step: 0.005, group: 'Applied field' },
        { key: 'Hz', label: 'Field H<sub>z</sub>', value: 0, min: -0.5, max: 0.5, step: 0.005, group: 'Applied field' },
        { key: 'damping', label: 'Damping α', value: 0.5, min: 0.02, max: 1, step: 0.01, group: 'Dynamics' },
        { key: 'dt', label: 'Time step Δt', value: 0.04, min: 0.005, max: 0.08, step: 0.005, group: 'Numerics' },
    ],
    init: [
        { key: 'start', label: 'Initial state', value: 0, type: 'select', options: ['Random (quench)', 'Up with noise', 'Helix along z'] },
    ],
    presets: [
        { label: 'Up/down domains', params: { Ku: 0.2, Dm: 0, Hz: 0 }, init: { start: 0 } },
        { label: 'Helices (DMI)', params: { Ku: 0, Dm: 0.3, Hz: 0 }, init: { start: 0 } },
        { label: 'Skyrmion tubes', params: { Ku: 0.02, Dm: 0.3, Hz: 0.045 }, init: { start: 0 } },
    ],
    layout: () => ({ fields: [{ name: 'mag', comps: 3 }], aux: [] }),
    passes: [
        {
            entry: 'llg',
            wgsl: `
    let m = vec3<f32>(mag(0u, i), mag(1u, i), mag(2u, i));
    let lapm = vec3<f32>(mag_lap(0u, p), mag_lap(1u, p), mag_lap(2u, p));
    let gx = mag_grad(0u, p);
    let gy = mag_grad(1u, p);
    let gz = mag_grad(2u, p);
    let curl = vec3<f32>(gz.y - gy.z, gx.z - gz.x, gy.x - gx.y);
    let H = Aex * lapm + vec3<f32>(0.0, 0.0, 2.0 * Ku * m.z) + vec3<f32>(Hx, Hy, Hz) - 2.0 * Dm * curl;
    let mxh = cross(m, H);
    let dm = -(mxh + damping * cross(m, mxh)) / (1.0 + damping * damping);
    var mn = m + dt * dm;
    let l = length(mn);
    mn = select(vec3<f32>(0.0, 0.0, 1.0), mn / l, l > 1e-6);
    mag_set(0u, i, mn.x);
    mag_set(1u, i, mn.y);
    mag_set(2u, i, mn.z);`,
        },
    ],
    views: [
        {
            id: 'mz', label: 'm<sub>z</sub>', source: ['mag', 2], range: [-1, 1], cmap: 'coolwarm', render: { mode: 'surface', iso: 0.5 },
            wgsl: `    let t = scalar01(mag(2u, i));
    return vec4<f32>(cmap(t), t);`,
        },
        {
            id: 'direction', label: 'Direction (hue = in-plane angle)', range: [0, 1], render: { mode: 'surface' },
            wgsl: `    let m = vec3<f32>(mag(0u, i), mag(1u, i), mag(2u, i));
    let h = atan2(m.y, m.x) / 6.2831853 + 0.5;
    let rgb = hsv2rgb(h, 1.0 - abs(m.z) * 0.85, 1.0);
    let c = select(mix(rgb, vec3<f32>(0.08), -m.z * 0.8), mix(rgb, vec3<f32>(1.0), m.z * 0.8), m.z > 0.0);
    return vec4<f32>(c, 0.5 * (m.z + 1.0));`,
        },
        {
            id: 'walls', label: 'Domain walls |∇m|', autoRange: false, range: [0, 0.6], cmap: 'magma', render: { mode: 'volume', opacity: 0.5, gamma: 1.3 },
            wgsl: `    let g = length(mag_grad(0u, p)) + length(mag_grad(1u, p)) + length(mag_grad(2u, p));
    let t = scalar01(g);
    return vec4<f32>(cmap(t), t);`,
        },
    ],
    initialState: (n, v, rand) => {
        const n3 = n ** 3;
        const mag = new Float32Array(3 * n3);
        for (let z = 0; z < n; z++) for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
            const i = (z * n + y) * n + x;
            let m;
            if (v.start === 1) {
                m = [0.2 * (rand() - 0.5), 0.2 * (rand() - 0.5), 1];
            } else if (v.start === 2) {
                const a = 2 * Math.PI * z / n * 4;
                m = [Math.cos(a), Math.sin(a), 0.05 * (rand() - 0.5)];
            } else {
                // Uniform random direction.
                const u = 2 * rand() - 1, t = 2 * Math.PI * rand(), s = Math.sqrt(1 - u * u);
                m = [s * Math.cos(t), s * Math.sin(t), u];
            }
            const l = Math.hypot(...m);
            mag[i] = m[0] / l; mag[n3 + i] = m[1] / l; mag[2 * n3 + i] = m[2] / l;
        }
        return { mag };
    },
};

