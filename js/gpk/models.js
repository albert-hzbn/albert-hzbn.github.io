// GPK Web: process models. Each model mirrors a ThermoModel in GPK's
// include/gpk_thermo.hpp and the matching driver in demos/gpk.cpp, with the
// reference input decks in inputs/ as defaults.
//
// Parameters are shown in convenient units; `scale` converts to SI.
// `reset: true` marks parameters that only take effect after a reset.

export function rng(seed) {
    let a = (seed >>> 0) || 1;
    const next = () => {
        a = (a + 0x6d2b79f5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    next.int = (lo, hi) => lo + Math.floor(next() * (hi - lo + 1));
    next.normal = () => Math.sqrt(-2 * Math.log(Math.max(next(), 1e-12))) * Math.cos(2 * Math.PI * next());
    return next;
}

const HALF_PI = Math.PI / 2;

// sigma = sigma0 (1 + eps4 cos 4 theta_c), theta_c measured in the crystal frame.
const cubicSigma = (sigma0, eps4) => (nx, ny, g) => {
    if (!eps4) return sigma0;
    const cg = Math.cos(g), sg = Math.sin(g);
    return sigma0 * (1 + eps4 * Math.cos(4 * Math.atan2(-nx * sg + ny * cg, nx * cg + ny * sg)));
};

// Read-Shockley misorientation-dependent mobility (ThermoModel::mobility).
const readShockley = (M0) => (gp, gq) => {
    let dg = Math.abs(gp - gq);
    dg = (dg + HALF_PI) % HALF_PI;
    if (dg > Math.PI / 4) dg = HALF_PI - dg;
    const thetaM = 15 * Math.PI / 180;
    if (dg < 1e-6) return 0;
    if (dg >= thetaM) return M0;
    const f = dg / thetaM;
    return M0 * f * (1 - Math.log(f));
};

// Periodic Voronoi tessellation; labels 1..n with orientations in [0, pi/2).
function voronoi(sim, n, rand) {
    const { Nx, Ny } = sim;
    const sx = [], sy = [], sg = [], used = new Set();
    for (let k = 0; k < n; k++) {
        let x, y;
        do { x = rand.int(0, Nx - 1); y = rand.int(0, Ny - 1); } while (used.has(y * Nx + x));
        used.add(y * Nx + x);
        sx.push(x); sy.push(y); sg.push(rand() * HALF_PI);
    }
    for (let y = 0; y < Ny; y++) for (let x = 0; x < Nx; x++) {
        let best = 0, bd = Infinity;
        for (let k = 0; k < n; k++) {
            let dx = Math.abs(x - sx[k]), dy = Math.abs(y - sy[k]);
            dx = Math.min(dx, Nx - dx); dy = Math.min(dy, Ny - dy);
            const d = dx * dx + dy * dy;
            if (d < bd) { bd = d; best = k; }
        }
        const i = y * Nx + x;
        sim.phi[i] = best + 1; sim.g[i] = sg[best]; sim.alpha[i] = 1;
    }
}

const eqCore = [
    String.raw`v_n = M\,\bigl[\Delta G - \sigma(\hat n)\,\kappa\bigr]`,
    String.raw`\hat n = -\frac{\nabla I_p}{|\nabla I_p|},\;\; \kappa = \nabla\!\cdot\hat n,\;\; I_p = G_{\ell}*\mathbf{1}[\phi=p]`,
];

const fmtNum = (v) => (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(4).replace(/\.?0+e/, 'e') : String(+v.toPrecision(6)));

function deckHeader(model, s) {
    return `# GPK input deck written by GPK Web (albertlinda.com/gpk.html)
# Run with:  ./gpk ${model}.ini   (from the GPK build folder)
# The web version uses its own random generator, so the initial
# microstructure for a given seed differs from the C++ code.

[run]
model = ${model}

[grid]
Nz = 1
Ny = ${s.n}
Nx = ${s.n}
dx = ${fmtNum(s.dx)}

[time]
t_end   = ${fmtNum(s.tEnd)}
n_snaps = 10
`;
}

const deckOutput = (model) => `
[output]
format  = vtk
dir     = output/${model}_web
prefix  = snap
compact = false
`;

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

export const MODELS = [
    {
        id: 'grain_growth',
        name: 'Grain growth',
        short: 'Grain growth',
        icon: '<path d="M3 9l6-5 7 3 5 6-4 7-8 1-6-5z"/><path d="M9 4l1 7 7-4M10 11l-5 5M10 11l4 9M17 7l4 6"/>',
        description: 'Curvature-driven coarsening of a Voronoi polycrystal. Every grain is a persistent label, so each crossing must pass the digital-topology certificate; grains below the Gaussian scale are completed and disappear. In 2D the mean grain area grows linearly in time (parabolic law).',
        equations: [...eqCore, String.raw`\Delta G = 0,\quad \sigma = \sigma_0\,[1+\varepsilon_4\cos 4\theta_c]`],
        n: 128, dx: 1, tEnd: 0.25,
        params: [
            { key: 'M0', label: 'Mobility M<sub>0</sub> (10⁻⁸ m³ J⁻¹ s⁻¹)', min: 0.1, max: 5, step: 0.1, value: 1, scale: 1e-8, group: 'Boundaries' },
            { key: 'sigma0', label: 'Boundary energy σ<sub>0</sub> (J m⁻²)', min: 0.1, max: 2, step: 0.05, value: 0.5, group: 'Boundaries' },
            { key: 'eps4', label: 'Four-fold anisotropy ε<sub>4</sub>', min: 0, max: 0.2, step: 0.01, value: 0, group: 'Boundaries' },
            { key: 'readShockley', label: 'Read–Shockley mobility (low-angle boundaries slow)', type: 'check', value: false, group: 'Boundaries' },
        ],
        init: [
            { key: 'nGrains', label: 'Initial grains', min: 2, max: 2000, step: 1, value: 80 },
        ],
        presets: [
            { label: 'Isotropic', params: { eps4: 0, readShockley: false }, init: { nGrains: 80 } },
            { label: 'Anisotropic', params: { eps4: 0.1, readShockley: false }, init: { nGrains: 80 } },
            { label: 'Read–Shockley', params: { eps4: 0, readShockley: true }, init: { nGrains: 120 } },
            { label: 'Fine grains', params: { eps4: 0 }, init: { nGrains: 400 }, n: 192 },
        ],
        views: ['grains', 'orient', 'kappa', 'vel', 'alpha'],
        thermo: (P) => ({
            grainNetwork: true,
            bulkDG: () => 0,
            sigma: cubicSigma(P.sigma0, P.eps4),
            mobility: P.readShockley ? readShockley(P.M0) : () => P.M0,
        }),
        setup(sim, P, I, rand) { voronoi(sim, Math.min(I.nGrains, sim.N), rand); },
        seed(sim, x, y, P, rand) {
            sim.paintDisc(x, y, 4, sim.maxLabel() + 1, rand() * HALF_PI);
            return 'Added a grain';
        },
        stats(sim) {
            const ng = sim.nLabels();
            const R = Math.sqrt(sim.N / ng / Math.PI) * sim.dx * 1e6;
            return { grains: ng, R, R2: R * R };
        },
        chart: { y: '⟨R⟩² (µm²)', series: [{ key: 'R2', label: '⟨R⟩²' }] },
        hud: (s) => [`${s.grains} grains`, `⟨R⟩ = ${s.R.toFixed(2)} µm`],
        deck: (s, P, I) => deckHeader('grain_growth', s) + `
[thermo]
M0     = ${fmtNum(P.M0)}
sigma0 = ${fmtNum(P.sigma0)}
eps4   = ${fmtNum(P.eps4)}
D      = 0.0
${P.readShockley ? '# Read-Shockley mobility is a GPK Web option; the C++ grain_growth model uses constant M0.\n' : ''}
[grain_growth]
n_grains = ${I.nGrains}
seed     = ${s.seed}
` + deckOutput('grain_growth'),
    },

    {
        id: 'recrystallisation',
        name: 'Recrystallisation',
        short: 'Recrystallisation',
        icon: '<circle cx="8" cy="9" r="3.2"/><circle cx="16" cy="15" r="4"/><circle cx="17" cy="6" r="1.6"/><path d="M3 20l3-2M20 21l1-2M3 3l2 2" />',
        description: 'Strain-free nuclei (labels ≥ 1) grow into deformed material (label 0, grey) driven by the stored dislocation energy E. A nucleus grows only above the critical radius R* = σ/E. After impingement the new grains coarsen by curvature alone. Click the map to add a nucleus.',
        equations: [...eqCore, String.raw`\Delta G_{0\to q} = -E_{\mathrm{stored}},\quad R^* = \sigma_0/E_{\mathrm{stored}}`],
        n: 128, dx: 1, tEnd: 200,
        params: [
            { key: 'M0', label: 'Mobility M<sub>0</sub> (10⁻¹² m³ J⁻¹ s⁻¹)', min: 0.1, max: 5, step: 0.1, value: 1, scale: 1e-12, group: 'Kinetics' },
            { key: 'sigma0', label: 'Boundary energy σ<sub>0</sub> (J m⁻²)', min: 0.05, max: 1, step: 0.05, value: 0.3, group: 'Kinetics' },
            { key: 'E', label: 'Stored energy E (MJ m⁻³)', min: 0.1, max: 5, step: 0.1, value: 1, scale: 1e6, group: 'Kinetics' },
            { key: 'Jdot', label: 'Continuous nucleation rate (10¹² m⁻³ s⁻¹)', min: 0, max: 100, step: 1, value: 0, scale: 1e12, group: 'Nucleation' },
        ],
        init: [
            { key: 'nNuclei', label: 'Initial nuclei', min: 0, max: 500, step: 1, value: 40 },
            { key: 'rNuc', label: 'Nucleus radius (cells)', min: 1, max: 10, step: 1, value: 3 },
        ],
        presets: [
            { label: 'Site saturated', params: { Jdot: 0, E: 1 }, init: { nNuclei: 40 } },
            { label: 'Continuous', params: { Jdot: 20, E: 1 }, init: { nNuclei: 5 } },
            { label: 'High energy', params: { Jdot: 0, E: 3 }, init: { nNuclei: 40 } },
        ],
        views: ['grains', 'orient', 'kappa', 'vel', 'alpha'],
        thermo: (P) => ({
            bulkDG: (p, q) => (p === 0 && q >= 1 ? -P.E : p >= 1 && q === 0 ? P.E : 0),
            sigma: () => P.sigma0,
            mobility: () => P.M0,
        }),
        setup(sim, P, I, rand) {
            for (let i = 0; i < sim.N; i++) sim.g[i] = rand() * HALF_PI;
            const r = I.rNuc;
            for (let n = 0; n < I.nNuclei; n++) {
                const cx = rand.int(r, sim.Nx - r - 1), cy = rand.int(r, sim.Ny - r - 1);
                const gid = n + 1, gv = rand() * HALF_PI;
                for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) {
                    if (dx * dx + dy * dy > r * r) continue;
                    const i = (cy + dy) * sim.Nx + cx + dx;
                    if (sim.phi[i] === 0) { sim.phi[i] = gid; sim.g[i] = gv; }
                }
            }
            sim.nextLabel = I.nNuclei;
        },
        // Poisson nucleation in still-deformed cells (rex_stochastic_nucleation).
        afterStep(sim, P, I, rand, dt) {
            if (!(P.Jdot > 0)) return;
            const prob = 1 - Math.exp(-P.Jdot * sim.dx ** 3 * dt);
            if (prob <= 0) return;
            const r = I.rNuc;
            // Expected count is tiny per step; sample the number of trials directly.
            const expected = prob * sim.N;
            let k = 0;
            for (let L = Math.exp(-expected), p = rand(); p > L; p *= rand()) k++;
            for (let t = 0; t < k; t++) {
                const cx = rand.int(r, sim.Nx - r - 1), cy = rand.int(r, sim.Ny - r - 1);
                if (sim.phi[cy * sim.Nx + cx] !== 0) continue;
                const gid = (sim.nextLabel = Math.max(sim.nextLabel || 0, sim.maxLabel()) + 1), gv = rand() * HALF_PI;
                for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) {
                    if (dx * dx + dy * dy > r * r) continue;
                    const i = (cy + dy) * sim.Nx + cx + dx;
                    if (sim.phi[i] === 0) { sim.phi[i] = gid; sim.g[i] = gv; sim.alpha[i] = 1; }
                }
            }
        },
        seed(sim, x, y, P, rand, I) {
            const gid = (sim.nextLabel = Math.max(sim.nextLabel || 0, sim.maxLabel()) + 1);
            sim.paintDisc(x, y, I.rNuc, gid, rand() * HALF_PI);
            return 'Added a nucleus';
        },
        stats(sim) {
            let rex = 0;
            for (let i = 0; i < sim.N; i++) if (sim.phi[i] !== 0) rex++;
            const X = rex / sim.N;
            const counts = sim.countLabels();
            let ng = 0;
            for (let l = 1; l < counts.length; l++) if (counts[l] > 0) ng++;
            return { X, grains: ng };
        },
        chart: { y: 'Recrystallised fraction X', series: [{ key: 'X', label: 'X' }], yMax: 1 },
        hud: (s) => [`X = ${(100 * s.X).toFixed(1)} %`, `${s.grains} new grains`],
        stop: (s) => s.X > 0.999,
        deck: (s, P, I) => deckHeader('recrystallisation', s) + `
[thermo]
M0       = ${fmtNum(P.M0)}
sigma0   = ${fmtNum(P.sigma0)}
E_stored = ${fmtNum(P.E)}

[nucleation]
n_nuclei  = ${I.nNuclei}
nucleus_r = ${I.rNuc}
Jdot      = ${fmtNum(P.Jdot)}
seed      = ${s.seed}
` + deckOutput('recrystallisation'),
    },

    {
        id: 'solidification',
        name: 'Dendritic solidification',
        short: 'Solidification',
        icon: '<path d="M12 3v18M3 12h18"/><path d="M12 7l-2.5-2M12 7l2.5-2M12 17l-2.5 2M12 17l2.5 2M7 12l-2-2.5M7 12l-2 2.5M17 12l2-2.5M17 12l2 2.5"/>',
        description: 'A solid seed (label 1) grows into an undercooled binary melt (label 0), driven by the undercooling below the local liquidus, ΔG = L ΔT / T<sub>m</sub>. The growing solid keeps k c<sub>L</sub> and rejects the rest into the liquid. A bump on the front reaches fresher, solute-poor melt and grows faster (the Mullins–Sekerka instability), so the front breaks into arms once the solute layer D/v spans several cells. Kinetic anisotropy ε<sub>k</sub> sets the ⟨10⟩ arm directions of each crystal, so a rotated seed grows rotated arms, although the square grid pulls off-axis arms part of the way back toward x and y. Click to add a seed with a random orientation.',
        equations: [...eqCore, String.raw`\Delta G_{L\to S} = -\frac{L}{T_m}\,(T_m + m_L c - T)`, String.raw`M(\hat n) = M_0\,[1+\varepsilon_k\cos 4\theta_c]`],
        n: 256, dx: 0.5, tEnd: 0.03,
        params: [
            { key: 'dT', label: 'Undercooling ΔT (K)', min: 1, max: 20, step: 0.5, value: 4, group: 'Thermodynamics' },
            { key: 'epsKin', label: 'Kinetic anisotropy ε<sub>k</sub>', min: 0, max: 0.8, step: 0.05, value: 0.5, group: 'Interface' },
            { key: 'eps4', label: 'Energy anisotropy ε<sub>4</sub>', min: 0, max: 0.15, step: 0.01, value: 0.04, group: 'Interface' },
            { key: 'sigma0', label: 'Interface energy σ<sub>0</sub> (J m⁻²)', min: 0.05, max: 0.6, step: 0.01, value: 0.2, group: 'Interface' },
            { key: 'M0', label: 'Mobility M<sub>0</sub> (10⁻⁸ m³ J⁻¹ s⁻¹)', min: 0.2, max: 10, step: 0.1, value: 1, scale: 1e-8, group: 'Interface' },
            { key: 'D', label: 'Solute diffusivity D (10⁻⁹ m² s⁻¹)', min: 0.1, max: 10, step: 0.1, value: 3, scale: 1e-9, group: 'Solute' },
            { key: 'k', label: 'Partition coefficient k', min: 0.05, max: 1, step: 0.05, value: 0.3, group: 'Solute' },
            { key: 'reject', label: 'Solute rejection at the growing solid', type: 'check', value: true, group: 'Solute' },
        ],
        init: [
            { key: 'c0', label: 'Alloy composition c₀', min: 0, max: 0.1, step: 0.005, value: 0.02 },
            { key: 'seedR', label: 'Seed radius (cells)', min: 2, max: 12, step: 1, value: 4 },
            { key: 'theta', label: 'Seed orientation (°)', min: 0, max: 45, step: 1, value: 0 },
            { key: 'nSeeds', label: 'Seeds', min: 1, max: 20, step: 1, value: 1 },
        ],
        fixed: { mLiq: -500, Tm: 933, L: 1e8 },
        presets: [
            { label: 'Dendrite', params: { dT: 4, epsKin: 0.5, eps4: 0.04, M0: 1, reject: true }, init: { c0: 0.02, theta: 0, nSeeds: 1 } },
            { label: 'Rotated 20°', params: { dT: 4, epsKin: 0.5, eps4: 0.04, M0: 1, reject: true }, init: { c0: 0.02, theta: 20, nSeeds: 1 } },
            { label: 'Equiaxed', params: { dT: 4, epsKin: 0.5, eps4: 0.04, M0: 1, reject: true }, init: { c0: 0.02, theta: 0, nSeeds: 6 } },
            { label: 'Isotropic', params: { dT: 4, epsKin: 0, eps4: 0, M0: 1, reject: true }, init: { c0: 0.02, theta: 0, nSeeds: 1 } },
            { label: 'Reference deck', params: { dT: 5, epsKin: 0, eps4: 0.04, M0: 5, reject: false }, init: { c0: 0.02, theta: 0, nSeeds: 1 }, n: 128, tEnd: 2e-3 },
        ],
        views: ['phase', 'conc', 'orient', 'kappa', 'vel', 'alpha'],
        thermo(P) {
            const { mLiq, Tm, L } = this.fixed;
            // Liquid cells carry g = 0, so the crystal frame is the solid side's.
            const frame = (gp, gq) => (gp !== 0 ? gp : gq);
            const theta = (nx, ny, g) => Math.atan2(-nx * Math.sin(g) + ny * Math.cos(g), nx * Math.cos(g) + ny * Math.sin(g));
            return {
                soluteRejection: !!P.reject, kPartition: P.k,
                bulkDG(p, q, c, T) {
                    const dG = (L / Tm) * (Tm + mLiq * c - T);
                    if (p === 0 && q === 1) return -dG;
                    if (p === 1 && q === 0) return dG;
                    return 0;
                },
                sigma: (nx, ny, gp, p, q, c, gq) => (P.eps4 ? P.sigma0 * (1 + P.eps4 * Math.cos(4 * theta(nx, ny, frame(gp, gq)))) : P.sigma0),
                // Kinetic anisotropy M = M0 (1 + eps_k cos 4 theta_c): a GPK Web option.
                mobility: (gp, gq, T, nx, ny) => {
                    if (!P.epsKin || (!nx && !ny)) return P.M0;
                    return P.M0 * (1 + P.epsKin * Math.cos(4 * theta(nx, ny, frame(gp, gq))));
                },
            };
        },
        diffusivity: (P) => P.D,
        temperature(P, I) { return this.fixed.Tm + this.fixed.mLiq * I.c0 - P.dT; },
        setup(sim, P, I, rand) {
            sim.c.fill(I.c0);
            sim.T.fill(this.temperature(P, I));
            const r = I.seedR;
            // g = 0 marks liquid, so a 0-degree crystal is stored as a tiny positive angle.
            const orient = (t) => Math.max(t, 1e-4);
            sim.paintDisc(sim.Nx >> 1, sim.Ny >> 1, r, 1, orient(I.theta * Math.PI / 180), I.c0 * P.k);
            for (let n = 1; n < I.nSeeds; n++) {
                const x = rand.int(0, sim.Nx - 1), y = rand.int(0, sim.Ny - 1);
                sim.paintDisc(x, y, r, 1, orient(rand() * HALF_PI), I.c0 * P.k);
            }
        },
        live(sim, P, I) { sim.T.fill(this.temperature(P, I)); },
        seed(sim, x, y, P, rand, I) {
            sim.paintDisc(x, y, Math.max(2, I.seedR), 1, Math.max(rand() * HALF_PI, 1e-4), I.c0 * P.k);
            sim.clickSeeded = true;
            return 'Added a seed';
        },
        stats(sim) {
            let cmax = 0;
            for (let i = 0; i < sim.N; i++) if (sim.c[i] > cmax) cmax = sim.c[i];
            const { Nx, Ny, phi } = sim;
            let edge = 0;
            for (let k = 0; k < Nx && !edge; k++) edge = phi[k] === 1 || phi[(Ny - 1) * Nx + k] === 1 ? 1 : 0;
            for (let k = 0; k < Ny && !edge; k++) edge = phi[k * Nx] === 1 || phi[k * Nx + Nx - 1] === 1 ? 1 : 0;
            return { fs: sim.phaseVolume(1) / sim.N, cmax, edge };
        },
        chart: { y: 'Solid fraction', series: [{ key: 'fs', label: 'f_s' }] },
        hud: (s) => [`f_s = ${(100 * s.fs).toFixed(1)} %`, `c_max = ${s.cmax.toFixed(4)}`],
        // A single dendrite pauses once, when an arm reaches the edge of the periodic box.
        stop(s, I, sim) {
            if (s.fs > 0.98) return true;
            if (!s.edge || I.nSeeds > 1 || sim.clickSeeded || sim.edgeStopped) return false;
            sim.edgeStopped = true;
            return 'The dendrite reached the edge of the periodic box. Press Run to continue.';
        },
        deck(s, P, I) {
            const { mLiq, Tm, L } = this.fixed;
            const notes = [];
            if (P.epsKin) notes.push(`# Kinetic anisotropy eps_k = ${fmtNum(P.epsKin)} is a GPK Web option; the C++ model has\n# isotropic mobility, so this deck grows without it.`);
            if (I.theta || I.nSeeds > 1) notes.push('# Seed orientation and extra seeds are GPK Web options; the C++ deck has one 0-degree seed.');
            return deckHeader('solidification', s) + `
[thermo]
M0     = ${fmtNum(P.M0)}
sigma0 = ${fmtNum(P.sigma0)}
eps4   = ${fmtNum(P.eps4)}
D      = ${fmtNum(P.D)}
m_liq  = ${mLiq}
k_part = ${fmtNum(P.k)}
T_melt = ${Tm}
L      = ${fmtNum(L)}
solute_rejection = ${P.reject ? 'true' : 'false'}
${notes.length ? notes.join('\n') + '\n' : ''}
[initial]
c0           = ${fmtNum(I.c0)}
undercooling = ${fmtNum(P.dT)}
seed_radius  = ${I.seedR}
` + deckOutput('solidification');
        },
    },

    {
        id: 'precipitation',
        name: 'Multiphase precipitation',
        short: 'Precipitation',
        icon: '<circle cx="7" cy="8" r="2.6"/><circle cx="15" cy="7" r="1.8"/><circle cx="16" cy="15.5" r="3.2"/><circle cx="7.5" cy="16.5" r="1.6"/>',
        description: 'Seeds of two product phases, β and γ, grow in a matrix α. Each phase has a polynomial Gibbs energy G<sub>q</sub>(c, T); at every interface cell the product is chosen by local Gibbs minimisation among the neighbouring phases. Both phases quickly consume the matrix; afterwards the more stable one grows into the other, at a rate set by their energy difference. Near 590 K the two are almost equal and only curvature acts. Click to add a seed.',
        equations: [...eqCore, String.raw`G_q = a_q + b_q c + c_q c^2 + d_q (T - T_{\mathrm{ref}})`, String.raw`q^\star = \arg\min_q\,[\,G_q - G_p\,]`],
        n: 128, dx: 0.5, tEnd: 5e-3,
        params: [
            { key: 'T', label: 'Temperature T (K)', min: 500, max: 1100, step: 5, value: 600, group: 'Thermodynamics' },
            { key: 'sigma0', label: 'Interface energy σ<sub>0</sub> (J m⁻²)', min: 0.02, max: 1, step: 0.01, value: 0.18, group: 'Interface' },
            { key: 'M0', label: 'Mobility M<sub>0</sub> (10⁻⁸ m³ J⁻¹ s⁻¹)', min: 0.1, max: 10, step: 0.1, value: 2, scale: 1e-8, group: 'Interface' },
        ],
        init: [
            { key: 'c0', label: 'Composition c₀', min: 0, max: 0.3, step: 0.005, value: 0.06 },
            { key: 'nSeeds', label: 'Seeds', min: 1, max: 300, step: 1, value: 60 },
            { key: 'seedR', label: 'Seed radius (cells)', min: 1, max: 8, step: 1, value: 2 },
        ],
        fixed: { Tref: 900, a: [0, -8e6, -7.5e6], b: [0, 2e7, -1e7], c2: [0, 2.5e8, 1.8e8], dT: [0, 1e4, 5e3] },
        presets: [
            { label: 'γ slowly wins', params: { T: 600 }, init: { c0: 0.06, nSeeds: 60 } },
            { label: 'β ≈ γ', params: { T: 590 }, init: { c0: 0.06, nSeeds: 60 } },
            { label: 'β favoured', params: { T: 880 }, init: { c0: 0.01, nSeeds: 40 } },
            { label: 'Reference deck', params: { T: 880 }, init: { c0: 0.06, nSeeds: 24 } },
        ],
        views: ['phase', 'kappa', 'vel', 'alpha'],
        diffusivity: () => 2e-9,
        G(q, c, T) {
            const f = this.fixed;
            return f.a[q] + f.b[q] * c + f.c2[q] * c * c + f.dT[q] * (T - f.Tref);
        },
        thermo(P) {
            const G = (q, c, T) => this.G(q, c, T);
            return {
                phaseSelection: true,
                bulkDG: (p, q, c, T) => (p < 0 || p > 2 || q < 0 || q > 2 ? 0 : G(q, c, T) - G(p, c, T)),
                selectProduct(p, c, T, cand, n) {
                    let best = 1e30, bq = -1;
                    for (let k = 0; k < n; k++) {
                        const q = cand[k];
                        if (q < 0 || q > 2 || q === p) continue;
                        const d = G(q, c, T) - G(p, c, T);
                        if (d < best) { best = d; bq = q; }
                    }
                    return bq;
                },
                sigma: () => P.sigma0,
                mobility: () => P.M0,
            };
        },
        setup(sim, P, I, rand) {
            sim.c.fill(I.c0);
            sim.T.fill(P.T);
            const r = I.seedR, { Nx, Ny } = sim;
            for (let n = 0; n < I.nSeeds; n++) {
                const cx = rand.int(r, Math.max(r, Nx - r - 1)), cy = rand.int(r, Math.max(r, Ny - r - 1));
                const p = rand.int(1, 2);
                let overlap = false;
                for (let dy = -r - 1; dy <= r + 1 && !overlap; dy++) for (let dx = -r - 1; dx <= r + 1 && !overlap; dx++) {
                    const x = cx + dx, y = cy + dy;
                    if (x >= 0 && x < Nx && y >= 0 && y < Ny && sim.phi[y * Nx + x] !== 0) overlap = true;
                }
                if (overlap) continue;
                for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) {
                    if (dx * dx + dy * dy > r * r) continue;
                    const x = cx + dx, y = cy + dy;
                    if (x >= 0 && x < Nx && y >= 0 && y < Ny) sim.phi[y * Nx + x] = p;
                }
            }
        },
        live(sim, P) { sim.T.fill(P.T); },
        seed(sim, x, y, P, rand, I) {
            const p = rand.int(1, 2);
            sim.paintDisc(x, y, Math.max(2, I.seedR), p);
            return `Added a ${p === 1 ? 'β' : 'γ'} seed`;
        },
        stats(sim) {
            const counts = sim.countLabels();
            return { f0: counts[0] / sim.N, f1: (counts[1] || 0) / sim.N, f2: (counts[2] || 0) / sim.N };
        },
        chart: { y: 'Phase fraction', series: [{ key: 'f0', label: 'α matrix', color: '#8a8f98' }, { key: 'f1', label: 'β', color: '#3373d1' }, { key: 'f2', label: 'γ', color: '#e0662f' }], yMax: 1 },
        hud: (s) => [`α ${(100 * s.f0).toFixed(1)} %`, `β ${(100 * s.f1).toFixed(1)} % · γ ${(100 * s.f2).toFixed(1)} %`],
        deck(s, P, I) {
            const f = this.fixed;
            return deckHeader('precipitation', s) + `
[thermo]
M0       = ${fmtNum(P.M0)}
sigma0   = ${fmtNum(P.sigma0)}
D        = 2.0e-9
n_phases = 3
T_ref    = ${f.Tref}
a        = ${f.a.map(fmtNum).join(', ')}
b        = ${f.b.map(fmtNum).join(', ')}
c2       = ${f.c2.map(fmtNum).join(', ')}
dT       = ${f.dT.map(fmtNum).join(', ')}

[initial]
c0          = ${fmtNum(I.c0)}
T0          = ${fmtNum(P.T)}
n_seeds     = ${I.nSeeds}
seed_radius = ${I.seedR}
seed        = ${s.seed}
` + deckOutput('precipitation');
        },
    },

    {
        id: 'wulff_shape',
        name: 'Wulff shape (tensor anisotropy)',
        short: 'Wulff shape',
        icon: '<ellipse cx="12" cy="12" rx="8" ry="5"/><path d="M12 3v3M12 18v3M2 12h2M20 12h2"/>',
        description: 'A single circular grain shrinks under capillarity with a gradient-energy tensor, σ(n̂) = σ<sub>0</sub>√(b<sub>xx</sub>n<sub>x</sub>² + b<sub>yy</sub>n<sub>y</sub>²). Boundaries with the higher energy retreat faster, so for b<sub>yy</sub> &gt; b<sub>xx</sub> the grain flattens and becomes wider than it is tall. This is the GPK wulff_shape verification case.',
        equations: [...eqCore, String.raw`\Delta G = 0,\quad \sigma(\hat n) = \sigma_0\sqrt{b_{xx} n_x^2 + b_{yy} n_y^2}`],
        n: 128, dx: 1, tEnd: 0.06,
        params: [
            { key: 'bxx', label: 'b<sub>xx</sub>', min: 0.25, max: 6, step: 0.25, value: 1, group: 'Anisotropy tensor' },
            { key: 'byy', label: 'b<sub>yy</sub>', min: 0.25, max: 6, step: 0.25, value: 4, group: 'Anisotropy tensor' },
            { key: 'sigma0', label: 'σ<sub>0</sub> (J m⁻²)', min: 0.1, max: 2, step: 0.05, value: 0.5, group: 'Kinetics' },
            { key: 'M0', label: 'Mobility M<sub>0</sub> (10⁻⁸ m³ J⁻¹ s⁻¹)', min: 0.1, max: 5, step: 0.1, value: 1, scale: 1e-8, group: 'Kinetics' },
        ],
        init: [
            { key: 'seedR', label: 'Grain radius (cells)', min: 5, max: 60, step: 1, value: 22 },
        ],
        presets: [
            { label: 'b_yy = 4', params: { bxx: 1, byy: 4 } },
            { label: 'Isotropic', params: { bxx: 1, byy: 1 } },
            { label: 'b_xx = 3', params: { bxx: 3, byy: 1 } },
        ],
        views: ['phase', 'kappa', 'vel', 'alpha'],
        thermo: (P) => ({
            bulkDG: () => 0,
            sigma: (nx, ny) => P.sigma0 * Math.sqrt(Math.max(P.bxx * nx * nx + P.byy * ny * ny, 1e-12)),
            mobility: () => P.M0,
        }),
        setup(sim, P, I) { sim.paintDisc(sim.Nx >> 1, sim.Ny >> 1, I.seedR, 1, 0); },
        stats(sim) {
            let x0 = sim.Nx, x1 = -1, y0 = sim.Ny, y1 = -1;
            for (let y = 0; y < sim.Ny; y++) for (let x = 0; x < sim.Nx; x++) {
                if (sim.phi[y * sim.Nx + x] !== 1) continue;
                if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
            }
            const um = sim.dx * 1e6;
            const hx = x1 >= x0 ? 0.5 * (x1 - x0 + 1) * um : 0, hy = y1 >= y0 ? 0.5 * (y1 - y0 + 1) * um : 0;
            return { hx, hy, ratio: hy > 0 ? hx / hy : 0, f: sim.phaseVolume(1) / sim.N };
        },
        chart: { y: 'Half-width (µm)', series: [{ key: 'hx', label: 'h_x' }, { key: 'hy', label: 'h_y' }] },
        hud: (s) => [`h_x / h_y = ${s.ratio.toFixed(3)}`, `h_x ${s.hx.toFixed(1)} µm · h_y ${s.hy.toFixed(1)} µm`],
        stop: (s) => s.f < 0.002,
        deck: (s, P, I) => deckHeader('wulff_shape', s) + `
[thermo]
M0     = ${fmtNum(P.M0)}
sigma0 = ${fmtNum(P.sigma0)}
bxx    = ${fmtNum(P.bxx)}
byy    = ${fmtNum(P.byy)}
bzz    = 1.0

[wulff]
seed_radius = ${I.seedR}
` + deckOutput('wulff_shape'),
    },
];

export const defaults = (list) => Object.fromEntries(list.map((p) => [p.key, p.value]));

// Convert UI values to SI.
export function resolveParams(model, ui) {
    const out = {};
    for (const p of model.params) out[p.key] = p.type === 'check' ? !!ui[p.key] : ui[p.key] * (p.scale || 1);
    return out;
}
