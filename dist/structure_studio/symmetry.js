// Space-group settings and symmetry operations.
// Data: all 530 Hall settings of the 230 space groups, extracted from spglib's
// database (BSD-3-Clause, https://github.com/spglib/spglib). Operations are kept
// in spglib's compact integer encoding and decoded here.
import { Structure, latticeFromParameters, vecmat } from './core.js';

let settings = null;

export async function loadSpaceGroups(url = new URL('./data/spacegroups.json', import.meta.url)) {
    if (settings) return settings;
    const res = await fetch(url);
    settings = await res.json();
    for (const s of settings) s.system = crystalSystem(s.n);
    return settings;
}

export function getSettings() { return settings; }

export function settingsForNumber(n) {
    return settings.filter((s) => s.n === n);
}

// Default setting per number: origin choice 2 and hexagonal axes where there is
// a choice, matching what most CIF files and databases use.
export function defaultSetting(n) {
    const list = settingsForNumber(n);
    return list.find((s) => s.choice === '2') || list.find((s) => s.choice === 'H') || list[0];
}

export function crystalSystem(n) {
    if (n <= 2) return 'triclinic';
    if (n <= 15) return 'monoclinic';
    if (n <= 74) return 'orthorhombic';
    if (n <= 142) return 'tetragonal';
    if (n <= 167) return 'trigonal';
    if (n <= 194) return 'hexagonal';
    return 'cubic';
}

// spglib encoding: rotation in base 3 (entries -1..1), translation in base 12.
export function decodeOp(code) {
    const rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    let r = code % 19683, digit = 6561;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
        rot[i][j] = Math.floor((r % (digit * 3)) / digit) - 1;
        digit /= 3;
    }
    const t = Math.floor(code / 19683);
    const trans = [0, 0, 0];
    let d = 144;
    for (let i = 0; i < 3; i++) {
        trans[i] = Math.floor((t % (d * 12)) / d) / 12;
        d /= 12;
    }
    return { rot, trans };
}

export function operations(setting) {
    if (!setting._ops) setting._ops = setting.ops.map(decodeOp);
    return setting._ops;
}

function applyOp(op, f) {
    const r = op.rot;
    return [
        r[0][0] * f[0] + r[0][1] * f[1] + r[0][2] * f[2] + op.trans[0],
        r[1][0] * f[0] + r[1][1] * f[1] + r[1][2] * f[2] + op.trans[1],
        r[2][0] * f[0] + r[2][1] * f[1] + r[2][2] * f[2] + op.trans[2],
    ];
}

const wrap01 = (v) => { let w = v - Math.floor(v); if (w > 1 - 1e-6) w = 0; return w; };

// Expand asymmetric-unit sites [{symbol, x, y, z}] with the setting's operations.
export function expandSites(setting, sites, cell, tol = 1e-3) {
    const ops = operations(setting);
    const out = [];
    const siteMultiplicities = [];
    for (const site of sites) {
        const orbit = [];
        for (const op of ops) {
            const f = applyOp(op, [site.x, site.y, site.z]).map(wrap01);
            // Compare in Cartesian space with minimum-image in fractional space.
            const dup = orbit.some((g) => {
                const d = [f[0] - g[0], f[1] - g[1], f[2] - g[2]].map((v) => v - Math.round(v));
                const c = vecmat(d, cell);
                return c[0] * c[0] + c[1] * c[1] + c[2] * c[2] < tol * tol * 100;
            });
            if (!dup) orbit.push(f);
        }
        siteMultiplicities.push(orbit.length);
        for (const f of orbit) out.push({ symbol: site.symbol, f });
    }
    return { atoms: out, multiplicities: siteMultiplicities };
}

export function buildCrystal({ setting, a, b, c, alpha, beta, gamma, sites, title }) {
    const cell = latticeFromParameters(a, b, c, alpha, beta, gamma);
    const { atoms, multiplicities } = expandSites(setting, sites, cell);
    const s = new Structure({ cell, title: title || setting.hm });
    for (const at of atoms) s.push(at.symbol, vecmat(at.f, cell));
    return { structure: s, multiplicities };
}

// Parse a CIF-style operation such as "-y+1/2, x, z+1/4".
export function parseXyzOp(str) {
    const parts = str.replace(/['"\s]/g, '').toLowerCase().split(',');
    if (parts.length !== 3) return null;
    const rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], trans = [0, 0, 0];
    parts.forEach((expr, i) => {
        const terms = expr.replace(/-/g, '+-').split('+').filter(Boolean);
        for (const term of terms) {
            const m = term.match(/^(-?)([\d.]*\/?[\d.]*)\*?([xyz]?)$/);
            if (!m) continue;
            const sign = m[1] ? -1 : 1;
            const axis = m[3];
            let coef = 1;
            if (m[2]) {
                if (m[2].includes('/')) { const [p, q] = m[2].split('/'); coef = parseFloat(p) / parseFloat(q); }
                else coef = parseFloat(m[2]);
            }
            if (axis) rot[i]['xyz'.indexOf(axis)] += sign * (m[2] ? coef : 1);
            else trans[i] += sign * coef;
        }
    });
    return { rot, trans };
}

export function expandWithOps(ops, sites, cell, tol = 1e-3) {
    const fake = { _ops: ops, ops: [] };
    return expandSites(fake, sites, cell, tol);
}

// Presets: common crystal structures with their setting, cell and Wyckoff sites.
export const PRESETS = [
    { id: 'fcc', label: 'FCC metal (Cu)', n: 225, a: 3.615, sites: [['Cu', 0, 0, 0]] },
    { id: 'bcc', label: 'BCC metal (Fe)', n: 229, a: 2.8665, sites: [['Fe', 0, 0, 0]] },
    { id: 'hcp', label: 'HCP metal (Mg)', n: 194, a: 3.209, c: 5.211, sites: [['Mg', 1 / 3, 2 / 3, 0.25]] },
    { id: 'sc', label: 'Simple cubic (Po)', n: 221, a: 3.359, sites: [['Po', 0, 0, 0]] },
    { id: 'diamond', label: 'Diamond (Si)', n: 227, a: 5.431, sites: [['Si', 0.125, 0.125, 0.125]] },
    { id: 'nacl', label: 'Rock salt (NaCl)', n: 225, a: 5.640, sites: [['Na', 0, 0, 0], ['Cl', 0.5, 0.5, 0.5]] },
    { id: 'cscl', label: 'B2 (CsCl)', n: 221, a: 4.123, sites: [['Cs', 0, 0, 0], ['Cl', 0.5, 0.5, 0.5]] },
    { id: 'zb', label: 'Zinc blende (GaAs)', n: 216, a: 5.653, sites: [['Ga', 0, 0, 0], ['As', 0.25, 0.25, 0.25]] },
    { id: 'wurtzite', label: 'Wurtzite (ZnO)', n: 186, a: 3.250, c: 5.207, sites: [['Zn', 1 / 3, 2 / 3, 0], ['O', 1 / 3, 2 / 3, 0.382]] },
    { id: 'perovskite', label: 'Perovskite (SrTiO₃)', n: 221, a: 3.905, sites: [['Sr', 0, 0, 0], ['Ti', 0.5, 0.5, 0.5], ['O', 0.5, 0.5, 0]] },
    { id: 'rutile', label: 'Rutile (TiO₂)', n: 136, a: 4.594, c: 2.959, sites: [['Ti', 0, 0, 0], ['O', 0.3049, 0.3049, 0]] },
    { id: 'fluorite', label: 'Fluorite (CaF₂)', n: 225, a: 5.463, sites: [['Ca', 0, 0, 0], ['F', 0.25, 0.25, 0.25]] },
    { id: 'l12', label: 'L1₂ (Ni₃Al)', n: 221, a: 3.572, sites: [['Al', 0, 0, 0], ['Ni', 0, 0.5, 0.5]] },
    { id: 'spinel', label: 'Spinel (MgAl₂O₄)', n: 227, a: 8.083, sites: [['Mg', 0.125, 0.125, 0.125], ['Al', 0.5, 0.5, 0.5], ['O', 0.2624, 0.2624, 0.2624]] },
    { id: 'corundum', label: 'Corundum (Al₂O₃)', n: 167, a: 4.759, c: 12.991, sites: [['Al', 0, 0, 0.35216], ['O', 0.30624, 0, 0.25]] },
    { id: 'graphite', label: 'Graphite (C)', n: 194, a: 2.464, c: 6.711, sites: [['C', 0, 0, 0.25], ['C', 1 / 3, 2 / 3, 0.25]] },
];

export function presetParameters(p) {
    const system = crystalSystem(p.n);
    const a = p.a, b = p.b ?? a, c = p.c ?? (system === 'cubic' ? a : p.c ?? a);
    const hex = system === 'hexagonal' || system === 'trigonal';
    return { a, b, c, alpha: 90, beta: 90, gamma: hex ? 120 : 90 };
}
