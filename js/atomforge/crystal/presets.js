// Common crystal structures: space group, cell and Wyckoff sites.
import { crystalSystem } from './spacegroups.js';

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
