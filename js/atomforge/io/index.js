// Structure file formats: format detection, readers and writers.
import { lines } from './text.js';
import { parseXYZ, writeXYZ } from './xyz.js';
import { parsePOSCAR, writePOSCAR } from './poscar.js';
import { parseCIF, writeCIF } from './cif.js';
import { parseLammpsData, parseLammpsDump, writeLammpsData } from './lammps.js';
import { parsePDB } from './pdb.js';

export { parseXYZ, writeXYZ, parsePOSCAR, writePOSCAR, parseCIF, writeCIF, parseLammpsData, parseLammpsDump, writeLammpsData, parsePDB };

export function detectFormat(name, text) {
    const n = (name || '').toLowerCase();
    if (n.endsWith('.cif')) return 'cif';
    if (n.endsWith('.xyz') || n.endsWith('.extxyz')) return 'xyz';
    if (n.endsWith('.pdb')) return 'pdb';
    if (n.includes('poscar') || n.includes('contcar') || n.endsWith('.vasp')) return 'poscar';
    if (n.endsWith('.data') || n.endsWith('.lmp') || n.includes('data.')) return 'lammps-data';
    if (n.includes('dump') || n.endsWith('.lammpstrj')) return 'lammps-dump';
    // Content sniffing
    if (/^\s*data_/m.test(text) && /_cell_length_a/.test(text)) return 'cif';
    if (/ITEM: TIMESTEP/.test(text)) return 'lammps-dump';
    if (/\batoms\b/.test(text) && /xlo xhi/.test(text)) return 'lammps-data';
    if (/^(ATOM|HETATM|CRYST1)/m.test(text)) return 'pdb';
    const L = lines(text);
    if (/^\s*\d+\s*$/.test(L[0])) return 'xyz';
    return 'poscar';
}

export function parseStructure(name, text, options = {}) {
    const fmt = options.format || detectFormat(name, text);
    const parsers = { xyz: parseXYZ, poscar: parsePOSCAR, cif: parseCIF, 'lammps-data': parseLammpsData, 'lammps-dump': parseLammpsDump, pdb: parsePDB };
    const s = parsers[fmt](text, options);
    s.title = s.title || name || fmt;
    return { structure: s, format: fmt };
}

export const WRITERS = {
    poscar: { label: 'VASP POSCAR', ext: 'vasp', fn: writePOSCAR },
    xyz: { label: 'Extended XYZ', ext: 'xyz', fn: writeXYZ },
    cif: { label: 'CIF (P1)', ext: 'cif', fn: writeCIF },
    lammps: { label: 'LAMMPS data', ext: 'data', fn: writeLammpsData },
};
