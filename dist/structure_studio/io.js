// Structure file readers and writers.
// Readers: XYZ / extended XYZ, VASP POSCAR/CONTCAR, CIF (with symmetry
// operations), LAMMPS data and dump, PDB. Writers: POSCAR, extended XYZ, CIF (P1),
// LAMMPS data.
import { Structure, vecmat, inv3, latticeFromParameters, cellParameters, dot, cross, removeDuplicates } from './core.js';
import { element, normalizeSymbol, ELEMENTS } from './elements.js';
import { parseXyzOp, expandWithOps, getSettings, defaultSetting, operations } from './symmetry.js';

const lines = (text) => text.replace(/\r/g, '').split('\n');
const toks = (line) => line.trim().split(/\s+/).filter(Boolean);

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

// ---------- XYZ / extended XYZ ----------

export function parseXYZ(text) {
    const L = lines(text);
    const n = parseInt(L[0], 10);
    if (!n) throw new Error('XYZ: first line must be the number of atoms.');
    const comment = L[1] || '';
    let cell = null;
    const lat = comment.match(/Lattice\s*=\s*"([^"]+)"/i);
    if (lat) {
        const v = lat[1].trim().split(/\s+/).map(Number);
        if (v.length === 9) cell = [v.slice(0, 3), v.slice(3, 6), v.slice(6, 9)];
    }
    // extxyz Properties: find species and pos columns
    let spCol = 0, posCol = 1;
    const props = comment.match(/Properties\s*=\s*([^\s]+)/i);
    if (props) {
        const parts = props[1].split(':');
        let col = 0;
        for (let i = 0; i + 2 < parts.length; i += 3) {
            const [nm, , cnt] = [parts[i], parts[i + 1], parseInt(parts[i + 2], 10)];
            if (/^species$/i.test(nm)) spCol = col;
            if (/^pos$/i.test(nm)) posCol = col;
            col += cnt;
        }
    }
    const s = new Structure({ cell, title: comment.replace(/Lattice\s*=\s*"[^"]*"|Properties\s*=\s*\S+|pbc\s*=\s*"[^"]*"/gi, '').trim() });
    for (let i = 0; i < n; i++) {
        const t = toks(L[i + 2] || '');
        if (t.length < 4) throw new Error(`XYZ: line ${i + 3} is incomplete.`);
        s.push(t[spCol], [+t[posCol], +t[posCol + 1], +t[posCol + 2]]);
    }
    return s;
}

export function writeXYZ(s) {
    let head = `${s.count}\n`;
    if (s.periodic) {
        head += `Lattice="${s.cell.flat().map((v) => v.toFixed(8)).join(' ')}" Properties=species:S:1:pos:R:3 pbc="T T T"`;
    } else head += 'Properties=species:S:1:pos:R:3';
    head += (s.title ? ` comment="${s.title.replace(/"/g, '')}"` : '') + '\n';
    return head + s.symbols.map((el, i) => `${el.padEnd(3)} ${s.positions[i].map((v) => v.toFixed(8).padStart(15)).join(' ')}`).join('\n') + '\n';
}

// ---------- VASP POSCAR ----------

export function parsePOSCAR(text) {
    const L = lines(text).map((l) => l.trim());
    const title = L[0];
    const sf = toks(L[1]).map(Number);
    let cell = [toks(L[2]), toks(L[3]), toks(L[4])].map((r) => r.slice(0, 3).map(Number));
    if (sf.length === 3) cell = cell.map((r) => r.map((v, k) => v * sf[k]));
    else if (sf[0] < 0) {
        const vol = Math.abs(dot(cell[0], cross(cell[1], cell[2])));
        const f = Math.cbrt(-sf[0] / vol);
        cell = cell.map((r) => r.map((v) => v * f));
    } else cell = cell.map((r) => r.map((v) => v * sf[0]));

    let li = 5;
    let species = toks(L[li]);
    let counts;
    if (species.every((t) => /^\d+$/.test(t))) {
        counts = species.map(Number);
        // VASP 4: species may be in the title line
        const guess = toks(title).filter((t) => /^[A-Z][a-z]?$/.test(t));
        species = guess.length === counts.length ? guess : counts.map((_, i) => ELEMENTS[i].symbol);
    } else {
        species = species.map((t) => t.split(/[_/]/)[0]);
        li++;
        counts = toks(L[li]).map(Number);
    }
    li++;
    if (/^s/i.test(L[li])) li++; // selective dynamics
    const cartesian = /^[ck]/i.test(L[li]);
    li++;
    const s = new Structure({ cell, title });
    let k = 0;
    species.forEach((sp, si) => {
        for (let n = 0; n < counts[si]; n++, k++) {
            const v = toks(L[li + k]).slice(0, 3).map(Number);
            const p = cartesian ? v.map((x) => x * (sf.length === 1 && sf[0] > 0 ? sf[0] : 1)) : vecmat(v, cell);
            s.push(sp, p);
        }
    });
    return s;
}

export function writePOSCAR(s, { direct = true } = {}) {
    if (!s.periodic) throw new Error('POSCAR needs a periodic cell. Add a box first (Modify → Cell).');
    const order = [...new Set(s.symbols)];
    const idx = order.flatMap((sp) => s.symbols.map((x, i) => (x === sp ? i : -1)).filter((i) => i >= 0));
    const inv = inv3(s.cell);
    const f = (v) => v.toFixed(10).padStart(16);
    let out = `${s.title || s.formula()}\n1.0\n`;
    out += s.cell.map((r) => r.map(f).join(' ')).join('\n') + '\n';
    out += order.map((x) => x.padStart(5)).join('') + '\n';
    out += order.map((sp) => String(s.symbols.filter((x) => x === sp).length).padStart(5)).join('') + '\n';
    out += direct ? 'Direct\n' : 'Cartesian\n';
    out += idx.map((i) => (direct ? vecmat(s.positions[i], inv) : s.positions[i]).map(f).join(' ')).join('\n') + '\n';
    return out;
}

// ---------- CIF ----------

function cifNumber(v) {
    if (v === undefined) return NaN;
    return parseFloat(String(v).replace(/\(.*\)/, ''));
}

// Minimal CIF tokenizer: handles loops, quoted strings and semicolon text fields.
function cifBlocks(text) {
    const L = lines(text);
    const items = {};
    const loops = [];
    let i = 0;
    const readValue = (rest) => {
        rest = rest.trim();
        if (rest) return rest.replace(/^['"]|['"]$/g, '');
        // value on following line(s)
        i++;
        if (L[i] && L[i].startsWith(';')) {
            const buf = [L[i].slice(1)];
            i++;
            while (i < L.length && !L[i].startsWith(';')) buf.push(L[i++]);
            return buf.join('\n').trim();
        }
        return (L[i] || '').trim().replace(/^['"]|['"]$/g, '');
    };
    const splitRow = (line) => {
        const out = [];
        const re = /'([^']*)'(?=\s|$)|"([^"]*)"(?=\s|$)|(\S+)/g;
        let m;
        while ((m = re.exec(line))) out.push(m[1] ?? m[2] ?? m[3]);
        return out;
    };
    while (i < L.length) {
        const line = L[i].trim();
        if (!line || line.startsWith('#')) { i++; continue; }
        if (/^loop_/i.test(line)) {
            const heads = [];
            i++;
            while (i < L.length && L[i].trim().startsWith('_')) heads.push(L[i++].trim().split(/\s+/)[0].toLowerCase());
            const vals = [];
            while (i < L.length) {
                const t = L[i].trim();
                if (!t || t.startsWith('_') || /^loop_/i.test(t) || /^data_/i.test(t)) break;
                if (t.startsWith('#')) { i++; continue; }
                if (t.startsWith(';')) {
                    const buf = [t.slice(1)]; i++;
                    while (i < L.length && !L[i].startsWith(';')) buf.push(L[i++]);
                    vals.push(buf.join('\n')); i++;
                    continue;
                }
                vals.push(...splitRow(t));
                i++;
            }
            const rows = [];
            for (let r = 0; r + heads.length <= vals.length; r += heads.length) {
                const row = {};
                heads.forEach((h, k) => { row[h] = vals[r + k]; });
                rows.push(row);
            }
            loops.push({ heads, rows });
            continue;
        }
        if (line.startsWith('_')) {
            const m = line.match(/^(\S+)\s*(.*)$/);
            items[m[1].toLowerCase()] = readValue(m[2]);
        }
        i++;
    }
    return { items, loops };
}

export function parseCIF(text) {
    const { items, loops } = cifBlocks(text);
    const a = cifNumber(items._cell_length_a), b = cifNumber(items._cell_length_b), c = cifNumber(items._cell_length_c);
    const al = cifNumber(items._cell_angle_alpha) || 90, be = cifNumber(items._cell_angle_beta) || 90, ga = cifNumber(items._cell_angle_gamma) || 90;
    if (!(a && b && c)) throw new Error('CIF: cell lengths missing.');
    const cell = latticeFromParameters(a, b, c, al, be, ga);

    const siteLoop = loops.find((l) => l.heads.includes('_atom_site_fract_x') || l.heads.includes('_atom_site_cartn_x'));
    if (!siteLoop) throw new Error('CIF: no _atom_site loop found.');
    const frac = siteLoop.heads.includes('_atom_site_fract_x');
    const sites = siteLoop.rows.map((r) => {
        const sym = r._atom_site_type_symbol || r._atom_site_label || 'X';
        let x = cifNumber(r[frac ? '_atom_site_fract_x' : '_atom_site_cartn_x']);
        let y = cifNumber(r[frac ? '_atom_site_fract_y' : '_atom_site_cartn_y']);
        let z = cifNumber(r[frac ? '_atom_site_fract_z' : '_atom_site_cartn_z']);
        if (!frac) [x, y, z] = vecmat([x, y, z], inv3(cell));
        return { symbol: normalizeSymbol(sym), x, y, z };
    });

    // Symmetry: explicit operations first, then Hall/HM/number lookup.
    const opLoop = loops.find((l) => l.heads.some((h) => /_symmetry_equiv_pos_as_xyz|_space_group_symop_operation_xyz/.test(h)));
    let ops = null;
    if (opLoop) {
        const key = opLoop.heads.find((h) => /_as_xyz|_operation_xyz/.test(h));
        ops = opLoop.rows.map((r) => parseXyzOp(r[key])).filter(Boolean);
    }
    if (!ops || !ops.length) {
        const settingsAll = getSettings();
        const hall = (items._symmetry_space_group_name_hall || items._space_group_name_hall || '').replace(/\s+/g, ' ').trim();
        const hm = (items['_symmetry_space_group_name_h-m'] || items['_space_group_name_h-m_alt'] || '').replace(/\s+/g, '').trim();
        const num = parseInt(items._symmetry_int_tables_number || items._space_group_it_number || '0', 10);
        let setting = null;
        if (settingsAll) {
            if (hall) setting = settingsAll.find((s) => s.hall.replace(/\s+/g, ' ') === hall);
            if (!setting && hm) setting = settingsAll.find((s) => s.hm.replace(/\s+/g, '') === hm || s.short === hm);
            if (!setting && num) setting = defaultSetting(num);
        }
        ops = setting ? operations(setting) : [{ rot: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], trans: [0, 0, 0] }];
    }
    const { atoms } = expandWithOps(ops, sites, cell);
    const s = new Structure({ cell, title: items._chemical_formula_sum || items._chemical_name_common || '' });
    for (const at of atoms) s.push(at.symbol, vecmat(at.f, cell));
    return removeDuplicates(s, 0.05);
}

export function writeCIF(s) {
    if (!s.periodic) throw new Error('CIF needs a periodic cell.');
    const p = cellParameters(s.cell);
    const inv = inv3(s.cell);
    const labels = {};
    let out = `data_${(s.formula() || 'structure').replace(/\s+/g, '')}\n`;
    out += `_chemical_formula_sum '${s.formula()}'\n`;
    out += `_cell_length_a ${p.a.toFixed(6)}\n_cell_length_b ${p.b.toFixed(6)}\n_cell_length_c ${p.c.toFixed(6)}\n`;
    out += `_cell_angle_alpha ${p.alpha.toFixed(4)}\n_cell_angle_beta ${p.beta.toFixed(4)}\n_cell_angle_gamma ${p.gamma.toFixed(4)}\n`;
    out += `_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number 1\n\nloop_\n_symmetry_equiv_pos_as_xyz\n  'x, y, z'\n\n`;
    out += 'loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n';
    s.symbols.forEach((el, i) => {
        labels[el] = (labels[el] || 0) + 1;
        const f = vecmat(s.positions[i], inv);
        out += `  ${el}${labels[el]} ${el} ${f.map((v) => v.toFixed(8)).join(' ')} 1.0\n`;
    });
    return out;
}

// ---------- LAMMPS ----------

// Map LAMMPS type masses back to elements (closest atomic mass).
function elementFromMass(m) {
    let best = ELEMENTS[0], d = Infinity;
    for (const e of ELEMENTS) { const dd = Math.abs(e.mass - m); if (dd < d) { d = dd; best = e; } }
    return d < 0.6 ? best.symbol : null;
}

function lammpsBoxToCell(xlo, xhi, ylo, yhi, zlo, zhi, xy = 0, xz = 0, yz = 0) {
    return { origin: [xlo, ylo, zlo], cell: [[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi - zlo]] };
}

export function parseLammpsData(text, { typeMap = null } = {}) {
    const L = lines(text);
    let xlo = 0, xhi = 0, ylo = 0, yhi = 0, zlo = 0, zhi = 0, xy = 0, xz = 0, yz = 0, nAtoms = 0;
    const masses = {};
    let style = 'atomic';
    let atomsStart = -1;
    for (let i = 0; i < L.length; i++) {
        const l = L[i].split('#')[0].trim();
        let m;
        if ((m = l.match(/^(\d+)\s+atoms$/))) nAtoms = +m[1];
        else if ((m = l.match(/^(\S+)\s+(\S+)\s+xlo\s+xhi/))) [xlo, xhi] = [+m[1], +m[2]];
        else if ((m = l.match(/^(\S+)\s+(\S+)\s+ylo\s+yhi/))) [ylo, yhi] = [+m[1], +m[2]];
        else if ((m = l.match(/^(\S+)\s+(\S+)\s+zlo\s+zhi/))) [zlo, zhi] = [+m[1], +m[2]];
        else if ((m = l.match(/^(\S+)\s+(\S+)\s+(\S+)\s+xy\s+xz\s+yz/))) [xy, xz, yz] = [+m[1], +m[2], +m[3]];
        else if (/^Masses/.test(l)) {
            let j = i + 1;
            while (j < L.length && !L[j].trim()) j++;
            while (j < L.length && L[j].trim()) {
                const t = toks(L[j].split('#')[0]);
                const comment = (L[j].split('#')[1] || '').trim();
                masses[t[0]] = comment && /^[A-Z][a-z]?$/.test(comment.split(/\s+/)[0]) ? comment.split(/\s+/)[0] : elementFromMass(+t[1]);
                j++;
            }
        } else if (/^Atoms/.test(l)) {
            const sm = L[i].match(/#\s*(\w+)/);
            if (sm) style = sm[1];
            atomsStart = i + 1;
            break;
        }
    }
    if (atomsStart < 0) throw new Error('LAMMPS data: no Atoms section.');
    const { origin, cell } = lammpsBoxToCell(xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz);
    const s = new Structure({ cell });
    let j = atomsStart;
    while (j < L.length && !L[j].trim()) j++;
    const rows = [];
    while (j < L.length && L[j].trim() && rows.length < nAtoms) rows.push(toks(L[j++].split('#')[0]));
    const typeCol = 1;
    const posCol = style === 'full' ? 4 : style === 'charge' ? 3 : style === 'molecular' || style === 'bond' ? 3 : 2;
    rows.sort((a, b) => +a[0] - +b[0]);
    for (const t of rows) {
        const type = t[style === 'full' || style === 'molecular' || style === 'bond' ? 2 : typeCol];
        const sym = (typeMap && typeMap[type]) || masses[type] || ELEMENTS[(+type - 1) % ELEMENTS.length].symbol;
        s.push(sym, [+t[posCol] - origin[0], +t[posCol + 1] - origin[1], +t[posCol + 2] - origin[2]], +type);
    }
    return s;
}

export function writeLammpsData(s, { style = 'atomic' } = {}) {
    if (!s.periodic) throw new Error('LAMMPS data needs a periodic cell.');
    // LAMMPS requires a lower-triangular cell: rotate into that frame.
    const p = cellParameters(s.cell);
    const cellL = latticeFromParameters(p.a, p.b, p.c, p.alpha, p.beta, p.gamma);
    const frac = s.fractionalPositions();
    const types = [...new Set(s.symbols)];
    const tilted = Math.abs(cellL[1][0]) > 1e-8 || Math.abs(cellL[2][0]) > 1e-8 || Math.abs(cellL[2][1]) > 1e-8;
    let out = `# ${s.title || s.formula()} (written by albertlinda.com structure builder)\n\n`;
    out += `${s.count} atoms\n${types.length} atom types\n\n`;
    out += `0.0 ${cellL[0][0].toFixed(10)} xlo xhi\n0.0 ${cellL[1][1].toFixed(10)} ylo yhi\n0.0 ${cellL[2][2].toFixed(10)} zlo zhi\n`;
    if (tilted) out += `${cellL[1][0].toFixed(10)} ${cellL[2][0].toFixed(10)} ${cellL[2][1].toFixed(10)} xy xz yz\n`;
    out += '\nMasses\n\n' + types.map((t, i) => `${i + 1} ${element(t).mass.toFixed(4)} # ${t}`).join('\n') + '\n';
    out += `\nAtoms # ${style}\n\n`;
    out += frac.map((f, i) => {
        const r = vecmat(f, cellL).map((v) => v.toFixed(10));
        const t = types.indexOf(s.symbols[i]) + 1;
        return style === 'charge' ? `${i + 1} ${t} 0.0 ${r.join(' ')}` : `${i + 1} ${t} ${r.join(' ')}`;
    }).join('\n') + '\n';
    return out;
}

export function parseLammpsDump(text, { typeMap = null, frame = -1 } = {}) {
    const L = lines(text);
    const starts = [];
    L.forEach((l, i) => { if (l.startsWith('ITEM: TIMESTEP')) starts.push(i); });
    if (!starts.length) throw new Error('LAMMPS dump: no ITEM: TIMESTEP found.');
    let i = starts[frame < 0 ? starts.length - 1 : Math.min(frame, starts.length - 1)];
    const timestep = +L[i + 1];
    const n = +L[i + 3];
    const boxHead = L[i + 4];
    const b = [L[i + 5], L[i + 6], L[i + 7]].map((l) => toks(l).map(Number));
    let xy = 0, xz = 0, yz = 0;
    let [xlo, xhi] = b[0], [ylo, yhi] = b[1], [zlo, zhi] = b[2];
    if (/xy xz yz/.test(boxHead)) {
        [xy, xz, yz] = [b[0][2], b[1][2], b[2][2]];
        xlo -= Math.min(0, xy, xz, xy + xz); xhi -= Math.max(0, xy, xz, xy + xz);
        ylo -= Math.min(0, yz); yhi -= Math.max(0, yz);
    }
    const { origin, cell } = lammpsBoxToCell(xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz);
    const cols = toks(L[i + 8]).slice(2);
    const c = (name) => cols.indexOf(name);
    let px = c('x'), py = c('y'), pz = c('z'), scaled = false;
    if (px < 0) { px = c('xu'); py = c('yu'); pz = c('zu'); }
    if (px < 0) { px = c('xs'); py = c('ys'); pz = c('zs'); scaled = true; }
    if (px < 0) throw new Error('LAMMPS dump: no x/y/z, xu/yu/zu or xs/ys/zs columns.');
    const ec = c('element'), tc = c('type'), idc = c('id');
    const rows = [];
    for (let k = 0; k < n; k++) rows.push(toks(L[i + 9 + k]));
    if (idc >= 0) rows.sort((a, b2) => +a[idc] - +b2[idc]);
    const s = new Structure({ cell, title: `LAMMPS dump, timestep ${timestep}` });
    for (const t of rows) {
        const type = tc >= 0 ? t[tc] : '1';
        const sym = ec >= 0 ? t[ec] : (typeMap && typeMap[type]) || ELEMENTS[(+type - 1) % ELEMENTS.length].symbol;
        const v = [+t[px], +t[py], +t[pz]];
        s.push(sym, scaled ? vecmat(v, cell) : [v[0] - origin[0], v[1] - origin[1], v[2] - origin[2]], +type || 0);
    }
    s.frames = starts.length;
    return s;
}

// ---------- PDB ----------

export function parsePDB(text) {
    const L = lines(text);
    let cell = null;
    const s = new Structure();
    for (const l of L) {
        if (l.startsWith('CRYST1')) {
            const v = [l.slice(6, 15), l.slice(15, 24), l.slice(24, 33), l.slice(33, 40), l.slice(40, 47), l.slice(47, 54)].map(Number);
            if (v[0] > 1.5) cell = latticeFromParameters(...v);
        } else if (l.startsWith('ATOM') || l.startsWith('HETATM')) {
            const sym = (l.slice(76, 78).trim() || l.slice(12, 16).trim());
            s.push(sym, [+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)]);
        }
    }
    s.cell = cell;
    return s;
}

export const WRITERS = {
    poscar: { label: 'VASP POSCAR', ext: 'POSCAR', fn: writePOSCAR },
    xyz: { label: 'Extended XYZ', ext: 'xyz', fn: writeXYZ },
    cif: { label: 'CIF (P1)', ext: 'cif', fn: writeCIF },
    lammps: { label: 'LAMMPS data', ext: 'data', fn: writeLammpsData },
};

