// CIF reader (loops, quoted values, text fields, symmetry operations or a
// space-group lookup) and a P1 writer.
import { Structure, vecmat, inv3, latticeFromParameters, cellParameters, removeDuplicates, normalizeSymbol } from '../core/index.js';
import { parseXyzOp, expandWithOps, getSettings, defaultSetting, operations } from '../crystal/spacegroups.js';
import { lines } from './text.js';

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
