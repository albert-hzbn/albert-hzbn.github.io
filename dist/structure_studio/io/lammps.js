// LAMMPS data files (atomic, charge, full, molecular) and dump files
// (orthogonal or triclinic boxes; scaled, wrapped or unwrapped coordinates).
import { Structure, vecmat, latticeFromParameters, cellParameters, element, ELEMENTS } from '../core/index.js';
import { lines, toks } from './text.js';

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
    let out = `# ${s.title || s.formula()} (written by AtomForge Web, albertlinda.com)\n\n`;
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
