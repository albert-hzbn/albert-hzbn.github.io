// XYZ and extended XYZ (Lattice / Properties comment keys).
import { Structure } from '../core/index.js';
import { lines, toks } from './text.js';

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
    const named = comment.match(/comment\s*=\s*"([^"]*)"/i);
    const title = named ? named[1] : comment.replace(/Lattice\s*=\s*"[^"]*"|Properties\s*=\s*\S+|pbc\s*=\s*"[^"]*"/gi, '');
    // Older saves nested the key (comment="comment=…"); drop the repeats.
    const s = new Structure({ cell, title: title.replace(/^\s*(comment\s*=\s*)+/i, '').trim() });
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
