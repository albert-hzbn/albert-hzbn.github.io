// Amorphous structure card (Random Sequential Addition). The preview packs a
// proportionally smaller system at the same density.
import { normalizeSymbol, BY_SYMBOL } from '../../core/index.js';
import { buildAmorphous, boxForDensity } from '../../builders/index.js';
import { state } from '../context.js';
import { $, num, int, segmented, removableRow, fmt } from '../dom.js';
import { setStructure, setSource } from '../history.js';
import { registerCard } from './card.js';

const PREVIEW_ATOMS = 1500;

function speciesRow(sym = 'Si', count = 64) {
    return removableRow(`<input value="${sym}" aria-label="Element"><input type="number" min="0" step="1" value="${count}" aria-label="Number of atoms">`, 'Remove element');
}

function pairRow(pair = 'Si-O', d = 1.45) {
    return removableRow(`<input value="${pair}" aria-label="Element pair, e.g. Si-O"><input type="number" step="0.05" value="${d}" aria-label="Minimum distance in Å">`, 'Remove pair');
}

const boxMode = () => $('am-mode').querySelector('[aria-pressed="true"]').dataset.v;

function readParams() {
    const composition = [...$('am-comp').querySelectorAll('.tr:not(.head)')].map((tr) => {
        const [e, n] = tr.querySelectorAll('input');
        return { symbol: normalizeSymbol(e.value), count: parseInt(n.value, 10) || 0 };
    }).filter((c) => c.count > 0);
    for (const c of composition) if (!BY_SYMBOL[c.symbol]) throw new Error(`Unknown element "${c.symbol}".`);
    const pairDistances = {};
    for (const tr of $('am-pairs').querySelectorAll('.tr:not(.head)')) {
        const [p, d] = tr.querySelectorAll('input');
        const parts = p.value.split(/[-–\s]+/).filter(Boolean);
        if (parts.length === 2 && parseFloat(d.value) > 0) pairDistances[parts.join('-')] = parseFloat(d.value);
    }
    return {
        composition,
        box: boxMode() === 'manual' ? [num('am-a'), num('am-b'), num('am-c')] : null,
        density: num('am-density'),
        scaleFactor: num('am-scale') || 1,
        tolerance: num('am-tol'),
        seed: int('am-seed'),
        maxAttempts: int('am-attempts') || 1000,
        pairDistances,
    };
}

function updateBoxHint(p) {
    if (boxMode() !== 'density' || !p.composition.length) { $('am-box-hint').textContent = ''; return; }
    const L = boxForDensity(p.composition, p.density) * p.scaleFactor;
    $('am-box-hint').textContent = `Cubic box of ${fmt(L, 2)} Å for ${p.composition.reduce((n, c) => n + c.count, 0)} atoms.`;
}

export function initAmorphousCard() {
    const comp = $('am-comp'), pairs = $('am-pairs');
    comp.innerHTML = '<div class="tr head"><span>Element</span><span>Atoms</span><span></span></div>';
    comp.append(speciesRow('Si', 64), speciesRow('O', 128));
    pairs.innerHTML = '<div class="tr head"><span>Pair</span><span>Min. distance (Å)</span><span></span></div>';
    pairs.append(pairRow('Si-O', 1.45), pairRow('O-O', 2.3), pairRow('Si-Si', 2.8));
    $('am-add-el').addEventListener('click', () => { comp.appendChild(speciesRow('Na', 16)); comp.dispatchEvent(new Event('change', { bubbles: true })); });
    $('am-add-pair').addEventListener('click', () => { pairs.appendChild(pairRow('Na-O', 2.1)); pairs.dispatchEvent(new Event('change', { bubbles: true })); });
    segmented($('am-mode'), (v) => {
        $('am-density-opts').hidden = v !== 'density';
        $('am-manual-opts').hidden = v !== 'manual';
        $('am-mode').dispatchEvent(new Event('change', { bubbles: true }));
    });

    registerCard('amorphous', {
        busy: 'Packing atoms…',
        preview: () => {
            const p = readParams();
            updateBoxHint(p);
            const total = p.composition.reduce((n, c) => n + c.count, 0);
            const k = total > PREVIEW_ATOMS ? PREVIEW_ATOMS / total : 1;
            const small = {
                ...p,
                composition: p.composition.map((c) => ({ ...c, count: Math.max(1, Math.round(c.count * k)) })),
                box: p.box ? p.box.map((v) => v * Math.cbrt(k)) : null,
            };
            const r = buildAmorphous(small);
            const note = k < 1 ? `<br><span>Preview packs ${r.requested} of ${total.toLocaleString()} atoms at the same density</span>` : '';
            const warn = r.skipped ? `<br><span class="bad">${r.skipped} atom(s) did not fit: loosen the constraints</span>` : '';
            return { structure: r.structure, caption: `<b>${r.structure.formula()}</b> · ${fmt(r.density, 3)} g/cm³ · box ${r.box.map((v) => fmt(v, 1)).join(' × ')} Å${note}${warn}` };
        },
        build: () => {
            const r = buildAmorphous(readParams());
            state.view.bonds = r.placed < 60000; $('v-bonds').checked = state.view.bonds;
            setStructure(r.structure, { message: r.message });
            setSource(r.structure);
        },
    });
}
