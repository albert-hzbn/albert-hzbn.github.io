// Random substitutional solid-solution card (acts on the current structure).
import { normalizeSymbol } from '../../core/index.js';
import { buildSolidSolution } from '../../builders/index.js';
import { state } from '../context.js';
import { $, int, removableRow } from '../dom.js';
import { setStructure } from '../history.js';
import { registerCard, PREVIEW_LIMIT } from './card.js';

function soluteRow(sym = 'Ni', pct = 25) {
    return removableRow(`<input value="${sym}" aria-label="Solute element"><input type="number" step="0.5" value="${pct}" aria-label="Percent of host sites">`, 'Remove solute');
}

function substitute() {
    const s = state.structure;
    if (!s.count) throw new Error('Build or open a structure first.');
    const solutes = [...$('ss-solutes').querySelectorAll('.tr:not(.head)')].map((tr) => {
        const [e, p] = tr.querySelectorAll('input');
        return { symbol: normalizeSymbol(e.value), fraction: parseFloat(p.value) / 100 };
    }).filter((x) => x.fraction > 0);
    if (solutes.reduce((a, b) => a + b.fraction, 0) > 1) throw new Error('Solute fractions add up to more than 100%.');
    return buildSolidSolution(s, $('ss-host').value, solutes, int('ss-seed'));
}

export function setFirstSolute(symbol) {
    const input = $('ss-solutes').querySelector('.tr:not(.head) input');
    if (input) input.value = symbol;
}

export function initSolidSolutionCard() {
    const table = $('ss-solutes');
    table.innerHTML = '<div class="tr head"><span>Solute</span><span>% of host sites</span><span></span></div>';
    table.appendChild(soluteRow('Ni', 25));
    $('ss-add').addEventListener('click', () => { table.appendChild(soluteRow('Al', 10)); table.dispatchEvent(new Event('change', { bubbles: true })); });
    registerCard('ss', {
        busy: 'Substituting…',
        preview: () => {
            if (state.structure.count > PREVIEW_LIMIT * 2) throw new Error(`Preview skipped for ${state.structure.count.toLocaleString()} atoms; Substitute still works.`);
            const { structure, report } = substitute();
            return { structure, caption: `<b>${structure.formula()}</b><br><span>${report}</span>` };
        },
        build: () => {
            const { structure, report } = substitute();
            setStructure(structure, { message: `Substituted ${report}`, fit: false });
        },
    });
}
