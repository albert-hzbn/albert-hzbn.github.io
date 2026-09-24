// FCC stacking-fault card: displacement slider with GSFE presets.
import { buildStackingFault } from '../../builders/index.js';
import { state } from '../context.js';
import { $, num, int } from '../dom.js';
import { setStructure, requireCubicSource } from '../history.js';
import { registerCard } from './card.js';

function isFcc(s) {
    const f = s.fractionalPositions().map((v) => v.map((x) => Math.round(x * 2) / 2 % 1));
    const keys = new Set(f.map((v) => v.join(',')));
    return s.count === 4 && ['0,0,0', '0,0.5,0.5', '0.5,0,0.5', '0.5,0.5,0'].every((k) => keys.has(k));
}

const LABELS = { 0: 'perfect crystal', 0.5: 'unstable stacking fault', 1: 'intrinsic stacking fault', 2: 'perfect crystal (full b)' };

function makeFault() {
    const s = requireCubicSource();
    if (!isFcc(s)) throw new Error('Needs a 4-atom FCC conventional cell as the source (e.g. the FCC preset).');
    return buildStackingFault(s, { reps: [int('sf-nx'), int('sf-ny'), int('sf-nz')], u: num('sf-u') });
}

export function initStackingFaultCard() {
    const sync = () => {
        $('sf-u-out').value = num('sf-u').toFixed(2);
        $('sf-presets').querySelectorAll('button').forEach((b) => b.setAttribute('aria-pressed', String(Math.abs(+b.dataset.u - num('sf-u')) < 1e-6)));
    };
    $('sf-u').addEventListener('input', sync);
    $('sf-presets').querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
        $('sf-u').value = b.dataset.u;
        sync();
        $('sf-u').dispatchEvent(new Event('change', { bubbles: true }));
    }));
    registerCard('sf', {
        busy: 'Building stacking fault…',
        preview: () => {
            const sf = makeFault();
            const u = num('sf-u');
            return { structure: sf, colorBy: 'tag', caption: `<b>u = ${u.toFixed(2)} b<sub>p</sub></b>${LABELS[u] ? ' · ' + LABELS[u] : ''} · ${sf.count} atoms<br><span>Upper half (colour 2) shifted along [11-2]</span>` };
        },
        build: () => {
            const sf = makeFault();
            state.view.color = 'tag'; $('v-color').value = 'tag';
            setStructure(sf, { message: `${sf.title}: ${sf.count} atoms` });
        },
    });
}
