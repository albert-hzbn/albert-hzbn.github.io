// Symmetric tilt grain-boundary card: pick the tilt axis, then a boundary
// from the Σ list; the preview shows the two grains in different colours.
import { listTiltBoundaries, buildTiltGB } from '../../builders/index.js';
import { state } from '../context.js';
import { $, num, int } from '../dom.js';
import { setStructure, requireCubicSource } from '../history.js';
import { registerCard, schedulePreview } from './card.js';

let choice = null;

export function setGBChoice(g) { choice = g; }

export function refreshGBList() {
    const axis = $('gb-axis').value.split(',').map(Number);
    const list = listTiltBoundaries(axis, 7, int('gb-maxsig') || 51);
    const box = $('gb-list');
    box.innerHTML = list.map((g, k) => `<button type="button" data-k="${k}" aria-pressed="${k === 0}"><span>Σ${g.sigma}</span><span>(${g.plane.join(' ')})</span><span>${g.angle.toFixed(2)}°</span></button>`).join('')
        || '<p class="note" style="padding:8px">No boundaries for this axis.</p>';
    choice = list[0] || null;
    box.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
        box.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', 'false'));
        b.setAttribute('aria-pressed', 'true');
        choice = list[+b.dataset.k];
        schedulePreview('gb');
    }));
}

function makeGB() {
    const s = requireCubicSource();
    if (!choice) throw new Error('Choose a boundary from the list.');
    const axis = $('gb-axis').value.split(',').map(Number);
    return buildTiltGB(s, axis, choice.plane, {
        reps: [int('gb-n1'), int('gb-n2'), int('gb-n3')], shift: [num('gb-sy') || 0, num('gb-sz') || 0], overlap: num('gb-ov'),
    });
}

export function initGBCard() {
    $('gb-axis').addEventListener('change', refreshGBList);
    $('gb-maxsig').addEventListener('change', refreshGBList);
    refreshGBList();
    registerCard('gb', {
        busy: 'Building bicrystal…',
        preview: () => {
            const gb = makeGB();
            return { structure: gb, colorBy: 'tag', caption: `<b>Σ${choice.sigma} (${choice.plane.join(' ')})</b> · ${choice.angle.toFixed(2)}° · ${gb.count} atoms<br><span>${gb.removed} overlapping atoms removed</span>` };
        },
        build: () => {
            const gb = makeGB();
            state.view.color = 'tag'; $('v-color').value = 'tag';
            setStructure(gb, { message: `${gb.title}: ${gb.count} atoms, ${gb.removed} overlapping atoms removed` });
        },
    });
}
