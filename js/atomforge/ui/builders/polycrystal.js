// Voronoi polycrystal card. The preview uses a smaller box with the same
// number of grains so it stays fast; grains are coloured individually.
import { cellVolume } from '../../core/index.js';
import { buildPolycrystal } from '../../builders/index.js';
import { state } from '../context.js';
import { $, num, int } from '../dom.js';
import { setStructure, requireSource } from '../history.js';
import { registerCard, previewScale } from './card.js';

const params = () => {
    const box = [num('pc-x'), num('pc-y'), num('pc-z')];
    if (box.some((v) => !(v > 2))) throw new Error('Box edges must be larger than 2 Å.');
    const grains = int('pc-n');
    if (!(grains >= 1)) throw new Error('Use at least one grain.');
    return { box, grains, seed: int('pc-seed'), overlap: num('pc-ov') };
};

export function initPolycrystalCard() {
    registerCard('poly', {
        busy: 'Building polycrystal…',
        preview: () => {
            const s = requireSource();
            const p = params();
            const estimate = p.box[0] * p.box[1] * p.box[2] / (cellVolume(s.cell) / s.count);
            const k = previewScale(estimate, 4000);
            const pc = buildPolycrystal(s, { ...p, box: p.box.map((v) => v * k) });
            const note = k < 1 ? `<br><span>Preview box scaled to ${(k * 100).toFixed(0)}%; about ${Math.round(estimate).toLocaleString()} atoms when built</span>` : '';
            return { structure: pc, colorBy: 'tag', caption: `<b>${p.grains} grains</b> · ${p.box.map((v) => v.toFixed(0)).join(' × ')} Å${note}` };
        },
        build: () => {
            const pc = buildPolycrystal(requireSource(), params());
            state.view.color = 'tag'; $('v-color').value = 'tag';
            state.view.bonds = pc.count < 60000; $('v-bonds').checked = state.view.bonds;
            setStructure(pc, { message: `${pc.title}: ${pc.count.toLocaleString()} atoms` });
        },
    });
}
