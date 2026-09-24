// Surface slab card.
import { supercell } from '../../core/index.js';
import { buildSlab } from '../../builders/index.js';
import { num, int } from '../dom.js';
import { setStructure, requireSource } from '../history.js';
import { registerCard } from './card.js';

function makeSlab(preview = false) {
    const s = requireSource();
    const hkl = [int('sl-h'), int('sl-k'), int('sl-l')];
    const layers = int('sl-layers');
    if (!(layers >= 1)) throw new Error('Layers must be at least 1.');
    let slab = buildSlab(s, hkl, layers, num('sl-vac'));
    const rep = [int('sl-ra') || 1, int('sl-rb') || 1];
    if (rep[0] > 1 || rep[1] > 1) {
        const scaled = preview ? rep.map((r) => Math.max(1, Math.min(r, Math.floor(Math.sqrt(4000 / slab.count)) || 1))) : rep;
        slab = supercell(slab, scaled[0], scaled[1], 1);
    }
    return { slab, hkl };
}

export function initSlabCard() {
    registerCard('slab', {
        busy: 'Cutting slab…',
        preview: () => {
            const { slab, hkl } = makeSlab(true);
            return { structure: slab, caption: `<b>(${hkl.join(' ')})</b> · ${slab.count} atoms · ${int('sl-layers')} layers, ${num('sl-vac')} Å vacuum` };
        },
        build: () => {
            const { slab, hkl } = makeSlab();
            setStructure(slab, { message: `(${hkl.join(' ')}) slab: ${slab.count} atoms` });
        },
    });
}
