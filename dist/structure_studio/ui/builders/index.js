// Build tab: builder cards, quick-access tiles and the source-crystal bar.
import { state, bus } from '../context.js';
import { $ } from '../dom.js';
import { setSource } from '../history.js';
import { initCrystalCard } from './crystal.js';
import { initTransformCard } from './transform.js';
import { initSlabCard } from './slab.js';
import { initGBCard } from './grain-boundary.js';
import { initNanoparticleCard } from './nanoparticle.js';
import { initPolycrystalCard } from './polycrystal.js';
import { initAmorphousCard } from './amorphous.js';
import { initStackingFaultCard } from './stacking-fault.js';
import { initSolidSolutionCard } from './solid-solution.js';
import { openCard, refreshOpenPreview } from './card.js';

export { loadPreset, buildFromPanel, applyLatticeConstraints, crystalFromPanel } from './crystal.js';
export { refreshGBList, setGBChoice } from './grain-boundary.js';
export { selectNanoMode } from './nanoparticle.js';
export { setFirstSolute } from './solid-solution.js';
export { openCard } from './card.js';

function refreshSourceBar() {
    const s = state.source;
    $('src-name').textContent = s ? state.sourceLabel : 'none yet';
    $('src-meta').textContent = s ? `${s.formula()} · ${s.count} atoms` : 'build or open a crystal';
    // "Use current" only makes sense when the current structure is a different periodic one.
    $('src-use').hidden = !state.structure.periodic || state.structure === state.source;
}

export function initBuilders() {
    initCrystalCard();
    initTransformCard();
    initSlabCard();
    initGBCard();
    initNanoparticleCard();
    initPolycrystalCard();
    initAmorphousCard();
    initStackingFaultCard();
    initSolidSolutionCard();

    document.querySelectorAll('.builder-tiles [data-open]').forEach((t) => t.addEventListener('click', () => openCard(t.dataset.open)));
    $('src-use').addEventListener('click', () => setSource(state.structure));
    bus.on('source', () => { refreshSourceBar(); refreshOpenPreview(); });
    // Cards that act on the current structure refresh their preview when it changes.
    bus.on('structure', () => { refreshSourceBar(); refreshOpenPreview(); });
    refreshSourceBar();
}
