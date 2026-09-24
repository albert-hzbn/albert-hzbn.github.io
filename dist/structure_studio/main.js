// AtomForge Web: entry point.
import { loadSpaceGroups } from './crystal/spacegroups.js';
import { parseStructure } from './io/index.js';
import { state, viewer, bus } from './ui/context.js';
import { $, toast } from './ui/dom.js';
import { setStructure, setSource } from './ui/history.js';
import { refreshAll, refreshPanels } from './ui/panels/overview.js';
import { initBuilders, loadPreset, buildFromPanel, crystalFromPanel } from './ui/builders/index.js';
import { initEdit } from './ui/panels/edit.js';
import { initView } from './ui/panels/view.js';
import { initAnalyze, refreshMeasurements } from './ui/panels/analyze.js';
import { initInteraction, setMode } from './ui/interaction.js';
import { initChrome } from './ui/chrome.js';

// Every structure change re-renders the scene and refreshes the panels.
bus.on('structure', ({ fit }) => {
    // Region colouring is meaningless without regions: fall back to elements.
    if (state.view.color === 'tag' && new Set(state.structure.tags).size < 2) {
        state.view.color = 'element';
        $('v-color').value = 'element';
    }
    refreshAll({ fit });
});

async function start() {
    viewer.setTheme(document.documentElement.dataset.theme === 'dark');
    $('hud-title').textContent = 'Loading space groups…';
    await loadSpaceGroups();
    initChrome();
    initBuilders();
    initEdit();
    initView();
    initAnalyze();
    initInteraction();
    refreshMeasurements();
    setMode('select');
    // Restore the last session if there is one, else start with FCC copper.
    let restored = false;
    try {
        const last = localStorage.getItem('studio-last');
        if (last) {
            const { structure } = parseStructure('last.xyz', last);
            if (structure.count) { setStructure(structure, { record: false, message: 'Restored your last structure' }); restored = true; }
        }
    } catch { /* ignore */ }
    loadPreset('fcc');
    $('cr-preset').value = 'fcc';
    if (restored) {
        // The restored structure may be a slab or cluster: use FCC Cu as the source.
        setSource(crystalFromPanel().structure, 'FCC metal (Cu)');
    } else {
        buildFromPanel({ quiet: true });
        state.undo = [];
        refreshPanels();
    }
}


start().catch((e) => { console.error(e); toast('Could not start: ' + e.message, true); });

// Expose for debugging in the console.
window.studio = { state, viewer };
