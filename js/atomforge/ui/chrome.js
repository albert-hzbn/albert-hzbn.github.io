// Application chrome: menus, examples, file open/export, keyboard, theme.
import { supercell } from '../core/index.js';
import { listTiltBoundaries } from '../builders/index.js';
import { parseStructure, WRITERS } from '../io/index.js';
import { state, viewer } from './context.js';
import { $, toast, withBusy, download } from './dom.js';
import { setStructure, setSource, undo, redo } from './history.js';
import { boundsOf } from './scene.js';
import { setSelection, selectAll, deleteSelection } from './panels/selection.js';
import { loadPreset, buildFromPanel, applyLatticeConstraints, refreshGBList, setGBChoice, selectNanoMode, setFirstSolute, openCard } from './builders/index.js';
import { setMode } from './interaction.js';
import { activateTab } from './tabs.js';

const EXAMPLES = [
    { label: 'Cu (FCC)', run: () => presetExample('fcc') },
    { label: 'Fe (BCC)', run: () => presetExample('bcc') },
    { label: 'Mg (HCP)', run: () => presetExample('hcp') },
    { label: 'Si (diamond)', run: () => presetExample('diamond') },
    { label: 'NaCl (rock salt)', run: () => presetExample('nacl') },
    { label: 'SrTiO₃ (perovskite)', run: () => presetExample('perovskite') },
    { label: 'Al₂O₃ (corundum)', run: () => presetExample('corundum') },
    { label: 'MgAl₂O₄ (spinel)', run: () => presetExample('spinel') },
    { sep: true },
    { label: 'Cu Σ5 (310)[001] grain boundary', run: () => { presetExample('fcc', true); openCard('gb'); $('gb-axis').value = '0,0,1'; refreshGBList(); setGBChoice(listTiltBoundaries([0, 0, 1], 7, 51).find((g) => g.sigma === 5 && g.angle < 40)); $('gb-n1').value = 4; $('gb-n3').value = 3; $('gb-build').click(); } },
    { label: 'Pt (111) slab', run: () => { presetExample('fcc', true, { Cu: 'Pt', a: 3.924 }); openCard('slab'); $('sl-h').value = 1; $('sl-k').value = 1; $('sl-l').value = 1; $('sl-layers').value = 4; $('sl-ra').value = 3; $('sl-rb').value = 3; $('sl-build').click(); } },
    { label: 'Au Wulff nanoparticle', run: () => { presetExample('fcc', true, { Cu: 'Au', a: 4.078 }); openCard('nano'); selectNanoMode('wulff'); $('np-build').click(); } },
    { label: 'Cu polycrystal (8 grains)', run: () => { presetExample('fcc', true); openCard('poly'); $('pc-x').value = $('pc-y').value = $('pc-z').value = 45; $('pc-build').click(); } },
    { label: 'Amorphous SiO₂', run: () => { openCard('amorphous'); $('am-build').click(); } },
    { label: 'Cu intrinsic stacking fault', run: () => { presetExample('fcc', true); openCard('sf'); $('sf-u').value = 1; $('sf-build').click(); } },
    { label: 'Cu₃Au-type random alloy', run: () => { presetExample('fcc', true); setStructure(supercell(state.structure, 4, 4, 4), { message: '' }); openCard('ss'); setFirstSolute('Au'); $('ss-build').click(); } },
];


function presetExample(id, quiet = false, subst = null) {
    $('cr-preset').value = id;
    loadPreset(id);
    if (subst) {
        const row = $('cr-sites').querySelector('.tr:not(.head) input');
        for (const [from, to] of Object.entries(subst)) if (from !== 'a' && row.value === from) row.value = to;
        if (subst.a) { $('cr-a').value = subst.a; applyLatticeConstraints(); }
        $('cr-preset').value = '';
    }
    buildFromPanel({ quiet });
    if (subst) state.structure.title = `${state.structure.formula().replace(/\d+/g, '')} (FCC)`;
}


export function initChrome() {
    document.querySelectorAll('.tabs [data-tab]').forEach((b) => b.addEventListener('click', () => activateTab(b.dataset.tab)));
    $('panel-toggle').addEventListener('click', () => document.body.classList.toggle('panel-open'));
    // Menus
    document.querySelectorAll('.menu').forEach((m) => {
        m.querySelector('[data-menu]').addEventListener('click', (e) => {
            e.stopPropagation();
            document.querySelectorAll('.menu.open').forEach((o) => { if (o !== m) o.classList.remove('open'); });
            m.classList.toggle('open');
        });
    });
    document.addEventListener('click', () => document.querySelectorAll('.menu.open').forEach((o) => o.classList.remove('open')));
    const ex = $('examples-list');
    for (const item of EXAMPLES) {
        if (item.sep) { ex.appendChild(document.createElement('hr')); continue; }
        const b = document.createElement('button');
        b.textContent = item.label;
        b.addEventListener('click', () => { try { item.run(); } catch (e) { toast(e.message, true); } });
        ex.appendChild(b);
    }
    document.querySelectorAll('[data-export]').forEach((b) => b.addEventListener('click', () => exportAs(b.dataset.export)));

    // Files
    $('btn-open').addEventListener('click', () => $('file-input').click());
    $('file-input').addEventListener('change', (e) => { if (e.target.files[0]) openFile(e.target.files[0]); e.target.value = ''; });
    const vp = $('viewport');
    let dragDepth = 0;
    vp.addEventListener('dragenter', (e) => { e.preventDefault(); dragDepth++; vp.classList.add('dragging'); });
    vp.addEventListener('dragover', (e) => e.preventDefault());
    vp.addEventListener('dragleave', () => { if (--dragDepth <= 0) { dragDepth = 0; vp.classList.remove('dragging'); } });
    vp.addEventListener('drop', (e) => {
        e.preventDefault(); dragDepth = 0; vp.classList.remove('dragging');
        if (e.dataTransfer.files[0]) openFile(e.dataTransfer.files[0]);
    });

    $('btn-undo').addEventListener('click', undo);
    $('btn-redo').addEventListener('click', redo);
    $('btn-theme').addEventListener('click', () => setTheme(document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark'));
    document.querySelectorAll('#modebar [data-mode]').forEach((b) => b.addEventListener('click', () => setMode(b.dataset.mode)));

    window.addEventListener('keydown', (e) => {
        const typing = /INPUT|SELECT|TEXTAREA/.test(document.activeElement?.tagName);
        const mod = e.ctrlKey || e.metaKey;
        if (mod && e.key.toLowerCase() === 'z' && !typing) { e.preventDefault(); e.shiftKey ? redo() : undo(); return; }
        if (mod && e.key.toLowerCase() === 'y' && !typing) { e.preventDefault(); redo(); return; }
        if (mod && e.key.toLowerCase() === 'o') { e.preventDefault(); $('file-input').click(); return; }
        if (typing) return;
        if (mod && e.key.toLowerCase() === 'a') { e.preventDefault(); selectAll(); return; }
        if (e.key === 'Delete' || e.key === 'Backspace') { if (state.selection.size) { e.preventDefault(); deleteSelection(); } return; }
        if (e.key === 'Escape') { setSelection([]); setMode('select'); document.body.classList.remove('panel-open'); return; }
        if (mod || e.altKey) return;
        const k = e.key.toLowerCase();
        if (k === 'f') viewer.fit(boundsOf(state.display));
        else if (k === 's') setMode('select');
        else if (k === 'd') setMode('distance');
        else if (k === 'a') setMode('angle');
        else if (k === '1' || k === '2' || k === '3') document.querySelector(`[data-view="${'abc'[+k - 1]}"]`).click();
    });
}


export function setTheme(t) {
    document.documentElement.dataset.theme = t;
    try { localStorage.setItem('studio-theme', t); } catch { /* ignore */ }
    viewer.setTheme(t === 'dark');
}


export async function openFile(file) {
    if (file.size > 80 * 1024 * 1024) return toast('File is larger than 80 MB.', true);
    const text = await file.text();
    await withBusy(`Reading ${file.name}…`, () => {
        const { structure, format } = parseStructure(file.name, text);
        if (!structure.count) throw new Error('No atoms found in the file.');
        state.view.color = 'element'; $('v-color').value = 'element';
        state.view.bonds = structure.count < 80000; $('v-bonds').checked = state.view.bonds;
        setStructure(structure, { message: `Opened ${file.name} (${format}, ${structure.count.toLocaleString()} atoms${structure.frames > 1 ? `, last of ${structure.frames} frames` : ''})` });
        if (structure.periodic) setSource(structure, file.name);
    });
}


function exportAs(kind) {
    const s = state.structure;
    const base = (s.formula() || 'structure').replace(/[^\w]+/g, '');
    try {
        if (kind === 'png' || kind === 'png-t') {
            const url = viewer.screenshot(kind === 'png' ? 2 : 4, kind === 'png-t');
            const a = document.createElement('a'); a.href = url; a.download = `${base}.png`; a.click();
            return;
        }
        const w = WRITERS[kind];
        const name = kind === 'poscar' ? `${base}.vasp` : `${base}.${w.ext}`;
        download(name, w.fn(s));
        toast(`Saved ${name}`);
    } catch (e) { toast(e.message, true); }
}

