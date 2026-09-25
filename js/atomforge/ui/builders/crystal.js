// Bulk crystal builder: crystal system → space group → setting, lattice
// constraints, asymmetric-unit sites and presets.
import { normalizeSymbol } from '../../core/index.js';
import { getSettings, settingsForNumber, defaultSetting, crystalSystem, buildCrystal } from '../../crystal/spacegroups.js';
import { PRESETS, presetParameters } from '../../crystal/presets.js';
import { $, num, int, toast, removableRow } from '../dom.js';
import { setStructure, setSource } from '../history.js';
import { registerCard } from './card.js';

export function fillPresetSelect() {
    $('cr-preset').innerHTML = '<option value="">Custom…</option>' + PRESETS.map((p) => `<option value="${p.id}">${p.label}</option>`).join('');
}


export function fillSpaceGroups(system) {
    const settings = getSettings();
    const nums = [...new Set(settings.filter((s) => s.system === system).map((s) => s.n))];
    $('cr-sg').innerHTML = nums.map((n) => {
        const st = defaultSetting(n);
        return `<option value="${n}">${n} · ${st.short || st.hm}</option>`;
    }).join('');
}


export function fillSettings(n) {
    const list = settingsForNumber(n);
    const def = defaultSetting(n);
    $('cr-setting').innerHTML = list.map((s) => `<option value="${s.h}" ${s === def ? 'selected' : ''}>${s.hm}${s.choice ? ` (${settingLabel(s.choice)})` : ''}</option>`).join('');
    $('cr-setting').disabled = list.length < 2;
}


function settingLabel(c) {
    if (c === '1' || c === '2') return `origin choice ${c}`;
    if (c === 'H') return 'hexagonal axes';
    if (c === 'R') return 'rhombohedral axes';
    return c;
}


export function currentSetting() {
    const h = int('cr-setting');
    return getSettings().find((s) => s.h === h);
}


// Enforce lattice constraints of the crystal system on the inputs.
export function applyLatticeConstraints() {
    const st = currentSetting();
    if (!st) return;
    const sys = st.system;
    const rh = st.choice === 'R';
    const set = (id, v, dis) => { if (v !== null) $(id).value = v; $(id).disabled = dis; };
    const a = num('cr-a');
    switch (sys) {
        case 'cubic': set('cr-b', a, true); set('cr-c', a, true); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 90, true); break;
        case 'tetragonal': set('cr-b', a, true); set('cr-c', null, false); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 90, true); break;
        case 'hexagonal':
        case 'trigonal':
            if (rh) { set('cr-b', a, true); set('cr-c', a, true); set('cr-al', null, false); set('cr-be', num('cr-al'), true); set('cr-ga', num('cr-al'), true); }
            else { set('cr-b', a, true); set('cr-c', null, false); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 120, true); }
            break;
        case 'orthorhombic': ['cr-b', 'cr-c'].forEach((id) => set(id, null, false)); set('cr-al', 90, true); set('cr-be', 90, true); set('cr-ga', 90, true); break;
        case 'monoclinic': ['cr-b', 'cr-c', 'cr-be'].forEach((id) => set(id, null, false)); set('cr-al', 90, true); set('cr-ga', 90, true); break;
        default: ['cr-b', 'cr-c', 'cr-al', 'cr-be', 'cr-ga'].forEach((id) => set(id, null, false));
    }
}


export function siteRow(sym = 'Cu', x = 0, y = 0, z = 0) {
    const f = (v) => +(+v).toFixed(5);
    return removableRow(`<input value="${sym}" aria-label="Element"><input type="number" step="0.01" value="${f(x)}" aria-label="x"><input type="number" step="0.01" value="${f(y)}" aria-label="y"><input type="number" step="0.01" value="${f(z)}" aria-label="z">`, 'Remove site');
}


export function setSites(sites) {
    const t = $('cr-sites');
    t.innerHTML = '<div class="tr head"><span>Element</span><span>x</span><span>y</span><span>z</span><span></span></div>';
    sites.forEach((s) => t.appendChild(siteRow(...s)));
}


export function readSites() {
    return [...$('cr-sites').querySelectorAll('.tr:not(.head)')].map((tr) => {
        const [e, x, y, z] = tr.querySelectorAll('input');
        const frac = (v) => {
            const s = String(v.value).trim();
            if (s.includes('/')) { const [p, q] = s.split('/'); return parseFloat(p) / parseFloat(q); }
            return parseFloat(s);
        };
        return { symbol: normalizeSymbol(e.value), x: frac(x), y: frac(y), z: frac(z) };
    }).filter((s) => [s.x, s.y, s.z].every(Number.isFinite));
}


export function loadPreset(id) {
    const p = PRESETS.find((x) => x.id === id);
    if (!p) return;
    const sys = crystalSystem(p.n);
    $('cr-system').value = sys;
    fillSpaceGroups(sys);
    $('cr-sg').value = p.n;
    fillSettings(p.n);
    const prm = presetParameters(p);
    $('cr-a').value = prm.a; $('cr-b').value = prm.b; $('cr-c').value = prm.c;
    $('cr-al').value = prm.alpha; $('cr-be').value = prm.beta; $('cr-ga').value = prm.gamma;
    applyLatticeConstraints();
    setSites(p.sites);
}


// Crystal described by the card's inputs (used by the preview and the build).
export function crystalFromPanel() {
    const setting = currentSetting();
    applyLatticeConstraints();
    const sites = readSites();
    if (!sites.length) throw new Error('Add at least one site.');
    const prm = { a: num('cr-a'), b: num('cr-b'), c: num('cr-c'), alpha: num('cr-al'), beta: num('cr-be'), gamma: num('cr-ga') };
    if (![prm.a, prm.b, prm.c].every((v) => v > 0)) throw new Error('Lattice lengths must be positive.');
    const preset = PRESETS.find((p) => p.id === $('cr-preset').value);
    const { structure, multiplicities } = buildCrystal({ setting, ...prm, sites });
    structure.title = preset ? preset.label : `${structure.formula()} (${setting.short || setting.hm})`;
    return { structure, setting, sites, multiplicities };
}

export function buildFromPanel({ quiet = false } = {}) {
    const { structure, setting } = crystalFromPanel();
    setStructure(structure, { message: quiet ? '' : `Built ${structure.formula()} in ${setting.hm}`, spaceGroup: `${setting.hm} (No. ${setting.n})` });
    setSource(structure);
}

export function initCrystalCard() {
    fillPresetSelect();
    $('cr-preset').addEventListener('change', (e) => { if (e.target.value) loadPreset(e.target.value); });
    $('cr-system').addEventListener('change', (e) => {
        fillSpaceGroups(e.target.value);
        fillSettings(int('cr-sg'));
        $('cr-preset').value = '';
        if (e.target.value === 'hexagonal' || e.target.value === 'trigonal') $('cr-ga').value = 120;
        applyLatticeConstraints();
    });
    $('cr-sg').addEventListener('change', () => { fillSettings(int('cr-sg')); $('cr-preset').value = ''; applyLatticeConstraints(); });
    $('cr-setting').addEventListener('change', applyLatticeConstraints);
    ['cr-a', 'cr-al'].forEach((id) => $(id).addEventListener('input', applyLatticeConstraints));
    $('cr-add-site').addEventListener('click', () => { $('cr-sites').appendChild(siteRow('O', 0, 0, 0)); $('cr-sites').dispatchEvent(new Event('change', { bubbles: true })); });
    registerCard('crystal', {
        busy: 'Building crystal…',
        preview: () => {
            const { structure, setting, sites, multiplicities } = crystalFromPanel();
            return {
                structure,
                caption: `<b>${structure.formula()}</b> · ${setting.hm} · ${structure.count} atoms<br><span>Multiplicities: ${sites.map((x, i) => `${x.symbol} ${multiplicities[i]}`).join(', ')}</span>`,
            };
        },
        build: () => { try { buildFromPanel(); } catch (e) { toast(e.message, true); } },
    });
}
