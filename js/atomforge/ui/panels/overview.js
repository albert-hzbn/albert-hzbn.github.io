// Panels that summarise the structure: HUD, status bar, legend, info,
// coordination, and the species / cell inputs that depend on it.
import { cellParameters, cellVolume, element } from '../../core/index.js';
import { state } from '../context.js';
import { $, fmt } from '../dom.js';
import { computeAnalysis, render, elementColor, TAG_COLORS, CN_COLORS, ramp } from '../scene.js';
import { refreshSelectionUI } from './selection.js';

export function refreshAll({ fit = false } = {}) {
    computeAnalysis();
    render({ fit });
    refreshPanels();
}


export function refreshPanels() {
    const s = state.structure;
    $('hud-title').textContent = s.title || s.formula() || 'Empty structure';
    const sub1 = [s.formula(), state.spaceGroupLabel, `${s.count.toLocaleString()} atoms`].filter(Boolean);
    $('hud-sub').textContent = sub1.join(' · ');
    $('st-atoms').textContent = s.count.toLocaleString();
    $('st-formula').textContent = s.formula();
    if (s.periodic) {
        const p = cellParameters(s.cell);
        $('st-cell').textContent = `a ${fmt(p.a)} · b ${fmt(p.b)} · c ${fmt(p.c)} Å · ${fmt(p.alpha, 1)}° ${fmt(p.beta, 1)}° ${fmt(p.gamma, 1)}°`;
    } else $('st-cell').textContent = 'no cell (cluster)';
    $('btn-undo').disabled = !state.undo.length;
    $('btn-redo').disabled = !state.redo.length;
    refreshSelectionUI();
    refreshLegend();
    refreshInfo();
    refreshSpeciesSelects();
    refreshCellInputs();
    refreshCoordination();
}


export function refreshLegend() {
    const s = state.structure, box = $('legend');
    box.innerHTML = '';
    const counts = s.speciesCounts();
    if (state.view.color === 'element') {
        for (const [sym, n] of Object.entries(counts)) {
            const li = document.createElement('div');
            li.className = 'li';
            li.innerHTML = `<input type="color" value="${elementColor(sym)}" aria-label="Colour of ${sym}"><span><b>${sym}</b> <span class="note">${element(sym).name}</span></span><span class="count">${n.toLocaleString()}</span><button class="eye" aria-pressed="${!state.hidden.has(sym)}" title="Show / hide ${sym}">${eyeIcon}</button>`;
            li.querySelector('input').addEventListener('input', (e) => { state.customColors[sym] = e.target.value; render(); });
            li.querySelector('.eye').addEventListener('click', () => {
                if (state.hidden.has(sym)) state.hidden.delete(sym); else state.hidden.add(sym);
                render(); refreshLegend();
            });
            box.appendChild(li);
        }
    } else if (state.view.color === 'tag') {
        const tags = [...new Set(s.tags)].sort((a, b) => a - b);
        box.innerHTML = tags.slice(0, 24).map((t) => `<div class="li"><span style="width:14px;height:14px;border-radius:50%;background:${TAG_COLORS[t % TAG_COLORS.length]}"></span><span>${t === 0 ? 'Untagged' : 'Region ' + t}</span><span class="count">${s.tags.filter((x) => x === t).length}</span><span></span></div>`).join('')
            + (tags.length > 24 ? `<p class="note">… ${tags.length - 24} more</p>` : '');
    } else if (state.view.color === 'cn') {
        const hist = {};
        state.cn.forEach((c) => { hist[c] = (hist[c] || 0) + 1; });
        box.innerHTML = Object.keys(hist).map(Number).sort((a, b) => a - b).map((c) => `<div class="li"><span style="width:14px;height:14px;border-radius:50%;background:${CN_COLORS[Math.min(c, CN_COLORS.length - 1)]}"></span><span>CN ${c}</span><span class="count">${hist[c]}</span><span></span></div>`).join('');
    } else {
        const { lo, hi } = state.heightRange;
        box.innerHTML = `<div style="height:12px;border-radius:6px;background:linear-gradient(90deg,${[0, 0.25, 0.5, 0.75, 1].map(ramp).join(',')})"></div><div class="row c2" style="font-family:var(--mono);font-size:11.5px"><span>${fmt(lo, 2)} Å</span><span style="text-align:right">${fmt(hi, 2)} Å</span></div>`;
    }
}


const eyeIcon = '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="1.8"><path d="M2 12s3.6-7 10-7 10 7 10 7-3.6 7-10 7S2 12 2 12z"/><circle cx="12" cy="12" r="3"/></svg>';

export function refreshInfo() {
    const s = state.structure;
    const rows = [['Formula', s.formula() || '—'], ['Atoms', s.count.toLocaleString()]];
    if (state.spaceGroupLabel) rows.push(['Space group', state.spaceGroupLabel]);
    if (s.periodic) {
        const p = cellParameters(s.cell);
        rows.push(['a, b, c (Å)', `${fmt(p.a)}, ${fmt(p.b)}, ${fmt(p.c)}`]);
        rows.push(['α, β, γ (°)', `${fmt(p.alpha, 2)}, ${fmt(p.beta, 2)}, ${fmt(p.gamma, 2)}`]);
        rows.push(['Volume (Å³)', fmt(cellVolume(s.cell), 3)]);
        rows.push(['Density (g/cm³)', fmt(s.density(), 4)]);
        rows.push(['Volume / atom (Å³)', s.count ? fmt(cellVolume(s.cell) / s.count, 3) : '—']);
    } else if (s.count) {
        const { lo, hi } = s.bounds();
        rows.push(['Extent (Å)', [0, 1, 2].map((k) => fmt(hi[k] - lo[k], 2)).join(' × ')]);
    }
    rows.push(['Mass (u)', fmt(s.mass(), 3)]);
    rows.push(['Bonds shown', state.bonds.length.toLocaleString()]);
    for (const [sym, n] of Object.entries(s.speciesCounts())) rows.push([`${sym} fraction`, `${(100 * n / s.count).toFixed(2)} %`]);
    $('info').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
}


export function refreshCoordination() {
    const s = state.structure;
    const per = {};
    s.symbols.forEach((sym, i) => { (per[sym] = per[sym] || []).push(state.cn[i] || 0); });
    $('cn-table').innerHTML = Object.entries(per).map(([sym, list]) => {
        const mean = list.reduce((a, b) => a + b, 0) / list.length;
        const hist = {};
        list.forEach((c) => { hist[c] = (hist[c] || 0) + 1; });
        const top = Object.entries(hist).sort((a, b) => b[1] - a[1]).slice(0, 3).map(([c, n]) => `${c}:${n}`).join(' ');
        return `<dt>${sym} (mean ${mean.toFixed(2)})</dt><dd>${top}</dd>`;
    }).join('') || '<dt>—</dt><dd></dd>';
}


export function refreshSpeciesSelects() {
    const species = Object.keys(state.structure.speciesCounts());
    const fill = (id, withAll = false) => {
        const el = $(id), cur = el.value;
        el.innerHTML = (withAll ? '<option value="">All</option>' : '') + species.map((x) => `<option>${x}</option>`).join('');
        if ([...el.options].some((o) => o.value === cur)) el.value = cur;
    };
    fill('sel-el'); fill('rn-from'); fill('ss-host'); fill('rdf-a', true); fill('rdf-b', true);
}


export function refreshCellInputs() {
    const s = state.structure;
    const ids = ['ce-a', 'ce-b', 'ce-c', 'ce-al', 'ce-be', 'ce-ga'];
    if (!s.periodic) { ids.forEach((id) => { $(id).value = ''; }); return; }
    const p = cellParameters(s.cell);
    [p.a, p.b, p.c, p.alpha, p.beta, p.gamma].forEach((v, k) => { $(ids[k]).value = +v.toFixed(k < 3 ? 5 : 4); });
    $('title-in').value = s.title || '';
}

