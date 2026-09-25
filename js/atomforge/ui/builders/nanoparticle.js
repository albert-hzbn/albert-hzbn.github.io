// Nanoparticle card: geometric shapes or a Wulff construction. The preview
// shows the particle (scaled down if it would exceed the preview budget).
import { cellVolume, norm } from '../../core/index.js';
import { SHAPES, buildNanoparticle, buildWulff, wulffShape, shapeOutline } from '../../builders/index.js';
import { FACET_COLORS } from '../../render/preview.js';
import { $, num, segmented, removableRow } from '../dom.js';
import { setStructure, requireSource } from '../history.js';
import { registerCard, previewScale } from './card.js';

// Approximate particle volumes (Å³), used to estimate the atom count.
const VOLUME = {
    sphere: (R) => (4 / 3) * Math.PI * R ** 3,
    cube: (R) => 8 * R ** 3,
    octahedron: (R) => (4 / 3) * R ** 3,
    cuboctahedron: (R) => (20 / 3) * R ** 3,
    truncoct: (R) => 4 * R ** 3,
    cylinder: (R, H) => Math.PI * R * R * H,
    ellipsoid: (R, H, E) => (4 / 3) * Math.PI * R * (E || R) * (H / 2 || R),
};

const mode = () => $('np-mode').querySelector('[aria-pressed="true"]').dataset.v;

function facetRow(h = 1, k = 1, l = 1, e = 1) {
    return removableRow(`<input type="number" value="${h}" aria-label="h"><input type="number" value="${k}" aria-label="k"><input type="number" value="${l}" aria-label="l"><input type="number" step="0.01" value="${e}" aria-label="Surface energy">`, 'Remove facet');
}

function readFacets() {
    const fs = [...$('wf-facets').querySelectorAll('.tr:not(.head)')].map((tr) => {
        const v = [...tr.querySelectorAll('input')].map((i) => parseFloat(i.value));
        return { hkl: v.slice(0, 3), energy: v[3] };
    }).filter((f) => f.energy > 0 && f.hkl.some((x) => x));
    if (!fs.length) throw new Error('Add at least one facet family with a positive energy.');
    return fs;
}

// Build with a linear scale factor (preview uses < 1 for large particles).
// The preview omits the vacuum box so the particle fills the view.
function makeParticle(scale = 1, preview = false) {
    const s = requireSource();
    const common = { centerOn: $('np-center').value, vacuum: preview ? 0 : num('np-vac') };
    if (mode() === 'wulff') {
        return buildWulff(s, readFacets(), { radius: num('wf-r') * scale, cubicSymmetry: $('wf-cubic').checked, ...common });
    }
    return buildNanoparticle(s, { shape: $('np-shape').value, radius: num('np-r') * scale, height: num('np-h') * scale, radiusY: num('np-ry') * scale, ...common });
}

function estimateAtoms() {
    const s = requireSource();
    const perAtom = cellVolume(s.cell) / s.count;
    const vol = mode() === 'wulff' ? (4 / 3) * Math.PI * (num('wf-r') * 1.15) ** 3 : VOLUME[$('np-shape').value](num('np-r'), num('np-h'), num('np-ry'));
    return vol / perAtom;
}

const viewMode = () => $('np-view').querySelector('[aria-pressed="true"]').dataset.v;

// Outline drawn in the preview, at the same scale as the previewed atoms.
function outline(scale) {
    if (mode() === 'wulff') {
        const w = wulffShape(requireSource(), readFacets(), { radius: num('wf-r') * scale, cubicSymmetry: $('wf-cubic').checked });
        if (!w.closed) throw new Error('These facets do not enclose a shape: add families or enable cubic symmetry.');
        return { kind: 'poly', ...w, extent: Math.max(...w.vertices.map(norm)) };
    }
    const o = shapeOutline($('np-shape').value, num('np-r') * scale, num('np-h') * scale, num('np-ry') * scale);
    o.extent = o.kind === 'poly' ? Math.max(...o.vertices.map(norm))
        : o.kind === 'cylinder' ? Math.hypot(o.radius, o.height / 2)
            : o.kind === 'ellipsoid' ? Math.max(...o.radii) : o.radius;
    return o;
}

// Colour key for the Wulff facet families, with each family's share of the surface.
function facetLegend(w) {
    return '<div class="facet-legend">' + readFacets().map((f, i) => {
        const pct = 100 * (w.fractions[i] || 0);
        const label = `{${f.hkl.join('')}}`;
        return pct > 0.05
            ? `<span><i style="background:${FACET_COLORS[i % FACET_COLORS.length]}"></i>${label} ${pct.toFixed(1)}%</span>`
            : `<span class="off"><i style="background:${FACET_COLORS[i % FACET_COLORS.length]};opacity:.35"></i>${label} not exposed</span>`;
    }).join('') + '</div>';
}

export function initNanoparticleCard() {
    $('np-shape').innerHTML = Object.entries(SHAPES).map(([k, v]) => `<option value="${k}">${v.label}</option>`).join('');
    segmented($('np-mode'), (v) => {
        $('np-shape-opts').hidden = v !== 'shape';
        $('np-wulff-opts').hidden = v !== 'wulff';
        $('np-mode').dispatchEvent(new Event('change', { bubbles: true }));
    });
    segmented($('np-view'), () => $('np-view').dispatchEvent(new Event('change', { bubbles: true })));
    const facets = $('wf-facets');
    facets.innerHTML = '<div class="tr head"><span>h</span><span>k</span><span>l</span><span>γ (J/m²)</span><span></span></div>';
    facets.append(facetRow(1, 1, 1, 1.0), facetRow(1, 0, 0, 1.15), facetRow(1, 1, 0, 1.25));
    $('wf-add').addEventListener('click', () => { facets.appendChild(facetRow(2, 1, 1, 1.3)); facets.dispatchEvent(new Event('change', { bubbles: true })); });
    registerCard('nano', {
        busy: 'Cutting nanoparticle…',
        preview: () => {
            const k = previewScale(estimateAtoms());
            const view = viewMode();
            const shape = outline(k);
            const np = makeParticle(k, true);
            const wulff = mode() === 'wulff';
            const label = wulff ? 'Wulff construction' : SHAPES[$('np-shape').value].label;
            const atoms = k < 1 ? ` · about ${Math.round(estimateAtoms()).toLocaleString()} atoms` : ` · ${np.count.toLocaleString()} atoms`;
            const note = k < 1 ? `<br><span>Preview at ${(k * 100).toFixed(0)}% size</span>` : '';
            return {
                structure: view === 'shape' ? null : np,
                shape,
                mode: view,
                caption: `<b>${label}</b>${atoms}${note}${wulff ? facetLegend(shape) : ''}`,
            };
        },
        build: () => {
            const np = makeParticle();
            setStructure(np, { message: `${np.title}: ${np.count.toLocaleString()} atoms` });
        },
    });
}

export function selectNanoMode(v) {
    $('np-mode').querySelector(`[data-v="${v}"]`).click();
}
