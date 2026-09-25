// View panel: style, colours, cell and images, planes and directions, camera.
import { vecmat, reciprocal } from '../../core/index.js';
import { state, viewer } from '../context.js';
import { $, num, int, toast, segmented } from '../dom.js';
import { computeAnalysis, render, renderOverlays, boundsOf } from '../scene.js';
import { refreshLegend, refreshCoordination, refreshInfo } from './overview.js';

export function initView() {
    segmented($('v-style'), (v) => { state.view.style = v; render(); });
    const bindRange = (id, key, digits = 2) => {
        const out = $(id + '-out');
        $(id).addEventListener('input', () => {
            state.view[key] = num(id);
            out.value = num(id).toFixed(digits);
            if (key === 'tolerance') computeAnalysis();
            render();
            if (key === 'tolerance') { refreshCoordination(); refreshInfo(); }
        });
    };
    bindRange('v-atom', 'atomScale');
    bindRange('v-bond', 'bondRadius');
    bindRange('v-tol', 'tolerance');
    $('v-bonds').addEventListener('change', (e) => { state.view.bonds = e.target.checked; render(); });
    $('v-labels').addEventListener('change', (e) => { state.view.labels = e.target.checked; renderOverlays(); });
    $('v-color').addEventListener('change', (e) => { state.view.color = e.target.value; render(); refreshLegend(); });
    $('v-cell').addEventListener('change', (e) => { state.view.cell = e.target.checked; render(); });
    $('v-boundary').addEventListener('change', (e) => { state.view.boundary = e.target.checked; render(); });
    ['v-ra', 'v-rb', 'v-rc'].forEach((id, k) => $(id).addEventListener('change', () => {
        const v = Math.max(1, Math.min(10, int(id) || 1));
        $(id).value = v;
        state.view.reps[k] = v;
        render({ fit: true });
    }));
    const updPlane = () => {
        state.view.plane = $('pl-on').checked ? { hkl: [int('pl-h'), int('pl-k'), int('pl-l')], d: num('pl-d') } : null;
        if (state.view.plane && state.view.plane.hkl.every((x) => !x)) { state.view.plane = null; toast('Miller indices cannot all be zero.', true); }
        renderOverlays();
    };
    ['pl-on', 'pl-h', 'pl-k', 'pl-l', 'pl-d'].forEach((id) => $(id).addEventListener('change', updPlane));
    const updDir = () => {
        state.view.direction = $('dir-on').checked ? [int('dir-u'), int('dir-v'), int('dir-w')] : null;
        if (state.view.direction && state.view.direction.every((x) => !x)) state.view.direction = null;
        renderOverlays();
    };
    ['dir-on', 'dir-u', 'dir-v', 'dir-w'].forEach((id) => $(id).addEventListener('change', updDir));
    segmented($('v-proj'), (v) => viewer.setProjection(v));
    $('v-rotate').addEventListener('change', (e) => { viewer.autoRotate = e.target.checked; });
    $('va-go').addEventListener('click', () => {
        const s = state.structure;
        const uvw = [int('va-u'), int('va-v'), int('va-w')];
        const dir = s.periodic ? vecmat(uvw, s.cell) : uvw;
        viewer.viewAlong(dir, s.periodic ? s.cell[1] : [0, 1, 0]);
    });
    document.querySelectorAll('[data-view]').forEach((b) => b.addEventListener('click', () => {
        const s = state.structure, v = b.dataset.view;
        if (v === 'fit') return viewer.fit(boundsOf(state.display));
        const idx = { a: 0, b: 1, c: 2 }[v];
        if (s.periodic) {
            // Look down the axis: along the reciprocal direction so the other two lie in the screen plane.
            const rec = reciprocal(s.cell);
            const up = idx === 2 ? s.cell[1] : s.cell[2];
            viewer.viewAlong(rec[idx], up);
        } else viewer.viewAlong([[1, 0, 0], [0, 1, 0], [0, 0, 1]][idx], idx === 2 ? [0, 1, 0] : [0, 0, 1]);
    }));
}

