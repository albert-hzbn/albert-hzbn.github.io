// Viewport interaction: hover card, click and box selection, measuring modes.
import { element, sub, norm, dot, cross, unit } from '../core/index.js';
import { state, viewer } from './context.js';
import { $, fmt, toast } from './dom.js';
import { renderHighlights, renderOverlays } from './scene.js';
import { setSelection } from './panels/selection.js';
import { refreshMeasurements } from './panels/analyze.js';
import { activateTab } from './tabs.js';

export function initInteraction() {
    const vp = $('viewport'), canvas = viewer.renderer.domElement;
    const card = $('hover-card'), rect = $('select-rect');
    let down = null, boxing = false, hoverRaf = 0, lastMove = null;

    canvas.addEventListener('pointerdown', (e) => {
        down = { x: e.clientX, y: e.clientY, shift: e.shiftKey };
        if (e.shiftKey && state.mode === 'select' && e.button === 0) {
            boxing = true;
            viewer.controls.enabled = false;
            rect.style.display = 'block';
            updateRect(e);
        }
    });
    const updateRect = (e) => {
        const r = vp.getBoundingClientRect();
        const x0 = Math.min(down.x, e.clientX) - r.left, y0 = Math.min(down.y, e.clientY) - r.top;
        Object.assign(rect.style, { left: x0 + 'px', top: y0 + 'px', width: Math.abs(e.clientX - down.x) + 'px', height: Math.abs(e.clientY - down.y) + 'px' });
    };
    window.addEventListener('pointermove', (e) => {
        if (boxing) { updateRect(e); return; }
        if (e.target !== canvas) { card.style.display = 'none'; return; }
        lastMove = e;
        if (!hoverRaf) hoverRaf = requestAnimationFrame(() => { hoverRaf = 0; hover(lastMove); });
    });
    window.addEventListener('pointerup', (e) => {
        if (boxing) {
            boxing = false;
            viewer.controls.enabled = true;
            rect.style.display = 'none';
            const ids = viewer.atomsInRect(down.x, down.y, e.clientX, e.clientY);
            const next = e.ctrlKey || e.metaKey ? new Set([...state.selection, ...ids]) : new Set(ids);
            setSelection([...next]);
            down = null;
            return;
        }
        if (!down || e.target !== canvas) { down = null; return; }
        const moved = Math.hypot(e.clientX - down.x, e.clientY - down.y);
        down = null;
        if (moved > 4 || e.button !== 0) return;
        click(e);
    });
    canvas.addEventListener('pointerleave', () => { card.style.display = 'none'; });

    function hover(e) {
        const hit = viewer.pick(e.clientX, e.clientY);
        if (!hit) { card.style.display = 'none'; canvas.style.cursor = ''; return; }
        canvas.style.cursor = 'pointer';
        const s = state.structure, i = hit.index;
        const el = element(s.symbols[i]);
        const f = s.periodic ? s.toFractional(s.positions[i]) : null;
        card.innerHTML = `<b>${el.symbol}</b> <span class="note">${el.name} · #${i + 1}${hit.image ? ' · image' : ''}</span><br>
            <span class="mono">xyz ${hit.position.map((v) => fmt(v)).join('  ')}</span><br>
            ${f ? `<span class="mono">frac ${f.map((v) => fmt(v, 4)).join('  ')}</span><br>` : ''}
            <span class="mono">CN ${state.cn[i] ?? '—'}${s.tags[i] ? ` · region ${s.tags[i]}` : ''}</span>`;
        const r = vp.getBoundingClientRect();
        let x = e.clientX - r.left + 16, y = e.clientY - r.top + 16;
        card.style.display = 'block';
        if (x + card.offsetWidth > r.width - 8) x = e.clientX - r.left - card.offsetWidth - 16;
        if (y + card.offsetHeight > r.height - 8) y = e.clientY - r.top - card.offsetHeight - 16;
        card.style.left = x + 'px'; card.style.top = y + 'px';
    }

    function click(e) {
        const hit = viewer.pick(e.clientX, e.clientY);
        if (state.mode === 'select') {
            if (!hit) { if (!(e.ctrlKey || e.metaKey)) setSelection([]); return; }
            const next = new Set(e.ctrlKey || e.metaKey ? state.selection : []);
            if ((e.ctrlKey || e.metaKey) && next.has(hit.index)) next.delete(hit.index); else next.add(hit.index);
            setSelection([...next]);
            if (next.size === 1) activateTab('edit', false);
            return;
        }
        if (!hit) return;
        const need = { distance: 2, angle: 3, dihedral: 4 }[state.mode];
        state.pending.push(hit);
        renderHighlights();
        if (state.pending.length === need) {
            const pts = state.pending.map((p) => p.position);
            const s = state.structure;
            const names = state.pending.map((p) => `${s.symbols[p.index]}${p.index + 1}`).join('–');
            let text;
            if (need === 2) text = `${fmt(norm(sub(pts[1], pts[0])))} Å`;
            else if (need === 3) {
                const u = sub(pts[0], pts[1]), v = sub(pts[2], pts[1]);
                text = `${fmt(Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (norm(u) * norm(v))))) * 180 / Math.PI, 2)}°`;
            } else {
                const b1 = sub(pts[1], pts[0]), b2 = sub(pts[2], pts[1]), b3 = sub(pts[3], pts[2]);
                const n1 = cross(b1, b2), n2 = cross(b2, b3);
                const m1 = cross(n1, unit(b2));
                text = `${fmt(Math.atan2(dot(m1, n2), dot(n1, n2)) * 180 / Math.PI, 2)}°`;
            }
            state.measurements.push({ points: pts, text, label: names });
            state.pending = [];
            renderOverlays(); renderHighlights(); refreshMeasurements();
            toast(`${names}: ${text}`);
        }
    }
}


export function setMode(mode) {
    state.mode = mode;
    state.pending = [];
    document.querySelectorAll('#modebar [data-mode]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.mode === mode)));
    const hints = {
        select: 'Click to select · Ctrl-click to add · Shift-drag for box selection · Delete removes the selection',
        distance: 'Click two atoms to measure their distance',
        angle: 'Click three atoms: the angle is measured at the second',
        dihedral: 'Click four atoms to measure the dihedral angle',
    };
    $('st-hint').textContent = hints[mode];
    renderHighlights();
}

