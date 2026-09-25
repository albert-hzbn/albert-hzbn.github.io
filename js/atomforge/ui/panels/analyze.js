// Analyze panel: radial distribution function and measurements.
import { radialDistribution } from '../../analysis/rdf.js';
import { state } from '../context.js';
import { $, num, int, fmt, withBusy } from '../dom.js';
import { renderOverlays, renderHighlights } from '../scene.js';

export function initAnalyze() {
    $('rdf-go').addEventListener('click', () => withBusy('Computing g(r)…', () => {
        const s = state.structure;
        if (s.count < 2) throw new Error('Need at least two atoms.');
        const res = radialDistribution(s, { rmax: num('rdf-r'), bins: int('rdf-bins'), a: $('rdf-a').value || null, b: $('rdf-b').value || null });
        drawPlot($('rdf-plot'), res.r, res.g, { xlabel: 'r (Å)', ylabel: 'g(r)' });
        $('rdf-summary').innerHTML = [
            ['First peak', res.firstPeak ? `${fmt(res.firstPeak)} Å` : '—'],
            ['First minimum', res.firstMin ? `${fmt(res.firstMin)} Å` : '—'],
            ['First-shell CN', res.firstShellCN ? fmt(res.firstShellCN, 2) : '—'],
            ['Boundary', res.periodic ? 'periodic' : 'open (cluster)'],
        ].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    }));
    $('measure-clear').addEventListener('click', () => { state.measurements = []; state.pending = []; refreshMeasurements(); renderOverlays(); renderHighlights(); });
}


export function drawPlot(canvas, xs, ys, { xlabel, ylabel }) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    const css = getComputedStyle(document.documentElement);
    const ink = css.getPropertyValue('--ink-2').trim(), line = css.getPropertyValue('--line').trim(), acc = css.getPropertyValue('--accent-2').trim();
    ctx.clearRect(0, 0, W, H);
    const m = { l: 58, r: 16, t: 16, b: 50 };
    const xmax = xs[xs.length - 1], ymax = Math.max(1.2, ...ys) * 1.08;
    const X = (x) => m.l + (x / xmax) * (W - m.l - m.r);
    const Y = (y) => H - m.b - (y / ymax) * (H - m.t - m.b);
    ctx.font = '22px Inter, sans-serif';
    ctx.fillStyle = ink; ctx.strokeStyle = line; ctx.lineWidth = 1.5;
    // grid + ticks
    for (let i = 0; i <= 4; i++) {
        const y = (ymax / 4) * i;
        ctx.beginPath(); ctx.moveTo(m.l, Y(y)); ctx.lineTo(W - m.r, Y(y)); ctx.stroke();
        ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(y.toFixed(y < 10 ? 1 : 0), m.l - 8, Y(y));
    }
    const step = xmax > 10 ? 2 : 1;
    for (let x = 0; x <= xmax + 1e-9; x += step) {
        ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(x.toString(), X(x), H - m.b + 8);
    }
    ctx.textAlign = 'center'; ctx.fillText(xlabel, (m.l + W - m.r) / 2, H - 24);
    ctx.save(); ctx.translate(16, (m.t + H - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(ylabel, 0, 0); ctx.restore();
    // g = 1 reference
    ctx.setLineDash([8, 6]); ctx.beginPath(); ctx.moveTo(m.l, Y(1)); ctx.lineTo(W - m.r, Y(1)); ctx.stroke(); ctx.setLineDash([]);
    ctx.strokeStyle = acc; ctx.lineWidth = 3; ctx.beginPath();
    xs.forEach((x, i) => { const px = X(x), py = Y(ys[i]); if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py); });
    ctx.stroke();
}


export function refreshMeasurements() {
    $('measure-list').innerHTML = state.measurements.map((m, k) => `<div><span>${k + 1}. ${m.label}</span><b>${m.text}</b></div>`).join('') || '<p class="note">No measurements yet.</p>';
}

