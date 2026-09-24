// Builder cards: each builder in the Build tab is a collapsible card with a
// live 3D preview that updates as its inputs change, and a Build button.
import { previewRenderer } from '../../render/preview.js';
import { withBusy } from '../dom.js';

const cards = new Map();
const PREVIEW_DELAY = 220;

// Preview atom budget; builders scale their preview down to stay under it.
export const PREVIEW_LIMIT = 5000;

// Factor (≤ 1) by which a linear size must shrink so `estimate` atoms fit the budget.
export function previewScale(estimate, limit = PREVIEW_LIMIT) {
    return estimate > limit ? Math.cbrt(limit / estimate) : 1;
}

/**
 * id: matches data-card="…" in the page.
 * preview(): returns { structure, caption, colorBy } for the current inputs.
 * build(): performs the build (runs behind the busy overlay).
 */
export function registerCard(id, { preview, build, busy = 'Building…' }) {
    const el = document.querySelector(`[data-card="${id}"]`);
    if (!el) throw new Error(`Missing builder card "${id}"`);
    const canvas = el.querySelector('.preview canvas');
    const caption = el.querySelector('.preview-cap');
    let timer = 0;

    const run = () => {
        if (!el.open || !canvas) return;
        try {
            const r = preview();
            previewRenderer().show(canvas, r.structure, { colorBy: r.colorBy || 'element' });
            caption.innerHTML = r.caption || '';
            caption.classList.remove('warn');
        } catch (e) {
            previewRenderer().show(canvas, null);
            caption.textContent = e.message;
            caption.classList.add('warn');
        }
    };
    const schedule = () => { clearTimeout(timer); timer = setTimeout(run, PREVIEW_DELAY); };

    el.addEventListener('toggle', () => {
        if (!el.open) return;
        // Accordion: one builder open at a time.
        for (const other of cards.values()) if (other.el !== el) other.el.open = false;
        document.querySelectorAll('.builder-tiles [data-open]').forEach((t) => t.setAttribute('aria-pressed', String(t.dataset.open === id)));
        run();
    });
    el.addEventListener('input', schedule);
    el.addEventListener('change', schedule);
    el.querySelector('[data-build]')?.addEventListener('click', () => withBusy(busy, build));
    cards.set(id, { el, run, schedule });
}

// Re-run the preview of the open card (e.g. after the source crystal changed).
export function refreshOpenPreview() {
    for (const c of cards.values()) if (c.el.open) c.schedule();
}

export function openCard(id) {
    const c = cards.get(id);
    if (!c) return;
    c.el.open = true;
    c.el.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
}

// Nudge a card's preview after programmatic input changes.
export function schedulePreview(id) {
    cards.get(id)?.schedule();
}
