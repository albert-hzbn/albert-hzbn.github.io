// DOM helpers: element lookup, number parsing, toasts, busy overlay, downloads.
export const $ = (id) => document.getElementById(id);
export const num = (id) => parseFloat($(id).value);
export const int = (id) => parseInt($(id).value, 10);
export const fmt = (v, d = 3) => (Math.abs(v) < 1e-9 ? 0 : v).toFixed(d);

let toastTimer;
export function toast(msg, isError = false) {
    const t = $('toast');
    t.textContent = msg;
    t.classList.toggle('error', isError);
    t.classList.add('show');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => t.classList.remove('show'), isError ? 5200 : 2600);
}


export function withBusy(label, fn) {
    $('busy-text').textContent = label;
    $('busy').classList.add('on');
    return new Promise((resolve) => {
        requestAnimationFrame(() => setTimeout(() => {
            try { resolve(fn()); } catch (e) { console.error(e); toast(e.message, true); resolve(null); }
            finally { $('busy').classList.remove('on'); }
        }, 20));
    });
}


export function download(name, content, type = 'text/plain') {
    const blob = content instanceof Blob ? content : new Blob([content], { type });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
}


export function segmented(el, onChange) {
    el.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
        el.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', String(x === b)));
        onChange(b.dataset.v, b);
    }));
}


// A table row of inputs with a remove button. Removing it fires a bubbling
// 'change' on the table so live previews update.
export function removableRow(innerHTML, removeLabel = 'Remove row') {
    const tr = document.createElement('div');
    tr.className = 'tr';
    tr.innerHTML = `${innerHTML}<button type="button" class="x-btn" aria-label="${removeLabel}" title="${removeLabel}">×</button>`;
    tr.querySelector('.x-btn').addEventListener('click', () => {
        const table = tr.parentElement;
        tr.remove();
        table?.dispatchEvent(new Event('change', { bubbles: true }));
    });
    return tr;
}
