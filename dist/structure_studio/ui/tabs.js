// Tool-panel tabs.
export function activateTab(name, openPanel = true) {
    document.querySelectorAll('.tabs [data-tab]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.tab === name)));
    document.querySelectorAll('[data-body]').forEach((b) => { b.hidden = b.dataset.body !== name; });
    if (openPanel && window.innerWidth <= 860) document.body.classList.add('panel-open');
}

