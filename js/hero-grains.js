// Live grain coarsening for the homepage header.
// A zero-temperature Monte Carlo Potts model on a square lattice with Moore
// neighbourhood, started from a Voronoi-like tessellation. Work per frame is
// capped by a time budget so the page stays smooth on slow devices, and the
// simulation restarts once the grains have coarsened.
// Clicking nucleates a recrystallised grain: flips from deformed grains into
// recrystallised ones gain a stored-energy term, so the nucleus grows.
(function () {
    const canvas = document.getElementById('grain-canvas');
    if (!canvas) return;

    const CELL = 2;               // screen pixels per lattice site
    const SEEDS_PER_SITE = 1 / 900;
    const RESTART_AFTER = 2400;   // frames
    const STORED_ENERGY = 2.5;    // driving force for recrystallised grains (in bond units)
    const NUCLEUS_RADIUS = 6;     // lattice sites
    const FRAME_BUDGET_MS = 5;    // Monte Carlo time per animation frame
    const CHUNK = 4000;           // flip attempts between budget checks
    const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

    const ctx = canvas.getContext('2d');
    const off = document.createElement('canvas');
    const octx = off.getContext('2d');

    let nx, ny, grid, image, palette, rx, frame, width = 0;
    let running = true, visible = true, ready = false;

    function isDark() {
        return getComputedStyle(document.documentElement).getPropertyValue('--grain-dark').trim() === '1';
    }

    function hslToRgb(h, s, l) {
        const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        const p = 2 * l - q;
        const f = (t) => {
            if (t < 0) t += 1;
            if (t > 1) t -= 1;
            if (t < 1 / 6) return p + (q - p) * 6 * t;
            if (t < 1 / 2) return q;
            if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
            return p;
        };
        return [Math.round(f(h + 1 / 3) * 255), Math.round(f(h) * 255), Math.round(f(h - 1 / 3) * 255)];
    }

    function makePalette(n, dark) {
        // Soft, low-contrast fills so text over the canvas stays readable.
        const p = [];
        for (let i = 0; i < n; i++) {
            const hue = 185 + Math.random() * 45;
            const sat = 18 + Math.random() * 22;
            const light = dark ? 13 + Math.random() * 9 : 88 + Math.random() * 8;
            p.push(hslToRgb(hue / 360, sat / 100, light / 100));
        }
        return p;
    }

    const DI = [-1, 0, 1, -1, 1, -1, 0, 1];
    const DJ = [-1, -1, -1, 0, 0, 1, 1, 1];

    function wrapX(i) { return i < 0 ? i + nx : (i >= nx ? i - nx : i); }
    function wrapY(j) { return j < 0 ? j + ny : (j >= ny ? j - ny : j); }

    function init() {
        const rect = canvas.getBoundingClientRect();
        const dpr = Math.min(window.devicePixelRatio || 1, 2);
        width = rect.width;
        canvas.width = Math.round(rect.width * dpr);
        canvas.height = Math.round(rect.height * dpr);
        nx = Math.max(40, Math.ceil(rect.width / CELL));
        ny = Math.max(20, Math.ceil(rect.height / CELL));
        off.width = nx;
        off.height = ny;
        image = octx.createImageData(nx, ny);
        grid = new Int32Array(nx * ny).fill(-1);

        // Tessellation by simultaneous flood fill from random seeds: linear in
        // the number of sites, and close enough to Voronoi once Potts smooths it.
        const nSeeds = Math.max(30, Math.round(nx * ny * SEEDS_PER_SITE));
        const queue = new Int32Array(nx * ny);
        let head = 0, tail = 0;
        for (let k = 0; k < nSeeds; k++) {
            const idx = ((Math.random() * ny) | 0) * nx + ((Math.random() * nx) | 0);
            if (grid[idx] === -1) { grid[idx] = k; queue[tail++] = idx; }
        }
        while (head < tail) {
            const idx = queue[head++];
            const i = idx % nx, j = (idx / nx) | 0, g = grid[idx];
            for (let n = 0; n < 8; n++) {
                const t = wrapY(j + DJ[n]) * nx + wrapX(i + DI[n]);
                if (grid[t] === -1) { grid[t] = g; queue[tail++] = t; }
            }
        }
        palette = makePalette(nSeeds, isDark());
        rx = new Uint8Array(nSeeds);
        frame = 0;
        ready = true;
    }

    function energy(i, j, s) {
        let e = 0;
        for (let n = 0; n < 8; n++) {
            if (grid[wrapY(j + DJ[n]) * nx + wrapX(i + DI[n])] !== s) e++;
        }
        return e;
    }

    function attempt(count) {
        for (let a = 0; a < count; a++) {
            const i = (Math.random() * nx) | 0, j = (Math.random() * ny) | 0;
            const idx = j * nx + i;
            const n = (Math.random() * 8) | 0;
            const s = grid[wrapY(j + DJ[n]) * nx + wrapX(i + DI[n])], cur = grid[idx];
            if (s === cur) continue;
            const drive = (rx[s] && !rx[cur]) ? STORED_ENERGY : 0;
            if (energy(i, j, s) - drive <= energy(i, j, cur)) grid[idx] = s;
        }
    }

    function step() {
        const t0 = performance.now();
        do { attempt(CHUNK); } while (performance.now() - t0 < FRAME_BUDGET_MS);
    }

    function draw() {
        const d = image.data;
        const edge = isDark() ? [70, 110, 125] : [120, 150, 160];
        for (let j = 0; j < ny; j++) {
            const row = j * nx, down = wrapY(j + 1) * nx;
            for (let i = 0; i < nx; i++) {
                const idx = row + i, s = grid[idx];
                const c = (s !== grid[row + wrapX(i + 1)] || s !== grid[down + i]) ? edge : palette[s];
                const o = idx * 4;
                d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = 255;
            }
        }
        octx.putImageData(image, 0, 0);
        ctx.imageSmoothingEnabled = true;
        ctx.drawImage(off, 0, 0, canvas.width, canvas.height);
    }

    function loop() {
        if (ready && running && visible) {
            step();
            draw();
            frame++;
            if (frame > RESTART_AFTER) init();
        }
        requestAnimationFrame(loop);
    }

    function settleStatic() {
        // Reduced motion: show a coarsened snapshot instead of an animation.
        for (let k = 0; k < 60; k++) attempt(nx * ny);
        draw();
    }

    function nucleate(clientX, clientY) {
        if (!ready) return;
        const rect = canvas.getBoundingClientRect();
        const ci = Math.floor((clientX - rect.left) / rect.width * nx);
        const cj = Math.floor((clientY - rect.top) / rect.height * ny);
        const id = palette.length;
        // Warm tint so recrystallised grains read as a different population.
        palette.push(hslToRgb((22 + Math.random() * 18) / 360, 0.45, isDark() ? 0.24 : 0.86));
        const grown = new Uint8Array(id + 1);
        grown.set(rx);
        grown[id] = 1;
        rx = grown;
        const r = NUCLEUS_RADIUS;
        for (let dj = -r; dj <= r; dj++) {
            for (let di = -r; di <= r; di++) {
                if (di * di + dj * dj <= r * r) grid[wrapY(cj + dj) * nx + wrapX(ci + di)] = id;
            }
        }
        frame = Math.min(frame, RESTART_AFTER - 900); // give the nucleus time to grow
        if (reduceMotion) { for (let k = 0; k < 30; k++) attempt(nx * ny); }
        draw();
    }

    canvas.addEventListener('pointerdown', (e) => nucleate(e.clientX, e.clientY));

    // Build the first microstructure after the page has painted, so the
    // entrance animation is not held up.
    requestAnimationFrame(() => setTimeout(() => {
        init();
        if (reduceMotion) settleStatic();
        else draw();
    }, 0));

    if (!reduceMotion) {
        requestAnimationFrame(loop);
        if ('IntersectionObserver' in window) {
            new IntersectionObserver((entries) => { visible = entries[0].isIntersecting; }).observe(canvas);
        }
        document.addEventListener('visibilitychange', () => { running = !document.hidden; });
    }

    // Only rebuild when the width changes; mobile URL-bar resizes change height only.
    let resizeTimer;
    window.addEventListener('resize', () => {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(() => {
            if (Math.abs(canvas.getBoundingClientRect().width - width) < 2) return;
            init();
            if (reduceMotion) settleStatic(); else draw();
        }, 200);
    });
})();
