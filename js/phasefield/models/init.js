// Initial-condition helpers shared by the models (run on the CPU once per reset).

export function rng(seed) {
    let a = (seed >>> 0) || 1;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

export const index = (n, x, y, z) => (z * n + y) * n + x;

// mean + amp·(uniform in [-1, 1]) everywhere.
export function noise(n, mean, amp, rand, comps = 1) {
    const out = new Float32Array(comps * n ** 3);
    for (let i = 0; i < out.length; i++) out[i] = mean + amp * (2 * rand() - 1);
    return out;
}

// Periodic Voronoi labels for `count` random seeds. Seeds are bucketed on a
// coarse grid so each cell only checks nearby seeds.
export function voronoi(n, count, rand) {
    const seeds = Array.from({ length: count }, () => [rand() * n, rand() * n, rand() * n]);
    const B = Math.max(1, Math.round(Math.cbrt(count / 2)));
    const bs = n / B;
    const buckets = Array.from({ length: B ** 3 }, () => []);
    seeds.forEach((s, k) => {
        const b = s.map((v) => Math.min(B - 1, Math.floor(v / bs)));
        buckets[(b[2] * B + b[1]) * B + b[0]].push(k);
    });
    const labels = new Int32Array(n ** 3);
    const reach = B >= 3 ? 1 : Math.ceil(B / 2);
    const md = (d) => d - n * Math.round(d / n);
    for (let z = 0; z < n; z++) for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
        const bx = Math.floor((x + 0.5) / bs), by = Math.floor((y + 0.5) / bs), bz = Math.floor((z + 0.5) / bs);
        let best = 0, bd = Infinity;
        for (let dz = -reach; dz <= reach; dz++) for (let dy = -reach; dy <= reach; dy++) for (let dx = -reach; dx <= reach; dx++) {
            const list = buckets[((((bz + dz) % B) + B) % B * B + (((by + dy) % B) + B) % B) * B + (((bx + dx) % B) + B) % B];
            for (const k of list) {
                const s = seeds[k];
                const d = md(x + 0.5 - s[0]) ** 2 + md(y + 0.5 - s[1]) ** 2 + md(z + 0.5 - s[2]) ** 2;
                if (d < bd) { bd = d; best = k; }
            }
        }
        labels[index(n, x, y, z)] = best;
    }
    return { labels, seeds };
}

// Assign each label one of `colors` so that touching labels differ where
// possible (greedy colouring, largest degree first).
export function colorLabels(n, labels, count, colors) {
    const adj = Array.from({ length: count }, () => new Set());
    for (let z = 0; z < n; z++) for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
        const a = labels[index(n, x, y, z)];
        for (const [u, v, w] of [[(x + 1) % n, y, z], [x, (y + 1) % n, z], [x, y, (z + 1) % n]]) {
            const b = labels[index(n, u, v, w)];
            if (a !== b) { adj[a].add(b); adj[b].add(a); }
        }
    }
    const order = [...Array(count).keys()].sort((a, b) => adj[b].size - adj[a].size);
    const color = new Int32Array(count).fill(-1);
    for (const g of order) {
        const used = new Array(colors).fill(0);
        for (const h of adj[g]) if (color[h] >= 0) used[color[h]]++;
        let best = 0;
        for (let c = 1; c < colors; c++) if (used[c] < used[best]) best = c;
        color[g] = best;
    }
    return color;
}

// Random non-overlapping sphere centres (periodic), as many as fit up to `count`.
export function spheres(n, count, radius, rand, minGap = 2) {
    const out = [];
    const md = (d) => d - n * Math.round(d / n);
    for (let tries = 0; out.length < count && tries < count * 60; tries++) {
        const c = [rand() * n, rand() * n, rand() * n];
        if (out.every((o) => Math.hypot(md(c[0] - o[0]), md(c[1] - o[1]), md(c[2] - o[2])) > 2 * radius + minGap)) out.push(c);
    }
    return out;
}

// Calls fn(i, r²) for every cell within `radius` of centre c (periodic).
export function forSphere(n, c, radius, fn) {
    const r = Math.ceil(radius);
    for (let dz = -r; dz <= r; dz++) for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) {
        const x = Math.floor(c[0]) + dx, y = Math.floor(c[1]) + dy, z = Math.floor(c[2]) + dz;
        const d2 = (x + 0.5 - c[0]) ** 2 + (y + 0.5 - c[1]) ** 2 + (z + 0.5 - c[2]) ** 2;
        if (d2 > radius * radius) continue;
        fn(index(n, ((x % n) + n) % n, ((y % n) + n) % n, ((z % n) + n) % n), d2);
    }
}
