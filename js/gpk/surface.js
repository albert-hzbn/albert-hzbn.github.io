// GPK Web: smoothed surface (naive surface nets) around the occupied cells
// of an n^3 volume, for the 3D view.

// Surface nets on a binary occupancy n^3 grid; cells outside the box are empty.
// colorOf(i, out) writes the RGB colour of occupied cell i into out[0..2].
// Each face takes the colour of the cell it bounds, so label boundaries stay
// sharp; vertex positions are smoothed. Returns indexed positions and one
// colour per triangle.
export function surfaceNets(n, occ, colorOf) {
    const V = n + 1, vid = new Int32Array(V * V * V).fill(-1);
    const pos = [], lo = [];
    const cell = (x, y, z) => (x < 0 || y < 0 || z < 0 || x >= n || y >= n || z >= n ? 0 : occ[(z * n + y) * n + x]);
    const corner = new Uint8Array(8);
    for (let vz = 0; vz <= n; vz++) for (let vy = 0; vy <= n; vy++) for (let vx = 0; vx <= n; vx++) {
        let k = 0, sum = 0;
        for (let c = 0; c < 2; c++) for (let b = 0; b < 2; b++) for (let a = 0; a < 2; a++) {
            const o = cell(vx - 1 + a, vy - 1 + b, vz - 1 + c);
            corner[k++] = o; sum += o;
        }
        if (sum === 0 || sum === 8) continue;
        // Average of the midpoints of the cube edges that cross the surface.
        let px = 0, py = 0, pz = 0, m = 0;
        for (let k1 = 0; k1 < 8; k1++) for (const bit of [1, 2, 4]) {
            const k2 = k1 | bit;
            if (k2 === k1 || corner[k1] === corner[k2]) continue;
            px += ((k1 & 1) + (k2 & 1)) / 2; py += (((k1 >> 1) & 1) + ((k2 >> 1) & 1)) / 2; pz += (((k1 >> 2) & 1) + ((k2 >> 2) & 1)) / 2;
            m++;
        }
        // Corner (vx, vy, vz) of the cell grid; cube corners are cell centres at vx - 0.5 + a.
        pos.push(vx - 0.5 + px / m, vy - 0.5 + py / m, vz - 0.5 + pz / m);
        lo.push(vx - 0.5, vy - 0.5, vz - 0.5);
        vid[(vz * V + vy) * V + vx] = pos.length / 3 - 1;
    }
    const index = [], tri = [], c3 = [0, 0, 0];
    // own: the occupied cell this face bounds; its colour fills both triangles.
    const quad = (a, b, c, d, flip, own) => {
        if (a < 0 || b < 0 || c < 0 || d < 0) return;
        if (flip) index.push(a, c, b, a, d, c); else index.push(a, b, c, a, c, d);
        colorOf(own, c3);
        tri.push(c3[0], c3[1], c3[2], c3[0], c3[1], c3[2]);
    };
    const v = (x, y, z) => vid[(z * V + y) * V + x];
    const at = (x, y, z) => (z * n + y) * n + x;
    for (let z = 0; z < n; z++) for (let y = 0; y < n; y++) for (let x = -1; x < n; x++) {
        const a = cell(x, y, z), b = cell(x + 1, y, z);
        if (a !== b) quad(v(x + 1, y, z), v(x + 1, y + 1, z), v(x + 1, y + 1, z + 1), v(x + 1, y, z + 1), !a, a ? at(x, y, z) : at(x + 1, y, z));
    }
    for (let z = 0; z < n; z++) for (let y = -1; y < n; y++) for (let x = 0; x < n; x++) {
        const a = cell(x, y, z), b = cell(x, y + 1, z);
        if (a !== b) quad(v(x, y + 1, z), v(x, y + 1, z + 1), v(x + 1, y + 1, z + 1), v(x + 1, y + 1, z), !a, a ? at(x, y, z) : at(x, y + 1, z));
    }
    for (let z = -1; z < n; z++) for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
        const a = cell(x, y, z), b = cell(x, y, z + 1);
        if (a !== b) quad(v(x, y, z + 1), v(x + 1, y, z + 1), v(x + 1, y + 1, z + 1), v(x, y + 1, z + 1), !a, a ? at(x, y, z) : at(x, y, z + 1));
    }

    // Two passes of Laplacian smoothing, each vertex kept inside its own cube.
    const nv = pos.length / 3, P = Float32Array.from(pos), Q = new Float32Array(P.length);
    const nbr = Array.from({ length: nv }, () => []);
    for (let t = 0; t < index.length; t += 3) {
        const [a, b, c] = [index[t], index[t + 1], index[t + 2]];
        nbr[a].push(b, c); nbr[b].push(a, c); nbr[c].push(a, b);
    }
    for (let it = 0; it < 2; it++) {
        for (let i = 0; i < nv; i++) {
            const list = nbr[i];
            for (let d = 0; d < 3; d++) {
                let s = 0;
                for (const j of list) s += P[3 * j + d];
                const target = list.length ? s / list.length : P[3 * i + d];
                Q[3 * i + d] = Math.min(lo[3 * i + d] + 1, Math.max(lo[3 * i + d], 0.5 * P[3 * i + d] + 0.5 * target));
            }
        }
        P.set(Q);
    }
    return { position: P, triColor: Float32Array.from(tri), index: nv > 65535 ? Uint32Array.from(index) : Uint16Array.from(index) };
}
