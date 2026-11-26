// GPK Web: three-dimensional port of the Gaussian-projected kinetics core
// (https://github.com/albert-hzbn/GPK, src/gpk_*.cpp with Nz > 1).
//
// The step sequence is the one in engine.js; the 3-D differences follow the
// C++ code: 6-connected topology certificates on the 3x3x3 neighbourhood,
// velocity and remap Gaussians from gpk_grid.cpp, two balanced remap sweeps,
// a spherical extinction threshold and a six-neighbour diffusion limit.
//
// The Gaussian-Hessian geometry is the same tensor-product stencil as the
// C++ local path, evaluated separably: for each resident label the x sums,
// then the y sums, are cached per cell and reused by every interface cell
// whose window covers them. The result is identical to the direct sum.
//
// Thermo callbacks receive the z components as trailing arguments:
//   sigma(nx, ny, g, p, q, c, gq, nz), mobility(gp, gq, T, nx, ny, nz)

const mod = (a, n) => ((a % n) + n) % n;
const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

function gaussWeights(sigma) {
    const s = Math.max(sigma, 0.5), r = Math.ceil(3.5 * s), w = new Float64Array(2 * r + 1);
    let sum = 0;
    for (let i = 0; i <= 2 * r; i++) { w[i] = Math.exp(-((i - r) ** 2) / (2 * s * s)); sum += w[i]; }
    for (let i = 0; i <= 2 * r; i++) w[i] /= sum;
    return { sigma: s, radius: r, w };
}

// 1-D factor of the cubic window exp(-|d|^2 / 2 sigma^2), |d|_inf <= ceil(3 sigma).
function window3(sigma) {
    const r = Math.max(1, Math.ceil(3 * sigma)), inv = 0.5 / (sigma * sigma);
    const w1 = new Float64Array(2 * r + 1);
    for (let i = -r; i <= r; i++) w1[i + r] = Math.exp(-i * i * inv);
    return { r, w1 };
}

// 6-adjacency inside a 3x3x3 (side 3) or 5x5x5 (side 5) box.
function boxAdjacency(side) {
    const n = side ** 3, adj = [];
    for (let s = 0; s < n; s++) {
        const z = (s / (side * side)) | 0, y = ((s / side) | 0) % side, x = s % side, list = [];
        if (x > 0) list.push(s - 1); if (x < side - 1) list.push(s + 1);
        if (y > 0) list.push(s - side); if (y < side - 1) list.push(s + side);
        if (z > 0) list.push(s - side * side); if (z < side - 1) list.push(s + side * side);
        adj.push(list);
    }
    return adj;
}
const ADJ3 = boxAdjacency(3), ADJ5 = boxAdjacency(5);

export class GPK3D {
    constructor(Nx, Ny, Nz, dx, thermo, D = 0, sigmaCells = 2) {
        this.Nx = Nx; this.Ny = Ny; this.Nz = Nz; this.N = Nx * Ny * Nz;
        this.dx = dx; this.th = thermo; this.D = D;
        this.time = 0; this.stepId = 0;
        const N = this.N;
        this.phi = new Int32Array(N);
        this.alpha = new Float32Array(N).fill(1);
        this.c = new Float32Array(N);
        this.T = new Float32Array(N).fill(800);
        this.g = new Float32Array(N);

        this.nxF = new Float32Array(N);
        this.nyF = new Float32Array(N);
        this.nzF = new Float32Array(N);
        this.kappa = new Float32Array(N);
        this.vn = new Float32Array(N);
        this.phiComp = new Int32Array(N);
        this.gComp = new Float32Array(N);
        this.mask = new Uint8Array(N);
        this.active = new Int32Array(N);
        this.nActive = 0;

        this.tmp = new Float32Array(N);
        this.cNew = new Float32Array(N);
        this.resScore = new Float32Array(N);
        this.bestScore = new Float32Array(N);
        this.solidOld = new Float32Array(N);
        this.solidNew = new Float32Array(N);
        this.visit = new Int32Array(N);
        this.visitStamp = 0;
        this.stack = new Int32Array(N);
        this.order = new Int32Array(N);

        // Separable geometry caches (x sums and y sums for the current label).
        this.X0 = new Float32Array(N); this.X1 = new Float32Array(N); this.X2 = new Float32Array(N);
        this.Y00 = new Float32Array(N); this.Y10 = new Float32Array(N); this.Y01 = new Float32Array(N);
        this.Y20 = new Float32Array(N); this.Y02 = new Float32Array(N); this.Y11 = new Float32Array(N);
        this.stX = new Int32Array(N); this.stY = new Int32Array(N);
        this.pass = 0;

        this.gauss = gaussWeights(sigmaCells);
        this.velocitySigma = clamp(0.375 * this.gauss.sigma, 0.85, 3.0);
        this.remapSigma = clamp(0.125 * this.gauss.sigma, 0.70, 1.0);
        this.buildGeometryKernels();
        this.velW = window3(this.velocitySigma);
        this.remapW = window3(this.remapSigma);

        const W = this.W = Math.max(this.geomR, this.velW.r, this.remapW.r, 3);
        const wrap = (n) => { const t = new Int32Array(n + 2 * W); for (let j = 0; j < t.length; j++) t[j] = mod(j - W, n); return t; };
        this.xw = wrap(Nx); this.yw = wrap(Ny); this.zw = wrap(Nz);

        this.counts = new Int32Array(16);
        this.stats = { accepted: 0, deferred: 0, remapped: 0 };
    }

    buildGeometryKernels() {
        const G = this.gauss, R = this.geomR = G.radius + 2, width = 2 * R + 1;
        this.w0 = new Float64Array(width);
        this.w1 = new Float64Array(width);
        this.w2 = new Float64Array(width);
        const g = (o) => (o < -G.radius || o > G.radius ? 0 : G.w[o + G.radius]);
        const o1 = [-2, -1, 1, 2], c1 = [1 / 12, -8 / 12, 8 / 12, -1 / 12];
        const o2 = [-2, -1, 0, 1, 2], c2 = [-1 / 12, 16 / 12, -30 / 12, 16 / 12, -1 / 12];
        for (let o = -R; o <= R; o++) {
            const j = o + R;
            this.w0[j] = g(o);
            for (let k = 0; k < 4; k++) this.w1[j] += c1[k] * g(o - o1[k]);
            for (let k = 0; k < 5; k++) this.w2[j] += c2[k] * g(o - o2[k]);
        }
    }

    idx(x, y, z = 0) { return (mod(z, this.Nz) * this.Ny + mod(y, this.Ny)) * this.Nx + mod(x, this.Nx); }

    // Six face neighbours of cell i: -x, +x, -y, +y, -z, +z.
    faces(i, out) {
        const { Nx, Ny, Nz } = this, x = i % Nx, t = (i / Nx) | 0, y = t % Ny, z = (t / Ny) | 0;
        const plane = Nx * Ny, row = i - x;
        out[0] = row + (x === 0 ? Nx - 1 : x - 1);
        out[1] = row + (x === Nx - 1 ? 0 : x + 1);
        out[2] = i + ((y === 0 ? Ny - 1 : y - 1) - y) * Nx;
        out[3] = i + ((y === Ny - 1 ? 0 : y + 1) - y) * Nx;
        out[4] = i + ((z === 0 ? Nz - 1 : z - 1) - z) * plane;
        out[5] = i + ((z === Nz - 1 ? 0 : z + 1) - z) * plane;
        return out;
    }

    maxLabel() {
        let m = 0;
        for (let i = 0; i < this.N; i++) if (this.phi[i] > m) m = this.phi[i];
        return m;
    }

    countLabels() {
        const m = this.maxLabel();
        if (this.counts.length <= m + 1) this.counts = new Int32Array(Math.max(m + 2, this.counts.length * 2));
        else this.counts.fill(0);
        for (let i = 0; i < this.N; i++) this.counts[this.phi[i]]++;
        return this.counts;
    }

    // ------------------------------------------------------------------
    // 1. Interface mask
    // ------------------------------------------------------------------
    stepInterfaceMask() {
        const { Nx, Ny, Nz, phi, mask, active } = this, plane = Nx * Ny;
        let n = 0;
        for (let z = 0; z < Nz; z++) {
            const zm = (z === 0 ? Nz - 1 : z - 1) * plane, zp = (z === Nz - 1 ? 0 : z + 1) * plane, zo = z * plane;
            for (let y = 0; y < Ny; y++) {
                const row = zo + y * Nx;
                const ym = (y === 0 ? Ny - 1 : y - 1) * Nx, yp = (y === Ny - 1 ? 0 : y + 1) * Nx, yo = y * Nx;
                for (let x = 0; x < Nx; x++) {
                    const i = row + x, v = phi[i];
                    const xm = x === 0 ? Nx - 1 : x - 1, xp = x === Nx - 1 ? 0 : x + 1;
                    const a = v !== phi[row + xm] || v !== phi[row + xp] || v !== phi[zo + ym + x] || v !== phi[zo + yp + x]
                        || v !== phi[zm + yo + x] || v !== phi[zp + yo + x];
                    mask[i] = a ? 1 : 0;
                    if (a) active[n++] = i;
                }
            }
        }
        this.nActive = n;
    }

    // ------------------------------------------------------------------
    // 2. Gaussian-Hessian geometry (separable evaluation per label)
    // ------------------------------------------------------------------
    stepGeometry() {
        const { Nx, Ny, phi, xw, yw, zw, W, w0, w1, w2, active, nActive, nxF, nyF, nzF, kappa, phiComp, order } = this;
        const { X0, X1, X2, Y00, Y10, Y01, Y20, Y02, Y11, stX, stY } = this;
        const R = this.geomR, dx = this.dx, invdx = 1 / dx, invdx2 = invdx * invdx;
        const klim = 1 / (this.gauss.sigma * dx);
        const plane = Nx * Ny;
        nxF.fill(0); nyF.fill(0); nzF.fill(0); kappa.fill(0);

        // Group interface cells by resident label (counting sort, stable in index).
        const counts = this.countLabels(), M = counts.length;
        const start = new Int32Array(M + 1);
        for (let a = 0; a < nActive; a++) start[phi[active[a]] + 1]++;
        for (let l = 0; l < M; l++) start[l + 1] += start[l];
        const fill = start.slice(0, M);
        for (let a = 0; a < nActive; a++) { const i = active[a]; order[fill[phi[i]]++] = i; }

        for (let p = 0; p < M; p++) {
            const a0 = start[p], a1 = start[p + 1];
            if (a0 === a1) continue;
            const pass = ++this.pass;
            for (let a = a0; a < a1; a++) {
                const i = order[a];
                const x0 = i % Nx, t = (i / Nx) | 0, y0 = t % Ny, z0 = (t / Ny) | 0;
                let Ix = 0, Iy = 0, Iz = 0, Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
                for (let dz = -R; dz <= R; dz++) {
                    const zz = zw[z0 + dz + W], jz = dz + R, zb = zz * plane;
                    const j = zb + y0 * Nx + x0;
                    if (stY[j] !== pass) {
                        // y sums of the x sums along this column.
                        let s00 = 0, s10 = 0, s01 = 0, s20 = 0, s02 = 0, s11 = 0;
                        for (let dy = -R; dy <= R; dy++) {
                            const rb = zb + yw[y0 + dy + W] * Nx, k = rb + x0, jy = dy + R;
                            if (stX[k] !== pass) {
                                let a0x = 0, a1x = 0, a2x = 0;
                                for (let ddx = -R; ddx <= R; ddx++) {
                                    if (phi[rb + xw[x0 + ddx + W]] !== p) continue;
                                    const jx = ddx + R;
                                    a0x += w0[jx]; a1x += w1[jx]; a2x += w2[jx];
                                }
                                X0[k] = a0x; X1[k] = a1x; X2[k] = a2x; stX[k] = pass;
                            }
                            const x0v = X0[k], x1v = X1[k], x2v = X2[k];
                            if (x0v === 0 && x1v === 0 && x2v === 0) continue;   // no cell of p in this row
                            const g0 = w0[jy], g1 = w1[jy], g2 = w2[jy];
                            s00 += x0v * g0; s10 += x1v * g0; s01 += x0v * g1;
                            s20 += x2v * g0; s02 += x0v * g2; s11 += x1v * g1;
                        }
                        Y00[j] = s00; Y10[j] = s10; Y01[j] = s01; Y20[j] = s20; Y02[j] = s02; Y11[j] = s11; stY[j] = pass;
                    }
                    const h0 = w0[jz], h1 = w1[jz], h2 = w2[jz];
                    Ix += Y10[j] * h0; Iy += Y01[j] * h0; Iz += Y00[j] * h1;
                    Ixx += Y20[j] * h0; Iyy += Y02[j] * h0; Izz += Y00[j] * h2;
                    Ixy += Y11[j] * h0; Ixz += Y10[j] * h1; Iyz += Y01[j] * h1;
                }
                Ix *= invdx; Iy *= invdx; Iz *= invdx;
                Ixx *= invdx2; Iyy *= invdx2; Izz *= invdx2; Ixy *= invdx2; Ixz *= invdx2; Iyz *= invdx2;
                const g2 = Ix * Ix + Iy * Iy + Iz * Iz;
                if (g2 * dx * dx > 1e-8) {
                    const gm = Math.sqrt(g2);
                    nxF[i] = -Ix / gm; nyF[i] = -Iy / gm; nzF[i] = -Iz / gm;
                    const gHg = Ix * Ix * Ixx + Iy * Iy * Iyy + Iz * Iz * Izz + 2 * (Ix * Iy * Ixy + Ix * Iz * Ixz + Iy * Iz * Iyz);
                    let k = (gHg - g2 * (Ixx + Iyy + Izz)) / (g2 * gm);
                    if (Math.abs(k * dx) < 1e-6) k = 0;
                    kappa[i] = clamp(k, -klim, klim);
                }
            }
        }

        // Competing label: best-supported unlike face neighbour, normal as tie-breaker.
        const nei = new Int32Array(6), score = new Float64Array(6);
        for (let a = 0; a < nActive; a++) {
            const i = active[a], x = i % Nx, t = (i / Nx) | 0, y = t % Ny, z = (t / Ny) | 0;
            this.faces(i, nei);
            score[0] = -nxF[i]; score[1] = nxF[i]; score[2] = -nyF[i]; score[3] = nyF[i]; score[4] = -nzF[i]; score[5] = nzF[i];
            let best = -1, bestSupport = -1, bestScore = -2;
            for (let k = 0; k < 6; k++) {
                const q = phi[nei[k]];
                if (q === phi[i]) continue;
                let support = 0;
                for (let dz = -1; dz <= 1; dz++) {
                    const zb = zw[z + dz + W] * plane;
                    for (let dy = -1; dy <= 1; dy++) {
                        const rb = zb + yw[y + dy + W] * Nx;
                        for (let ddx = -1; ddx <= 1; ddx++) {
                            if (!ddx && !dy && !dz) continue;
                            if (phi[rb + xw[x + ddx + W]] === q) support++;
                        }
                    }
                }
                if (support > bestSupport || (support === bestSupport && score[k] > bestScore)) {
                    best = k; bestSupport = support; bestScore = score[k];
                }
            }
            phiComp[i] = best >= 0 ? phi[nei[best]] : phi[i];
        }
    }

    // ------------------------------------------------------------------
    // 3. Constitutive interface velocity
    // ------------------------------------------------------------------
    stepVelocities() {
        const { phi, c, T, g, active, nActive, nxF, nyF, nzF, kappa, vn, phiComp, gComp, th } = this;
        vn.fill(0);
        const cand = new Int32Array(6), nei = new Int32Array(6);
        for (let a = 0; a < nActive; a++) {
            const i = active[a], p = phi[i];
            this.faces(i, nei);
            if (th.phaseSelection) {
                let n = 0;
                for (let k = 0; k < 6; k++) {
                    const q = phi[nei[k]];
                    if (q === p) continue;
                    let seen = false;
                    for (let m = 0; m < n; m++) if (cand[m] === q) { seen = true; break; }
                    if (!seen) cand[n++] = q;
                }
                const sel = th.selectProduct(p, c[i], T[i], cand, n);
                if (sel >= 0) phiComp[i] = sel;
            }
            const q = phiComp[i], gp = g[i];
            let gq = gp, cq = c[i];
            for (let k = 0; k < 6; k++) if (phi[nei[k]] === q) { gq = g[nei[k]]; cq = c[nei[k]]; break; }
            gComp[i] = gq;
            const ci = th.interfaceComp ? th.interfaceComp(p, q, c[i], cq) : c[i];
            const dG = th.bulkDG(p, q, ci, T[i]);
            const sig = th.sigma(nxF[i], nyF[i], gp, p, q, c[i], gq, nzF[i]);
            const mob = th.mobility(gp, gq, T[i], nxF[i], nyF[i], nzF[i]);
            vn[i] = mob * (dG - sig * kappa[i]);
        }
    }

    // ------------------------------------------------------------------
    // Separable masked Gaussian sums. For request r at cell i with key k:
    //   outS[r] = sum_d w(d) [key(i+d) = k],  outV[r] = sum_d w(d) [key(i+d) = k] val(i+d)
    // with w(d) = w1(dx) w1(dy) w1(dz) over |d|_inf <= rad. x and y partial
    // sums are cached per cell for the key being processed.
    // ------------------------------------------------------------------
    maskedSums(key, val, w1, rad, cells, keys, nReq, outS, outV) {
        const { Nx, Ny, xw, yw, zw, W, X0, X1, Y00, Y10, stX, stY, order } = this, plane = Nx * Ny;
        let K = 0;
        for (let r = 0; r < nReq; r++) if (keys[r] + 1 > K) K = keys[r] + 1;
        const start = new Int32Array(K + 1);
        for (let r = 0; r < nReq; r++) start[keys[r] + 1]++;
        for (let k = 0; k < K; k++) start[k + 1] += start[k];
        const fill = start.slice(0, K);
        for (let r = 0; r < nReq; r++) order[fill[keys[r]]++] = r;
        for (let k = 0; k < K; k++) {
            if (start[k] === start[k + 1]) continue;
            const pass = ++this.pass;
            for (let a = start[k]; a < start[k + 1]; a++) {
                const r = order[a], i = cells[r];
                const x0 = i % Nx, t = (i / Nx) | 0, y0 = t % Ny, z0 = (t / Ny) | 0;
                let s = 0, v = 0;
                for (let dz = -rad; dz <= rad; dz++) {
                    const zb = zw[z0 + dz + W] * plane, j = zb + y0 * Nx + x0;
                    if (stY[j] !== pass) {
                        let ys = 0, yv = 0;
                        for (let dy = -rad; dy <= rad; dy++) {
                            const rb = zb + yw[y0 + dy + W] * Nx, m = rb + x0;
                            if (stX[m] !== pass) {
                                let xs = 0, xv = 0;
                                for (let dx = -rad; dx <= rad; dx++) {
                                    const c = rb + xw[x0 + dx + W];
                                    if (key[c] !== k) continue;
                                    xs += w1[dx + rad];
                                    if (val) xv += w1[dx + rad] * val[c];
                                }
                                X0[m] = xs; X1[m] = xv; stX[m] = pass;
                            }
                            ys += w1[dy + rad] * X0[m]; yv += w1[dy + rad] * X1[m];
                        }
                        Y00[j] = ys; Y10[j] = yv; stY[j] = pass;
                    }
                    s += w1[dz + rad] * Y00[j]; v += w1[dz + rad] * Y10[j];
                }
                outS[r] = s;
                if (outV) outV[r] = v;
            }
        }
    }

    // ------------------------------------------------------------------
    // 4. Pair-wise Gaussian projection of v_n
    // ------------------------------------------------------------------
    stepVelocityProjection() {
        const { phi, phiComp, active, nActive, vn, tmp } = this;
        const pid = this.pairKey || (this.pairKey = new Int32Array(this.N));
        const sv = this.signedV || (this.signedV = new Float32Array(this.N));
        const cells = this.reqCells || (this.reqCells = new Int32Array(this.N));
        const keys = this.reqKeys || (this.reqKeys = new Int32Array(this.N));
        const outS = this.reqS || (this.reqS = new Float64Array(this.N));
        const outV = this.reqV || (this.reqV = new Float64Array(this.N));
        pid.fill(-1);
        // Integer id for each unordered label pair present at the interface.
        const ids = new Map();
        let nReq = 0;
        for (let a = 0; a < nActive; a++) {
            const i = active[a], p = phi[i], q = phiComp[i];
            if (p === q) continue;
            const low = p < q ? p : q, high = p < q ? q : p, key = low * 1048576 + high;
            let id = ids.get(key);
            if (id === undefined) ids.set(key, id = ids.size);
            pid[i] = id;
            sv[i] = (p === low ? 1 : -1) * vn[i];
            cells[nReq] = i; keys[nReq] = id; nReq++;
        }
        const { r, w1 } = this.velW;
        this.maskedSums(pid, sv, w1, r, cells, keys, nReq, outS, outV);
        tmp.set(vn);
        for (let k = 0; k < nReq; k++) {
            if (!(outS[k] > 0)) continue;
            const i = cells[k];
            tmp[i] = (phi[i] < phiComp[i] ? 1 : -1) * outV[k] / outS[k];
        }
        this.tmp = vn; this.vn = tmp;
    }

    maxDt() {
        let vmax = 1e-30;
        for (let a = 0; a < this.nActive; a++) {
            const v = Math.abs(this.vn[this.active[a]]);
            if (v > vmax) vmax = v;
        }
        const dtIface = 0.125 * this.dx / vmax;
        const dtDiff = this.D > 0 ? 0.4 * this.dx * this.dx / (6 * this.D) : 1e30;
        return Math.min(dtIface, dtDiff);
    }

    // ------------------------------------------------------------------
    // Digital-topology certificate for one label change (3-D)
    // ------------------------------------------------------------------
    topologySimpleMove(cell, newLabel, oldCount, newCount, localOnly = false) {
        const { Nx, Ny, phi, xw, yw, zw, W } = this, plane = Nx * Ny;
        const old = phi[cell];
        if (old === newLabel) return false;
        const x0 = cell % Nx, t = (cell / Nx) | 0, y0 = t % Ny, z0 = (t / Ny) | 0;
        const occ = this._occ || (this._occ = new Uint8Array(125));
        const seen = this._seen || (this._seen = new Uint8Array(125));
        const queue = this._queue || (this._queue = new Int16Array(125));

        // Face-connected components of `label` in the 3x3x3 box (centre removed)
        // that touch the centre through a face.
        const topologicalNumber = (label) => {
            for (let lz = 0; lz < 3; lz++) {
                const zb = zw[z0 + lz - 1 + W] * plane;
                for (let ly = 0; ly < 3; ly++) {
                    const rb = zb + yw[y0 + ly - 1 + W] * Nx;
                    for (let lx = 0; lx < 3; lx++) occ[lz * 9 + ly * 3 + lx] = phi[rb + xw[x0 + lx - 1 + W]] === label ? 1 : 0;
                }
            }
            occ[13] = 0;
            seen.fill(0, 0, 27);
            let components = 0;
            for (let s = 0; s < 27; s++) {
                if (!occ[s] || seen[s]) continue;
                let head = 0, tail = 0, touches = false;
                queue[tail++] = s; seen[s] = 1;
                while (head < tail) {
                    const cur = queue[head++];
                    if (cur === 4 || cur === 10 || cur === 12 || cur === 14 || cur === 16 || cur === 22) touches = true;
                    for (const nx of ADJ3[cur]) if (occ[nx] && !seen[nx]) { seen[nx] = 1; queue[tail++] = nx; }
                }
                if (touches) components++;
            }
            return components;
        };

        const nb = new Int32Array(6);
        const globalBranchesConnected = (label, excludeCentre) => {
            this.faces(cell, nb);
            const branches = [];
            for (let k = 0; k < 6; k++) if (phi[nb[k]] === label) branches.push(nb[k]);
            if (branches.length <= 1) return true;
            const stamp = ++this.visitStamp, visit = this.visit, stack = this.stack;
            const targets = new Set(branches.slice(1));
            let top = 0;
            stack[top++] = branches[0]; visit[branches[0]] = stamp;
            const nn = new Int32Array(6);
            while (top > 0) {
                const cur = stack[--top];
                this.faces(cur, nn);
                for (let k = 0; k < 6; k++) {
                    const nxt = nn[k];
                    if (visit[nxt] === stamp || phi[nxt] !== label || (excludeCentre && nxt === cell)) continue;
                    visit[nxt] = stamp;
                    if (targets.delete(nxt) && targets.size === 0) return true;
                    stack[top++] = nxt;
                }
            }
            return targets.size === 0;
        };

        // Constant-cost 5x5x5 proof used by the geometry-only remap.
        const boundedBranchesConnected = (label, excludeCentre) => {
            const side = 5, rad = 2, centre = 62;
            for (let lz = 0; lz < side; lz++) {
                const zb = zw[z0 + lz - rad + W] * plane;
                for (let ly = 0; ly < side; ly++) {
                    const rb = zb + yw[y0 + ly - rad + W] * Nx;
                    for (let lx = 0; lx < side; lx++) occ[(lz * side + ly) * side + lx] = phi[rb + xw[x0 + lx - rad + W]] === label ? 1 : 0;
                }
            }
            if (excludeCentre) occ[centre] = 0;
            const br = [centre - 1, centre + 1, centre - 5, centre + 5, centre - 25, centre + 25].filter((k) => occ[k]);
            if (br.length <= 1) return true;
            seen.fill(0);
            let head = 0, tail = 0;
            queue[tail++] = br[0]; seen[br[0]] = 1;
            while (head < tail) {
                const cur = queue[head++];
                for (const nx of ADJ5[cur]) if (occ[nx] && !seen[nx]) { seen[nx] = 1; queue[tail++] = nx; }
            }
            for (let k = 1; k < br.length; k++) if (!seen[br[k]]) return false;
            return true;
        };

        const oldObject = topologicalNumber(old);
        const oldExtinction = oldCount <= 1;
        if (oldObject === 0 && !oldExtinction) return false;
        if (oldObject > 1 && !(localOnly ? boundedBranchesConnected(old, true) : globalBranchesConnected(old, true))) return false;
        const newObject = topologicalNumber(newLabel);
        if (newObject === 0) return false;
        if (newObject > 1 && !(localOnly ? boundedBranchesConnected(newLabel, false) : globalBranchesConnected(newLabel, false))) return false;
        return true;
    }

    // ------------------------------------------------------------------
    // 6. Integrate alpha and certify label crossings
    // ------------------------------------------------------------------
    stepAlpha(dt) {
        const { phi, alpha, g, vn, phiComp, gComp, active, nActive, th } = this;
        const idt = dt / this.dx;
        const counts = this.countLabels();
        let accepted = 0, deferred = 0;
        for (let a = 0; a < nActive; a++) {
            const i = active[a];
            const an = Math.min(alpha[i] + vn[i] * idt, 1);
            if (!(an <= 1e-6 && phiComp[i] !== phi[i])) { alpha[i] = an; continue; }
            const oldL = phi[i], newL = phiComp[i];
            if (!th.grainNetwork || this.topologySimpleMove(i, newL, counts[oldL], counts[newL])) {
                phi[i] = newL; g[i] = gComp[i]; alpha[i] = 1;
                counts[oldL]--; counts[newL]++;
                accepted++;
            } else {
                alpha[i] = Math.max(an, 1e-6);
                deferred++;
            }
        }
        this.stats.accepted = accepted;
        this.stats.deferred = deferred;
        this.stats.remapped = 0;
        if (accepted > 0) {
            // A second sweep removes orthogonal one-cell ledges left by the first.
            const first = this.stepBalancedRemap();
            if (first > 0) this.stats.remapped = first + this.stepBalancedRemap();
            else this.stats.remapped = first;
        }
        this.stepCompleteUnresolved();
    }

    // ------------------------------------------------------------------
    // 7. Balanced Gaussian remap
    // ------------------------------------------------------------------
    stepBalancedRemap() {
        this.stepInterfaceMask();
        const { phi, g, alpha, phiComp, active, nActive, resScore, bestScore } = this;
        resScore.fill(-Infinity); bestScore.fill(-Infinity);
        phiComp.set(phi);
        // Requests: the resident label and each distinct unlike face label of every interface cell.
        const cap = 7 * nActive;
        if (!this.rmCells || this.rmCells.length < cap) {
            this.rmCells = new Int32Array(cap); this.rmKeys = new Int32Array(cap); this.rmS = new Float64Array(cap);
        }
        const cells = this.rmCells, keys = this.rmKeys, outS = this.rmS;
        const face = new Int32Array(6), first = new Int32Array(nActive + 1);
        let nReq = 0;
        for (let a = 0; a < nActive; a++) {
            const i = active[a], p = phi[i];
            first[a] = nReq;
            cells[nReq] = i; keys[nReq] = p; nReq++;
            this.faces(i, face);
            for (let k = 0; k < 6; k++) {
                const v = phi[face[k]];
                if (v === p) continue;
                let known = false;
                for (let r = first[a] + 1; r < nReq; r++) if (keys[r] === v) { known = true; break; }
                if (!known) { cells[nReq] = i; keys[nReq] = v; nReq++; }
            }
        }
        first[nActive] = nReq;
        if (this.order.length < nReq) this.order = new Int32Array(nReq);
        const { r, w1 } = this.remapW;
        this.maskedSums(phi, null, w1, r, cells, keys, nReq, outS, null);
        for (let a = 0; a < nActive; a++) {
            const i = active[a], r0 = first[a];
            let best = -1, bestLabel = phi[i];
            for (let m = r0 + 1; m < first[a + 1]; m++) if (outS[m] > best) { best = outS[m]; bestLabel = keys[m]; }
            resScore[i] = outS[r0];
            bestScore[i] = best;
            phiComp[i] = bestLabel;
        }

        const groups = new Map();
        for (let a = 0; a < nActive; a++) {
            const i = active[a], p = phi[i], q = phiComp[i];
            const gain = bestScore[i] - resScore[i];
            if (p === q || !(gain > 1e-7)) continue;
            const low = Math.min(p, q), high = Math.max(p, q), key = low * 1048576 + high;
            let gr = groups.get(key);
            if (!gr) groups.set(key, gr = { low, high, fwd: [], rev: [] });
            (p === low ? gr.fwd : gr.rev).push(i);
        }
        if (!groups.size) { this.stepInterfaceMask(); return 0; }

        const gainOf = (i) => bestScore[i] - resScore[i];
        const rank = (a, b) => (gainOf(b) - gainOf(a)) || (a - b);
        const counts = this.countLabels();
        const fa = new Int32Array(6);
        const targetOrientation = (cell, target) => {
            this.faces(cell, fa);
            for (let k = 0; k < 6; k++) if (phi[fa[k]] === target) return g[fa[k]];
            return g[cell];
        };
        const mismatch = (cell, proposed, paired, pairedProposed) => {
            this.faces(cell, fa);
            let m = 0;
            for (let k = 0; k < 6; k++) m += ((fa[k] === paired ? pairedProposed : phi[fa[k]]) !== proposed) ? 1 : 0;
            return m;
        };
        let exchanges = 0;
        const groupKeys = [...groups.keys()].sort((a, b) => a - b);
        for (const key of groupKeys) {
            const { low, high, fwd, rev } = groups.get(key);
            fwd.sort(rank); rev.sort(rank);
            const n = Math.min(fwd.length, rev.length);
            for (let k = 0; k < n; k++) {
                const a = fwd[k], b = rev[k];
                if (phi[a] !== low || phi[b] !== high) continue;
                const ga = targetOrientation(a, high), gb = targetOrientation(b, low);
                const before = mismatch(a, low, b, high) + mismatch(b, high, a, low);
                const after = mismatch(a, high, b, low) + mismatch(b, low, a, high);
                if (after >= before) continue;
                if (!this.topologySimpleMove(a, high, counts[low], counts[high], true)) continue;
                phi[a] = high;
                if (!this.topologySimpleMove(b, low, counts[high] + 1, counts[low] - 1, true)) { phi[a] = low; continue; }
                phi[b] = low; g[a] = ga; g[b] = gb; alpha[a] = alpha[b] = 1;
                exchanges += 2;
            }
        }
        this.stepInterfaceMask();
        return exchanges;
    }

    // ------------------------------------------------------------------
    // 8. Complete the extinction of labels below the Gaussian scale
    // ------------------------------------------------------------------
    stepCompleteUnresolved() {
        if (!this.th.grainNetwork || (this.stepId & 7) !== 0) return;
        const { Nx, Ny, N, phi, g, alpha, xw, yw, zw, W } = this, plane = Nx * Ny;
        const counts = this.countLabels();
        const sigma = Math.max(this.gauss.sigma, 1);
        const minCells = Math.max(2, Math.ceil(4 * Math.PI * sigma ** 3 / 3));
        const flag = new Uint8Array(counts.length);
        let nUnresolved = 0, nLabels = 0;
        for (let l = 0; l < counts.length; l++) {
            if (counts[l] <= 0) continue;
            nLabels++;
            if (counts[l] < minCells) { flag[l] = 1; nUnresolved++; }
        }
        if (!nUnresolved || nUnresolved === nLabels) return;
        const face = new Int32Array(6);
        for (;;) {
            const cells = [], targets = [], tg = [];
            for (let i = 0; i < N; i++) {
                if (!flag[phi[i]]) continue;
                const x = i % Nx, t = (i / Nx) | 0, y = t % Ny, z = (t / Ny) | 0;
                this.faces(i, face);
                let best = -1, bestSupport = -1;
                for (let k = 0; k < 6; k++) {
                    const q = phi[face[k]];
                    if (flag[q]) continue;
                    let support = 0;
                    for (let dz = -1; dz <= 1; dz++) {
                        const zb = zw[z + dz + W] * plane;
                        for (let dy = -1; dy <= 1; dy++) {
                            const rb = zb + yw[y + dy + W] * Nx;
                            for (let ddx = -1; ddx <= 1; ddx++) if (phi[rb + xw[x + ddx + W]] === q) support++;
                        }
                    }
                    if (support > bestSupport) { best = k; bestSupport = support; }
                }
                if (best >= 0) { cells.push(i); targets.push(phi[face[best]]); tg.push(g[face[best]]); }
            }
            if (!cells.length) break;
            for (let k = 0; k < cells.length; k++) { phi[cells[k]] = targets[k]; g[cells[k]] = tg[k]; alpha[cells[k]] = 1; }
        }
        this.stepInterfaceMask();
    }

    // ------------------------------------------------------------------
    // 9. Transport
    // ------------------------------------------------------------------
    stepDiffusion(dt) {
        if (this.D <= 0) return;
        const { N, c, cNew, mask, alpha } = this;
        const coef = dt * this.D / (this.dx * this.dx);
        const nb = new Int32Array(6);
        for (let i = 0; i < N; i++) {
            this.faces(i, nb);
            const lap = c[nb[0]] + c[nb[1]] + c[nb[2]] + c[nb[3]] + c[nb[4]] + c[nb[5]] - 6 * c[i];
            cNew[i] = c[i] + coef * (mask[i] ? alpha[i] : 1) * lap;
        }
        this.cNew = c; this.c = cNew;
    }

    cellSolidFraction(s) {
        for (let i = 0; i < this.N; i++) s[i] = this.phi[i] !== 0 ? this.alpha[i] : 1 - this.alpha[i];
    }

    stepSoluteRejection() {
        const k = this.th.kPartition;
        if (!(k < 1)) return;
        const { N, phi, c, solidOld } = this;
        const sNew = this.solidNew;
        this.cellSolidFraction(sNew);
        const delta = this.tmp;
        delta.fill(0);
        const nb = new Int32Array(6), w = new Float64Array(6);
        for (let i = 0; i < N; i++) {
            const ds = sNew[i] - solidOld[i];
            if (ds <= 0) continue;
            const s0 = Math.min(1, Math.max(0, solidOld[i]));
            const cL = c[i] / Math.max(1 - s0 * (1 - k), 1e-6);
            let q = (1 - k) * cL * Math.min(ds, 1);
            q = Math.min(q, c[i]);
            if (q <= 0) continue;
            this.faces(i, nb);
            let wsum = 0;
            for (let m = 0; m < 6; m++) { w[m] = phi[nb[m]] === 0 ? Math.max(1 - sNew[nb[m]], 0) : 0; wsum += w[m]; }
            if (wsum <= 0) continue;
            delta[i] -= q;
            for (let m = 0; m < 6; m++) if (w[m] > 0) delta[nb[m]] += q * w[m] / wsum;
        }
        for (let i = 0; i < N; i++) c[i] += delta[i];
    }

    stepDiffusionPartitioned(dt) {
        if (this.D <= 0) return;
        const k = this.th.kPartition;
        const { N, c, cNew } = this;
        const s = this.solidNew, cL = this.tmp;
        this.cellSolidFraction(s);
        for (let i = 0; i < N; i++) cL[i] = c[i] / Math.max(1 - s[i] * (1 - k), 1e-6);
        const coef = dt * this.D / (this.dx * this.dx);
        const nb = new Int32Array(6);
        for (let i = 0; i < N; i++) {
            const li = 1 - s[i];
            let flux = 0;
            this.faces(i, nb);
            for (let m = 0; m < 6; m++) {
                const j = nb[m], lf = Math.min(li, 1 - s[j]);
                if (lf > 0) flux += lf * (cL[j] - cL[i]);
            }
            cNew[i] = c[i] + coef * flux;
        }
        this.cNew = c; this.c = cNew;
    }

    // ------------------------------------------------------------------
    // One complete time step
    // ------------------------------------------------------------------
    step(dtCap = 1e30) {
        this.stepInterfaceMask();
        this.stepGeometry();
        this.stepVelocities();
        this.stepVelocityProjection();
        const dt = Math.min(this.maxDt(), dtCap);
        const partition = !!this.th.soluteRejection;
        if (partition) this.cellSolidFraction(this.solidOld);
        this.stepAlpha(dt);
        if (partition) { this.stepSoluteRejection(); this.stepDiffusionPartitioned(dt); }
        else this.stepDiffusion(dt);
        this.time += dt;
        this.stepId++;
        return dt;
    }

    // ------------------------------------------------------------------
    // Diagnostics and painting
    // ------------------------------------------------------------------
    phaseVolume(label) {
        let v = 0;
        for (let i = 0; i < this.N; i++) {
            if (this.phi[i] === label) v += this.alpha[i];
            if (this.mask[i] && this.phiComp[i] === label && this.phiComp[i] !== this.phi[i]) v += 1 - this.alpha[i];
        }
        return v;
    }

    nLabels() {
        const counts = this.countLabels();
        let n = 0;
        for (let l = 0; l < counts.length; l++) if (counts[l] > 0) n++;
        return n;
    }

    // Paint a ball of `label` (periodic) centred on (cx, cy, cz).
    paintDisc(cx, cy, r, label, gval, cval, cz = this.Nz >> 1) {
        for (let dz = -r; dz <= r; dz++) for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) {
            if (dx * dx + dy * dy + dz * dz > r * r) continue;
            const i = this.idx(cx + dx, cy + dy, cz + dz);
            this.phi[i] = label; this.alpha[i] = 1;
            if (gval !== undefined) this.g[i] = gval;
            if (cval !== undefined) this.c[i] = cval;
        }
    }
}
