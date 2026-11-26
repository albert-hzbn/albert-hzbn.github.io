// GPK Web: two-dimensional port of the Gaussian-projected kinetics core
// (https://github.com/albert-hzbn/GPK, src/gpk_*.cpp with Nz = 1).
//
// State per cell: an integer label phi (grain or phase), the resident
// interface position alpha in [0, 1], composition c, temperature T and an
// orientation angle g. One step follows GPKGrid::step():
//
//   1. interface mask        cells with an unlike face neighbour
//   2. geometry              Gaussian-projected normal and curvature of the
//                            resident label, evaluated only at interface cells
//   3. velocities            v_n = M (dG - sigma kappa) from the thermo model
//   4. pair projection       Gaussian average of v_n over the same label pair
//   5. stable time step      dt = min(dx / 8 vmax, 0.4 dx^2 / 4 D)
//   6. alpha and crossings   alpha += v_n dt / dx; a crossing that passes the
//                            digital-topology certificate switches the label
//   7. balanced remap        volume-preserving anti-alias exchanges
//   8. extinction            labels below the Gaussian scale are completed
//   9. transport             composition diffusion (optionally partitioned)
//
// The thermo model is a plain object (see models.js):
//   { grainNetwork, phaseSelection, soluteRejection, kPartition,
//     bulkDG(p, q, c, T), interfaceComp?(p, q, cp, cq),
//     sigma(nx, ny, g, p, q, c), mobility(gp, gq, T),
//     selectProduct?(p, c, T, candidates, n) }

const mod = (a, n) => ((a % n) + n) % n;

class GaussKernel {
    constructor(sigma = 2) {
        this.sigma = Math.max(sigma, 0.5);
        this.radius = Math.ceil(3.5 * this.sigma);
        const size = 2 * this.radius + 1;
        this.w = new Float64Array(size);
        let sum = 0;
        for (let i = 0; i < size; i++) {
            const x = i - this.radius;
            this.w[i] = Math.exp(-x * x / (2 * this.sigma * this.sigma));
            sum += this.w[i];
        }
        for (let i = 0; i < size; i++) this.w[i] /= sum;
    }
}

export class GPK2D {
    constructor(Nx, Ny, dx, thermo, D = 0, sigmaCells = 2) {
        this.Nx = Nx; this.Ny = Ny; this.N = Nx * Ny;
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

        this.gauss = new GaussKernel(sigmaCells);
        this.velocitySigma = 0.85;   // 2-D values of velocity_gauss_ / remap_gauss_
        this.remapSigma = 0.70;
        this.buildGeometryKernels();
        this.velW = GPK2D.windowWeights(this.velocitySigma);
        this.remapW = GPK2D.windowWeights(this.remapSigma);

        // Periodic index tables: xw[x + d + W] = (x + d) mod Nx for |d| <= W.
        const W = this.W = Math.max(this.geomR, 3);
        this.xw = new Int32Array(Nx + 2 * W);
        this.yw = new Int32Array(Ny + 2 * W);
        for (let j = 0; j < this.xw.length; j++) this.xw[j] = mod(j - W, Nx);
        for (let j = 0; j < this.yw.length; j++) this.yw[j] = mod(j - W, Ny);

        this.counts = new Int32Array(16);
        this.stats = { accepted: 0, deferred: 0, remapped: 0 };
    }

    static windowWeights(sigma) {
        const r = Math.max(1, Math.ceil(3 * sigma));
        const s = 2 * r + 1;
        const w = new Float64Array(s * s);
        const inv = 0.5 / (sigma * sigma);
        for (let dy = -r; dy <= r; dy++)
            for (let dx = -r; dx <= r; dx++) w[(dy + r) * s + dx + r] = Math.exp(-(dx * dx + dy * dy) * inv);
        return { r, s, w };
    }

    // Exact convolution of the Gaussian with fourth-order centred stencils.
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

    idx(x, y) { return mod(y, this.Ny) * this.Nx + mod(x, this.Nx); }

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
    // 1. Interface mask and the sorted compact interface set
    // ------------------------------------------------------------------
    stepInterfaceMask() {
        const { Nx, Ny, phi, mask, active } = this;
        let n = 0;
        for (let y = 0; y < Ny; y++) {
            const row = y * Nx;
            const up = (y === 0 ? Ny - 1 : y - 1) * Nx, dn = (y === Ny - 1 ? 0 : y + 1) * Nx;
            for (let x = 0; x < Nx; x++) {
                const i = row + x, v = phi[i];
                const xm = x === 0 ? Nx - 1 : x - 1, xp = x === Nx - 1 ? 0 : x + 1;
                const a = v !== phi[row + xm] || v !== phi[row + xp] || v !== phi[up + x] || v !== phi[dn + x];
                mask[i] = a ? 1 : 0;
                if (a) active[n++] = i;
            }
        }
        this.nActive = n;
    }

    // ------------------------------------------------------------------
    // 2. Gaussian-Hessian geometry: indicator -> normal -> curvature
    // ------------------------------------------------------------------
    stepGeometry() {
        const { Nx, phi, xw, yw, W, w0, w1, w2, active, nActive, nxF, nyF, kappa, phiComp } = this;
        const R = this.geomR, dx = this.dx, invdx = 1 / dx, invdx2 = invdx * invdx;
        const klim = 1 / (this.gauss.sigma * dx);
        nxF.fill(0); nyF.fill(0); kappa.fill(0);
        for (let a = 0; a < nActive; a++) {
            const i = active[a], res = phi[i];
            const x0 = i % Nx, y0 = (i / Nx) | 0;
            let Ix = 0, Iy = 0, Ixx = 0, Iyy = 0, Ixy = 0;
            for (let dy = -R; dy <= R; dy++) {
                const row = yw[y0 + dy + W] * Nx, jy = dy + R;
                let s0 = 0, s1 = 0, s2 = 0;
                for (let dxs = -R; dxs <= R; dxs++) {
                    if (phi[row + xw[x0 + dxs + W]] !== res) continue;
                    const jx = dxs + R;
                    s0 += w0[jx]; s1 += w1[jx]; s2 += w2[jx];
                }
                Ix += s1 * w0[jy]; Iy += s0 * w1[jy];
                Ixx += s2 * w0[jy]; Iyy += s0 * w2[jy]; Ixy += s1 * w1[jy];
            }
            Ix *= invdx; Iy *= invdx; Ixx *= invdx2; Iyy *= invdx2; Ixy *= invdx2;
            const g2 = Ix * Ix + Iy * Iy;
            if (g2 * dx * dx > 1e-8) {
                const gm = Math.sqrt(g2);
                nxF[i] = -Ix / gm;
                nyF[i] = -Iy / gm;
                const gHg = Ix * Ix * Ixx + Iy * Iy * Iyy + 2 * Ix * Iy * Ixy;
                let k = (gHg - g2 * (Ixx + Iyy)) / (g2 * gm);
                if (Math.abs(k * dx) < 1e-6) k = 0;
                kappa[i] = Math.min(klim, Math.max(-klim, k));
            }
        }

        // Competing label: the best locally supported unlike face neighbour,
        // with outward-normal alignment as the tie-breaker.
        for (let a = 0; a < nActive; a++) {
            const i = active[a], x = i % Nx, y = (i / Nx) | 0;
            const nei = [y * Nx + xw[x - 1 + W], y * Nx + xw[x + 1 + W], yw[y - 1 + W] * Nx + x, yw[y + 1 + W] * Nx + x];
            const score = [-nxF[i], nxF[i], -nyF[i], nyF[i]];
            let best = -1, bestSupport = -1, bestScore = -2;
            for (let k = 0; k < 4; k++) {
                const q = phi[nei[k]];
                if (q === phi[i]) continue;
                let support = 0;
                for (let dy = -1; dy <= 1; dy++) {
                    const row = yw[y + dy + W] * Nx;
                    for (let dxs = -1; dxs <= 1; dxs++) {
                        if (!dxs && !dy) continue;
                        if (phi[row + xw[x + dxs + W]] === q) support++;
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
        const { Nx, phi, c, T, g, xw, yw, W, active, nActive, nxF, nyF, kappa, vn, phiComp, gComp, th } = this;
        vn.fill(0);
        const cand = new Int32Array(4);
        for (let a = 0; a < nActive; a++) {
            const i = active[a], x = i % Nx, y = (i / Nx) | 0, p = phi[i];
            const nei = [y * Nx + xw[x - 1 + W], y * Nx + xw[x + 1 + W], yw[y - 1 + W] * Nx + x, yw[y + 1 + W] * Nx + x];
            if (th.phaseSelection) {
                let n = 0;
                for (let k = 0; k < 4; k++) {
                    const q = phi[nei[k]];
                    if (q === p) continue;
                    let seen = false;
                    for (let m = 0; m < n; m++) if (cand[m] === q) { seen = true; break; }
                    if (!seen) cand[n++] = q;
                }
                const sel = th.selectProduct(p, c[i], T[i], cand, n);
                if (sel >= 0) phiComp[i] = sel;
            }
            const q = phiComp[i];
            const gp = g[i];
            let gq = gp, cq = c[i];
            for (let k = 0; k < 4; k++) if (phi[nei[k]] === q) { gq = g[nei[k]]; cq = c[nei[k]]; break; }
            gComp[i] = gq;
            const ci = th.interfaceComp ? th.interfaceComp(p, q, c[i], cq) : c[i];
            const dG = th.bulkDG(p, q, ci, T[i]);
            const sig = th.sigma(nxF[i], nyF[i], gp, p, q, c[i], gq);
            const mob = th.mobility(gp, gq, T[i], nxF[i], nyF[i]);
            vn[i] = mob * (dG - sig * kappa[i]);
        }
    }

    // ------------------------------------------------------------------
    // 4. Project v_n among interface cells of the same label pair
    // ------------------------------------------------------------------
    stepVelocityProjection() {
        const { Nx, phi, phiComp, mask, active, nActive, vn, tmp, xw, yw, W } = this;
        const { r, s, w } = this.velW;
        tmp.set(vn);
        for (let a = 0; a < nActive; a++) {
            const i = active[a], p = phi[i], q = phiComp[i];
            if (p === q) continue;
            const low = Math.min(p, q), high = Math.max(p, q);
            const oi = p === low ? 1 : -1;
            const x0 = i % Nx, y0 = (i / Nx) | 0;
            let num = 0, den = 0;
            for (let dy = -r; dy <= r; dy++) {
                const row = yw[y0 + dy + W] * Nx;
                for (let dxs = -r; dxs <= r; dxs++) {
                    const j = row + xw[x0 + dxs + W];
                    if (!mask[j]) continue;
                    const pj = phi[j], qj = phiComp[j];
                    if (pj === qj) continue;
                    if ((pj < qj ? pj : qj) !== low || (pj > qj ? pj : qj) !== high) continue;
                    const wt = w[(dy + r) * s + dxs + r];
                    num += wt * (pj === low ? 1 : -1) * vn[j];
                    den += wt;
                }
            }
            if (den > 0) tmp[i] = oi * num / den;
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
        const dtDiff = this.D > 0 ? 0.4 * this.dx * this.dx / (4 * this.D) : 1e30;
        return Math.min(dtIface, dtDiff);
    }

    // ------------------------------------------------------------------
    // Digital-topology certificate for one label change (2-D)
    // ------------------------------------------------------------------
    topologySimpleMove(cell, newLabel, oldCount, newCount, localOnly = false) {
        const { Nx, phi, xw, yw, W } = this;
        const old = phi[cell];
        if (old === newLabel) return false;
        const x0 = cell % Nx, y0 = (cell / Nx) | 0;
        const occ = new Uint8Array(9), seen = new Uint8Array(9), queue = new Int8Array(9);

        // Number of face-connected components of `label` in the 3x3
        // neighbourhood (centre removed) that touch the centre by a face.
        const topologicalNumber = (label) => {
            for (let ly = 0; ly < 3; ly++) {
                const row = yw[y0 + ly - 1 + W] * Nx;
                for (let lx = 0; lx < 3; lx++) occ[ly * 3 + lx] = phi[row + xw[x0 + lx - 1 + W]] === label ? 1 : 0;
            }
            occ[4] = 0;
            seen.fill(0);
            let components = 0;
            for (let s = 0; s < 9; s++) {
                if (!occ[s] || seen[s]) continue;
                let head = 0, tail = 0, touches = false;
                queue[tail++] = s; seen[s] = 1;
                while (head < tail) {
                    const cur = queue[head++], cy = (cur / 3) | 0, cx = cur % 3;
                    if (Math.abs(cy - 1) + Math.abs(cx - 1) === 1) touches = true;
                    const nb = [[cy, cx - 1], [cy, cx + 1], [cy - 1, cx], [cy + 1, cx]];
                    for (const [ny, nx] of nb) {
                        if (ny < 0 || ny > 2 || nx < 0 || nx > 2) continue;
                        const nxt = ny * 3 + nx;
                        if (occ[nxt] && !seen[nxt]) { seen[nxt] = 1; queue[tail++] = nxt; }
                    }
                }
                if (touches) components++;
            }
            return components;
        };

        const faces = () => [y0 * Nx + xw[x0 - 1 + W], y0 * Nx + xw[x0 + 1 + W], yw[y0 - 1 + W] * Nx + x0, yw[y0 + 1 + W] * Nx + x0];

        // Early-exit flood fill over the whole periodic grid.
        const globalBranchesConnected = (label, excludeCentre) => {
            const branches = faces().filter((j) => phi[j] === label);
            if (branches.length <= 1) return true;
            const stamp = ++this.visitStamp;
            const visit = this.visit, stack = this.stack;
            const targets = new Set(branches.slice(1));
            let top = 0;
            stack[top++] = branches[0]; visit[branches[0]] = stamp;
            while (top > 0) {
                const cur = stack[--top], cx = cur % Nx, cy = (cur / Nx) | 0;
                const nexts = [cy * Nx + xw[cx - 1 + W], cy * Nx + xw[cx + 1 + W], yw[cy - 1 + W] * Nx + cx, yw[cy + 1 + W] * Nx + cx];
                for (const nxt of nexts) {
                    if (visit[nxt] === stamp || phi[nxt] !== label || (excludeCentre && nxt === cell)) continue;
                    visit[nxt] = stamp;
                    if (targets.delete(nxt) && targets.size === 0) return true;
                    stack[top++] = nxt;
                }
            }
            return targets.size === 0;
        };

        // Constant-cost 5x5 proof used by the geometry-only remap.
        const boundedBranchesConnected = (label, excludeCentre) => {
            const side = 5, rad = 2, centre = rad * side + rad;
            const o = new Uint8Array(25), sn = new Uint8Array(25), q = new Int8Array(25);
            for (let ly = 0; ly < side; ly++) {
                const row = yw[y0 + ly - rad + W] * Nx;
                for (let lx = 0; lx < side; lx++) o[ly * side + lx] = phi[row + xw[x0 + lx - rad + W]] === label ? 1 : 0;
            }
            if (excludeCentre) o[centre] = 0;
            const br = [centre - 1, centre + 1, centre - side, centre + side].filter((k) => o[k]);
            if (br.length <= 1) return true;
            let head = 0, tail = 0;
            q[tail++] = br[0]; sn[br[0]] = 1;
            while (head < tail) {
                const cur = q[head++], cy = (cur / side) | 0, cx = cur % side;
                for (const [ny, nx] of [[cy, cx - 1], [cy, cx + 1], [cy - 1, cx], [cy + 1, cx]]) {
                    if (ny < 0 || ny >= side || nx < 0 || nx >= side) continue;
                    const nxt = ny * side + nx;
                    if (o[nxt] && !sn[nxt]) { sn[nxt] = 1; q[tail++] = nxt; }
                }
            }
            for (let k = 1; k < br.length; k++) if (!sn[br[k]]) return false;
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
        if (accepted > 0) this.stepBalancedRemap();
        this.stepCompleteUnresolved();
    }

    // ------------------------------------------------------------------
    // 7. Balanced Gaussian remap: volume-preserving anti-alias exchanges
    // ------------------------------------------------------------------
    stepBalancedRemap() {
        this.stepInterfaceMask();
        const { Nx, phi, g, alpha, phiComp, active, nActive, xw, yw, W, resScore, bestScore } = this;
        const { r, s, w } = this.remapW;
        resScore.fill(-Infinity); bestScore.fill(-Infinity);
        phiComp.set(phi);
        const cand = [0, 0, 0, 0], support = [0, 0, 0, 0];
        for (let a = 0; a < nActive; a++) {
            const i = active[a], x0 = i % Nx, y0 = (i / Nx) | 0, p = phi[i];
            const face = [y0 * Nx + xw[x0 - 1 + W], y0 * Nx + xw[x0 + 1 + W], yw[y0 - 1 + W] * Nx + x0, yw[y0 + 1 + W] * Nx + x0];
            let nc = 0;
            for (let k = 0; k < 4; k++) {
                const v = phi[face[k]];
                if (v === p) continue;
                let known = false;
                for (let m = 0; m < nc; m++) if (cand[m] === v) known = true;
                if (!known) cand[nc++] = v;
            }
            let resident = 0, best = -1, bestLabel = p;
            support.fill(0);
            for (let dy = -r; dy <= r; dy++) {
                const row = yw[y0 + dy + W] * Nx;
                for (let dxs = -r; dxs <= r; dxs++) {
                    const v = phi[row + xw[x0 + dxs + W]], wt = w[(dy + r) * s + dxs + r];
                    if (v === p) resident += wt;
                    for (let m = 0; m < nc; m++) if (v === cand[m]) support[m] += wt;
                }
            }
            resScore[i] = resident;
            for (let m = 0; m < nc; m++) if (support[m] > best) { best = support[m]; bestLabel = cand[m]; }
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
        if (!groups.size) { this.stepInterfaceMask(); return; }

        const gainOf = (i) => bestScore[i] - resScore[i];
        const rank = (a, b) => (gainOf(b) - gainOf(a)) || (a - b);
        const counts = this.countLabels();
        const faces = (cell) => {
            const x = cell % Nx, y = (cell / Nx) | 0;
            return [y * Nx + xw[x - 1 + W], y * Nx + xw[x + 1 + W], yw[y - 1 + W] * Nx + x, yw[y + 1 + W] * Nx + x];
        };
        const targetOrientation = (cell, target) => {
            for (const j of faces(cell)) if (phi[j] === target) return g[j];
            return g[cell];
        };
        const mismatch = (cell, proposed, paired, pairedProposed) => {
            let m = 0;
            for (const j of faces(cell)) m += ((j === paired ? pairedProposed : phi[j]) !== proposed) ? 1 : 0;
            return m;
        };
        let exchanges = 0;
        const keys = [...groups.keys()].sort((a, b) => a - b);   // std::map order
        for (const key of keys) {
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
        this.stats.remapped = exchanges;
        this.stepInterfaceMask();
    }

    // ------------------------------------------------------------------
    // 8. Complete the extinction of labels below the Gaussian scale
    // ------------------------------------------------------------------
    stepCompleteUnresolved() {
        if (!this.th.grainNetwork || (this.stepId & 7) !== 0) return;
        const { Nx, N, phi, g, alpha, xw, yw, W } = this;
        const counts = this.countLabels();
        const sigma = Math.max(this.gauss.sigma, 1);
        const minCells = Math.max(2, Math.ceil(Math.PI * sigma * sigma));
        const flag = new Uint8Array(counts.length);
        let nUnresolved = 0, nLabels = 0;
        for (let l = 0; l < counts.length; l++) {
            if (counts[l] <= 0) continue;
            nLabels++;
            if (counts[l] < minCells) { flag[l] = 1; nUnresolved++; }
        }
        if (!nUnresolved || nUnresolved === nLabels) return;
        for (;;) {
            const cells = [], targets = [], tg = [];
            for (let i = 0; i < N; i++) {
                if (!flag[phi[i]]) continue;
                const x = i % Nx, y = (i / Nx) | 0;
                const face = [y * Nx + xw[x - 1 + W], y * Nx + xw[x + 1 + W], yw[y - 1 + W] * Nx + x, yw[y + 1 + W] * Nx + x];
                let best = -1, bestSupport = -1;
                for (let k = 0; k < 4; k++) {
                    const q = phi[face[k]];
                    if (flag[q]) continue;
                    let support = 0;
                    for (let dy = -1; dy <= 1; dy++) {
                        const row = yw[y + dy + W] * Nx;
                        for (let dxs = -1; dxs <= 1; dxs++) if (phi[row + xw[x + dxs + W]] === q) support++;
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
        const { Nx, Ny, c, cNew, mask, alpha } = this;
        const coef = dt * this.D / (this.dx * this.dx);
        for (let y = 0; y < Ny; y++) {
            const row = y * Nx, up = (y === 0 ? Ny - 1 : y - 1) * Nx, dn = (y === Ny - 1 ? 0 : y + 1) * Nx;
            for (let x = 0; x < Nx; x++) {
                const i = row + x, xm = x === 0 ? Nx - 1 : x - 1, xp = x === Nx - 1 ? 0 : x + 1;
                const lap = c[row + xm] + c[row + xp] + c[up + x] + c[dn + x] - 4 * c[i];
                cNew[i] = c[i] + coef * (mask[i] ? alpha[i] : 1) * lap;
            }
        }
        this.cNew = c; this.c = cNew;
    }

    // Solid fraction of a cell: alpha for a solid resident, 1 - alpha for liquid (label 0).
    cellSolidFraction(s) {
        for (let i = 0; i < this.N; i++) s[i] = this.phi[i] !== 0 ? this.alpha[i] : 1 - this.alpha[i];
    }

    neighbours4(i) {
        const { Nx, Ny } = this, x = i % Nx, y = (i / Nx) | 0;
        return [y * Nx + (x === Nx - 1 ? 0 : x + 1), y * Nx + (x === 0 ? Nx - 1 : x - 1),
            (y === Ny - 1 ? 0 : y + 1) * Nx + x, (y === 0 ? Ny - 1 : y - 1) * Nx + x];
    }

    // Growing solid keeps k c_L and passes (1 - k) c_L ds to liquid neighbours.
    stepSoluteRejection() {
        const k = this.th.kPartition;
        if (!(k < 1)) return;
        const { N, phi, c, solidOld } = this;
        const sNew = this.solidNew;
        this.cellSolidFraction(sNew);
        const delta = this.tmp;
        delta.fill(0);
        for (let i = 0; i < N; i++) {
            const ds = sNew[i] - solidOld[i];
            if (ds <= 0) continue;
            const s0 = Math.min(1, Math.max(0, solidOld[i]));
            const cL = c[i] / Math.max(1 - s0 * (1 - k), 1e-6);
            let q = (1 - k) * cL * Math.min(ds, 1);
            q = Math.min(q, c[i]);
            if (q <= 0) continue;
            const nb = this.neighbours4(i);
            let wsum = 0;
            const w = nb.map((j) => (phi[j] === 0 ? Math.max(1 - sNew[j], 0) : 0));
            for (const v of w) wsum += v;
            if (wsum <= 0) continue;
            delta[i] -= q;
            for (let m = 0; m < 4; m++) if (w[m] > 0) delta[nb[m]] += q * w[m] / wsum;
        }
        for (let i = 0; i < N; i++) c[i] += delta[i];
    }

    // Conservative face-flux diffusion of the liquid composition.
    stepDiffusionPartitioned(dt) {
        if (this.D <= 0) return;
        const k = this.th.kPartition;
        const { N, c, cNew } = this;
        const s = this.solidNew, cL = this.tmp;
        this.cellSolidFraction(s);
        for (let i = 0; i < N; i++) cL[i] = c[i] / Math.max(1 - s[i] * (1 - k), 1e-6);
        const coef = dt * this.D / (this.dx * this.dx);
        for (let i = 0; i < N; i++) {
            const li = 1 - s[i];
            let flux = 0;
            for (const j of this.neighbours4(i)) {
                const lf = Math.min(li, 1 - s[j]);
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
    // Diagnostics
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

    // Paint a disc of `label` (periodic), used for nuclei and click seeding.
    paintDisc(cx, cy, r, label, gval, cval) {
        for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) {
            if (dx * dx + dy * dy > r * r) continue;
            const i = this.idx(cx + dx, cy + dy);
            this.phi[i] = label; this.alpha[i] = 1;
            if (gval !== undefined) this.g[i] = gval;
            if (cval !== undefined) this.c[i] = cval;
        }
    }
}
