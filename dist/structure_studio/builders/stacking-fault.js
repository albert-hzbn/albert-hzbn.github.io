// FCC stacking faults along the generalised stacking-fault path (from
// AtomForge's StackingFaultBuilder, tilted-cell set-up).
import { add } from '../core/index.js';
import { species, orientedBox } from './common.js';

// Displacement u is in units of the Shockley partial b_p = a/√6 along [11-2].
// The cell's c vector is tilted by the same displacement, so the cell holds a
// single fault plane and stays periodic (the usual GSFE set-up).
export function buildStackingFault(fccCrystal, { reps = [2, 1, 4], u = 1 } = {}) {
    const box = orientedBox(fccCrystal, [[1, -1, 0], [1, 1, -2], [1, 1, 1]], reps);
    const a = fccCrystal.cell[0][0];
    const bp = a / Math.sqrt(6);
    const zs = [...new Set(box.positions.map((p) => +p[2].toFixed(4)))].sort((x, y) => x - y);
    const mid = zs[Math.floor(zs.length / 2)] - 1e-3;
    const shift = [0, u * bp, 0];
    box.positions = box.positions.map((p) => (p[2] >= mid ? add(p, shift) : p));
    box.tags = box.positions.map((p) => (p[2] >= mid ? 2 : 1));
    box.cell[2] = add(box.cell[2], shift);
    box.title = `${species(fccCrystal)} stacking fault, u = ${u.toFixed(2)} bₚ`;
    return box.wrap();
}
