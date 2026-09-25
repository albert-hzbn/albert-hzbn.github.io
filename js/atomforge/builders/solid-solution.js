// Random substitutional solid solutions (from AtomForge's
// SubstitutionalSolidSolutionBuilder).
import { mulberry32 } from '../core/index.js';

// solutes: [{ symbol, fraction }] replacing randomly chosen host atoms.
export function buildSolidSolution(s, host, solutes, seed = 1) {
    const out = s.clone();
    const rand = mulberry32(seed);
    const sites = out.symbols.map((x, i) => (x === host ? i : -1)).filter((i) => i >= 0);
    if (!sites.length) throw new Error(`No ${host} atoms in the structure.`);
    // Fisher–Yates shuffle
    for (let i = sites.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [sites[i], sites[j]] = [sites[j], sites[i]]; }
    let k = 0;
    const report = [];
    for (const sol of solutes) {
        const n = Math.round(sol.fraction * sites.length);
        for (let m = 0; m < n && k < sites.length; m++, k++) out.symbols[sites[k]] = sol.symbol;
        report.push(`${n} ${sol.symbol}`);
    }
    out.title = `${out.formula()} solid solution`;
    return { structure: out, report: report.join(', ') + ` on ${sites.length} ${host} sites` };
}
