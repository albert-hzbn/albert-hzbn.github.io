// GPK Web: simulation worker. Runs the engine off the main thread and posts
// copies of the fields after each batch of steps.
import { GPK2D } from './engine.js';
import { MODELS, resolveParams, rng } from './models.js';

let sim = null, model = null, P = null, I = null, rand = null, cfg = null;
let running = false, timer = 0, stepsPerBatch = 0, steps = 0;
let history = [], pending = [], lastRecord = -Infinity;

function build(msg) {
    model = MODELS.find((m) => m.id === msg.model);
    cfg = msg;
    P = resolveParams(model, msg.params);
    I = msg.init;
    rand = rng(msg.seed);
    sim = new GPK2D(msg.n, msg.n, msg.dx * 1e-6, model.thermo(P), model.diffusivity ? model.diffusivity(P) : 0);
    model.setup(sim, P, I, rand);
    steps = 0;
    history = []; pending = []; lastRecord = -Infinity;
    // Geometry for the first frame (normals, curvature, competitors).
    sim.stepInterfaceMask();
    sim.stepGeometry();
    sim.stepVelocities();
    sim.stepVelocityProjection();
    record();
    post(true);
}

// Statistics are sampled by simulated time (about 1500 points per run), so
// the chart has the same resolution however fast the steps are.
function record(force = true, s = null) {
    if (!force && sim.time - lastRecord < cfg.tEnd / 1500) return;
    const p = { t: sim.time, ...(s || model.stats(sim)) };
    lastRecord = sim.time;
    history.push(p);
    pending.push(p);
    if (history.length > 4000) history = history.filter((_, k) => k % 2 === 0);
}

function post(full = false) {
    const f = {
        phi: sim.phi.slice(), alpha: sim.alpha.slice(), comp: sim.phiComp.slice(), mask: sim.mask.slice(),
        c: sim.c.slice(), g: sim.g.slice(), kappa: sim.kappa.slice(), vn: sim.vn.slice(),
        nx: sim.nxF.slice(), ny: sim.nyF.slice(),
    };
    const last = history[history.length - 1];
    self.postMessage({
        type: 'frame', n: sim.Nx, step: steps, time: sim.time, stats: last, running,
        history: full ? history : null, points: full ? null : pending, maxLabel: sim.maxLabel(), ...f,
    }, Object.values(f).map((a) => a.buffer));
    pending = [];
}

function advance(maxSteps, budgetMs) {
    const t0 = performance.now();
    let k = 0;
    while (k < maxSteps && performance.now() - t0 < budgetMs) {
        const remaining = cfg.tEnd - sim.time;
        if (remaining <= 1e-15 * cfg.tEnd) return 'end';
        const dt = sim.step(remaining);
        if (model.afterStep) model.afterStep(sim, P, I, rand, dt);
        steps++; k++;
        const due = sim.time - lastRecord >= cfg.tEnd / 1500;
        const s = model.stop || due ? model.stats(sim) : null;
        if (due) record(true, s);
        const why = model.stop && model.stop(s, I, sim);
        if (why) return typeof why === 'string' ? why : 'done';
    }
    return null;
}

function loop() {
    timer = 0;
    if (!running || !sim) return;
    const end = advance(stepsPerBatch || Infinity, stepsPerBatch ? Infinity : 30);
    record();
    if (end) {
        running = false;
        self.postMessage({ type: 'stopped', reason: end });
    }
    post();
    if (running) timer = setTimeout(loop, 0);
}

self.onmessage = (e) => {
    const m = e.data;
    try {
        if (m.type === 'init') {
            running = false;
            clearTimeout(timer); timer = 0;
            build(m);
        } else if (!sim) {
            return;
        } else if (m.type === 'run') {
            running = m.on;
            if (running && !timer) timer = setTimeout(loop, 0);
        } else if (m.type === 'speed') {
            stepsPerBatch = m.steps;
        } else if (m.type === 'step') {
            advance(m.steps, Infinity);
            record();
            post();
        } else if (m.type === 'params') {
            P = resolveParams(model, m.params);
            sim.th = model.thermo(P);
            sim.D = model.diffusivity ? model.diffusivity(P) : 0;
            if (model.live) model.live(sim, P, I);
            if (!running) { sim.stepInterfaceMask(); sim.stepGeometry(); sim.stepVelocities(); sim.stepVelocityProjection(); post(); }
        } else if (m.type === 'tEnd') {
            cfg.tEnd = m.tEnd;
        } else if (m.type === 'seed') {
            if (!model.seed) return;
            const note = model.seed(sim, m.x, m.y, P, rand, I);
            sim.stepInterfaceMask(); sim.stepGeometry(); sim.stepVelocities(); sim.stepVelocityProjection();
            post();
            self.postMessage({ type: 'note', text: note });
        }
    } catch (err) {
        running = false;
        self.postMessage({ type: 'error', message: err.message });
    }
};
