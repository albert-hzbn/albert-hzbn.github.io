// Runs a phase-field model on the GPU: double-buffered state fields, auxiliary
// fields, compute pipelines for each pass, and a display pass that writes a
// colour-mapped rgba8 3D texture (rgb = colour, a = scalar used for volume and
// isosurface rendering).
import { stepModule, displayModule, paramIndex, PARAM_SLOTS, WORKGROUP } from './wgsl.js';

const MAX_SUBSTEPS = 256;          // steps per frame (one uniform slot each)
const SLOT = 256;                  // dynamic-offset alignment

export class Simulation {
    constructor(device) {
        this.device = device;
        this.model = null;
        this.n = 0;
        this.parity = 0;           // which buffer of each pair holds the current state
        this.step = 0;             // total steps taken
        this.time = 0;             // simulated time (steps × dt)
        this.seed = 1;
        this.paramsBuf = device.createBuffer({ size: 256, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.gridBuf = device.createBuffer({ size: SLOT * MAX_SUBSTEPS, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.gridStatic = device.createBuffer({ size: 16, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.displayBuf = device.createBuffer({ size: 32, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.gridData = new Uint32Array(SLOT / 4 * MAX_SUBSTEPS);
        this.displayPipelines = new Map();
        this.displayModules = new Map();
    }

    // Bytes of GPU memory the model needs at grid size n.
    static memoryEstimate(model, n) {
        const comps = model.fields.reduce((a, f) => a + 2 * f.comps, 0) + (model.aux || []).reduce((a, f) => a + f.comps, 0);
        return comps * n ** 3 * 4 + n ** 3 * 4;
    }

    // Build buffers and pipelines for `model` at grid size n.
    load(model, n) {
        const d = this.device;
        this.destroyBuffers();
        this.model = model;
        this.n = n;
        const n3 = n ** 3;
        const usage = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC;
        const maxBinding = d.limits.maxStorageBufferBindingSize;
        this.fields = model.fields.map((f) => {
            const size = f.comps * n3 * 4;
            if (size > maxBinding) throw new Error(`${f.name} needs ${(size / 2 ** 20).toFixed(0)} MB, above this GPU's ${(maxBinding / 2 ** 20).toFixed(0)} MB buffer limit. Use a smaller grid.`);
            return { ...f, bufs: [d.createBuffer({ size, usage }), d.createBuffer({ size, usage })] };
        });
        this.aux = (model.aux || []).map((a) => ({ ...a, buf: d.createBuffer({ size: a.comps * n3 * 4, usage }) }));
        this.texture = d.createTexture({
            size: [n, n, n], dimension: '3d', format: 'rgba8unorm',
            usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC,
        });
        d.queue.writeBuffer(this.gridStatic, 0, new Uint32Array([n, 0, this.seed, 0]));

        // Step pipelines -----------------------------------------------------
        const entries = [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform', hasDynamicOffset: true } },
        ];
        let b = 2;
        for (let j = 0; j < this.fields.length; j++) {
            entries.push({ binding: b++, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } });
            entries.push({ binding: b++, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } });
        }
        for (let j = 0; j < this.aux.length; j++) entries.push({ binding: b++, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } });
        this.stepLayout = d.createBindGroupLayout({ entries });
        const module = d.createShaderModule({ code: stepModule(model), label: `${model.id} step` });
        this.stepModule = module;
        const layout = d.createPipelineLayout({ bindGroupLayouts: [this.stepLayout] });
        this.passes = model.passes.map((ps) => d.createComputePipeline({ layout, compute: { module, entryPoint: ps.entry }, label: ps.entry }));
        this.stepGroups = [0, 1].map((par) => d.createBindGroup({
            layout: this.stepLayout,
            entries: [
                { binding: 0, resource: { buffer: this.paramsBuf } },
                { binding: 1, resource: { buffer: this.gridBuf, size: 16 } },
                ...this.fields.flatMap((f, j) => [
                    { binding: 2 + 2 * j, resource: { buffer: f.bufs[par] } },
                    { binding: 3 + 2 * j, resource: { buffer: f.bufs[1 - par] } },
                ]),
                ...this.aux.map((a, j) => ({ binding: 2 + 2 * this.fields.length + j, resource: { buffer: a.buf } })),
            ],
        }));

        // Display pipelines (compiled lazily per view) -----------------------
        const dEntries = [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        ];
        b = 3;
        for (let j = 0; j < this.fields.length + this.aux.length; j++) dEntries.push({ binding: b++, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } });
        dEntries.push({ binding: b, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'rgba8unorm', viewDimension: '3d' } });
        this.displayLayout = d.createBindGroupLayout({ entries: dEntries });
        this.displayPipelines = new Map();
        this.displayModules = new Map();
        const texView = this.texture.createView();
        this.displayGroups = [0, 1].map((par) => d.createBindGroup({
            layout: this.displayLayout,
            entries: [
                { binding: 0, resource: { buffer: this.paramsBuf } },
                { binding: 1, resource: { buffer: this.gridStatic } },
                { binding: 2, resource: { buffer: this.displayBuf } },
                ...this.fields.map((f, j) => ({ binding: 3 + j, resource: { buffer: f.bufs[par] } })),
                ...this.aux.map((a, j) => ({ binding: 3 + this.fields.length + j, resource: { buffer: a.buf } })),
                { binding: 3 + this.fields.length + this.aux.length, resource: texView },
            ],
        }));
        this.parity = 0;
        this.step = 0;
        this.time = 0;
    }

    // Compilation errors of the step shaders and of every display view.
    async compileErrors() {
        const modules = [this.stepModule];
        for (const v of this.model.views) {
            this.displayPipeline(v);
            modules.push(this.displayModules.get(v.id));
        }
        const out = [];
        for (const m of modules) {
            const info = await m.getCompilationInfo();
            for (const msg of info.messages) if (msg.type === 'error') out.push(`${m.label} ${msg.lineNum}:${msg.linePos} ${msg.message}`);
        }
        return out;
    }

    displayPipeline(view) {
        let p = this.displayPipelines.get(view.id);
        if (!p) {
            const module = this.device.createShaderModule({ code: displayModule(this.model, view), label: `${this.model.id} view ${view.id}` });
            this.displayModules.set(view.id, module);
            p = this.device.createComputePipeline({ layout: this.device.createPipelineLayout({ bindGroupLayouts: [this.displayLayout] }), compute: { module, entryPoint: 'display' } });
            this.displayPipelines.set(view.id, p);
        }
        return p;
    }

    setParams(values) {
        const m = this.model;
        const data = new Float32Array(64);
        m.params.forEach((p, j) => { data[j] = values[p.key] ?? p.value; });
        const dx = paramIndex(m, 'dx') >= 0 ? values.dx : 1;
        data[63] = 1 / (dx * dx);
        if (m.params.length > PARAM_SLOTS) throw new Error('Too many parameters');
        this.device.queue.writeBuffer(this.paramsBuf, 0, data);
        this.dt = values.dt ?? 1;
    }

    // Upload the initial state: { fieldName: Float32Array(comps·n³) }.
    upload(state, seed = 1) {
        this.seed = seed;
        this.parity = 0;
        this.step = 0;
        this.time = 0;
        for (const f of this.fields) {
            this.device.queue.writeBuffer(f.bufs[0], 0, state[f.name]);
            this.device.queue.writeBuffer(f.bufs[1], 0, state[f.name]);
        }
        for (const a of this.aux) this.device.queue.writeBuffer(a.buf, 0, state[a.name] || new Float32Array(a.comps * this.n ** 3));
        this.device.queue.writeBuffer(this.gridStatic, 0, new Uint32Array([this.n, 0, seed, 0]));
    }

    // Advance `count` time steps in one submission.
    advance(count) {
        count = Math.max(1, Math.min(MAX_SUBSTEPS, count | 0));
        const d = this.device;
        const g = this.gridData;
        for (let s = 0; s < count; s++) {
            const o = s * (SLOT / 4);
            g[o] = this.n; g[o + 1] = this.step + s; g[o + 2] = this.seed; g[o + 3] = 0;
        }
        d.queue.writeBuffer(this.gridBuf, 0, g, 0, count * (SLOT / 4));
        const wg = Math.ceil(this.n / WORKGROUP);
        const enc = d.createCommandEncoder();
        const pass = enc.beginComputePass();
        for (let s = 0; s < count; s++) {
            const group = this.stepGroups[this.parity];
            for (const pipe of this.passes) {
                pass.setPipeline(pipe);
                pass.setBindGroup(0, group, [s * SLOT]);
                pass.dispatchWorkgroups(wg, wg, wg);
            }
            this.parity ^= 1;
        }
        pass.end();
        d.queue.submit([enc.finish()]);
        this.step += count;
        this.time += count * this.dt;
    }

    // Fill the display texture from the current state.
    display(view, { lo = 0, hi = 1, cmap = 0 } = {}) {
        const d = this.device;
        d.queue.writeBuffer(this.displayBuf, 0, new Float32Array([lo, hi, cmap, 0, 0, 0, 0, 0]));
        const enc = d.createCommandEncoder();
        const pass = enc.beginComputePass();
        pass.setPipeline(this.displayPipeline(view));
        pass.setBindGroup(0, this.displayGroups[this.parity]);
        const wg = Math.ceil(this.n / WORKGROUP);
        pass.dispatchWorkgroups(wg, wg, wg);
        pass.end();
        d.queue.submit([enc.finish()]);
    }

    // Copy the current state of a field back to the CPU.
    async readField(name) {
        const f = this.fields.find((x) => x.name === name);
        const buf = f ? f.bufs[this.parity] : this.aux.find((x) => x.name === name).buf;
        const size = buf.size;
        const staging = this.device.createBuffer({ size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
        const enc = this.device.createCommandEncoder();
        enc.copyBufferToBuffer(buf, 0, staging, 0, size);
        this.device.queue.submit([enc.finish()]);
        await staging.mapAsync(GPUMapMode.READ);
        const out = new Float32Array(staging.getMappedRange().slice(0));
        staging.unmap();
        staging.destroy();
        return out;
    }

    destroyBuffers() {
        for (const f of this.fields || []) f.bufs.forEach((b) => b.destroy());
        for (const a of this.aux || []) a.buf.destroy();
        this.texture?.destroy();
        this.fields = [];
        this.aux = [];
    }
}
