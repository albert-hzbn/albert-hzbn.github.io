// Live previews for the builder cards. One off-screen WebGL renderer is shared
// by every preview canvas: it renders the active preview and copies the frame
// into that card's 2D canvas, so any number of previews cost one GL context.
import * as THREE from 'three';
import { element } from '../core/index.js';

// Colours for facet families of a Wulff shape (also used by the legend).
export const FACET_COLORS = ['#e8a33d', '#4e8fd6', '#5bb07a', '#c7648b', '#8f79c9', '#d0765a', '#56b4c9', '#b0a35c'];

const TAG_COLORS = ['#9aa5b1', '#4e79a7', '#f28e2b', '#59a14f', '#e15759', '#76b7b2', '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#86bcb6', '#d37295', '#a0cbe8', '#ffbe7d', '#8cd17d'];

class PreviewRenderer {
    constructor() {
        this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
        this.renderer.setClearColor(0x000000, 0);
        this.renderer.outputColorSpace = THREE.SRGBColorSpace;
        this.scene = new THREE.Scene();
        this.camera = new THREE.OrthographicCamera(-10, 10, 10, -10, -1000, 1000);
        this.scene.add(new THREE.HemisphereLight(0xffffff, 0x80858f, 1.4));
        const key = new THREE.DirectionalLight(0xffffff, 1.7);
        key.position.set(0.5, 0.8, 1);
        this.camera.add(key);
        this.scene.add(this.camera);
        this.group = new THREE.Group();
        this.scene.add(this.group);
        this.active = null;           // { canvas, ctx }
        this.spin = { x: -0.45, y: 0.6 };
        this.drag = null;
        this.autoRotate = !window.matchMedia('(prefers-reduced-motion: reduce)').matches;
        this.tick = this.tick.bind(this);
        requestAnimationFrame(this.tick);
    }

    // Attach drag-to-rotate to a preview canvas (once per canvas).
    bind(canvas) {
        if (canvas.dataset.bound) return;
        canvas.dataset.bound = '1';
        canvas.addEventListener('pointerdown', (e) => {
            this.drag = { x: e.clientX, y: e.clientY, spin: { ...this.spin } };
            canvas.setPointerCapture(e.pointerId);
        });
        canvas.addEventListener('pointermove', (e) => {
            if (!this.drag) return;
            this.spin.y = this.drag.spin.y + (e.clientX - this.drag.x) * 0.01;
            this.spin.x = this.drag.spin.x + (e.clientY - this.drag.y) * 0.01;
        });
        const end = () => { this.drag = null; };
        canvas.addEventListener('pointerup', end);
        canvas.addEventListener('pointercancel', end);
    }

    clear() {
        for (const child of [...this.group.children]) {
            this.group.remove(child);
            child.geometry?.dispose();
            child.material?.dispose();
        }
    }

    // Show `s` in `canvas`. colorBy: 'element' | 'tag'.
    // shape: optional outline (see builders/nanoparticle shapeOutline / wulffShape),
    // drawn around the origin; mode: 'atoms' | 'shape' | 'both'.
    show(canvas, s, { colorBy = 'element', shape = null, mode = 'atoms' } = {}) {
        this.bind(canvas);
        this.active = { canvas, ctx: canvas.getContext('2d') };
        this.clear();
        const drawAtoms = !shape || mode !== 'shape';
        const drawShape = shape && mode !== 'atoms';
        if ((!s || !s.count) && !drawShape) { this.drawEmpty(); return; }
        // Centre on the origin when a shape is shown (particles are built around it),
        // otherwise on the structure's bounds.
        let center = new THREE.Vector3();
        let radius = 2;
        if (s && s.count) {
            const { lo, hi } = s.bounds();
            if (!shape) center.set((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2);
            radius = Math.max(radius, Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2 + 1);
        }
        if (drawShape) radius = Math.max(shape.extent ? shape.extent * 1.08 : 0, drawAtoms ? radius : 0, 2);
        this.radius = radius;
        if (drawAtoms && s && s.count) this.addAtoms(s, center, colorBy);
        if (drawShape) this.addShape(shape, drawAtoms ? 0.28 : 1);
        if (drawAtoms && s && s.periodic) this.addCell(s, center);
    }

    addAtoms(s, center, colorBy) {
        const n = s.count;
        const seg = n < 1500 ? [20, 14] : n < 8000 ? [12, 8] : [8, 6];
        const mesh = new THREE.InstancedMesh(new THREE.SphereGeometry(1, seg[0], seg[1]), new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.05 }), n);
        const m = new THREE.Matrix4(), q = new THREE.Quaternion(), c = new THREE.Color();
        const k = n > 3000 ? 0.36 / 0.42 : 1;
        for (let i = 0; i < n; i++) {
            const p = s.positions[i];
            const r = (0.3 + element(s.symbols[i]).radius * 0.28) * k;
            m.compose(new THREE.Vector3(p[0] - center.x, p[1] - center.y, p[2] - center.z), q, new THREE.Vector3(r, r, r));
            mesh.setMatrixAt(i, m);
            mesh.setColorAt(i, c.set(colorBy === 'tag' ? TAG_COLORS[(s.tags[i] || 0) % TAG_COLORS.length] : element(s.symbols[i]).color));
        }
        this.group.add(mesh);
    }

    addCell(s, center) {
        const o = new THREE.Vector3(-center.x, -center.y, -center.z);
        const [a, b, cc] = s.cell.map((v) => new THREE.Vector3(...v));
        const P = [o, o.clone().add(a), o.clone().add(b), o.clone().add(cc)];
        P.push(P[1].clone().add(b), P[1].clone().add(cc), P[2].clone().add(cc), P[1].clone().add(b).add(cc));
        const E = [[0, 1], [0, 2], [0, 3], [1, 4], [1, 5], [2, 4], [2, 6], [3, 5], [3, 6], [4, 7], [5, 7], [6, 7]];
        const dark = document.documentElement.dataset.theme === 'dark';
        const geo = new THREE.BufferGeometry().setFromPoints(E.flatMap(([i, j]) => [P[i], P[j]]));
        this.group.add(new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color: dark ? '#8a94a0' : '#4a525c', transparent: true, opacity: 0.7 })));
    }

    // Facets coloured by family, with dark edges. opacity < 1 when atoms show through.
    addShape(shape, opacity) {
        const translucent = opacity < 1;
        // Without a colour the per-vertex (per-facet) colours are used.
        const material = (color) => new THREE.MeshStandardMaterial({
            ...(color ? { color } : { vertexColors: true }),
            roughness: 0.55, metalness: 0.05, flatShading: true, side: THREE.DoubleSide,
            transparent: translucent, opacity, depthWrite: !translucent,
        });
        const dark = document.documentElement.dataset.theme === 'dark';
        const edgeMat = new THREE.LineBasicMaterial({ color: dark ? '#dde3ea' : '#2b3138', transparent: true, opacity: translucent ? 0.5 : 0.85 });
        if (shape.kind === 'poly') {
            const pos = [], col = [], edges = [];
            const c = new THREE.Color();
            for (const f of shape.faces) {
                c.set(FACET_COLORS[(f.plane.family || 0) % FACET_COLORS.length]);
                for (let t = 1; t + 1 < f.points.length; t++) {
                    for (const v of [f.points[0], f.points[t], f.points[t + 1]]) { pos.push(...v); col.push(c.r, c.g, c.b); }
                }
                f.points.forEach((v, i) => { edges.push(...v, ...f.points[(i + 1) % f.points.length]); });
            }
            const geo = new THREE.BufferGeometry();
            geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
            geo.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
            geo.computeVertexNormals();
            this.group.add(new THREE.Mesh(geo, material(null)));
            const eg = new THREE.BufferGeometry();
            eg.setAttribute('position', new THREE.Float32BufferAttribute(edges, 3));
            this.group.add(new THREE.LineSegments(eg, edgeMat));
            return;
        }
        let geo;
        if (shape.kind === 'cylinder') {
            geo = new THREE.CylinderGeometry(shape.radius, shape.radius, shape.height, 48, 1);
            geo.rotateX(Math.PI / 2);                      // axis along z
        } else if (shape.kind === 'ellipsoid') {
            geo = new THREE.SphereGeometry(1, 48, 32);
            geo.scale(...shape.radii);
        } else geo = new THREE.SphereGeometry(shape.radius, 48, 32);
        const mat = material(FACET_COLORS[1]);
        mat.flatShading = false;
        this.group.add(new THREE.Mesh(geo, mat));
        this.group.add(new THREE.LineSegments(new THREE.EdgesGeometry(geo, 30), edgeMat));
    }

    drawEmpty() {
        const { canvas, ctx } = this.active;
        ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    tick() {
        requestAnimationFrame(this.tick);
        const a = this.active;
        if (!a || !a.canvas.isConnected || !a.canvas.offsetParent || !this.group.children.length) return;
        const w = a.canvas.clientWidth, h = a.canvas.clientHeight;
        if (!w || !h) return;
        const dpr = Math.min(window.devicePixelRatio || 1, 2);
        if (a.canvas.width !== Math.round(w * dpr) || a.canvas.height !== Math.round(h * dpr)) {
            a.canvas.width = Math.round(w * dpr);
            a.canvas.height = Math.round(h * dpr);
        }
        if (this.autoRotate && !this.drag) this.spin.y += 0.006;
        this.group.rotation.set(this.spin.x, this.spin.y, 0);
        const aspect = w / h, R = this.radius * 1.05;
        Object.assign(this.camera, { left: -R * aspect, right: R * aspect, top: R, bottom: -R });
        this.camera.position.set(0, 0, 100);
        this.camera.lookAt(0, 0, 0);
        this.camera.updateProjectionMatrix();
        this.renderer.setPixelRatio(1);
        this.renderer.setSize(a.canvas.width, a.canvas.height, false);
        this.renderer.render(this.scene, this.camera);
        a.ctx.clearRect(0, 0, a.canvas.width, a.canvas.height);
        a.ctx.drawImage(this.renderer.domElement, 0, 0);
    }
}

let instance = null;
export function previewRenderer() {
    if (!instance) instance = new PreviewRenderer();
    return instance;
}
