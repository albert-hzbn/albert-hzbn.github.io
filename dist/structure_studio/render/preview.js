// Live previews for the builder cards. One off-screen WebGL renderer is shared
// by every preview canvas: it renders the active preview and copies the frame
// into that card's 2D canvas, so any number of previews cost one GL context.
import * as THREE from 'three';
import { element } from '../core/index.js';

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
    show(canvas, s, { colorBy = 'element' } = {}) {
        this.bind(canvas);
        this.active = { canvas, ctx: canvas.getContext('2d') };
        this.clear();
        if (!s || !s.count) { this.drawEmpty(); return; }
        // Centre on the cell (or the atoms) and size the camera to fit.
        const { lo, hi } = s.bounds();
        const center = new THREE.Vector3((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2);
        const radius = Math.max(2, Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2 + 1);
        this.radius = radius;
        const n = s.count;
        const seg = n < 1500 ? [20, 14] : n < 8000 ? [12, 8] : [8, 6];
        const mesh = new THREE.InstancedMesh(new THREE.SphereGeometry(1, seg[0], seg[1]), new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.05 }), n);
        const m = new THREE.Matrix4(), q = new THREE.Quaternion(), c = new THREE.Color();
        const scaleR = n > 3000 ? 0.36 : 0.42;
        for (let i = 0; i < n; i++) {
            const p = s.positions[i];
            const r = (0.3 + element(s.symbols[i]).radius * 0.28) * (scaleR / 0.42);
            m.compose(new THREE.Vector3(p[0] - center.x, p[1] - center.y, p[2] - center.z), q, new THREE.Vector3(r, r, r));
            mesh.setMatrixAt(i, m);
            mesh.setColorAt(i, c.set(colorBy === 'tag' ? TAG_COLORS[(s.tags[i] || 0) % TAG_COLORS.length] : element(s.symbols[i]).color));
        }
        this.group.add(mesh);
        if (s.periodic) {
            const o = new THREE.Vector3(-center.x, -center.y, -center.z);
            const [a, b, cc] = s.cell.map((v) => new THREE.Vector3(...v));
            const pts = [o, o.clone().add(a), o.clone().add(b), o.clone().add(cc)];
            const P = [pts[0], pts[1], pts[2], pts[3], pts[1].clone().add(b), pts[1].clone().add(cc), pts[2].clone().add(cc), pts[1].clone().add(b).add(cc)];
            const E = [[0, 1], [0, 2], [0, 3], [1, 4], [1, 5], [2, 4], [2, 6], [3, 5], [3, 6], [4, 7], [5, 7], [6, 7]];
            const dark = document.documentElement.dataset.theme === 'dark';
            const geo = new THREE.BufferGeometry().setFromPoints(E.flatMap(([i, j]) => [P[i], P[j]]));
            this.group.add(new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color: dark ? '#8a94a0' : '#4a525c', transparent: true, opacity: 0.7 })));
        }
    }

    drawEmpty() {
        const { canvas, ctx } = this.active;
        ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    tick(t) {
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
        void t;
    }
}

let instance = null;
export function previewRenderer() {
    if (!instance) instance = new PreviewRenderer();
    return instance;
}
