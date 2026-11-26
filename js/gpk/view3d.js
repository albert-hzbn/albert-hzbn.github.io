// GPK Web: 3D view. Builds a smoothed surface (naive surface nets) around
// the visible cells of the volume and draws it with three.js. Loaded only
// when a 3D run is shown, so 2D use never fetches three.js.
import * as THREE from 'three';
import { OrbitControls } from 'jsm/controls/OrbitControls.js';
export { surfaceNets } from './surface.js';

export class View3D {
    constructor(canvas) {
        this.canvas = canvas;
        this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
        this.renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(32, 1, 0.1, 5000);
        this.camera.up.set(0, 0, 1);
        this.controls = new OrbitControls(this.camera, canvas);
        this.controls.enableDamping = false;
        this.controls.addEventListener('change', () => this.render());
        this.scene.add(new THREE.HemisphereLight(0xffffff, 0x606878, 1.1));
        const key = new THREE.DirectionalLight(0xffffff, 1.5);
        key.position.set(1.2, -0.8, 1.6);
        this.camera.add(key);
        this.scene.add(this.camera);
        this.material = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.7, metalness: 0.0, side: THREE.DoubleSide });
        this.mesh = null;
        this.box = null;
        this.n = 0;
    }

    fit(n) {
        if (this.n === n) return;
        this.n = n;
        if (this.box) { this.scene.remove(this.box); this.box.geometry.dispose(); }
        const g = new THREE.EdgesGeometry(new THREE.BoxGeometry(n, n, n));
        this.box = new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: 0x8a8f98, transparent: true, opacity: 0.7 }));
        this.scene.add(this.box);
        const d = 2.45 * n;
        this.camera.position.set(0.62 * d, -0.95 * d, 0.62 * d);
        this.controls.target.set(0, 0, 0);
        this.controls.update();
    }

    setMesh(n, surf) {
        this.fit(n);
        if (this.mesh) { this.scene.remove(this.mesh); this.mesh.geometry.dispose(); this.mesh = null; }
        if (!surf.index.length) { this.render(); return; }
        const g = new THREE.BufferGeometry();
        // Cell centres sit at integer coordinates, so the box spans [-0.5, n - 0.5].
        const p = surf.position, h = n / 2 - 0.5;
        for (let i = 0; i < p.length; i += 3) { p[i] -= h; p[i + 1] -= h; p[i + 2] -= h; }
        g.setAttribute('position', new THREE.BufferAttribute(p, 3));
        g.setIndex(new THREE.BufferAttribute(surf.index, 1));
        // Smooth normals from the shared vertices, then one colour per triangle.
        g.computeVertexNormals();
        const flat = g.toNonIndexed();
        g.dispose();
        const tc = surf.triColor, nt = tc.length / 3, col = new Float32Array(9 * nt);
        // The page colours are sRGB; three.js expects linear vertex colours.
        const lin = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
        for (let t = 0; t < nt; t++) {
            const r = lin(tc[3 * t]), gg = lin(tc[3 * t + 1]), b = lin(tc[3 * t + 2]);
            for (let k = 0; k < 3; k++) { col[9 * t + 3 * k] = r; col[9 * t + 3 * k + 1] = gg; col[9 * t + 3 * k + 2] = b; }
        }
        flat.setAttribute('color', new THREE.BufferAttribute(col, 3));
        this.mesh = new THREE.Mesh(flat, this.material);
        this.scene.add(this.mesh);
        this.render();
    }

    setTheme(bg, line) {
        this.renderer.setClearColor(new THREE.Color(bg));
        if (this.box) this.box.material.color.set(line);
        this.render();
    }

    resize() {
        const r = this.canvas.getBoundingClientRect();
        const w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));
        const s = this.renderer.getSize(new THREE.Vector2());
        if (s.x === w && s.y === h) return false;
        this.renderer.setSize(w, h, false);
        this.camera.aspect = w / h;
        this.camera.updateProjectionMatrix();
        this.render();
        return true;
    }

    resetCamera() { const n = this.n; this.n = 0; this.fit(n); this.render(); }

    render() { this.renderer.render(this.scene, this.camera); }
}
