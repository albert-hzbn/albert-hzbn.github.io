// three.js viewer: instanced atoms and bonds, cell, overlays, picking, gizmo.
import * as THREE from 'three';
import { TrackballControls } from 'jsm/controls/TrackballControls.js';

const UP = new THREE.Vector3(0, 1, 0);

export class Viewer {
    constructor(container) {
        this.container = container;
        this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true, alpha: false });
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        this.renderer.outputColorSpace = THREE.SRGBColorSpace;
        container.appendChild(this.renderer.domElement);
        this.renderer.domElement.classList.add('viewport-canvas');

        this.scene = new THREE.Scene();
        this.background = new THREE.Color('#f4f5f7');
        this.scene.background = this.background;

        this.ortho = new THREE.OrthographicCamera(-10, 10, 10, -10, -1000, 1000);
        this.persp = new THREE.PerspectiveCamera(35, 1, 0.1, 5000);
        this.camera = this.ortho;
        this.camera.position.set(0, 0, 50);

        // Lights: soft sky/ground fill plus a headlight that follows the camera.
        this.scene.add(new THREE.HemisphereLight(0xffffff, 0x8a8f99, 1.35));
        this.headlight = new THREE.DirectionalLight(0xffffff, 1.6);
        this.headlight.position.set(0.4, 0.6, 1);
        this.camera.add(this.headlight);
        this.scene.add(this.camera);

        this.controls = this.makeControls(this.camera);

        this.root = new THREE.Group();
        this.scene.add(this.root);
        this.overlays = new THREE.Group();
        this.scene.add(this.overlays);

        this.atomMeshes = [];      // [{ mesh, indices }]
        this.atomData = [];        // per displayed instance: { index, position, radius }
        this.labels = [];          // HTML labels anchored in 3D
        this.labelLayer = document.createElement('div');
        this.labelLayer.className = 'label-layer';
        container.appendChild(this.labelLayer);

        this.gizmo = document.createElement('canvas');
        this.gizmo.className = 'axes-gizmo';
        this.gizmo.width = this.gizmo.height = 180;
        container.appendChild(this.gizmo);
        this.gizmoAxes = null;     // { vectors: [[x,y,z]...], labels: [...] }

        this.raycaster = new THREE.Raycaster();
        this.pointer = new THREE.Vector2();
        this.autoRotate = false;

        this.resizeObserver = new ResizeObserver(() => this.resize());
        this.resizeObserver.observe(container);
        this.resize();
        this.animate = this.animate.bind(this);
        requestAnimationFrame(this.animate);
    }

    makeControls(camera) {
        const c = new TrackballControls(camera, this.renderer.domElement);
        c.rotateSpeed = 3.2;
        c.zoomSpeed = 1.4;
        c.panSpeed = 0.9;
        c.dynamicDampingFactor = 0.18;
        c.staticMoving = false;
        return c;
    }

    setTheme(dark) {
        this.dark = dark;
        this.background.set(dark ? '#16191d' : '#f4f5f7');
        if (this.cellLines) this.cellLines.material.color.set(dark ? '#9aa4ae' : '#39404a');
    }

    resize() {
        const w = this.container.clientWidth || 1, h = this.container.clientHeight || 1;
        this.renderer.setSize(w, h, false);
        this.aspect = w / h;
        this.updateOrtho();
        this.persp.aspect = this.aspect;
        this.persp.updateProjectionMatrix();
        this.controls.handleResize();
    }

    updateOrtho() {
        const half = this.viewHalf || 10;
        this.ortho.left = -half * this.aspect;
        this.ortho.right = half * this.aspect;
        this.ortho.top = half;
        this.ortho.bottom = -half;
        this.ortho.updateProjectionMatrix();
    }

    setProjection(kind) {
        const from = this.camera;
        const to = kind === 'perspective' ? this.persp : this.ortho;
        if (from === to) return;
        to.position.copy(from.position);
        to.up.copy(from.up);
        to.quaternion.copy(from.quaternion);
        from.remove(this.headlight);
        to.add(this.headlight);
        this.scene.remove(from);
        this.scene.add(to);
        const target = this.controls.target.clone();
        this.controls.dispose();
        this.camera = to;
        this.controls = this.makeControls(to);
        this.controls.target.copy(target);
        if (to === this.persp) {
            const dist = (this.viewHalf || 10) / Math.tan(THREE.MathUtils.degToRad(this.persp.fov / 2));
            const dir = new THREE.Vector3().subVectors(to.position, target).normalize();
            to.position.copy(target).addScaledVector(dir, dist);
        } else {
            to.zoom = 1;
            to.updateProjectionMatrix();
        }
    }

    // ---------- Scene content ----------

    clearGroup(g) {
        for (const child of [...g.children]) {
            g.remove(child);
            child.traverse?.((o) => {
                o.geometry?.dispose();
                if (o.material) (Array.isArray(o.material) ? o.material : [o.material]).forEach((m) => m.dispose());
            });
        }
    }

    // atoms: [{ index, position:[x,y,z], radius, color:'#rrggbb', ghost }]
    // bonds: [{ a:[x,y,z], b:[x,y,z], colorA, colorB, radius }]
    setContent({ atoms, bonds, cell, bondRadius = 0.12, quality = 'auto' }) {
        this.clearGroup(this.root);
        this.atomMeshes = [];
        this.atomData = atoms;
        const n = atoms.length;
        const seg = quality === 'auto' ? (n < 3000 ? [36, 22] : n < 30000 ? [18, 12] : [10, 7]) : [24, 16];
        const sphere = new THREE.SphereGeometry(1, seg[0], seg[1]);
        const material = new THREE.MeshStandardMaterial({ roughness: 0.42, metalness: 0.08 });
        if (n) {
            const mesh = new THREE.InstancedMesh(sphere, material, n);
            const m = new THREE.Matrix4(), q = new THREE.Quaternion(), c = new THREE.Color();
            atoms.forEach((at, k) => {
                m.compose(new THREE.Vector3(...at.position), q, new THREE.Vector3(at.radius, at.radius, at.radius));
                mesh.setMatrixAt(k, m);
                mesh.setColorAt(k, c.set(at.color));
            });
            mesh.instanceMatrix.needsUpdate = true;
            mesh.instanceColor.needsUpdate = true;
            mesh.userData.kind = 'atoms';
            this.root.add(mesh);
            this.atomMeshes.push(mesh);
        }
        if (bonds && bonds.length) {
            const cyl = new THREE.CylinderGeometry(1, 1, 1, n < 30000 ? 14 : 8, 1, true);
            const bmat = new THREE.MeshStandardMaterial({ roughness: 0.5, metalness: 0.05 });
            const mesh = new THREE.InstancedMesh(cyl, bmat, bonds.length * 2);
            const m = new THREE.Matrix4(), q = new THREE.Quaternion(), c = new THREE.Color();
            const A = new THREE.Vector3(), B = new THREE.Vector3(), M = new THREE.Vector3(), D = new THREE.Vector3();
            let k = 0;
            for (const bd of bonds) {
                A.set(...bd.a); B.set(...bd.b);
                M.addVectors(A, B).multiplyScalar(0.5);
                // Two half-bonds, each coloured by its atom.
                for (const [P, col] of [[A, bd.colorA], [B, bd.colorB]]) {
                    if (col === null) { m.makeScale(0, 0, 0); mesh.setMatrixAt(k++, m); continue; }
                    D.subVectors(M, P);
                    const len = D.length();
                    q.setFromUnitVectors(UP, D.clone().normalize());
                    const mid = P.clone().addScaledVector(D, 0.5);
                    const r = bd.radius ?? bondRadius;
                    m.compose(mid, q, new THREE.Vector3(r, len, r));
                    mesh.setMatrixAt(k, m);
                    mesh.setColorAt(k, c.set(col));
                    k++;
                }
            }
            mesh.count = k;
            mesh.instanceMatrix.needsUpdate = true;
            if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
            mesh.userData.kind = 'bonds';
            this.root.add(mesh);
        }
        this.cellLines = null;
        if (cell) {
            const o = new THREE.Vector3();
            const [a, b, c] = cell.map((v) => new THREE.Vector3(...v));
            const corners = [o, a, b, c, a.clone().add(b), a.clone().add(c), b.clone().add(c), a.clone().add(b).add(c)];
            const E = [[0, 1], [0, 2], [0, 3], [1, 4], [1, 5], [2, 4], [2, 6], [3, 5], [3, 6], [4, 7], [5, 7], [6, 7]];
            const pts = [];
            for (const [i, j] of E) pts.push(corners[i], corners[j]);
            const geo = new THREE.BufferGeometry().setFromPoints(pts);
            this.cellLines = new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color: this.dark ? '#9aa4ae' : '#39404a' }));
            this.root.add(this.cellLines);
        }
    }

    // Highlight spheres for selected / hovered atoms.
    setHighlights(items) {
        if (this.highlightMesh) {
            this.root.remove(this.highlightMesh);
            this.highlightMesh.geometry.dispose();
            this.highlightMesh.material.dispose();
            this.highlightMesh = null;
        }
        if (!items.length) return;
        const geo = new THREE.SphereGeometry(1, 24, 16);
        const mat = new THREE.MeshBasicMaterial({ color: '#2f7de1', transparent: true, opacity: 0.38, depthWrite: false });
        const mesh = new THREE.InstancedMesh(geo, mat, items.length);
        const m = new THREE.Matrix4(), q = new THREE.Quaternion(), c = new THREE.Color();
        items.forEach((it, k) => {
            const r = it.radius * 1.28 + 0.08;
            m.compose(new THREE.Vector3(...it.position), q, new THREE.Vector3(r, r, r));
            mesh.setMatrixAt(k, m);
            mesh.setColorAt(k, c.set(it.color || '#2f7de1'));
        });
        mesh.renderOrder = 2;
        this.highlightMesh = mesh;
        this.root.add(mesh);
    }

    // ---------- Overlays: measurement lines, planes, arrows, labels ----------

    clearOverlays() {
        this.clearGroup(this.overlays);
        this.labels.forEach((l) => l.el.remove());
        this.labels = [];
    }

    addLine(points, color = '#e0662f', dashed = true) {
        const geo = new THREE.BufferGeometry().setFromPoints(points.map((p) => new THREE.Vector3(...p)));
        const mat = dashed
            ? new THREE.LineDashedMaterial({ color, dashSize: 0.25, gapSize: 0.15, depthTest: false })
            : new THREE.LineBasicMaterial({ color, depthTest: false });
        const line = new THREE.Line(geo, mat);
        if (dashed) line.computeLineDistances();
        line.renderOrder = 5;
        this.overlays.add(line);
        return line;
    }

    addPolygon(points, color = '#2f7de1', opacity = 0.28) {
        if (points.length < 3) return;
        const pos = [];
        const p0 = points[0];
        for (let i = 1; i + 1 < points.length; i++) pos.push(...p0, ...points[i], ...points[i + 1]);
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
        geo.computeVertexNormals();
        const mesh = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color, transparent: true, opacity, side: THREE.DoubleSide, depthWrite: false }));
        mesh.renderOrder = 3;
        this.overlays.add(mesh);
        this.addLine([...points, points[0]], color, false);
    }

    addArrow(origin, dir, length, color = '#c2410c') {
        const arrow = new THREE.ArrowHelper(new THREE.Vector3(...dir).normalize(), new THREE.Vector3(...origin), length, color, Math.min(1.2, length * 0.18), Math.min(0.6, length * 0.08));
        this.overlays.add(arrow);
    }

    addLabel(position, text, cls = '') {
        const el = document.createElement('div');
        el.className = 'vlabel ' + cls;
        el.textContent = text;
        this.labelLayer.appendChild(el);
        this.labels.push({ el, position: new THREE.Vector3(...position) });
    }

    updateLabels() {
        const w = this.container.clientWidth, h = this.container.clientHeight;
        const v = new THREE.Vector3();
        for (const l of this.labels) {
            v.copy(l.position).project(this.camera);
            const visible = v.z < 1 && v.z > -1;
            l.el.style.display = visible ? '' : 'none';
            l.el.style.transform = `translate(${(v.x * 0.5 + 0.5) * w}px, ${(-v.y * 0.5 + 0.5) * h}px) translate(-50%, -130%)`;
        }
    }

    // ---------- Picking ----------

    pick(clientX, clientY) {
        const rect = this.renderer.domElement.getBoundingClientRect();
        this.pointer.set(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
        this.raycaster.setFromCamera(this.pointer, this.camera);
        const hits = this.raycaster.intersectObjects(this.atomMeshes, false);
        if (!hits.length) return null;
        const inst = hits[0].instanceId;
        return this.atomData[inst] ? { ...this.atomData[inst], instance: inst } : null;
    }

    // Screen-space rectangle selection: returns displayed atom indices inside.
    atomsInRect(x0, y0, x1, y1) {
        const rect = this.renderer.domElement.getBoundingClientRect();
        const [lx, hx] = [Math.min(x0, x1) - rect.left, Math.max(x0, x1) - rect.left];
        const [ly, hy] = [Math.min(y0, y1) - rect.top, Math.max(y0, y1) - rect.top];
        const v = new THREE.Vector3();
        const out = [];
        for (const at of this.atomData) {
            if (at.ghost) continue;
            v.set(...at.position).project(this.camera);
            const sx = (v.x * 0.5 + 0.5) * rect.width, sy = (-v.y * 0.5 + 0.5) * rect.height;
            if (sx >= lx && sx <= hx && sy >= ly && sy <= hy) out.push(at.index);
        }
        return out;
    }

    // ---------- Camera helpers ----------

    fit(bounds, keepDirection = true) {
        const lo = new THREE.Vector3(...bounds.lo), hi = new THREE.Vector3(...bounds.hi);
        const center = lo.clone().add(hi).multiplyScalar(0.5);
        const radius = Math.max(2, lo.distanceTo(hi) / 2 + 1.5);
        const dir = keepDirection
            ? new THREE.Vector3().subVectors(this.camera.position, this.controls.target).normalize()
            : new THREE.Vector3(0.9, 0.55, 1.2).normalize();
        if (!isFinite(dir.x) || dir.lengthSq() < 0.5) dir.set(0, 0, 1);
        this.controls.target.copy(center);
        this.viewHalf = radius * 1.08;
        this.updateOrtho();
        this.ortho.zoom = 1;
        this.ortho.updateProjectionMatrix();
        const dist = this.camera === this.persp ? this.viewHalf / Math.tan(THREE.MathUtils.degToRad(this.persp.fov / 2)) : radius * 4;
        this.camera.position.copy(center).addScaledVector(dir, dist);
        this.camera.lookAt(center);
        this.controls.update();
    }

    viewAlong(direction, up = [0, 1, 0]) {
        const d = new THREE.Vector3(...direction).normalize();
        let u = new THREE.Vector3(...up);
        if (Math.abs(u.clone().normalize().dot(d)) > 0.99) u = Math.abs(d.z) < 0.9 ? new THREE.Vector3(0, 0, 1) : new THREE.Vector3(0, 1, 0);
        const target = this.controls.target.clone();
        const dist = this.camera.position.distanceTo(target) || 50;
        this.camera.position.copy(target).addScaledVector(d, dist);
        this.camera.up.copy(u);
        this.camera.lookAt(target);
        this.controls.update();
    }

    setGizmoAxes(vectors, labels, colors) {
        this.gizmoAxes = { vectors, labels, colors };
    }

    drawGizmo() {
        const g = this.gizmo, ctx = g.getContext('2d');
        const size = g.width;
        ctx.clearRect(0, 0, size, size);
        if (!this.gizmoAxes) return;
        const q = this.camera.quaternion.clone().invert();
        const cx = size / 2, cy = size / 2, L = size * 0.32;
        const items = this.gizmoAxes.vectors.map((vec, i) => {
            const v = new THREE.Vector3(...vec).normalize().applyQuaternion(q);
            return { v, label: this.gizmoAxes.labels[i], color: this.gizmoAxes.colors[i] };
        }).sort((a, b) => a.v.z - b.v.z);
        ctx.lineCap = 'round';
        ctx.font = `600 ${size * 0.12}px Inter, system-ui, sans-serif`;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        for (const it of items) {
            const x = cx + it.v.x * L, y = cy - it.v.y * L;
            ctx.globalAlpha = it.v.z < -0.2 ? 0.45 : 1;
            ctx.strokeStyle = it.color;
            ctx.lineWidth = size * 0.028;
            ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(x, y); ctx.stroke();
            ctx.fillStyle = it.color;
            ctx.beginPath(); ctx.arc(x, y, size * 0.075, 0, Math.PI * 2); ctx.fill();
            ctx.fillStyle = '#fff';
            ctx.fillText(it.label, x, y + 1);
        }
        ctx.globalAlpha = 1;
    }

    screenshot(scale = 2, transparent = false) {
        const r = this.renderer;
        const size = r.getSize(new THREE.Vector2());
        const pr = r.getPixelRatio();
        const bg = this.scene.background;
        if (transparent) { this.scene.background = null; r.setClearColor(0x000000, 0); }
        r.setPixelRatio(scale);
        r.setSize(size.x, size.y, false);
        r.render(this.scene, this.camera);
        const url = r.domElement.toDataURL('image/png');
        r.setPixelRatio(pr);
        r.setSize(size.x, size.y, false);
        this.scene.background = bg;
        r.setClearColor(0x000000, 1);
        return url;
    }

    animate() {
        requestAnimationFrame(this.animate);
        if (this.autoRotate) {
            const axis = this.camera.up.clone().normalize();
            const offset = new THREE.Vector3().subVectors(this.camera.position, this.controls.target);
            offset.applyAxisAngle(axis, 0.004);
            this.camera.position.copy(this.controls.target).add(offset);
            this.camera.lookAt(this.controls.target);
        }
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
        this.updateLabels();
        this.drawGizmo();
    }
}
