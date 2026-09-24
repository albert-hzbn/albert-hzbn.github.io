// Small vector / 3x3 matrix helpers and seeded random numbers.
// Matrices are arrays of rows; cells are stored as rows [a, b, c].

export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const norm = (a) => Math.hypot(a[0], a[1], a[2]);
export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
export const unit = (a) => { const n = norm(a) || 1; return [a[0] / n, a[1] / n, a[2] / n]; };

export function det3(m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

export function inv3(m) {
    const d = det3(m);
    if (Math.abs(d) < 1e-12) throw new Error('Singular matrix');
    const r = [
        [m[1][1] * m[2][2] - m[1][2] * m[2][1], m[0][2] * m[2][1] - m[0][1] * m[2][2], m[0][1] * m[1][2] - m[0][2] * m[1][1]],
        [m[1][2] * m[2][0] - m[1][0] * m[2][2], m[0][0] * m[2][2] - m[0][2] * m[2][0], m[0][2] * m[1][0] - m[0][0] * m[1][2]],
        [m[1][0] * m[2][1] - m[1][1] * m[2][0], m[0][1] * m[2][0] - m[0][0] * m[2][1], m[0][0] * m[1][1] - m[0][1] * m[1][0]],
    ];
    return r.map((row) => row.map((v) => v / d));
}

export function matmul(a, b) {
    return a.map((row) => [0, 1, 2].map((j) => row[0] * b[0][j] + row[1] * b[1][j] + row[2] * b[2][j]));
}

export const transpose = (m) => [0, 1, 2].map((i) => [m[0][i], m[1][i], m[2][i]]);

// Row vector times matrix: v · M
export const vecmat = (v, m) => [
    v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
    v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
    v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
];

// Matrix times column vector: M · v
export const matvec = (m, v) => [dot(m[0], v), dot(m[1], v), dot(m[2], v)];

export function gcd(a, b) {
    a = Math.abs(a); b = Math.abs(b);
    while (b) [a, b] = [b, a % b];
    return a;
}

// Rotation matrix (acting on column vectors) about unit axis by angle (rad).
export function rotationMatrix(axis, angle) {
    const [x, y, z] = unit(axis), c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
    return [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ];
}

// Uniformly distributed random rotation (Shoemake's method).
export function randomRotation(rand = Math.random) {
    const u1 = rand(), u2 = rand() * 2 * Math.PI, u3 = rand() * 2 * Math.PI;
    const a = Math.sqrt(1 - u1), b = Math.sqrt(u1);
    const [w, x, y, z] = [a * Math.sin(u2), a * Math.cos(u2), b * Math.sin(u3), b * Math.cos(u3)];
    return [
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ];
}

// Small seeded PRNG so builders are reproducible.
export function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}
