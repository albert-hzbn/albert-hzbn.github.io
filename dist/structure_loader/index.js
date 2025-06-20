import * as THREE from "three";
import { OrbitControls } from "jsm/controls/OrbitControls.js";
import { CSS2DRenderer, CSS2DObject } from "jsm/renderers/CSS2DRenderer.js";

const atomsInfo = {
  H: { radius: 0.31, color: 0xffffff },
  He: { radius: 0.28, color: 0x00ffff },
  Li: { radius: 1.28, color: 0xcc80ff },
  Be: { radius: 0.96, color: 0x2a802a },
  B: { radius: 0.84, color: 0x00ff00 },
  C: { radius: 0.76, color: 0xc8c8c8 },
  N: { radius: 0.71, color: 0x8f8fff },
  O: { radius: 0.66, color: 0xf00000 },
  F: { radius: 0.57, color: 0x00ff00 },
  Ne: { radius: 0.58, color: 0x00ffff },

  Na: { radius: 1.66, color: 0x0000ff },
  Mg: { radius: 1.41, color: 0x2a802a },
  Al: { radius: 1.21, color: 0xc8c8c8 },
  Si: { radius: 1.11, color: 0xf0c8a0 },
  P: { radius: 1.07, color: 0xffa500 },
  S: { radius: 1.05, color: 0xffc832 },
  Cl: { radius: 1.02, color: 0x00ff00 },
  Ar: { radius: 1.06, color: 0x00ffff },

  K: { radius: 2.03, color: 0xcc80ff },
  Ca: { radius: 1.76, color: 0x808090 },
  Sc: { radius: 1.7, color: 0xc8c8c8 },
  Ti: { radius: 1.6, color: 0xc8c8c8 },
  V: { radius: 1.53, color: 0xc8c8c8 },
  Cr: { radius: 1.39, color: 0xc8c8c8 },
  Mn: { radius: 1.5, color: 0xc8c8c8 },
  Fe: { radius: 1.42, color: 0xffa500 },
  Co: { radius: 1.38, color: 0xc8c8c8 },
  Ni: { radius: 1.24, color: 0xc8c8c8 },
  Cu: { radius: 1.32, color: 0xca7833 },
  Zn: { radius: 1.22, color: 0xc8c8c8 },
  Ga: { radius: 1.22, color: 0xc8c8c8 },
  Ge: { radius: 1.2, color: 0xc8c8c8 },
  As: { radius: 1.19, color: 0xc8c8c8 },
  Se: { radius: 1.2, color: 0xc8c8c8 },
  Br: { radius: 1.2, color: 0x8f40d4 },
  Kr: { radius: 1.16, color: 0x00ffff },

  Rb: { radius: 2.2, color: 0xcc80ff },
  Sr: { radius: 1.95, color: 0x2a802a },
  Y: { radius: 1.9, color: 0xc8c8c8 },
  Zr: { radius: 1.75, color: 0xc8c8c8 },
  Nb: { radius: 1.64, color: 0xc8c8c8 },
  Mo: { radius: 1.54, color: 0xc8c8c8 },
  Tc: { radius: 1.47, color: 0xc8c8c8 },
  Ru: { radius: 1.46, color: 0xc8c8c8 },
  Rh: { radius: 1.42, color: 0xc8c8c8 },
  Pd: { radius: 1.39, color: 0xc8c8c8 },
  Ag: { radius: 1.45, color: 0xc8c8c8 },
  Cd: { radius: 1.44, color: 0xc8c8c8 },
  In: { radius: 1.42, color: 0xc8c8c8 },
  Sn: { radius: 1.39, color: 0xc8c8c8 },
  Sb: { radius: 1.39, color: 0xc8c8c8 },
  Te: { radius: 1.38, color: 0xc8c8c8 },
  I: { radius: 1.39, color: 0x940094 },
  Xe: { radius: 1.4, color: 0x00ffff },

  Cs: { radius: 2.44, color: 0xcc80ff },
  Ba: { radius: 2.15, color: 0x2a802a },

  La: { radius: 2.07, color: 0xc8c8c8 },
  Ce: { radius: 2.04, color: 0xc8c8c8 },
  Pr: { radius: 2.03, color: 0xc8c8c8 },
  Nd: { radius: 2.01, color: 0xc8c8c8 },
  Pm: { radius: 1.99, color: 0xc8c8c8 },
  Sm: { radius: 1.98, color: 0xc8c8c8 },
  Eu: { radius: 1.98, color: 0xc8c8c8 },
  Gd: { radius: 1.96, color: 0xc8c8c8 },
  Tb: { radius: 1.94, color: 0xc8c8c8 },
  Dy: { radius: 1.92, color: 0xc8c8c8 },
  Ho: { radius: 1.92, color: 0xc8c8c8 },
  Er: { radius: 1.89, color: 0xc8c8c8 },
  Tm: { radius: 1.9, color: 0xc8c8c8 },
  Yb: { radius: 1.87, color: 0xc8c8c8 },
  Lu: { radius: 1.87, color: 0xc8c8c8 },

  Hf: { radius: 1.75, color: 0xc8c8c8 },
  Ta: { radius: 1.7, color: 0xc8c8c8 },
  W: { radius: 1.62, color: 0xc8c8c8 },
  Re: { radius: 1.51, color: 0xc8c8c8 },
  Os: { radius: 1.44, color: 0xc8c8c8 },
  Ir: { radius: 1.41, color: 0xc8c8c8 },
  Pt: { radius: 1.36, color: 0xc8c8c8 },
  Au: { radius: 1.36, color: 0xc8c8c8 },
  Hg: { radius: 1.32, color: 0xc8c8c8 },

  Tl: { radius: 1.45, color: 0xc8c8c8 },
  Pb: { radius: 1.46, color: 0xc8c8c8 },
  Bi: { radius: 1.48, color: 0xc8c8c8 },
  Po: { radius: 1.4, color: 0xc8c8c8 },
  At: { radius: 1.5, color: 0x8f40d4 },
  Rn: { radius: 1.5, color: 0x00ffff },

  Fr: { radius: 2.6, color: 0xcc80ff },
  Ra: { radius: 2.21, color: 0x2a802a },

  Ac: { radius: 2.15, color: 0xc8c8c8 },
  Th: { radius: 2.06, color: 0xc8c8c8 },
  Pa: { radius: 2.0, color: 0xc8c8c8 },
  U: { radius: 1.96, color: 0xc8c8c8 },
  Np: { radius: 1.9, color: 0xc8c8c8 },
  Pu: { radius: 1.87, color: 0xc8c8c8 },
  Am: { radius: 1.8, color: 0xc8c8c8 },
  Cm: { radius: 1.76, color: 0xc8c8c8 },
  Bk: { radius: 1.68, color: 0xc8c8c8 },
  Cf: { radius: 1.68, color: 0xc8c8c8 },
  Es: { radius: 1.65, color: 0xc8c8c8 },
  Fm: { radius: 1.67, color: 0xc8c8c8 },
  Md: { radius: 1.73, color: 0xc8c8c8 },
  No: { radius: 1.76, color: 0xc8c8c8 },
  Lr: { radius: 1.61, color: 0xc8c8c8 },
  Rf: { radius: 1.57, color: 0xc8c8c8 },
  Db: { radius: 1.49, color: 0xc8c8c8 },
  Sg: { radius: 1.43, color: 0xc8c8c8 },
  Bh: { radius: 1.41, color: 0xc8c8c8 },
  Hs: { radius: 1.34, color: 0xc8c8c8 },
  Mt: { radius: 1.29, color: 0xc8c8c8 },
  Ds: { radius: 1.28, color: 0xc8c8c8 },
  Rg: { radius: 1.21, color: 0xc8c8c8 },
  Cn: { radius: 1.22, color: 0xc8c8c8 },
  Nh: { radius: NaN, color: 0xc8c8c8 },
  Fl: { radius: NaN, color: 0xc8c8c8 },
  Mc: { radius: NaN, color: 0xc8c8c8 },
  Lv: { radius: NaN, color: 0xc8c8c8 },
  Ts: { radius: NaN, color: 0xc8c8c8 },
  Og: { radius: NaN, color: 0xc8c8c8 },
};

const bondInfo = {
  O: ["Al", "Ca", "Na", "Si"],
  Al: ["O"],
  Ca: ["0"],
  Na: ["0", "Cl"],
  Si: ["O"],
  C: ["C", "O", "H"],
  O: ["C", "O", "O"],
  H: ["C", "O", "H"],
  Na: ["Cl"],
};

const cutoffDistance = 1.6;

const elementIdToElement = {
  1: "Ca",
  2: "O",
  3: "Al",
  4: "Na",
  5: "Si",
};
let latticeVectors = [];
let allAtoms = [];

const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2(1, 1);
const color = new THREE.Color();
let atomsAdded = false;

// Scene

const scene = new THREE.Scene();

const canvas = document.getElementById("webgl-canvas");
const renderer = new THREE.WebGLRenderer({ canvas });
renderer.setSize(canvas.clientWidth, canvas.clientHeight);
renderer.setPixelRatio(window.devicePixelRatio);
const aspect = canvas.clientWidth / canvas.clientHeight;
const frustumSize = 20; // increase this from 10 to make a bigger view box

const camera = new THREE.OrthographicCamera(
  (-frustumSize * aspect) / 2,
  (frustumSize * aspect) / 2,
  frustumSize / 2,
  -frustumSize / 2,
  -100, // set near to a more negative value
  1000 // and far to a larger value
);
camera.position.z = 100;

// Lights
const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
directionalLight.position.set(10, 10, 10);
scene.add(directionalLight);

const ambient = new THREE.AmbientLight(0xffffff, 1); // Soft fill
scene.add(ambient);

// Controls
const controls = new OrbitControls(camera, renderer.domElement);

let atoms = [];

function parseXYZ(data) {
  const lines = data.trim().split("\n");
  const numAtoms = parseInt(lines[0]);
  const parsedAtoms = [];

  for (let i = 0; i < numAtoms; i++) {
    const line = lines[i + 2]; // Skip count and comment line
    const [element, x, y, z] = line.trim().split(/\s+/);
    parsedAtoms.push({
      element,
      x: parseFloat(x),
      y: parseFloat(y),
      z: parseFloat(z),
    });
  }

  return parsedAtoms;
}

function addOverlayScene() {
  // Overlay Scene and Camera
  const overlayScene = new THREE.Scene();
  const overlayCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 1, 10);
  overlayCamera.position.z = 5;

  // Create Axes using ArrowHelpers
  const axisLength = 0.5;
  const origin = new THREE.Vector3(0, 0, 0);

  const arrowX = new THREE.ArrowHelper(
    new THREE.Vector3(1, 0, 0),
    origin,
    axisLength,
    0xff0000
  );
  const arrowY = new THREE.ArrowHelper(
    new THREE.Vector3(0, 1, 0),
    origin,
    axisLength,
    0x00ff00
  );
  const arrowZ = new THREE.ArrowHelper(
    new THREE.Vector3(0, 0, 1),
    origin,
    axisLength,
    0x0000ff
  );

  overlayScene.add(arrowX, arrowY, arrowZ);

  // Create a small viewport in the corner (optional)
  const overlayViewport = {
    x: 10,
    y: 10,
    width: 150,
    height: 150,
  };
}
addOverlayScene();

function parseLammpsDump(data) {
  const lines = data.trim().split("\n");
  const numAtoms = parseInt(lines[3]);
  const parsedAtoms = [];

  const [xLow, xHigh] = lines[5].trim().split(/\s+/);
  const [yLow, yHigh] = lines[6].trim().split(/\s+/);
  const [zLow, zHigh] = lines[7].trim().split(/\s+/);
  const latticeVectors = [
    [xLow, xHigh],
    [yLow, yHigh],
    [zLow, zHigh],
  ];

  for (let i = 0; i < numAtoms; i++) {
    const line = lines[i + 9]; // Skip count and comment line
    const [atomId, elementId, x, y, z, vx, vy, vz] = line.trim().split(/\s+/);
    const element = elementIdToElement[elementId];
    parsedAtoms.push({
      element,
      x: parseFloat(x),
      y: parseFloat(y),
      z: parseFloat(z),
    });
  }
  return [latticeVectors, parsedAtoms];
}

function parsePOSCAR(content) {
  const lines = content
    .trim()
    .split(/\r?\n/)
    .map((line) => line.trim());

  const scale = parseFloat(lines[1]);
  const latticeVectors = [
    lines[2].split(/\s+/).map(Number),
    lines[3].split(/\s+/).map(Number),
    lines[4].split(/\s+/).map(Number),
  ].map((vec) => vec.map((val) => val * scale));

  // Detect elements and counts
  let elementLineIndex = 5;
  let elements = lines[elementLineIndex].split(/\s+/);
  let counts = lines[elementLineIndex + 1].split(/\s+/).map(Number);
  if (counts.some(isNaN)) {
    // No element names, only counts
    elements = [];
    counts = lines[elementLineIndex].split(/\s+/).map(Number);
    elementLineIndex--;
  }

  const totalAtoms = counts.reduce((a, b) => a + b, 0);
  const coordTypeLine = lines[elementLineIndex + 2].toLowerCase();
  const isCartesian = coordTypeLine.startsWith("c");
  const atomStartLine = elementLineIndex + 3;

  const atoms = [];
  let atomIndex = 0;

  for (let i = 0; i < elements.length; i++) {
    const element = elements[i];
    const count = counts[i];
    for (let j = 0; j < count; j++) {
      const line = lines[atomStartLine + atomIndex];
      const [x, y, z] = line.split(/\s+/).slice(0, 3).map(Number);

      let cartesian;
      if (isCartesian) {
        cartesian = [x * scale, y * scale, z * scale];
      } else {
        // Direct to cartesian using scaled lattice vectors
        cartesian = [
          x * latticeVectors[0][0] +
            y * latticeVectors[1][0] +
            z * latticeVectors[2][0],
          x * latticeVectors[0][1] +
            y * latticeVectors[1][1] +
            z * latticeVectors[2][1],
          x * latticeVectors[0][2] +
            y * latticeVectors[1][2] +
            z * latticeVectors[2][2],
        ];
      }

      atoms.push({
        element: element || "X",
        x: cartesian[0],
        y: cartesian[1],
        z: cartesian[2],
      });

      atomIndex++;
    }
  }

  return [latticeVectors, atoms];
}

function addPeriodicImages(atoms, latticeVectors) {
  // Start with an empty array for the result
  const periodicAtoms = [];

  // We want a 3x3x3 grid of cells centered on the original cell
  // This creates 27 cells - the original cell plus 26 surrounding cells
  for (let na = -1; na <= 1; na++) {
    for (let nb = -1; nb <= 1; nb++) {
      for (let nc = -1; nc <= 1; nc++) {
        // Calculate the translation vector for this cell
        const dx =
          na * latticeVectors[0][0] +
          nb * latticeVectors[1][0] +
          nc * latticeVectors[2][0];
        const dy =
          na * latticeVectors[0][1] +
          nb * latticeVectors[1][1] +
          nc * latticeVectors[2][1];
        const dz =
          na * latticeVectors[0][2] +
          nb * latticeVectors[1][2] +
          nc * latticeVectors[2][2];

        // Add all atoms with this translation
        for (const atom of atoms) {
          periodicAtoms.push({
            element: atom.element,
            x: atom.x + dx,
            y: atom.y + dy,
            z: atom.z + dz,
          });
        }
      }
    }
  }

  // Filter out atoms outside the box
  // The simulation box is defined by the lattice vectors
  const filteredAtoms = filterAtomsInBox(periodicAtoms, latticeVectors);

  return filteredAtoms;
}

// Function to check if an atom is inside the simulation box
function filterAtomsInBox(atoms, latticeVectors) {
  // Calculate the inverse of the lattice matrix to convert Cartesian coordinates to fractional coordinates
  const latticeMatrix = [
    [latticeVectors[0][0], latticeVectors[1][0], latticeVectors[2][0]],
    [latticeVectors[0][1], latticeVectors[1][1], latticeVectors[2][1]],
    [latticeVectors[0][2], latticeVectors[1][2], latticeVectors[2][2]],
  ];

  const inverseMatrix = calculateInverseMatrix(latticeMatrix);

  // Filter atoms to keep only those inside the box (fractional coordinates between 0 and 1)
  return atoms.filter((atom) => {
    // Convert to fractional coordinates
    const fracCoords = [
      inverseMatrix[0][0] * atom.x +
        inverseMatrix[0][1] * atom.y +
        inverseMatrix[0][2] * atom.z,
      inverseMatrix[1][0] * atom.x +
        inverseMatrix[1][1] * atom.y +
        inverseMatrix[1][2] * atom.z,
      inverseMatrix[2][0] * atom.x +
        inverseMatrix[2][1] * atom.y +
        inverseMatrix[2][2] * atom.z,
    ];

    // Check if the atom is inside the box (fractional coordinates between 0 and 1)
    // Using a small epsilon value to include atoms exactly on the boundary
    const epsilon = 1e-10;
    return (
      fracCoords[0] >= -epsilon &&
      fracCoords[0] < 1 + epsilon &&
      fracCoords[1] >= -epsilon &&
      fracCoords[1] < 1 + epsilon &&
      fracCoords[2] >= -epsilon &&
      fracCoords[2] < 1 + epsilon
    );
  });
}

// Function to calculate the inverse of a 3x3 matrix
function calculateInverseMatrix(matrix) {
  // Calculate the determinant
  const det =
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
    matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
    matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

  // Check if the determinant is zero (or very close to zero)
  if (Math.abs(det) < 1e-10) {
    throw new Error("Matrix is singular and cannot be inverted");
  }

  // Calculate the adjugate matrix
  const adjugate = [
    [
      matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1],
      matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2],
      matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1],
    ],
    [
      matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2],
      matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0],
      matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2],
    ],
    [
      matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0],
      matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1],
      matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0],
    ],
  ];

  // Calculate the inverse matrix
  const inverse = adjugate.map((row) => row.map((value) => value / det));

  return inverse;
}

function createCylinderBetweenPoints(
  point1,
  point2,
  radius = 0.1,
  material = new THREE.MeshNormalMaterial()
) {
  const start = new THREE.Vector3(...point1);
  const end = new THREE.Vector3(...point2);

  // Calculate direction vector
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();

  // Create cylinder geometry aligned along Y-axis by default
  const geometry = new THREE.CylinderGeometry(radius, radius, length, 6);

  // Create mesh
  const cylinder = new THREE.Mesh(geometry, material);

  // Compute midpoint and position the cylinder
  const midpoint = new THREE.Vector3()
    .addVectors(start, end)
    .multiplyScalar(0.5);
  cylinder.position.copy(midpoint);

  // Orient the cylinder to align with direction
  const up = new THREE.Vector3(0, 1, 0); // default cylinder axis in Three.js
  const quaternion = new THREE.Quaternion().setFromUnitVectors(
    up,
    direction.clone().normalize()
  );
  cylinder.setRotationFromQuaternion(quaternion);

  return cylinder;
}

function getPairs(atoms) {
  const numAtoms = atoms.length;
  const bondPairs = [];
  const cutoffSq = cutoffDistance * cutoffDistance;

  for (let ind1 = 0; ind1 < numAtoms; ind1++) {
    const atom1 = atoms[ind1];
    const { element: element, x: x1, y: y1, z: z1 } = atom1;
    const elem1 = element;
    const validBonds1 = bondInfo[elem1];

    if (!validBonds1) continue;

    for (let ind2 = ind1 + 1; ind2 < numAtoms; ind2++) {
      const atom2 = atoms[ind2];
      const { element: element, x: x2, y: y2, z: z2 } = atom2;
      const elem2 = element;
      const validBonds2 = bondInfo[elem2];
      // console.log(elem1, elem2);
      // if (!validBonds2) continue;

      // Check if bonding is allowed
      if (validBonds1.includes(elem2) || validBonds2.includes(elem1)) {
        // Compute squared distance
        const dx = x2 - x1;
        const dy = y2 - y1;
        const dz = z2 - z1;
        const distSq = dx * dx + dy * dy + dz * dz;

        if (distSq < cutoffSq) {
          console.log(distSq);
          bondPairs.push([
            [x1, y1, z1],
            [x2, y2, z2],
          ]);
        }
      }
    }
  }

  return bondPairs;
}

function addBonds(atoms) {
  const bondPairs = getPairs(atoms);
  const numPairs = bondPairs.length;
  console.log(bondPairs);
  for (let ind = 0; ind < numPairs; ind++) {
    const point1 = bondPairs[ind][0];
    const point2 = bondPairs[ind][1];
    const material = new THREE.MeshStandardMaterial({ color: 0xff0000 });
    const cylinder = createCylinderBetweenPoints(
      point1,
      point2,
      0.05,
      material
    );
    scene.add(cylinder);
  }
}

function addBoundingBox(latticeVectors) {
  const material = new THREE.LineBasicMaterial({
    color: 0x000000,
    linewidth: 1.5,
  });

  // latticeVectors = [a, b, c] where each is [x, y, z]
  const a = latticeVectors[0];
  const b = latticeVectors[1];
  const c = latticeVectors[2];

  // Generate all 8 corners of the unit cell
  const v000 = new THREE.Vector3(0, 0, 0);
  const v100 = new THREE.Vector3(...a);
  const v010 = new THREE.Vector3(...b);
  const v110 = new THREE.Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
  const v001 = new THREE.Vector3(...c);
  const v101 = new THREE.Vector3(a[0] + c[0], a[1] + c[1], a[2] + c[2]);
  const v011 = new THREE.Vector3(b[0] + c[0], b[1] + c[1], b[2] + c[2]);
  const v111 = new THREE.Vector3(
    a[0] + b[0] + c[0],
    a[1] + b[1] + c[1],
    a[2] + b[2] + c[2]
  );

  const cubeVertices = [v000, v100, v010, v110, v001, v101, v011, v111];

  // Define edges between vertex indices (generalized box)
  const cubeEdges = [
    [0, 1],
    [0, 2],
    [0, 4],
    [1, 3],
    [1, 5],
    [2, 3],
    [2, 6],
    [3, 7],
    [4, 5],
    [4, 6],
    [5, 7],
    [6, 7],
  ];

  // Draw the edges
  for (const [start, end] of cubeEdges) {
    const geometry = new THREE.BufferGeometry().setFromPoints([
      cubeVertices[start],
      cubeVertices[end],
    ]);
    const line = new THREE.Line(geometry, material);
    scene.add(line);
  }

  createLatticeVectorLines(latticeVectors);
}

const labelRenderer = new CSS2DRenderer();
labelRenderer.setSize(canvas.clientWidth, canvas.clientHeight);
labelRenderer.domElement.style.position = "absolute";
labelRenderer.domElement.style.top = "0px";
labelRenderer.domElement.style.pointerEvents = "none";
document.body.appendChild(labelRenderer.domElement);
canvas.labelRenderer = labelRenderer;
canvas.labelRendererInitialized = true;

function createLatticeVectorLines(latticeVectors) {
  const colors = [0xff0000, 0x3cb371, 0x0000ff];
  const textContents = ["a", "b", "c"];

  for (let ind = 0; ind < latticeVectors.length; ind++) {
    const latticeVector = new THREE.Vector3(...latticeVectors[ind]);

    // Create the line
    const points = [new THREE.Vector3(0, 0, 0), latticeVector];
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({
      color: colors[ind],
      linewidth: 2,
    });
    const line = new THREE.Line(geometry, material);
    scene.add(line);

    // Create the label
    const div = document.createElement("div");
    div.className = "label";
    div.textContent = textContents[ind];
    div.style.marginTop = "-1em";
    div.style.color = "#" + colors[ind].toString(16).padStart(6, "0");
    div.style.fontSize = "20px";
    div.style.padding = "2px 2px";
    div.style.borderRadius = "1px";

    const label = new CSS2DObject(div);
    label.position.copy(latticeVector.clone().multiplyScalar(0.2));
    scene.add(label);
  }
}

function addBoundingBoxLammps(latticeVectors) {
  const material = new THREE.LineBasicMaterial({
    color: 0x000000,
  });

  // Create the 8 vertices of the cube using lattice vectors
  // assuming the xlow, ylow and zlow to be origin
  const x1 = latticeVectors[0][1]; // Assuming [0, x1]
  const y1 = latticeVectors[1][1]; // Assuming [0, y1]
  const z1 = latticeVectors[2][1]; // Assuming [0, z1]

  const cubeVertices = [
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(x1, 0, 0),
    new THREE.Vector3(x1, y1, 0),
    new THREE.Vector3(0, y1, 0),
    new THREE.Vector3(0, 0, z1),
    new THREE.Vector3(x1, 0, z1),
    new THREE.Vector3(x1, y1, z1),
    new THREE.Vector3(0, y1, z1),
  ];

  // Define the 12 edges by connecting the vertices
  const cubeEdges = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0], // Bottom face
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4], // Top face
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7], // Vertical edges
  ];

  // Add lines for each edge
  for (const [start, end] of cubeEdges) {
    const edgePoints = [cubeVertices[start], cubeVertices[end]];
    const geometry = new THREE.BufferGeometry().setFromPoints(edgePoints);
    const line = new THREE.Line(geometry, material);
    scene.add(line);
  }
}

function addAtoms(atoms) {
  const numAtoms = atoms.length;
  for (let ind = 0; ind < numAtoms; ind += 1) {
    const element = atoms[ind].element;
    const radius = atomsInfo[element].radius / 2;
    const color = atomsInfo[element].color;
    const geometry = new THREE.SphereGeometry(radius, 16, 16);
    const material = new THREE.MeshPhongMaterial({
      color: color,
      shininess: 200,
    });
    const atomMesh = new THREE.Mesh(geometry, material);
    atomMesh.name = element.concat("_", ind);
    atomMesh.position.set(atoms[ind].x, atoms[ind].y, atoms[ind].z);
    scene.add(atomMesh);
  }
  atomsAdded = true;
}

function removeObjects() {
  for (var i = scene.children.length - 1; i >= 0; i--) {
    const obj = scene.children[i];
    if (!(obj instanceof THREE.Light)) {
      scene.remove(obj);
    }
  }
}

// Make it globally available
window.transformStructure = transformStructure;

function transformStructure() {
  // Obtain transformation matrix from textbox
  const ids = [
    ["a1", "b1", "c1"],
    ["a2", "b2", "c2"],
    ["a3", "b3", "c3"],
  ];
  const transformationMatrix = ids.map((row) =>
    row.map((id) => parseFloat(document.getElementById(id).value) || 0)
  );
  console.log("Transformation matrix:", transformationMatrix);

  // Helper: matrix multiplication for 3D vectors
  const applyMatrixToVector = (vec, matrix) => {
    return {
      x: vec.x * matrix[0][0] + vec.y * matrix[0][1] + vec.z * matrix[0][2],
      y: vec.x * matrix[1][0] + vec.y * matrix[1][1] + vec.z * matrix[1][2],
      z: vec.x * matrix[2][0] + vec.y * matrix[2][1] + vec.z * matrix[2][2],
    };
  };

  // Helper: matrix transpose
  const transpose = (M) => M[0].map((_, i) => M.map((row) => row[i]));

  // Transpose the transformation matrix (consistent with right multiplication)
  const T = transpose(transformationMatrix);

  // Create a deep copy of original lattice vectors to ensure we don't modify the orginal variable
  const originalLatticeVectors = JSON.parse(JSON.stringify(latticeVectors));

  // Now create the transformed lattice vectors
  const transformedLatticeVectors = originalLatticeVectors.map((vec) => {
    const transformed = applyMatrixToVector(
      { x: vec[0], y: vec[1], z: vec[2] },
      T
    );
    return [transformed.x, transformed.y, transformed.z];
  });

  // Replicate the atoms within a range from -range to + range
  const range = 5; // Try replicating in a 5x5x5 block
  let replicatedAtoms = [];
  for (let i = -range; i <= range; i++) {
    for (let j = -range; j <= range; j++) {
      for (let k = -range; k <= range; k++) {
        allAtoms.forEach((atom) => {
          const pos = {
            x:
              atom.x +
              (i * latticeVectors[0][0] +
                j * latticeVectors[1][0] +
                k * latticeVectors[2][0]),
            y:
              atom.y +
              (i * latticeVectors[0][1] +
                j * latticeVectors[1][1] +
                k * latticeVectors[2][1]),
            z:
              atom.z +
              (i * latticeVectors[0][2] +
                j * latticeVectors[1][2] +
                k * latticeVectors[2][2]),
          };
          replicatedAtoms.push({
            element: atom.element,
            x: pos.x,
            y: pos.y,
            z: pos.z,
          });
        });
      }
    }
  }

  // Clean up existing objects
  removeObjects();

  // Combine original atoms with replicated ones
  const atomsWithReplicated = atoms.concat(replicatedAtoms);

  // Filter atoms to keep only those within the transformed unit cell
  const filteredAtoms = filterAtomsInBox(
    atomsWithReplicated,
    transformedLatticeVectors
  );

  // Render the new structure
  addBoundingBox(transformedLatticeVectors);
  addAtoms(filteredAtoms);
}

document
  .getElementById("fileInput")
  .addEventListener("change", function (event) {
    console.log("test");
    const file = event.target.files[0];
    const reader = new FileReader();

    // Get file name and extension
    const fileName = file.name;
    const fileExtension = fileName.includes(".")
      ? fileName.split(".").pop().toLowerCase()
      : "";
    const isPOSCAR = fileName.startsWith("POSCAR") && !fileName.includes(".");

    console.log("File name:", fileName);
    console.log("File extension:", fileExtension);
    console.log("Is POSCAR format:", isPOSCAR);

    reader.onload = function (e) {
      removeObjects();
      const content = e.target.result;
      let atoms = [];

      // Check the type of file loaded
      if (isPOSCAR) {
        [latticeVectors, atoms] = parsePOSCAR(content);
        const extraAtoms = addPeriodicImages(atoms, latticeVectors);
        allAtoms = atoms.concat(extraAtoms);
        console.log(allAtoms);
        addBoundingBox(latticeVectors);
      } else if (fileExtension == "dump") {
        [latticeVectors, atoms] = parseLammpsDump(content);
        const extraAtoms = addPeriodicImages(atoms, latticeVectors);
        allAtoms = atoms.concat(extraAtoms);
        addBoundingBoxLammps(latticeVectors);
      } else if (fileExtension == "xyz") {
        atoms = parseXYZ(content);
        allAtoms = atoms;
      }

      console.log(allAtoms);
      console.log(latticeVectors);
      addAtoms(allAtoms);
      addBonds(allAtoms);
    };

    reader.readAsText(file);
  });

function onMouseMove(event) {
  event.preventDefault();
  mouse.x = (event.clientX / canvas.clientWidth) * 2 - 1;
  mouse.y = -(event.clientY / canvas.clientHeight) * 2 + 1;
}

document.addEventListener("mousemove", onMouseMove);

let previouslyHovered = null;

function hoverAtom() {
  if (!atomsAdded) return;

  raycaster.setFromCamera(mouse, camera);

  const meshObjects = [];
  scene.traverse((obj) => {
    if (obj instanceof THREE.Mesh) {
      meshObjects.push(obj);
    }
  });

  const intersections = raycaster.intersectObjects(meshObjects, true);
  const infoBox = document.getElementById("atom-info");

  if (intersections.length > 0) {
    const hovered = intersections[0].object;

    if (previouslyHovered !== hovered) {
      const pos = hovered.position;
      const name = hovered.name || hovered.userData?.element || "Unknown";
      const element = name.split("_")[0];
      infoBox.innerHTML = `${element}<br>X: ${pos.x.toFixed(
        5
      )}<br>Y: ${pos.y.toFixed(5)}<br>Z: ${pos.z.toFixed(5)}`;
      infoBox.style.display = "block";

      previouslyHovered = hovered;
    }
  } else {
    if (infoBox) infoBox.style.display = "none";
    previouslyHovered = null;
  }
}

function highlightAtom(mesh) {
  const baseColor = mesh.material.color.clone(); // get the current color
  mesh.material.emissive = baseColor; // use it as the glow color
  mesh.material.emissiveIntensity = 1.5;
  mesh.material.needsUpdate = true;
}

let enableBondDistance = false;
let bondLengthElements = [];

document.addEventListener("click", measureBondLengthClick);

function measureBondLengthClick(event) {
  if (!enableBondDistance) return;
  event.preventDefault();
  mouse.x = (event.clientX / canvas.clientWidth) * 2 - 1;
  mouse.y = -(event.clientY / canvas.clientHeight) * 2 + 1;
  console.log(mouse.x, mouse.y);

  raycaster.setFromCamera(mouse, camera);

  const meshObjects = [];
  scene.traverse((obj) => {
    if (obj instanceof THREE.Mesh) {
      meshObjects.push(obj);
    }
  });

  const intersections = raycaster.intersectObjects(meshObjects, true);
  const infoBox = document.getElementById("atom-info");

  if (intersections.length > 0) {
    const clickedAtom = intersections[0].object;
    const pos = clickedAtom.position;
    const name = clickedAtom.name || clickedAtom.userData?.element || "Unknown";
    const element = name.split("_")[0];
    infoBox.innerHTML = `${element}<br>X: ${pos.x.toFixed(
      5
    )}<br>Y: ${pos.y.toFixed(5)}<br>Z: ${pos.z.toFixed(5)}`;
    infoBox.style.display = "block";

    if (bondLengthElements.length < 2) {
      highlightAtom(clickedAtom);
      bondLengthElements.push(clickedAtom);
    }
    if (bondLengthElements.length == 2) {
      enableBondDistance = false;
      const pos1 = bondLengthElements[0].position;
      const pos2 = bondLengthElements[1].position;

      const bondDistance = pos1.distanceTo(pos2);

      // Create dotted material
      const material = new THREE.LineDashedMaterial({
        color: 0x000000,
        dashSize: 0.2,
        gapSize: 0.2,
        scale: 1, // This won't work on most platforms, but fine to leave
      });

      // Create geometry
      const geometry = new THREE.BufferGeometry().setFromPoints([pos1, pos2]);
      // geometry.computeLineDistances(); // essential for lines to appear dotted
      const dottedLine = new THREE.Line(geometry, material);
      dottedLine.computeLineDistances();
      dottedLine.name = "bond_length_line";
      scene.add(dottedLine);

      // Show bond info
      const bondLengthInfo = document.getElementById("bond-length-info");
      const elem1 = bondLengthElements[0].name.split("_")[0];
      const elem2 = bondLengthElements[1].name.split("_")[0];
      bondLengthInfo.innerHTML = `${elem1}-${elem2} bond distance = ${bondDistance.toFixed(
        5
      )} <span>&#8491;</span>`;
    }
  }
}

// Make it globally available
window.clearBondLengthInfo = clearBondLengthInfo;

function clearBondLengthInfo() {
  // Remove the bond length line
  const lineToRemove = scene.getObjectByName("bond_length_line");
  if (lineToRemove) {
    scene.remove(lineToRemove);
    lineToRemove.geometry.dispose(); // free memory
    lineToRemove.material.dispose();
    console.log("line removed");
  }

  // Remove glowing from mesh
  bondLengthElements.forEach((mesh) => {
    const baseColor = mesh.material.color.clone();
    mesh.material.emissive = baseColor;
    mesh.material.emissiveIntensity = 0;
    mesh.material.needsUpdate = true;
  });

  bondLengthElements = []; // Reset bond length array
  const bondLengthInfo = document.getElementById("bond-length-info");
  bondLengthInfo.innerHTML = "N/A";
  document.getElementById("bond-length-btn").disabled = false;
}

// Make it globally available
window.measureBondLength = measureBondLength;

function measureBondLength() {
  clearAll();
  enableBondDistance = true;
  document.getElementById("bond-length-btn").disabled = true;
}

let enableBondAngle = false;
let bondAngleElements = [];

document.addEventListener("click", measureBondAngleClick);

function measureBondAngleClick(event) {
  if (!enableBondAngle) return;
  event.preventDefault();
  mouse.x = (event.clientX / canvas.clientWidth) * 2 - 1;
  mouse.y = -(event.clientY / canvas.clientHeight) * 2 + 1;
  console.log(mouse.x, mouse.y);

  raycaster.setFromCamera(mouse, camera);

  const meshObjects = [];
  scene.traverse((obj) => {
    if (obj instanceof THREE.Mesh) {
      meshObjects.push(obj);
    }
  });

  const intersections = raycaster.intersectObjects(meshObjects, true);
  const infoBox = document.getElementById("atom-info");

  if (intersections.length > 0) {
    const clickedAtom = intersections[0].object;
    const pos = clickedAtom.position;
    const name = clickedAtom.name || clickedAtom.userData?.element || "Unknown";
    const element = name.split("_")[0];
    infoBox.innerHTML = `${element}<br>X: ${pos.x.toFixed(
      5
    )}<br>Y: ${pos.y.toFixed(5)}<br>Z: ${pos.z.toFixed(5)}`;
    infoBox.style.display = "block";

    if (bondAngleElements.length < 3) {
      highlightAtom(clickedAtom);
      bondAngleElements.push(clickedAtom);
    }
    if (bondAngleElements.length == 3) {
      enableBondAngle = false;
      const A = bondAngleElements[0].position;
      const B = bondAngleElements[1].position;
      const C = bondAngleElements[2].position;

      const AB = new THREE.Vector3().subVectors(A, B); // vector from B to A
      const CB = new THREE.Vector3().subVectors(C, B); // vector from B to C

      const angleRad = AB.angleTo(CB); // in radians
      const angleDeg = THREE.MathUtils.radToDeg(angleRad); // convert to degrees

      // Create dotted material
      const material = new THREE.LineDashedMaterial({
        color: 0x000000,
        dashSize: 0.2,
        gapSize: 0.2,
        scale: 1,
      });

      // Create geometry
      let geometry = new THREE.BufferGeometry().setFromPoints([A, B]);
      let dottedLine = new THREE.Line(geometry, material);
      dottedLine.computeLineDistances();
      dottedLine.name = "bond_angle_line";
      scene.add(dottedLine);

      geometry = new THREE.BufferGeometry().setFromPoints([B, C]);
      dottedLine = new THREE.Line(geometry, material);
      dottedLine.computeLineDistances();
      dottedLine.name = "bond_angle_line";
      scene.add(dottedLine);

      // Show bond info
      const bondLengthInfo = document.getElementById("bond-angle-info");
      const elem1 = bondAngleElements[0].name.split("_")[0];
      const elem2 = bondAngleElements[1].name.split("_")[0];
      const elem3 = bondAngleElements[2].name.split("_")[0];
      bondLengthInfo.innerHTML = `${elem1}-${elem2}-${elem3} bond angle = ${angleDeg.toFixed(
        3
      )}<span>&deg;</span>`;
    }
  }
}

// Make it globally available
window.clearBondAngleInfo = clearBondAngleInfo;

function clearBondAngleInfo() {
  // Remove both bond angle lines
  const meshesToRemove = [];
  scene.traverse((obj) => {
    if (obj.isLine && obj.name === "bond_angle_line") {
      meshesToRemove.push(obj);
    }
  });
  meshesToRemove.forEach((mesh) => {
    scene.remove(mesh);
    mesh.geometry.dispose();
    mesh.material.dispose();
    console.log("Removed:", mesh.name);
  });

  // Remove glowing from mesh
  bondAngleElements.forEach((mesh) => {
    const baseColor = mesh.material.color.clone();
    mesh.material.emissive = baseColor;
    mesh.material.emissiveIntensity = 0;
    mesh.material.needsUpdate = true;
  });

  bondAngleElements = []; // Reset bond length array
  const bondAngleInfo = document.getElementById("bond-angle-info");
  bondAngleInfo.innerHTML = "N/A";
  document.getElementById("bond-angle-btn").disabled = false;
}

// Make it globally available
window.measureBondAngle = measureBondAngle;

function measureBondAngle() {
  clearAll();
  enableBondAngle = true;
  document.getElementById("bond-angle-btn").disabled = true;
}

function clearAll() {
  clearBondLengthInfo();
  clearBondAngleInfo();
}

function resizeRendererToDisplaySize() {
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;

  if (canvas.width !== width || canvas.height !== height) {
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
  }
}

// Animate
function animate() {
  requestAnimationFrame(animate);
  controls.update(); // Optional but good for damping
  resizeRendererToDisplaySize();
  renderer.render(scene, camera);
  renderer.setClearColor(0xffffff, 0);
  canvas.labelRenderer.render(scene, camera);
  hoverAtom();
}
animate();

window.addEventListener("resize", () => {
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;

  camera.aspect = width / height;
  camera.updateProjectionMatrix();

  renderer.setSize(width, height);
  labelRenderer.setSize(width, height);
});
