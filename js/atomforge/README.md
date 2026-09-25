# AtomForge Web

Browser port of the structure-building and visualisation tools in
[AtomForge](https://github.com/albert-hzbn/AtomForge). Served by `atomforge.html`;
plain ES modules, no build step. three.js is loaded from jsDelivr via an import map.

## Layout

| Folder | Contents |
| --- | --- |
| `core/` | Maths, lattice parameters, the `Structure` model, neighbour search, supercells and cell transformations, element data. `core/index.js` re-exports everything. |
| `crystal/` | Space-group settings and symmetry expansion (`spacegroups.js`) and crystal presets. |
| `io/` | One module per file format (XYZ, POSCAR, CIF, LAMMPS, PDB); `io/index.js` does format detection and exposes the writers. |
| `builders/` | Pure structure builders: slab, grain boundary, nanoparticle, polycrystal, amorphous, stacking fault, solid solution. No DOM access. |
| `analysis/` | Bonds / coordination and the radial distribution function. |
| `render/` | The main three.js `Viewer`, the shared preview renderer and overlay geometry. |
| `ui/` | Application code: shared context and event bus, undo history, scene building, one module per builder card (`ui/builders/`), the other panels (`ui/panels/`), viewport interaction and the app chrome. |
| `data/` | `spacegroups.json`: all 530 Hall settings, from the spglib database (BSD-3-Clause). |

`main.js` is the entry point. Structure changes go through `ui/history.js`
(`setStructure`, undo/redo), which emits a `structure` event on the bus; the
panels and the scene re-render in response, so modules do not call each other
in cycles.

Each builder card registers `{ preview, build }` with `ui/builders/card.js`.
The preview runs (debounced) whenever the card's inputs change and is drawn by a
single shared WebGL renderer, scaled down automatically for large systems.
