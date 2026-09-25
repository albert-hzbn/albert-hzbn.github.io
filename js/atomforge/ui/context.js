// Shared application context: state, the 3D viewer and a tiny event bus.
import { Structure } from '../core/index.js';
import { Viewer } from '../render/viewer.js';

export const state = {
    structure: new Structure(),
    undo: [],
    redo: [],
    selection: new Set(),
    mode: 'select',
    pending: [],            // picked instances for the current measurement
    measurements: [],
    customColors: {},
    hidden: new Set(),
    bonds: [],
    cn: [],
    display: [],
    spaceGroupLabel: '',
    source: null,           // crystal the derived builders start from
    sourceLabel: '',
    view: {
        style: 'ball', atomScale: 1, bondRadius: 0.14, tolerance: 1.15, bonds: true, labels: false,
        color: 'element', cell: true, boundary: true, reps: [1, 1, 1],
        plane: null, direction: null,
    },
};


export const viewer = new Viewer(document.getElementById('viewport'));

// Minimal publish/subscribe used to decouple modules (e.g. history → panels).
const listeners = {};
export const bus = {
    on(event, fn) { (listeners[event] = listeners[event] || []).push(fn); },
    emit(event, payload) { (listeners[event] || []).forEach((fn) => fn(payload)); },
};
