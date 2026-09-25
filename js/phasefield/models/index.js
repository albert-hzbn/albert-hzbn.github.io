// Registry of phase-field models (order = order of the tiles in the panel).
import spinodal from './spinodal.js';
import copolymer from './copolymer.js';
import grainGrowth from './grain-growth.js';
import precipitation from './precipitation.js';
import multiphase from './multiphase.js';
import dendrite from './dendrite.js';
import magnetic from './magnetic.js';
import swiftHohenberg from './swift-hohenberg.js';
import allenCahn from './allen-cahn.js';
import nucleation from './nucleation.js';

export const MODELS = [spinodal, grainGrowth, precipitation, multiphase, magnetic, dendrite, copolymer, nucleation, allenCahn, swiftHohenberg];

// Concrete model for given initial-condition values: resolves the field
// layout (which may depend on e.g. the number of order parameters).
export function resolveModel(model, initValues) {
    const layout = model.layout(initValues);
    return {
        ...model,
        fields: layout.fields,
        aux: layout.aux || [],
        wgslCommon: (layout.common || '') + (model.wgslCommon || ''),
    };
}

export function defaults(list) {
    return Object.fromEntries(list.map((p) => [p.key, p.value]));
}
