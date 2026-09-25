// Line and token helpers shared by the readers.
export const lines = (text) => text.replace(/\r/g, '').split('\n');
export const toks = (line) => line.trim().split(/\s+/).filter(Boolean);
