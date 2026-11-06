#!/usr/bin/env node
// Build wesbite/world.js: the world outline and the centre of every country,
// used by the "world map" pop-ups. Only needs to be re-run if you want a
// different map; it does not depend on the plasmid data.
//
//   npm install world-atlas@2.0.2 topojson-client@3.1.0
//   node scripts/make_world_js.js > wesbite/world.js
//
// Outline: Natural Earth 1:110m (public domain, via world-atlas), Antarctica
// dropped, equirectangular projection in a 1000 x 500 frame.
// Centres: centroid of each country's LARGEST landmass at 1:50m, so France sits
// in Europe rather than being pulled towards French Guiana.

const topo = require("topojson-client");
const c50 = require("world-atlas/countries-50m.json");
const c110 = require("world-atlas/countries-110m.json");

const W = 1000, H = 500;
const P = ([lon, lat]) => [(lon + 180) / 360 * W, (90 - lat) / 180 * H];
const polys = g => g.type === "Polygon" ? [g.coordinates] : g.coordinates;

function ringArea(r) { let a = 0; for (let i = 0; i < r.length - 1; i++) a += r[i][0] * r[i + 1][1] - r[i + 1][0] * r[i][1]; return a / 2; }
function ringCentroid(r) {
  let a = 0, x = 0, y = 0;
  for (let i = 0; i < r.length - 1; i++) {
    const k = r[i][0] * r[i + 1][1] - r[i + 1][0] * r[i][1];
    a += k; x += (r[i][0] + r[i + 1][0]) * k; y += (r[i][1] + r[i + 1][1]) * k;
  }
  return [x / (3 * a), y / (3 * a)];
}

// ---- country centres (Natural Earth names) ----------------------------------
const centres = {};
for (const f of topo.feature(c50, c50.objects.countries).features) {
  const best = polys(f.geometry).reduce((b, p) => Math.abs(ringArea(p[0])) > Math.abs(ringArea(b[0])) ? p : b);
  centres[f.properties.name] = P(ringCentroid(best[0])).map(v => Math.round(v * 10) / 10);
}
// PLSDB spellings -> Natural Earth names
const ALIASES = {
  "USA": "United States of America", "United States": "United States of America",
  "Czech Republic": "Czechia", "Viet Nam": "Vietnam", "Republic of Korea": "South Korea",
  "Korea": "South Korea", "Russian Federation": "Russia", "UK": "United Kingdom",
  "Democratic Republic of the Congo": "Dem. Rep. Congo", "Dominican Republic": "Dominican Rep.",
  "Ivory Coast": "Côte d'Ivoire", "Cote d'Ivoire": "Côte d'Ivoire",
  "Bosnia and Herzegovina": "Bosnia and Herz.", "North Macedonia": "Macedonia",
  "Eswatini": "eSwatini", "Hong Kong": "China", "Iran, Islamic Republic of": "Iran",
  "Turkiye": "Turkey", "Türkiye": "Turkey", "Lao PDR": "Laos",
};
for (const [k, v] of Object.entries(ALIASES)) if (centres[v]) centres[k] = centres[v];

// ---- outline path -----------------------------------------------------------
const path = topo.feature(c110, c110.objects.countries).features
  .filter(f => f.properties.name !== "Antarctica")
  .map(f => polys(f.geometry).map(p => p.map(ring => {
    let s = "", prev = null;
    for (const pt of ring) {
      const q = P(pt);
      // a jump across the antimeridian starts a new sub-path instead of a line across the map
      s += (prev && Math.abs(q[0] - prev[0]) > W / 2) ? "M" : (s ? "L" : "M");
      s += q[0].toFixed(1) + "," + q[1].toFixed(1); prev = q;
    }
    return s + "Z";
  }).join("")).join("")).join("");

process.stdout.write(
  `// World outline (Natural Earth 1:110m via world-atlas, public domain) in an equirectangular ${W}x${H} frame,\n` +
  `// and the centre of each country (centroid of its largest landmass, from 1:50m).\n` +
  `// Built by scripts/make_world_js.js.\n` +
  `window.WORLD_MAP = ${JSON.stringify({ w: W, h: H, path, centres })};\n`);
