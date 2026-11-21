// The story's engine: loads data/story.json (written by build_data.jl), fills the
// numbers and tables, draws the charts, and changes each docked chart as its steps
// scroll into view.

"use strict";

const NS = "http://www.w3.org/2000/svg";
const $ = (id) => document.getElementById(id);
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
let D = null;

// ------------------------------------------------------------------ numbers

function get(path) {
  return path.split(".").reduce((o, k) => (o == null ? undefined : o[k]), D);
}

function fmt(v, f) {
  if (v === null || v === undefined || (typeof v === "number" && !Number.isFinite(v))) return "–";
  if (typeof v !== "number") return String(v);
  switch (f) {
    case "int": return Math.round(v).toLocaleString("en-GB");
    case "f2": return v.toFixed(2);
    case "f3": return v.toFixed(3);
    case "x3": return "×" + v.toFixed(3);
    case "p1": return v.toFixed(1) + "%";
    case "ord": {
      const n = Math.round(v), s = ["th", "st", "nd", "rd"], m = n % 100;
      return n + (s[(m - 20) % 10] || s[m] || s[0]);
    }
    default: return Math.abs(v) >= 1000 ? v.toLocaleString("en-GB") : String(+v.toPrecision(4));
  }
}

function sig(v, n = 3) {
  if (v === null || !Number.isFinite(v)) return "–";
  const a = Math.abs(v);
  if (a >= 100) return Math.round(v).toLocaleString("en-GB");
  return String(+v.toPrecision(n));
}

function bindNumbers() {
  for (const el of document.querySelectorAll("[data-v]")) el.textContent = fmt(get(el.dataset.v), el.dataset.f);
}

// ------------------------------------------------------------------ tables

function table(id, head, rows) {
  const t = $(id);
  const cap = t.querySelector("caption");
  t.innerHTML = "";
  if (cap) t.appendChild(cap);
  const th = t.createTHead().insertRow();
  for (const h of head) { const c = document.createElement("th"); c.textContent = h; th.appendChild(c); }
  const tb = t.createTBody();
  for (const r of rows) {
    const tr = tb.insertRow();
    for (const cell of r) {
      const td = tr.insertCell();
      if (cell && typeof cell === "object") { td.innerHTML = cell.html; if (cell.cls) td.className = cell.cls; }
      else td.textContent = cell;
    }
  }
}

const iqr = (q) => ({ html: `${sig(q[0])}<span class="q">${sig(q[1])}–${sig(q[2])}</span>` });

function buildTables() {
  const strata = $("strata");
  strata.innerHTML = "";
  for (const s of D.design.strata) {
    const d = document.createElement("div");
    d.className = "stratum";
    d.innerHTML = `<b>${s.L.toLocaleString("en-GB")} loci</b><span>${s.u} crossovers · ${s.n.toLocaleString("en-GB")} populations · ${s.sample} rerun</span>`;
    strata.appendChild(d);
  }
  table("params-table", ["Parameter", "5%", "median", "95%"], D.design.parameters.map((p) => [p.name, ...p.q.map((v) => sig(v))]));
  table("paper-table", ["Quantity", "Paper", "Recomputed"], D.paper.map((r) => [r.quantity, r.paper, iqr(r.table)]));
  table("survivors-table", ["Survivors'", "Paper", "All published survivors", "Published runs, sample", "QTLSim reruns"],
    D.survivors.stats.map((r) => [r.quantity, r.paper, iqr(r.all), iqr(r.published), iqr(r.reruns)]));
  table("summary-table", ["Result", "Paper", "Paper's table", "Published runs, sample", "QTLSim reruns"],
    D.summary.map((r) => [r.quantity, r.paper || "", sig(r.table), sig(r.published), sig(r.reruns)]));
}

// ------------------------------------------------------------------ an SVG chart

function el(name, attrs = {}, parent) {
  const e = document.createElementNS(NS, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}

function niceLin(lo, hi, n = 5) {
  if (!(hi > lo)) { const d = Math.abs(lo) || 1; lo -= d / 2; hi += d / 2; }
  const raw = (hi - lo) / n, p = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw);
  const t = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) t.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return t;
}
function niceLog(lo, hi) {
  const t = [], a = Math.floor(Math.log10(lo)), b = Math.ceil(Math.log10(hi));
  const mult = b - a <= 3 ? [1, 2, 5] : [1];
  for (let e = a; e <= b; e++) for (const m of mult) { const v = m * 10 ** e; if (v >= lo * 0.999 && v <= hi * 1.001) t.push(v); }
  return t;
}
const tick = (v) => {
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(0);
  return String(+v.toPrecision(4));
};

// spec: {x:{label,log,min,max}, y:{...}, series:[{kind, data, color, dash, r, open, label}], hover}
function draw(svg, spec) {
  svg.innerHTML = "";
  const W = svg.clientWidth || 600, H = svg.clientHeight || 400;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const m = { l: 62, r: 16, t: 10, b: 46 };
  const all = (k) => spec.series.flatMap((s) => (s.kind === "vline" ? (k === "x" ? [s.at] : []) : s.kind === "hline" ? (k === "y" ? [s.at] : []) :
    s.kind === "interval" ? s.data.flatMap((d) => (k === "x" ? [d.lo, d.hi, d.est] : [d.y])) : s.data.map((p) => (k === "x" ? p[0] : p[1]))))
    .filter((v) => v !== null && Number.isFinite(v) && (!spec[k].log || v > 0));
  const dom = (k) => {
    const v = all(k);
    let lo = spec[k].min ?? Math.min(...v), hi = spec[k].max ?? Math.max(...v);
    if (spec[k].log) {
      if (spec[k].min === undefined) lo /= 1.15;
      if (spec[k].max === undefined) hi *= 1.15;
    } else if (spec[k].min === undefined || spec[k].max === undefined) {
      const pad = (hi - lo) * 0.05 || Math.abs(hi) * 0.1 || 1;
      if (spec[k].min === undefined) lo -= pad;
      if (spec[k].max === undefined) hi += pad;
    }
    return [lo, hi];
  };
  const [x0, x1] = dom("x"), [y0, y1] = dom("y");
  const sx = spec.x.log ? (v) => m.l + (Math.log(v / x0) / Math.log(x1 / x0)) * (W - m.l - m.r) : (v) => m.l + ((v - x0) / (x1 - x0)) * (W - m.l - m.r);
  const sy = spec.y.log ? (v) => H - m.b - (Math.log(v / y0) / Math.log(y1 / y0)) * (H - m.t - m.b) : (v) => H - m.b - ((v - y0) / (y1 - y0)) * (H - m.t - m.b);
  const ok = (p) => p && p[0] !== null && p[1] !== null && Number.isFinite(p[0]) && Number.isFinite(p[1]) && (!spec.x.log || p[0] > 0) && (!spec.y.log || p[1] > 0);

  const g = el("g", { class: "grid" }, svg);
  for (const v of spec.y.ticks || (spec.y.log ? niceLog(y0, y1) : niceLin(y0, y1, 5))) {
    if (v < y0 || v > y1) continue;
    el("line", { x1: m.l, x2: W - m.r, y1: sy(v), y2: sy(v) }, g);
    el("text", { x: m.l - 8, y: sy(v) + 4, "text-anchor": "end" }, svg).textContent = spec.y.fmt ? spec.y.fmt(v) : tick(v);
  }
  for (const v of spec.x.ticks || (spec.x.log ? niceLog(x0, x1) : niceLin(x0, x1, 6))) {
    if (v < x0 || v > x1) continue;
    el("text", { x: sx(v), y: H - m.b + 18, "text-anchor": "middle" }, svg).textContent = spec.x.fmt ? spec.x.fmt(v) : tick(v);
  }
  el("text", { class: "axis-title", x: (m.l + W - m.r) / 2, y: H - 8, "text-anchor": "middle" }, svg).textContent = spec.x.label;
  el("text", { class: "axis-title", x: 14, y: (m.t + H - m.b) / 2, "text-anchor": "middle", transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})` }, svg).textContent = spec.y.label;

  for (const s of spec.series) {
    const color = s.color.startsWith("--") ? css(s.color) : s.color;
    if (s.kind === "line" || s.kind === "linepoints") {
      let d = "", pen = false;
      for (const p of s.data) {
        if (!ok(p)) { pen = false; continue; }
        d += (pen ? "L" : "M") + sx(p[0]).toFixed(1) + "," + sy(p[1]).toFixed(1);
        pen = true;
      }
      el("path", { d, fill: "none", stroke: color, "stroke-width": s.width || 2, "stroke-dasharray": s.dash || "", "stroke-linejoin": "round" }, svg);
    }
    if (s.kind === "points" || s.kind === "linepoints") {
      for (const p of s.data) {
        if (!ok(p)) continue;
        el("circle", { cx: sx(p[0]), cy: sy(p[1]), r: s.r || 3.5, fill: s.open ? "none" : color, stroke: color, "stroke-width": 1.3, "fill-opacity": s.alpha ?? 1 }, svg);
      }
    }
    if (s.kind === "vline" && s.at >= x0 && s.at <= x1) {
      el("line", { x1: sx(s.at), x2: sx(s.at), y1: m.t, y2: H - m.b, stroke: color, "stroke-dasharray": s.dash || "4 4" }, svg);
      if (s.label) el("text", { class: "note", x: sx(s.at) + 6, y: m.t + 14 }, svg).textContent = s.label;
    }
    if (s.kind === "hline" && s.at >= y0 && s.at <= y1) {
      el("line", { x1: m.l, x2: W - m.r, y1: sy(s.at), y2: sy(s.at), stroke: color, "stroke-dasharray": s.dash || "4 4" }, svg);
    }
    if (s.kind === "interval") {
      for (const d of s.data) {
        const yy = sy(d.y);
        el("line", { x1: sx(d.lo), x2: sx(d.hi), y1: yy, y2: yy, stroke: d.color ? css(d.color) : color, "stroke-width": 3, "stroke-linecap": "round" }, svg);
        el("circle", { cx: sx(d.est), cy: yy, r: 6, fill: d.color ? css(d.color) : color }, svg);
        el("text", { x: sx(d.hi) + 10, y: yy + 4, fill: css("--ink") }, svg).textContent = d.text || "";
      }
    }
  }
  return { sx, sy, m, W, H };
}

function legend(id, items) {
  $(id).innerHTML = items.map(([c, t, dash]) =>
    `<span><i style="background:${dash ? "transparent" : css(c)};${dash ? `border-top:2px dashed ${css(c)};height:0` : ""}"></i>${t}</span>`).join("");
}

// ------------------------------------------------------------------ the docked charts

const zipXY = (xs, ys) => xs.map((x, i) => [x, ys[i]]);
const after = (e, from) => e.gen.map((g, i) => i).filter((i) => e.gen[i] >= from);

const SCENES = {
  model: [
    () => {
      const e = D.examples[0];
      return { cap: "Genetic variance builds up from nothing", legend: [["--c-pub", "genetic variance"]],
        spec: { x: { label: "generation", min: 0, max: 10000 }, y: { label: "variance", min: 0 },
          series: [{ kind: "line", data: zipXY(e.gen, e.vg).filter((p) => p[0] <= 10000), color: "--c-pub" }] } };
    },
    () => {
      const e = D.examples[0];
      const pre = (a) => zipXY(e.gen, a).filter((p) => p[0] <= 10000);
      return { cap: "Mutation, selection and drift in balance", legend: [["--c-all", "phenotypic variance"], ["--c-pub", "genetic variance"]],
        spec: { x: { label: "generation", min: 0, max: 10000 }, y: { label: "variance", min: 0 },
          series: [{ kind: "line", data: pre(e.vp), color: "--c-all" }, { kind: "line", data: pre(e.vg), color: "--c-pub" }] } };
    },
    () => {
      const e = D.examples[0], i = after(e, 9600);
      return { cap: `The optimum moves by ${e.d} per generation: the population keeps up`, legend: [["--c-rerun", "optimum", true], ["--c-pub", "mean phenotype"]],
        spec: { x: { label: "generation" }, y: { label: "phenotype" },
          series: [{ kind: "line", data: i.map((j) => [e.gen[j], e.theta[j]]), color: "--c-rerun", dash: "6 4" },
                   { kind: "line", data: i.map((j) => [e.gen[j], e.z[j]]), color: "--c-pub" },
                   { kind: "vline", at: 10000, color: "--muted", label: "the move" }] } };
    },
    () => {
      const e = D.examples[1], i = after(e, 9600);
      return { cap: `At ${e.d} per generation it falls behind and dies`, legend: [["--c-rerun", "optimum", true], ["--c-pub", "mean phenotype"]],
        spec: { x: { label: "generation" }, y: { label: "phenotype" },
          series: [{ kind: "line", data: i.map((j) => [e.gen[j], e.theta[j]]), color: "--c-rerun", dash: "6 4" },
                   { kind: "line", data: i.map((j) => [e.gen[j], e.z[j]]), color: "--c-pub" },
                   { kind: "vline", at: 10000, color: "--muted", label: "the move" }] } };
    },
    () => {
      const lag = (e) => after(e, 10000).map((j) => [e.gen[j] - 10000, (e.theta[j] - e.z[j]) / Math.sqrt(e.vp[j])]);
      const [a, b] = D.examples;
      return { cap: "The lag, in phenotypic standard deviations", legend: [["--c-pub", `d = ${a.d}`], ["--c-rerun", `d = ${b.d}`]],
        spec: { x: { label: "generations since the move", min: 0 }, y: { label: "lag (θ − z̄) / SD", min: 0 },
          series: [{ kind: "line", data: lag(a), color: "--c-pub" }, { kind: "line", data: lag(b), color: "--c-rerun" }] } };
    },
  ],
  rerun: [
    () => ({ cap: "Alive 2000 generations after the move", legend: [["--c-all", "all published"], ["--c-pub", "published runs, sample"], ["--c-rerun", "QTLSim reruns"]],
      spec: { x: { label: "standardised rate of environmental change x", log: true }, y: { label: "share alive", min: 0, max: 1 },
        series: [{ kind: "line", data: D.survival.all, color: "--c-all", width: 3 },
                 { kind: "linepoints", data: D.survival.published, color: "--c-pub" },
                 { kind: "linepoints", data: D.survival.reruns, color: "--c-rerun" }] } }),
    () => ({ cap: "Generations to extinction (Figure 1a)", legend: [["--c-all", "paper"], ["--c-pub", "published runs, sample"], ["--c-rerun", "QTLSim reruns"]],
      spec: { x: { label: "standardised rate of environmental change x", log: true }, y: { label: "median generations to extinction", log: true },
        series: [{ kind: "line", data: D.fig1.a.paper, color: "--c-all", width: 3 },
                 { kind: "linepoints", data: D.fig1.a.published, color: "--c-pub" },
                 { kind: "linepoints", data: D.fig1.a.reruns, color: "--c-rerun" },
                 { kind: "hline", at: 2000, color: "--muted" }] } }),
    () => ({ cap: "Standardised lag (Figure 1b)", legend: [["--c-all", "paper"], ["--c-pub", "published runs, sample"], ["--c-rerun", "QTLSim reruns"]],
      spec: { x: { label: "standardised rate of environmental change x", log: true }, y: { label: "median standardised lag", log: true },
        series: [{ kind: "line", data: D.fig1.b.paper, color: "--c-all", width: 3 },
                 { kind: "linepoints", data: D.fig1.b.published, color: "--c-pub" },
                 { kind: "linepoints", data: D.fig1.b.reruns, color: "--c-rerun" }] } }),
    () => ({ cap: "Mortality 45 generations into the move (Figure 1c)", legend: [["--c-all", "paper"], ["--c-pub", "published runs, sample"], ["--c-rerun", "QTLSim reruns"]],
      spec: { x: { label: "median standardised lag", log: true }, y: { label: "share dying in a generation", min: 0, max: 1 },
        series: [{ kind: "line", data: D.fig1.c.paper, color: "--c-all", width: 3 },
                 { kind: "linepoints", data: D.fig1.c.published, color: "--c-pub" },
                 { kind: "linepoints", data: D.fig1.c.reruns, color: "--c-rerun" }] } }),
    () => {
      const p = D.reruns.filter((r) => r.VG > 0 && r.shc > 0).map((r) => [r.shc, r.VG]);
      const lo = Math.min(...p.flat()), hi = Math.max(...p.flat());
      return { cap: "Genetic variance before the move, against the house of cards", legend: [["--c-rerun", "QTLSim reruns"], ["--muted", "equality", true]],
        spec: { x: { label: "stochastic house-of-cards approximation", log: true }, y: { label: "simulated genetic variance", log: true },
          series: [{ kind: "line", data: [[lo, lo], [hi, hi]], color: "--muted", dash: "4 4", width: 1.5 },
                   { kind: "points", data: p, color: "--c-rerun", r: 2.6, alpha: 0.45 }] } };
    },
    () => {
      const a = D.agreement;
      return { cap: "QTLSim against the published runs of the same populations", legend: [["--c-rerun", "ratio, with 95% interval"]],
        spec: { x: { label: "ratio, QTLSim rerun / published run (geometric mean)", min: 0.95, max: 1.05 }, y: { label: "", min: 0, max: 3, ticks: [], fmt: () => "" },
          series: [{ kind: "vline", at: 1, color: "--muted" },
                   { kind: "interval", color: "--c-rerun", data: [
                     { y: 2, est: a.vg[0], lo: a.vg[1], hi: a.vg[2], text: `genetic variance ×${a.vg[0].toFixed(3)}` },
                     { y: 1, est: a.time[0], lo: a.time[1], hi: a.time[2], text: `time to extinction ×${a.time[0].toFixed(3)}` }] }] } };
    },
  ],
  moves: [
    () => {
      const c = D.crossover;
      return { cap: "QTLSim against the published runs: two crossover models", legend: [["--c-pub", "QTSim's (lowest point only)"], ["--c-rerun", "every crossover point"]],
        spec: { x: { label: "ratio to the published runs (geometric mean, 95% interval)", min: 0.95, max: 1.25 }, y: { label: "", min: 0, max: 5, ticks: [], fmt: () => "" },
          series: [{ kind: "vline", at: 1, color: "--muted" },
                   { kind: "interval", color: "--c-pub", data: [
                     { y: 4, est: c.qtsim.time[0], lo: c.qtsim.time[1], hi: c.qtsim.time[2], text: "time to extinction", color: "--c-pub" },
                     { y: 3, est: c.qtsim.vg[0], lo: c.qtsim.vg[1], hi: c.qtsim.vg[2], text: "genetic variance", color: "--c-pub" },
                     { y: 2, est: c.all.time[0], lo: c.all.time[1], hi: c.all.time[2], text: "time to extinction", color: "--c-rerun" },
                     { y: 1, est: c.all.vg[0], lo: c.all.vg[1], hi: c.all.vg[2], text: "genetic variance", color: "--c-rerun" }] }] } };
    },
    () => {
      const cols = ["--c-pub", "--c-4", "--c-rerun", "--c-5"];
      return { cap: "Alive 2000 generations after the move, by demography", legend: D.demography.map((s, i) => [cols[i], s.name]),
        spec: { x: { label: "standardised rate of environmental change x", log: true }, y: { label: "share alive", min: 0, max: 1 },
          series: D.demography.flatMap((s, i) => [{ kind: "line", data: s.curve, color: cols[i] }, { kind: "points", data: s.bins, color: cols[i] }]) } };
    },
    () => {
      const cols = ["--c-pub", "--c-4", "--c-rerun", "--c-5"];
      return { cap: "The rate at which half survive 2000 generations", legend: D.demography.map((s, i) => [cols[i], s.name]),
        spec: { x: { label: "standardised rate x at which half are alive", min: 0.04, max: 0.085 }, y: { label: "", min: 0, max: 5, ticks: [], fmt: () => "" },
          series: [{ kind: "interval", color: "--c-pub", data: D.demography.map((s, i) => ({
            y: 4 - i, est: s.x50, lo: s.x50, hi: s.x50, color: cols[i],
            text: `${s.name}: ${s.x50.toFixed(3)}` + (i ? `, time to extinction ×${s.time_ratio.toFixed(2)}` : "") })) }] } };
    },
  ],
};

const current = {};
function show(chart, step) {
  if (current[chart] === step) return;
  current[chart] = step;
  const scene = SCENES[chart][step]();
  $(`${chart}-cap`).textContent = scene.cap;
  legend(`${chart}-legend`, scene.legend || []);
  draw($(`${chart}-chart`), scene.spec);
}

function observeSteps() {
  const steps = [...document.querySelectorAll(".step")];
  const io = new IntersectionObserver((entries) => {
    for (const e of entries) {
      if (!e.isIntersecting) continue;
      const s = e.target;
      for (const o of steps) if (o.dataset.chart === s.dataset.chart) o.classList.toggle("active", o === s);
      show(s.dataset.chart, Number(s.dataset.step));
    }
  }, { rootMargin: "-45% 0px -45% 0px" });
  steps.forEach((s) => io.observe(s));
  for (const c of Object.keys(SCENES)) show(c, 0);
}

function observeSections() {
  const links = [...document.querySelectorAll(".sections a")];
  const io = new IntersectionObserver((entries) => {
    for (const e of entries) {
      if (!e.isIntersecting) continue;
      for (const a of links) a.classList.toggle("active", a.getAttribute("href") === "#" + e.target.id);
    }
  }, { rootMargin: "-50% 0px -50% 0px" });
  for (const a of links) { const t = document.querySelector(a.getAttribute("href")); if (t) io.observe(t); }
}

// ------------------------------------------------------------------ explore

const EX_Y = {
  T: { label: "generations to extinction (2000 = alive)", log: true, get: (r) => Math.min(r.T, 2000) },
  lag: { label: "standardised lag", log: true, get: (r) => r.lag },
  rate: { label: "rate of phenotypic change", log: true, get: (r) => r.rate },
  vgratio: { label: "genetic variance / house of cards", log: true, get: (r) => (r.shc > 0 ? r.VG / r.shc : null) },
};

function explorer() {
  const colors = { 100: "--c-pub", 500: "--c-4", 1000: "--c-rerun", 5000: "--c-5" };
  const loci = [...new Set(D.reruns.map((r) => r.L))].sort((a, b) => a - b);
  const us = [...new Set(D.reruns.map((r) => r.u))].sort((a, b) => a - b);
  const box = (id, vals, label) => {
    const f = $(id);
    for (const v of vals) {
      const l = document.createElement("label");
      l.className = "chk";
      l.innerHTML = `<input type="checkbox" value="${v}" checked> ${label(v)}`;
      f.appendChild(l);
    }
    f.addEventListener("change", render);
  };
  const swatch = (l) => `<i style="display:inline-block;width:10px;height:10px;border-radius:50%;background:${css(colors[l])}"></i>`;
  box("ex-loci", loci, (v) => swatch(v) + " " + v.toLocaleString("en-GB"));
  box("ex-u", us, (v) => String(v));
  $("ex-y").addEventListener("change", render);

  function render() {
    const L = new Set([...$("ex-loci").querySelectorAll("input:checked")].map((i) => Number(i.value)));
    const U = new Set([...$("ex-u").querySelectorAll("input:checked")].map((i) => Number(i.value)));
    const y = EX_Y[$("ex-y").value];
    const rs = D.reruns.filter((r) => L.has(r.L) && U.has(r.u) && r.x > 0 && y.get(r) > 0);
    $("ex-count").textContent = `${rs.length} populations; ${rs.filter((r) => !r.extinct).length} alive after 2000 generations`;
    const series = loci.filter((l) => L.has(l)).flatMap((l) => [
      { kind: "points", data: rs.filter((r) => r.L === l && r.extinct).map((r) => [r.x, y.get(r)]), color: colors[l], r: 3.4, alpha: 0.6 },
      { kind: "points", data: rs.filter((r) => r.L === l && !r.extinct).map((r) => [r.x, y.get(r)]), color: colors[l], r: 3.4, open: true }]);
    if (!series.some((s) => s.data.length)) { $("ex-chart").innerHTML = ""; return; }
    const g = draw($("ex-chart"), { x: { label: "standardised rate of environmental change x", log: true }, y: { label: y.label, log: y.log }, series });
    const svg = $("ex-chart"), tip = $("ex-tip");
    svg.onmousemove = (ev) => {
      const r0 = svg.getBoundingClientRect();
      const mx = ((ev.clientX - r0.left) / r0.width) * g.W, my = ((ev.clientY - r0.top) / r0.height) * g.H;
      let best = null, bd = 144;
      for (const r of rs) {
        const dx = g.sx(r.x) - mx, dy = g.sy(y.get(r)) - my, d = dx * dx + dy * dy;
        if (d < bd) { bd = d; best = r; }
      }
      if (!best) { tip.hidden = true; return; }
      tip.hidden = false;
      tip.innerHTML = `${best.L.toLocaleString("en-GB")} loci, ${best.u} crossovers<br>x = ${sig(best.x)}<br>` +
        (best.extinct ? `died after ${best.T.toLocaleString("en-GB")} generations` : "alive after 2000 generations") +
        `<br>lag ${sig(best.lag)} · rate of change ${sig(best.rate)}<br>V<sub>G</sub> ${sig(best.VG)} (house of cards ${sig(best.shc)})`;
      tip.style.left = Math.min(ev.clientX - r0.left + 14, r0.width - 230) + "px";
      tip.style.top = ev.clientY - r0.top + 14 + "px";
    };
    svg.onmouseleave = () => { tip.hidden = true; };
  }
  render();
  return render;
}

// ------------------------------------------------------------------ start

async function start() {
  try {
    const r = await fetch("data/story.json");
    if (!r.ok) throw new Error(`data/story.json: ${r.status}`);
    D = await r.json();
  } catch (e) {
    $("err").hidden = false;
    $("err").textContent = "Could not load the data (open the page through a web server, not from disk): " + e.message;
    return;
  }
  for (const e of D.examples) e.died = e.extinct ? e.last + 1 - D.meta.move : null;
  bindNumbers();
  buildTables();
  observeSteps();
  observeSections();
  const redrawExplorer = explorer();
  let t = 0;
  const redraw = () => {
    clearTimeout(t);
    t = setTimeout(() => {
      for (const c of Object.keys(current)) { const s = current[c]; current[c] = -1; show(c, s); }
      redrawExplorer();
    }, 150);
  };
  window.addEventListener("resize", redraw);
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redraw);
}

start();
