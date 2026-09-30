// Numerical methods for polynomial interpolation and Newton–Cotes quadrature.

const NA = (() => {
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => (n === 1 ? a : a + ((b - a) * i) / (n - 1)));

  // ---------- Nodes ----------
  function uniformNodes(a, b, N) {
    return linspace(a, b, N + 1);
  }
  // Chebyshev points (roots of T_{N+1}) mapped to [a,b], in increasing order.
  function chebyshevNodes(a, b, N) {
    const xs = [];
    for (let i = 0; i <= N; i++) xs.push((a + b) / 2 + ((b - a) / 2) * Math.cos(((2 * i + 1) * Math.PI) / (2 * N + 2)));
    return xs.sort((p, q) => p - q);
  }
  function checkDistinct(xs, label = 'Nodes') {
    const s = [...xs].sort((p, q) => p - q);
    for (let i = 1; i < s.length; i++) if (s[i] === s[i - 1]) throw new Error(`${label} must be distinct (${s[i]} appears twice).`);
    return s;
  }

  // ---------- Linear interpolation ----------
  function linearInterpolant(x0, x1, f0, f1) {
    const slope = (f1 - f0) / (x1 - x0);
    return { slope, eval: (x) => f0 + slope * (x - x0) };
  }

  // ---------- Piecewise-linear interpolation ----------
  function piecewiseLinear(ys, fs) {
    const J = ys.length - 1;
    return (x) => {
      let j = 1;
      if (x <= ys[0]) j = 1;
      else if (x >= ys[J]) j = J;
      else { let lo = 0, hi = J; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (ys[m] <= x) lo = m; else hi = m; } j = hi; }
      const h = ys[j] - ys[j - 1];
      return fs[j - 1] + ((fs[j] - fs[j - 1]) / h) * (x - ys[j - 1]);
    };
  }

  // ---------- Lagrange form ----------
  function lagrangeBasis(xs, j) {
    let denom = 1;
    for (let i = 0; i < xs.length; i++) if (i !== j) denom *= xs[j] - xs[i];
    return (x) => { let p = 1; for (let i = 0; i < xs.length; i++) if (i !== j) p *= x - xs[i]; return p / denom; };
  }
  function lagrangeInterpolant(xs, fs) {
    const L = xs.map((_, j) => lagrangeBasis(xs, j));
    return { basis: L, eval: (x) => { let s = 0; for (let j = 0; j < xs.length; j++) s += fs[j] * L[j](x); return s; } };
  }

  // ---------- Newton's divided differences ----------
  // table[k][i] = f[x_i, ..., x_{i+k}]
  function dividedDifferences(xs, fs) {
    const n = xs.length, table = [fs.slice()];
    for (let k = 1; k < n; k++) {
      const prev = table[k - 1], col = [];
      for (let i = 0; i + k < n; i++) col.push((prev[i + 1] - prev[i]) / (xs[i + k] - xs[i]));
      table.push(col);
    }
    return table;
  }
  function newtonInterpolant(xs, fs) {
    const table = dividedDifferences(xs, fs);
    const coeffs = table.map((col) => col[0]); // f[x0], f[x0,x1], ...
    const evalN = (x) => { // nested (Horner-like) evaluation
      let p = coeffs[coeffs.length - 1];
      for (let k = coeffs.length - 2; k >= 0; k--) p = coeffs[k] + (x - xs[k]) * p;
      return p;
    };
    return { table, coeffs, eval: evalN };
  }

  // Monomial coefficients c[0] + c[1] x + ... of the interpolant (via Newton form).
  function monomialCoefficients(xs, coeffs) {
    let poly = [coeffs[coeffs.length - 1]];
    for (let k = coeffs.length - 2; k >= 0; k--) {
      // poly = coeffs[k] + (x - xs[k]) * poly
      const next = new Array(poly.length + 1).fill(0);
      for (let i = 0; i < poly.length; i++) { next[i + 1] += poly[i]; next[i] -= xs[k] * poly[i]; }
      next[0] += coeffs[k];
      poly = next;
    }
    return poly;
  }

  // w_{N+1}(x) = (x - x0)...(x - xN)
  const nodePolynomial = (xs) => (x) => xs.reduce((p, xi) => p * (x - xi), 1);

  function maxAbsOn(g, a, b, samples = 4001) {
    let best = 0, at = a;
    for (const x of linspace(a, b, samples)) { const v = Math.abs(g(x)); if (v > best || !isFinite(v)) { best = v; at = x; } }
    return { value: best, at };
  }

  // ---------- Newton-Cotes rules on [0,1] (exact rational weights) ----------
  // w_i = ∫_0^1 L_i(x) dx with x_i = i/N.
  const ncCache = {};
  function newtonCotes(N) {
    if (ncCache[N]) return ncCache[N];
    const F = (n, d = 1) => math.fraction(n, d);
    const nodes = [], weights = [];
    for (let i = 0; i <= N; i++) {
      // Build L_i as a polynomial with fraction coefficients.
      let poly = [F(1)];
      let denom = F(1);
      for (let k = 0; k <= N; k++) {
        if (k === i) continue;
        const xk = F(k, N);
        const next = new Array(poly.length + 1).fill(null).map(() => F(0));
        for (let m = 0; m < poly.length; m++) {
          next[m + 1] = next[m + 1].add(poly[m]);
          next[m] = next[m].sub(poly[m].mul(xk));
        }
        poly = next;
        denom = denom.mul(F(i, N).sub(xk));
      }
      let integral = F(0);
      for (let m = 0; m < poly.length; m++) integral = integral.add(poly[m].div(F(m + 1)));
      const w = N === 0 ? F(1) : integral.div(denom);
      nodes.push(N === 0 ? F(1, 2) : F(i, N));
      weights.push(w);
    }
    const rule = {
      N,
      nodesFrac: nodes,
      weightsFrac: weights,
      nodes: nodes.map((q) => q.valueOf()),
      weights: weights.map((q) => q.valueOf()),
    };
    Object.assign(rule, degreeOfPrecision(rule));
    ncCache[N] = rule;
    return rule;
  }

  // Degree of precision d with E(x^r) = 0 for r <= d, E(x^{d+1}) != 0. Exact (fractions).
  function degreeOfPrecision(rule) {
    const rows = [];
    let d = -1;
    for (let r = 0; r <= rule.N + 3; r++) {
      const I = math.fraction(1, r + 1);
      let Q = math.fraction(0);
      rule.nodesFrac.forEach((x, i) => { Q = Q.add(rule.weightsFrac[i].mul(r === 0 ? math.fraction(1) : x.pow(r))); });
      const E = I.sub(Q);
      rows.push({ r, I, Q, E });
      if (E.valueOf() === 0 && d === r - 1) d = r;
      if (E.valueOf() !== 0) break;
    }
    // Error constant K with E(f) = K f^{(d+1)}(ξ), K = E(x^{d+1})/(d+1)!
    const Efirst = rows[rows.length - 1].E;
    const K = Efirst.div(math.fraction(Expr.FACT[d + 1]));
    return { dop: d, dopRows: rows, errorConst: K, errorConstValue: K.valueOf() };
  }

  // Rule on [a,b]: Q^{[a,b]}(f) = (b-a) Σ w_i f(a + (b-a) x_i)
  function applyRule(rule, f, a, b) {
    let s = 0;
    for (let i = 0; i < rule.nodes.length; i++) s += rule.weights[i] * f(a + (b - a) * rule.nodes[i]);
    return (b - a) * s;
  }

  // Composite rule on mesh ys
  function applyComposite(rule, f, ys) {
    let s = 0;
    for (let j = 1; j < ys.length; j++) s += applyRule(rule, f, ys[j - 1], ys[j]);
    return s;
  }

  // ---------- Reference integral: adaptive Gauss–Kronrod (7,15) ----------
  const GK_X = [0.991455371120812639, 0.949107912342758525, 0.864864423359769073, 0.741531185599394440,
    0.586087235467691130, 0.405845151377397167, 0.207784955007898468, 0.0];
  const GK_WK = [0.022935322010529225, 0.063092092629978553, 0.104790010322250184, 0.140653259715525919,
    0.169004726639267903, 0.190350578064785410, 0.204432940075298892, 0.209482141084727828];
  const GK_WG = [0, 0.129484966168869693, 0, 0.279705391489276668, 0, 0.381830050505118945, 0, 0.417959183673469388];

  function gk15(f, a, b) {
    const c = (a + b) / 2, h = (b - a) / 2;
    let k = 0, g = 0;
    for (let i = 0; i < 8; i++) {
      if (i === 7) { const v = f(c); k += GK_WK[i] * v; g += GK_WG[i] * v; continue; }
      const v = f(c - h * GK_X[i]) + f(c + h * GK_X[i]);
      k += GK_WK[i] * v; g += GK_WG[i] * v;
    }
    return { value: k * h, err: Math.abs((k - g) * h) };
  }

  function referenceIntegral(f, a, b, tol = 1e-13) {
    let total = 0, errTotal = 0, evals = 0;
    const stack = [[a, b, 0]];
    while (stack.length) {
      const [l, r, depth] = stack.pop();
      const res = gk15(f, l, r); evals += 15;
      if (!isFinite(res.value)) throw new Error(`f is not finite on [${l}, ${r}] — cannot compute a reference value.`);
      if (res.err <= Math.max(tol * Math.abs(res.value), 1e-15) * ((r - l) / (b - a)) || depth > 50 || evals > 2e6) {
        total += res.value; errTotal += res.err;
      } else {
        const m = (l + r) / 2;
        stack.push([l, m, depth + 1], [m, r, depth + 1]);
      }
    }
    return { value: total, errEstimate: errTotal };
  }

  return {
    linspace, uniformNodes, chebyshevNodes, checkDistinct,
    linearInterpolant, piecewiseLinear, lagrangeBasis, lagrangeInterpolant,
    dividedDifferences, newtonInterpolant, monomialCoefficients, nodePolynomial, maxAbsOn,
    newtonCotes, applyRule, applyComposite, referenceIntegral,
  };
})();
