// Parsing user-entered functions f(x) and computing high-order derivatives.
//
// Derivatives are computed with truncated Taylor-series arithmetic (automatic
// differentiation) over the math.js expression tree, so f^(k)(x) is exact up to
// rounding even for large k, e.g. the 10th derivative of 1/(25x^2+1).

const Expr = (() => {
  function normalise(src) {
    return String(src)
      .replace(/\*\*/g, '^')
      .replace(/\bln\s*\(/g, 'log(');
  }

  const CONSTANTS = ['pi', 'e', 'E', 'PI'];

  // Parse an expression in the given variables (default x). Returns { src, node, f, vars }
  // where f(v1, v2, ...) is a plain JS function taking the variables in order.
  function compile(src, vars = ['x'], label = 'f(x)') {
    const text = normalise(src).trim();
    if (!text) throw new Error(`Please enter ${label}.`);
    let node;
    try {
      node = math.parse(text);
    } catch (e) {
      throw new Error(`Could not read ${label}: ${e.message}`);
    }
    const unknown = [];
    node.traverse((n, path, parent) => {
      if (n.isSymbolNode && !(parent && parent.isFunctionNode && path === 'fn')) {
        if (!vars.includes(n.name) && !CONSTANTS.includes(n.name)) unknown.push(n.name);
      }
    });
    if (unknown.length) {
      const list = vars.length === 1 ? vars[0] : vars.slice(0, -1).join(', ') + ' and ' + vars[vars.length - 1];
      throw new Error(`Unknown variable "${unknown[0]}" in ${label} — the variable${vars.length > 1 ? 's are' : ' is'} ${list}.`);
    }
    const code = node.compile();
    const scope = {};
    const f = (...args) => {
      for (let i = 0; i < vars.length; i++) scope[vars[i]] = args[i];
      const v = code.evaluate(scope);
      return typeof v === 'number' ? v : (v && typeof v.re === 'number' && Math.abs(v.im) < 1e-14 ? v.re : NaN);
    };
    // Sanity check that it evaluates to a number somewhere.
    f(...vars.map(() => 0.5));
    return { src: text, node, f, vars, tex: safeTex(node) };
  }

  function safeTex(node) {
    try { return node.toTex({ parenthesis: 'auto', implicit: 'hide' }); } catch { return null; }
  }

  // Evaluate a numeric expression such as "1/4" or "pi/2".
  function number(src, label) {
    const text = normalise(src).trim();
    if (!text) throw new Error(`Please enter a value for ${label}.`);
    let v;
    try { v = math.evaluate(text); } catch { throw new Error(`Could not read ${label}: "${src}"`); }
    if (typeof v !== 'number' || !isFinite(v)) throw new Error(`${label} must be a finite number.`);
    return v;
  }

  function numberList(src, label) {
    const parts = String(src).split(/[,;\s]+/).filter(Boolean);
    if (!parts.length) throw new Error(`Please enter at least one value for ${label}.`);
    return parts.map((p, i) => number(p, `${label} (entry ${i + 1})`));
  }

  // ---------- Taylor-series arithmetic ----------
  // A series is a Float64Array c of length K+1 with c[k] = f^(k)(x0)/k!.

  const S = {
    constant(v, K) { const c = new Float64Array(K + 1); c[0] = v; return c; },
    variable(x0, K) { const c = new Float64Array(K + 1); c[0] = x0; if (K >= 1) c[1] = 1; return c; },
    add(a, b) { return a.map((v, i) => v + b[i]); },
    sub(a, b) { return a.map((v, i) => v - b[i]); },
    neg(a) { return a.map((v) => -v); },
    scale(a, s) { return a.map((v) => v * s); },
    mul(a, b) {
      const K = a.length - 1, c = new Float64Array(K + 1);
      for (let k = 0; k <= K; k++) { let s = 0; for (let i = 0; i <= k; i++) s += a[i] * b[k - i]; c[k] = s; }
      return c;
    },
    div(a, b) {
      const K = a.length - 1, c = new Float64Array(K + 1);
      for (let k = 0; k <= K; k++) {
        let s = a[k];
        for (let i = 1; i <= k; i++) s -= b[i] * c[k - i];
        c[k] = s / b[0];
      }
      return c;
    },
    exp(a) {
      const K = a.length - 1, e = new Float64Array(K + 1);
      e[0] = Math.exp(a[0]);
      for (let k = 1; k <= K; k++) { let s = 0; for (let j = 1; j <= k; j++) s += j * a[j] * e[k - j]; e[k] = s / k; }
      return e;
    },
    log(a) {
      const K = a.length - 1, l = new Float64Array(K + 1);
      l[0] = Math.log(a[0]);
      for (let k = 1; k <= K; k++) {
        let s = a[k];
        for (let j = 1; j < k; j++) s -= (j / k) * l[j] * a[k - j];
        l[k] = s / a[0];
      }
      return l;
    },
    powConst(a, p) {
      if (Number.isInteger(p) && p >= 0 && p <= 64) {
        let result = S.constant(1, a.length - 1), base = a, n = p;
        while (n > 0) { if (n & 1) result = S.mul(result, base); base = S.mul(base, base); n >>= 1; }
        return result;
      }
      if (Number.isInteger(p) && p < 0 && p >= -64) return S.div(S.constant(1, a.length - 1), S.powConst(a, -p));
      const K = a.length - 1, y = new Float64Array(K + 1);
      y[0] = Math.pow(a[0], p);
      for (let k = 1; k <= K; k++) {
        let s = 0;
        for (let j = 1; j <= k; j++) s += ((p + 1) * j - k) * a[j] * y[k - j];
        y[k] = s / (k * a[0]);
      }
      return y;
    },
    sincos(a, hyperbolic) {
      const K = a.length - 1, s = new Float64Array(K + 1), c = new Float64Array(K + 1);
      s[0] = hyperbolic ? Math.sinh(a[0]) : Math.sin(a[0]);
      c[0] = hyperbolic ? Math.cosh(a[0]) : Math.cos(a[0]);
      const sign = hyperbolic ? 1 : -1;
      for (let k = 1; k <= K; k++) {
        let ss = 0, cc = 0;
        for (let j = 1; j <= k; j++) { ss += j * a[j] * c[k - j]; cc += j * a[j] * s[k - j]; }
        s[k] = ss / k; c[k] = (sign * cc) / k;
      }
      return [s, c];
    },
    // Series for y where y' = q (q given as series of a', i.e. integrate), y0 given.
    integrate(q, y0) {
      const K = q.length - 1, y = new Float64Array(K + 1);
      y[0] = y0;
      for (let k = 1; k <= K; k++) y[k] = q[k - 1] / k;
      return y;
    },
    deriv(a) {
      const K = a.length - 1, d = new Float64Array(K + 1);
      for (let k = 0; k < K; k++) d[k] = (k + 1) * a[k + 1];
      return d;
    },
  };

  const FUNCS = {
    sin: (a) => S.sincos(a, false)[0],
    cos: (a) => S.sincos(a, false)[1],
    tan: (a) => { const [s, c] = S.sincos(a, false); return S.div(s, c); },
    sec: (a) => S.div(S.constant(1, a.length - 1), S.sincos(a, false)[1]),
    csc: (a) => S.div(S.constant(1, a.length - 1), S.sincos(a, false)[0]),
    cot: (a) => { const [s, c] = S.sincos(a, false); return S.div(c, s); },
    sinh: (a) => S.sincos(a, true)[0],
    cosh: (a) => S.sincos(a, true)[1],
    tanh: (a) => { const [s, c] = S.sincos(a, true); return S.div(s, c); },
    exp: (a) => S.exp(a),
    log: (a, b) => (b ? S.div(S.log(a), S.log(b)) : S.log(a)),
    log10: (a) => S.scale(S.log(a), 1 / Math.LN10),
    log2: (a) => S.scale(S.log(a), 1 / Math.LN2),
    sqrt: (a) => S.powConst(a, 0.5),
    cbrt: (a) => (a[0] < 0 ? S.neg(S.powConst(S.neg(a), 1 / 3)) : S.powConst(a, 1 / 3)),
    abs: (a) => (a[0] > 0 ? a : a[0] < 0 ? S.neg(a) : a.map((v, i) => (i === 0 ? 0 : NaN))),
    atan: (a) => {
      const one = S.constant(1, a.length - 1);
      return S.integrate(S.div(S.deriv(a), S.add(one, S.mul(a, a))), Math.atan(a[0]));
    },
    asin: (a) => {
      const one = S.constant(1, a.length - 1);
      return S.integrate(S.div(S.deriv(a), S.powConst(S.sub(one, S.mul(a, a)), 0.5)), Math.asin(a[0]));
    },
    acos: (a) => {
      const one = S.constant(1, a.length - 1);
      return S.integrate(S.neg(S.div(S.deriv(a), S.powConst(S.sub(one, S.mul(a, a)), 0.5))), Math.acos(a[0]));
    },
  };

  // env maps variable names to series (Float64Array of length K+1).
  function series(node, env, K) {
    switch (node.type) {
      case 'ConstantNode': return S.constant(Number(node.value), K);
      case 'SymbolNode':
        if (env[node.name]) return env[node.name];
        if (node.name === 'pi' || node.name === 'PI') return S.constant(Math.PI, K);
        if (node.name === 'e' || node.name === 'E') return S.constant(Math.E, K);
        throw new Error('Unknown symbol ' + node.name);
      case 'ParenthesisNode': return series(node.content, env, K);
      case 'OperatorNode': {
        const args = node.args;
        if (args.length === 1) {
          const a = series(args[0], env, K);
          if (node.fn === 'unaryMinus') return S.neg(a);
          if (node.fn === 'unaryPlus') return a;
          break;
        }
        const a = series(args[0], env, K), b = series(args[1], env, K);
        switch (node.fn) {
          case 'add': return S.add(a, b);
          case 'subtract': return S.sub(a, b);
          case 'multiply': return S.mul(a, b);
          case 'divide': return S.div(a, b);
          case 'pow': {
            if (b.slice(1).every((v) => v === 0)) return S.powConst(a, b[0]);
            return S.exp(S.mul(b, S.log(a)));
          }
        }
        break;
      }
      case 'FunctionNode': {
        const name = node.fn.name;
        const impl = FUNCS[name];
        if (!impl) throw new Error(`Derivatives of "${name}" are not supported.`);
        const args = node.args.map((n) => series(n, env, K));
        return impl(...args);
      }
    }
    throw new Error('Unsupported expression for differentiation: ' + node.toString());
  }

  const FACT = [1];
  for (let i = 1; i <= 170; i++) FACT[i] = FACT[i - 1] * i;

  // f^(k)(x0) for a function of one variable.
  function derivative(parsed, k, x0) {
    const c = series(parsed.node, { [parsed.vars[0]]: S.variable(x0, k) }, k);
    return c[k] * FACT[k];
  }

  // Partial derivative ∂f/∂name at the point given by values (in parsed.vars order).
  function partial(parsed, name, values) {
    const env = {};
    parsed.vars.forEach((v, i) => { env[v] = v === name ? S.variable(values[i], 1) : S.constant(values[i], 1); });
    return series(parsed.node, env, 1)[1];
  }

  // For the scalar ODE dy/dt = f(y, t) with y(t0) = y0, return [y, y', y'', ..., y^(K)] at t0,
  // computed from the Taylor coefficients of the solution.
  function odeDerivatives(parsed, y0, t0, K) {
    const Y = new Float64Array(K + 1);
    Y[0] = y0;
    const T = S.variable(t0, K);
    for (let k = 0; k < K; k++) {
      const F = series(parsed.node, { y: Y, t: T }, K);
      Y[k + 1] = F[k] / (k + 1);
    }
    return Array.from(Y, (c, k) => c * FACT[k]);
  }

  // Estimate ||f^(k)||_{inf,[a,b]} by sampling. Returns { value, at, finite }.
  function supNormDerivative(parsed, k, a, b, samples = 2001) {
    let best = 0, at = a, finite = true;
    for (let i = 0; i < samples; i++) {
      const x = a + ((b - a) * i) / (samples - 1);
      const v = Math.abs(derivative(parsed, k, x));
      if (!isFinite(v)) { finite = false; at = x; best = Infinity; break; }
      if (v > best) { best = v; at = x; }
    }
    return { value: best, at, finite };
  }

  return { compile, number, numberList, derivative, partial, odeDerivatives, supNormDerivative, FACT };
})();
