// Root-finding methods and one-step methods for initial-value problems.

const Solvers = (() => {
  // ---------- Root finding ----------

  // Bisection on [a, b] with f(a) f(b) < 0.
  function bisection(f, a, b, tol, maxIter) {
    let fa = f(a), fb = f(b);
    if (!isFinite(fa) || !isFinite(fb)) throw new Error('f must be defined at both a and b.');
    if (fa === 0) return { steps: [{ n: 0, a, b, m: a, fm: 0, width: b - a }], converged: true, root: a };
    if (fb === 0) return { steps: [{ n: 0, a, b, m: b, fm: 0, width: b - a }], converged: true, root: b };
    if (fa * fb > 0) throw new Error(`f(a) = ${fa.toPrecision(4)} and f(b) = ${fb.toPrecision(4)} have the same sign — bisection needs a sign change on [a, b].`);
    const steps = [];
    let converged = false;
    for (let n = 0; n < maxIter; n++) {
      const m = (a + b) / 2, fm = f(m);
      steps.push({ n, a, b, m, fm, width: b - a });
      if (fm === 0 || (b - a) / 2 < tol) { converged = true; break; }
      if (fa * fm < 0) { b = m; fb = fm; } else { a = m; fa = fm; }
    }
    return { steps, converged, root: steps[steps.length - 1].m };
  }

  // Fixed-point iteration x_{n+1} = g(x_n).
  function fixedPoint(g, x0, tol, maxIter) {
    const xs = [x0];
    let converged = false, diverged = false;
    for (let n = 0; n < maxIter; n++) {
      const x = g(xs[n]);
      if (!isFinite(x) || Math.abs(x) > 1e12) { diverged = true; break; }
      xs.push(x);
      if (Math.abs(x - xs[n]) < tol) { converged = true; break; }
    }
    return { xs, converged, diverged };
  }

  // Newton's method x_{n+1} = x_n - f(x_n)/f'(x_n).
  function newton(f, df, x0, tol, maxIter) {
    const xs = [x0];
    let converged = false, diverged = false, stalled = false;
    for (let n = 0; n < maxIter; n++) {
      const x = xs[n], d = df(x);
      if (d === 0 || !isFinite(d)) { stalled = true; break; }
      const next = x - f(x) / d;
      if (!isFinite(next) || Math.abs(next) > 1e12) { diverged = true; break; }
      xs.push(next);
      if (Math.abs(next - x) < tol) { converged = true; break; }
    }
    return { xs, converged, diverged, stalled };
  }

  // Polish a root estimate with a few Newton steps on h (h' given), keeping the best.
  function polish(h, dh, x) {
    let best = x, bestVal = Math.abs(h(x));
    for (let i = 0; i < 8; i++) {
      const d = dh(x);
      if (!d || !isFinite(d)) break;
      x = x - h(x) / d;
      const v = Math.abs(h(x));
      if (!isFinite(v)) break;
      if (v <= bestVal) { best = x; bestVal = v; }
    }
    return best;
  }

  // ---------- One-step methods for y' = f(y, t) (y a vector) ----------

  const add = (a, b, s = 1) => a.map((v, i) => v + s * b[i]);
  const norm = (a) => Math.max(...a.map(Math.abs));

  // Solve G(Y) = Y - (h/2) f(Y, t) - rhs = 0 by Newton's method with a finite-difference Jacobian.
  function solveImplicit(F, rhs, t, h, guess) {
    let Y = guess.slice();
    const n = Y.length;
    const G = (Z) => { const fz = F(Z, t); return Z.map((z, i) => z - (h / 2) * fz[i] - rhs[i]); };
    for (let it = 0; it < 50; it++) {
      const g = G(Y);
      if (norm(g) < 1e-14 * Math.max(1, norm(Y))) return Y;
      // Jacobian
      const Jm = [];
      for (let k = 0; k < n; k++) {
        const eps = 1e-7 * Math.max(1, Math.abs(Y[k]));
        const Z = Y.slice(); Z[k] += eps;
        const gz = G(Z);
        Jm.push(gz.map((v, i) => (v - g[i]) / eps)); // column k
      }
      let delta;
      if (n === 1) delta = [g[0] / Jm[0][0]];
      else {
        // 2x2 solve: J[i][k] = Jm[k][i]
        const a = Jm[0][0], b = Jm[1][0], c = Jm[0][1], d = Jm[1][1];
        const det = a * d - b * c;
        delta = [(d * g[0] - b * g[1]) / det, (-c * g[0] + a * g[1]) / det];
      }
      if (!delta.every(isFinite)) break;
      Y = add(Y, delta, -1);
      if (norm(delta) < 1e-15 * Math.max(1, norm(Y))) return Y;
    }
    const g = G(Y);
    if (norm(g) < 1e-9 * Math.max(1, norm(Y))) return Y;
    throw new Error(`Crank–Nicolson: the implicit equation could not be solved at t = ${t.toPrecision(4)} (try a smaller h).`);
  }

  const METHODS = {
    euler: {
      name: "Euler's method", short: 'Euler', order: 1, implicit: false,
      step: (F, Y, t, h) => add(Y, F(Y, t), h),
    },
    cn: {
      name: 'Crank–Nicolson (trapezoidal) method', short: 'Crank–Nicolson', order: 2, implicit: true,
      step: (F, Y, t, h) => {
        const fY = F(Y, t);
        const rhs = add(Y, fY, h / 2);
        return solveImplicit(F, rhs, t + h, h, add(Y, fY, h));
      },
    },
    improved: {
      name: 'Improved Euler method', short: 'Improved Euler', order: 2, implicit: false,
      step: (F, Y, t, h) => {
        const k1 = F(Y, t);
        const pred = add(Y, k1, h);
        const k2 = F(pred, t + h);
        return Y.map((y, i) => y + (h / 2) * (k1[i] + k2[i]));
      },
    },
    rk4: {
      name: 'Fourth-order Runge–Kutta (RK4)', short: 'RK4', order: 4, implicit: false,
      step: (F, Y, t, h) => {
        const k1 = F(Y, t);
        const k2 = F(add(Y, k1, h / 2), t + h / 2);
        const k3 = F(add(Y, k2, h / 2), t + h / 2);
        const k4 = F(add(Y, k3, h), t + h);
        return Y.map((y, i) => y + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
      },
    },
  };

  // Integrate from t=0 to T with n steps. Returns { ts, Ys, blowup }.
  function solveIVP(method, F, y0, T, n) {
    const h = T / n, m = METHODS[method];
    const ts = [0], Ys = [y0.slice()];
    let blowup = null;
    for (let j = 1; j <= n; j++) {
      const Y = m.step(F, Ys[j - 1], ts[j - 1], h);
      if (!Y.every((v) => isFinite(v) && Math.abs(v) < 1e150)) { blowup = ts[j - 1]; break; }
      ts.push(j * h); Ys.push(Y);
    }
    return { ts, Ys, blowup, h };
  }

  // One step of the method applied to the exact solution, used for the local truncation error:
  // τ_j = (y(t_j) - y(t_{j-1}))/h - F_h(y(t_{j-1}))  (explicit)
  // τ_j = (y(t_j) - y(t_{j-1}))/h - ½(f(y(t_{j-1})) + f(y(t_j)))  (Crank–Nicolson)
  function truncationErrors(method, F, ys, ts, h) {
    const taus = [];
    for (let j = 1; j < ys.length; j++) {
      const diff = ys[j].map((v, i) => (v - ys[j - 1][i]) / h);
      let incr;
      if (method === 'cn') {
        const a = F(ys[j - 1], ts[j - 1]), b = F(ys[j], ts[j]);
        incr = a.map((v, i) => (v + b[i]) / 2);
      } else {
        const next = METHODS[method].step(F, ys[j - 1], ts[j - 1], h);
        incr = next.map((v, i) => (v - ys[j - 1][i]) / h);
      }
      taus.push(norm(diff.map((v, i) => v - incr[i])));
    }
    return taus;
  }

  return { bisection, fixedPoint, newton, polish, METHODS, solveIVP, truncationErrors, norm };
})();
