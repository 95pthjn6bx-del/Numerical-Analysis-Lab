// Initial-value problem window: Euler, Crank–Nicolson, improved Euler and RK4.

const IVP = (() => {
  const { $, fmt, fmtE, texNum, T: Tx, table, card, plot } = UI;
  const root = () => $('#tab-ivp');
  const form = () => $('#ivp-form');
  const METHOD_KEYS = ['euler', 'cn', 'improved', 'rk4'];
  const COLORS = { euler: '--c-m1', cn: '--c-m2', improved: '--c-m3', rk4: '--c-m4' };

  const PRESETS = {
    ysq: { dim: '1', f: 'y^2', y0: '1', exact: '1/(1-t)', T: '1/2', h: '1/8', method: 'euler' },
    decay: { dim: '1', f: '-y', y0: '1', exact: 'exp(-t)', T: '5', h: '0.25', method: 'euler' },
    logistic: { dim: '1', f: 'y*(1-y)', y0: '0.1', exact: '1/(1+9*exp(-t))', T: '10', h: '0.5', method: 'improved' },
    gauss: { dim: '1', f: '-2*t*y', y0: '1', exact: 'exp(-t^2)', T: '3', h: '0.2', method: 'cn' },
    shm: { dim: '2', f1: 'y2', f2: '-y1', y10: '1', y20: '0', exact1: 'cos(t)', exact2: '-sin(t)', T: '20', h: '0.1', method: 'euler' },
    pendulum: { dim: '2', f1: 'y2', f2: '-sin(y1)', y10: '1', y20: '0', exact1: '', exact2: '', T: '20', h: '0.1', method: 'rk4' },
  };

  function showError(msg) {
    const box = $('#ivp-error');
    box.hidden = !msg;
    box.textContent = msg || '';
  }

  // ---------- Problem set-up ----------
  function readProblem(o) {
    if (o.dim === '1') {
      const fy = Expr.compile(o.f, ['y', 't'], 'f(y, t)');
      const y0 = [Expr.number(o.y0, 'y₀')];
      let exact = null;
      if (String(o.exact).trim()) { const ex = Expr.compile(o.exact, ['t'], 'the exact solution y(t)'); exact = (t) => [ex.f(t)]; }
      return { dim: 1, fy, F: (Y, t) => [fy.f(Y[0], t)], y0, exact };
    }
    const f1 = Expr.compile(o.f1, ['y1', 'y2', 't'], 'f₁(y₁, y₂, t)');
    const f2 = Expr.compile(o.f2, ['y1', 'y2', 't'], 'f₂(y₁, y₂, t)');
    const y0 = [Expr.number(o.y10, 'y₁(0)'), Expr.number(o.y20, 'y₂(0)')];
    const e1 = String(o.exact1).trim(), e2 = String(o.exact2).trim();
    if (!!e1 !== !!e2) throw new Error('Give both exact solution components y₁(t) and y₂(t), or neither.');
    let exact = null;
    if (e1) {
      const x1 = Expr.compile(e1, ['t'], 'the exact y₁(t)'), x2 = Expr.compile(e2, ['t'], 'the exact y₂(t)');
      exact = (t) => [x1.f(t), x2.f(t)];
    }
    return { dim: 2, F: (Y, t) => [f1.f(Y[0], Y[1], t), f2.f(Y[0], Y[1], t)], y0, exact };
  }

  // Reference values of y at the grid points (exact solution, or RK4 with many sub-steps).
  function referenceOnGrid(P, T, n) {
    if (P.exact) return { ys: NA.linspace(0, T, n + 1).map(P.exact), kind: 'exact' };
    const m = Math.max(1, Math.min(64, Math.floor(200000 / n)));
    const ref = Solvers.solveIVP('rk4', P.F, P.y0, T, n * m);
    const ys = [];
    for (let j = 0; j <= n && j * m < ref.Ys.length; j++) ys.push(ref.Ys[j * m]);
    return { ys, kind: 'reference', fine: ref, m };
  }

  function referenceAtT(P, T, nMax) {
    if (P.exact) return P.exact(T);
    const nRef = Math.min(400000, Math.max(20000, 4 * nMax));
    const ref = Solvers.solveIVP('rk4', P.F, P.y0, T, nRef);
    return ref.blowup != null ? null : ref.Ys[ref.Ys.length - 1];
  }

  const errNorm = (a, b) => Solvers.norm(a.map((v, i) => v - b[i]));

  // Lipschitz constant estimated as max ||∂f/∂y|| along the reference solution.
  function lipschitz(P, ys, ts) {
    let L = 0;
    for (let j = 0; j < ys.length; j++) {
      const Y = ys[j], t = ts[j];
      if (P.dim === 1) {
        let d;
        try { d = Expr.partial(P.fy, 'y', [Y[0], t]); } catch { d = NaN; }
        if (!isFinite(d)) { const eps = 1e-7 * Math.max(1, Math.abs(Y[0])); d = (P.F([Y[0] + eps], t)[0] - P.F([Y[0] - eps], t)[0]) / (2 * eps); }
        L = Math.max(L, Math.abs(d));
      } else {
        const rows = [0, 0];
        for (let k = 0; k < 2; k++) {
          const eps = 1e-7 * Math.max(1, Math.abs(Y[k]));
          const Yp = Y.slice(), Ym = Y.slice(); Yp[k] += eps; Ym[k] -= eps;
          const fp = P.F(Yp, t), fm = P.F(Ym, t);
          rows[0] += Math.abs((fp[0] - fm[0]) / (2 * eps)); rows[1] += Math.abs((fp[1] - fm[1]) / (2 * eps));
        }
        L = Math.max(L, rows[0], rows[1]);
      }
    }
    return L;
  }

  function update() {
    const o = UI.readForm(form());
    UI.applyVisibility(root());
    try {
      showError('');
      const P = readProblem(o);
      const T = Expr.number(o.T, 'T');
      if (!(T > 0)) throw new Error('The final time T must be positive.');
      const hIn = Expr.number(o.h, 'h');
      if (!(hIn > 0)) throw new Error('The step size h must be positive.');
      const n = Math.max(1, Math.round(T / hIn));
      if (n > 100000) throw new Error('That is more than 100 000 steps — please use a larger h.');
      const h = T / n;
      const method = o.method;

      const sol = Solvers.solveIVP(method, P.F, P.y0, T, n);
      const ref = referenceOnGrid(P, T, n);
      const nCommon = Math.min(sol.Ys.length, ref.ys.length);
      const errs = [];
      for (let j = 0; j < nCommon; j++) errs.push(errNorm(ref.ys[j], sol.Ys[j]));
      const refFinite = ref.ys.every((Y) => Y.every(isFinite));
      const taus = refFinite ? Solvers.truncationErrors(method, P.F, ref.ys, NA.linspace(0, T, n + 1), h) : [];

      const all = {};
      if (o.compare) for (const k of METHOD_KEYS) { try { all[k] = k === method ? sol : Solvers.solveIVP(k, P.F, P.y0, T, n); } catch (e) { all[k] = null; } }

      drawSolution(P, T, sol, ref, method, all);
      drawErrors(P, sol, ref, errs, taus, method, all);
      renderCards(P, o, T, n, h, hIn, method, sol, ref, errs, taus);
    } catch (e) {
      showError(e.message);
      console.error(e);
    }
  }

  // ---------- Plots ----------
  function refCurve(P, T, ref) {
    if (P.exact) {
      const ts = NA.linspace(0, T, 600);
      return { ts, Ys: ts.map(P.exact) };
    }
    const f = ref.fine, stride = Math.max(1, Math.floor(f.Ys.length / 1500));
    const ts = [], Ys = [];
    for (let i = 0; i < f.Ys.length; i += stride) { ts.push(f.ts[i]); Ys.push(f.Ys[i]); }
    return { ts, Ys };
  }

  function drawSolution(P, T, sol, ref, method, all) {
    const rc = refCurve(P, T, ref);
    const refName = P.exact ? 'exact' : 'reference (RK4, small h)';
    const traces = [];
    const comps = P.dim === 1 ? [0] : [0, 1];
    const compName = (i) => (P.dim === 1 ? 'y' : `y${i === 0 ? '₁' : '₂'}`);
    const finiteY = (v) => (isFinite(v) && Math.abs(v) < 1e8 ? v : null);
    const compare = Object.keys(all).length > 0;
    // When comparing methods the reference is drawn in a neutral colour so the method colours stand out.
    const refColor = (i) => UI.cssVar(compare ? (i === 0 ? '--c-node' : '--axis') : i === 0 ? '--c-f' : '--c-w');
    comps.forEach((i) => traces.push({ x: rc.ts, y: rc.Ys.map((Y) => finiteY(Y[i])), name: `${compName(i)}(t) ${refName}`, mode: 'lines',
      line: { color: refColor(i), width: 2.5 } }));
    if (compare) {
      for (const k of METHOD_KEYS) {
        const s = all[k]; if (!s) continue;
        traces.push({ x: s.ts, y: s.Ys.map((Y) => finiteY(Y[0])), name: `${Solvers.METHODS[k].short}${P.dim === 2 ? ' (y₁)' : ''}`, mode: s.ts.length <= 200 ? 'lines+markers' : 'lines',
          line: { color: UI.cssVar(COLORS[k]), width: 1.6, dash: 'dash' }, marker: { size: 5 } });
      }
    } else {
      comps.forEach((i) => traces.push({ x: sol.ts, y: sol.Ys.map((Y) => finiteY(Y[i])), name: `${compName(i)} — ${Solvers.METHODS[method].short}`, mode: sol.ts.length <= 200 ? 'lines+markers' : 'lines',
        line: { color: UI.cssVar(i === 0 ? '--c-p' : '--c-err'), width: 1.8, dash: 'dash' }, marker: { size: 6 } }));
    }
    plot($('#ivp-plot'), traces, {
      title: { text: `Solution on [0, ${fmt(T, 6)}] with h = ${fmt(sol.h, 6)}`, font: { size: 14 } },
      xaxis: { title: { text: 't' } },
    });

    const phaseCard = $('#ivp-phase-card');
    phaseCard.hidden = P.dim !== 2;
    if (P.dim === 2) {
      const pt = [{ x: rc.Ys.map((Y) => finiteY(Y[0])), y: rc.Ys.map((Y) => finiteY(Y[1])), name: refName, mode: 'lines', line: { color: refColor(0), width: 2.5 } }];
      const list = compare ? METHOD_KEYS.filter((k) => all[k]).map((k) => [k, all[k]]) : [[method, sol]];
      for (const [k, s] of list) pt.push({ x: s.Ys.map((Y) => finiteY(Y[0])), y: s.Ys.map((Y) => finiteY(Y[1])), name: Solvers.METHODS[k].short, mode: 'lines',
        line: { color: UI.cssVar(compare ? COLORS[k] : '--c-p'), width: 1.6, dash: 'dash' } });
      pt.push({ x: [sol.Ys[0][0]], y: [sol.Ys[0][1]], name: 'start', mode: 'markers', marker: { color: UI.cssVar('--c-node'), size: 9 } });
      plot($('#ivp-phase-plot'), pt, {
        title: { text: 'Phase plane (y₁, y₂)', font: { size: 14 } },
        xaxis: { title: { text: 'y₁' } }, yaxis: { title: { text: 'y₂' }, scaleanchor: 'x', scaleratio: 1 },
      });
    }
  }

  function drawErrors(P, sol, ref, errs, taus, method, all) {
    const pos = (v) => (v > 0 && isFinite(v) ? v : null);
    const traces = [];
    const compare = Object.keys(all).length > 0;
    if (compare) {
      for (const k of METHOD_KEYS) {
        const s = all[k]; if (!s) continue;
        const e = [];
        for (let j = 0; j < Math.min(s.Ys.length, ref.ys.length); j++) e.push(errNorm(ref.ys[j], s.Ys[j]));
        traces.push({ x: s.ts.slice(0, e.length), y: e.map(pos), name: `|eⱼ| ${Solvers.METHODS[k].short}`, mode: 'lines', line: { color: UI.cssVar(COLORS[k]), width: 2 } });
      }
    } else {
      traces.push({ x: sol.ts.slice(0, errs.length), y: errs.map(pos), name: 'global error |eⱼ| = |y(tⱼ) − Yⱼ|', mode: 'lines+markers', line: { color: UI.cssVar('--c-err'), width: 2 }, marker: { size: 4 } });
      traces.push({ x: sol.ts.slice(1, taus.length + 1), y: taus.map(pos), name: 'local truncation error |τⱼ|', mode: 'lines+markers', line: { color: UI.cssVar('--c-tau'), width: 1.6, dash: 'dot' }, marker: { size: 4 } });
    }
    plot($('#ivp-err-plot'), traces, {
      title: { text: 'Errors (log scale)', font: { size: 14 } },
      xaxis: { title: { text: 't' } },
      yaxis: { type: 'log', exponentformat: 'power' },
    });
  }

  // ---------- Cards ----------
  const FORMULAS = {
    euler: `Y_j = Y_{j-1} + h\\,f(Y_{j-1}, t_{j-1})`,
    cn: `Y_j = Y_{j-1} + \\frac{h}{2}\\big(f(Y_{j-1}, t_{j-1}) + f(Y_j, t_j)\\big)`,
    improved: `\\begin{aligned}\\hat Y_j &= Y_{j-1} + h\\,f(Y_{j-1}, t_{j-1}) &&\\text{(prediction)}\\\\ Y_j &= Y_{j-1} + \\tfrac{h}{2}\\big(f(Y_{j-1}, t_{j-1}) + f(\\hat Y_j, t_j)\\big) &&\\text{(correction)}\\end{aligned}`,
    rk4: `\\begin{aligned}K_1 &= f(Y_{j-1}, t_{j-1}), & K_2 &= f\\big(Y_{j-1} + \\tfrac{h}{2}K_1,\\ t_{j-1} + \\tfrac{h}{2}\\big),\\\\ K_3 &= f\\big(Y_{j-1} + \\tfrac{h}{2}K_2,\\ t_{j-1} + \\tfrac{h}{2}\\big), & K_4 &= f\\big(Y_{j-1} + hK_3,\\ t_j\\big),\\\\ Y_j &= Y_{j-1} + \\tfrac{h}{6}\\big(K_1 + 2K_2 + 2K_3 + K_4\\big) \\end{aligned}`,
  };
  const NOTES = {
    euler: 'Explicit, one evaluation of f per step. Found by replacing the integral of f over [t<sub>j−1</sub>, t<sub>j</sub>] by h times its value at the left end-point.',
    cn: 'Implicit: each step needs the solution of a (possibly nonlinear) equation for Y<sub>j</sub> — here solved by Newton\'s method. Uses the trapezium rule for the integral of f; implicit methods have good stability properties.',
    improved: 'Explicit, two evaluations of f per step: an Euler prediction followed by a trapezium-rule correction.',
    rk4: 'Explicit, four evaluations of f per step. The classical fourth-order Runge–Kutta method.',
  };
  const TAU_DEF = {
    euler: `\\tau_j = \\frac{y(t_j) - y(t_{j-1})}{h} - f(y(t_{j-1}), t_{j-1})`,
    cn: `\\tau_j = \\frac{y(t_j) - y(t_{j-1})}{h} - \\tfrac{1}{2}\\big(f(y(t_{j-1}), t_{j-1}) + f(y(t_j), t_j)\\big)`,
    improved: `\\tau_j = \\frac{y(t_j) - y(t_{j-1})}{h} - F_h(y(t_{j-1}))`,
    rk4: `\\tau_j = \\frac{y(t_j) - y(t_{j-1})}{h} - F_h(y(t_{j-1}))`,
  };

  function vecFmt(Y, sig = 10) { return Y.length === 1 ? fmt(Y[0], sig) : `(${Y.map((v) => fmt(v, sig)).join(', ')})`; }

  function renderCards(P, o, T, n, h, hIn, method, sol, ref, errs, taus) {
    const M = Solvers.METHODS[method];
    const last = sol.Ys[sol.Ys.length - 1];
    const refT = ref.ys[ref.ys.length - 1];
    const blew = sol.blowup != null;
    const eT = !blew && ref.ys.length === n + 1 ? errNorm(refT, last) : null;
    const hNote = Math.abs(h - hIn) > 1e-12 * hIn ? `<p class="note">h was adjusted to T/n = ${fmt(h, 8)} so that n = ${n} steps land exactly on T.</p>` : '';

    const result = card('Result', `
      ${blew ? `<p class="bad">The numerical solution blew up near t = ${fmt(sol.blowup, 5)} — try a smaller h or a shorter T.</p>` : ''}
      <div class="kv big">
        <div>Steps ${Tx('n = T/h')}</div><div>${n} &nbsp;(h = ${fmt(h, 8)})</div>
        <div>${Tx(`Y_n \\approx y(${texNum(T, 5)})`)}</div><div><b>${blew ? '—' : vecFmt(last, 12)}</b></div>
        <div>${Tx(`y(${texNum(T, 5)})`)} (${ref.kind})</div><div>${ref.ys.length === n + 1 ? vecFmt(refT, 12) : '—'}</div>
        <div>Error ${Tx(`|e_n| = |y(T) - Y_n|`)}</div><div><b>${eT != null ? fmtE(eT, 3) : '—'}</b></div>
        <div>Max error ${Tx('\\max_j |e_j|')}</div><div>${errs.length ? fmtE(Math.max(...errs), 3) : '—'}</div>
      </div>
      ${hNote}
      ${ref.kind === 'reference' ? '<p class="note">No exact solution given: the reference is RK4 with a much smaller step. For systems, |·| is the maximum norm.</p>' : P.dim === 2 ? '<p class="note">For systems, |·| is the maximum norm.</p>' : ''}`);

    const methodCard = card(`The method: ${M.name}`, `
      <p>${Tx(FORMULAS[method], true)}</p>
      <p class="note">${NOTES[method]}</p>
      <p>Order ${M.order}: the local truncation error is ${Tx(`\\tau_j = O(h^{${M.order}})`)} and the global error is ${Tx(`|e_n| = O(h^{${M.order}})`)}.</p>
      <p>${Tx(TAU_DEF[method], true)}</p>`);

    // ---------- Error bound ----------
    let boundCard = '';
    const gridTs = NA.linspace(0, T, n + 1);
    if (ref.ys.length === n + 1 && ref.ys.every((Y) => Y.every(isFinite)) && taus.length) {
      const L = lipschitz(P, ref.ys, gridTs);
      const tauMax = Math.max(...taus);
      const growth = (Lf, tt) => (Lf > 1e-12 ? (Math.exp(tt * Lf) - 1) / Lf : tt);
      let Lf, bound, boundTex, condNote = '';
      if (method === 'cn') {
        Lf = L; bound = growth(2 * L, T) * 2; // (e^{2TL}-1)/L = 2 * (e^{2TL}-1)/(2L)
        boundTex = `|e_n| \\le \\frac{e^{2TL} - 1}{L}\\,\\max_{j}|\\tau_j| \\quad (hL \\le 1)`;
        condNote = h * L <= 1 ? `<span class="ok">hL = ${fmt(h * L, 4)} ≤ 1 ✓</span>` : `<span class="bad">hL = ${fmt(h * L, 4)} > 1 — the bound does not apply</span>`;
      } else {
        Lf = method === 'euler' ? L : method === 'improved' ? L + (h * L * L) / 2 : L * (1 + (h * L) / 2 + (h * L) ** 2 / 6 + (h * L) ** 3 / 24);
        bound = growth(Lf, T);
        boundTex = `|e_n| \\le \\frac{e^{TL_F} - 1}{L_F}\\,\\max_{j}|\\tau_j|`;
      }
      bound *= tauMax;
      const LfTex = { euler: 'L_F = L', improved: 'L_F = L + \\tfrac{1}{2}hL^2', rk4: 'L_F = L\\big(1 + \\tfrac{hL}{2} + \\tfrac{(hL)^2}{6} + \\tfrac{(hL)^3}{24}\\big)', cn: '' }[method];

      // A priori estimate of τ from derivatives of the solution (single equations only)
      let apriori = '';
      if (P.dim === 1 && (method === 'euler' || method === 'cn')) {
        let d2 = 0, d3 = 0, ok = true;
        try {
          for (let j = 0; j <= n; j++) {
            const ds = Expr.odeDerivatives(P.fy, ref.ys[j][0], gridTs[j], 3);
            d2 = Math.max(d2, Math.abs(ds[2])); d3 = Math.max(d3, Math.abs(ds[3]));
          }
        } catch { ok = false; }
        if (ok && isFinite(d2) && isFinite(d3)) {
          apriori = method === 'euler'
            ? `<div>Taylor estimate ${Tx("|\\tau_j| \\le \\tfrac{h}{2}\\max|y''|")}</div><div>${fmtE((h / 2) * d2, 3)} &nbsp;<span class="dim">(max|y″| = ${fmt(d2, 5)})</span></div>`
            : `<div>Taylor estimate ${Tx("|\\tau_j| \\approx \\tfrac{h^2}{12}\\max|y'''|")}</div><div>${fmtE((h * h / 12) * d3, 3)} &nbsp;<span class="dim">(max|y‴| = ${fmt(d3, 5)})</span></div>`;
        }
      }

      boundCard = card('Error analysis', `
        <p>If the increment function ${Tx('F_h')} is Lipschitz with constant ${Tx('L_F')}, the global error is controlled by the local truncation errors:</p>
        <p>${Tx(`|e_j| \\le (1 + hL_F)|e_{j-1}| + h|\\tau_j|`, true)}</p>
        <p>${Tx(boundTex, true)}</p>
        <div class="kv">
          <div>Lipschitz constant ${Tx("L \\approx \\max|\\partial f/\\partial y|")}</div><div>${fmt(L, 6)}</div>
          ${LfTex ? `<div>${Tx(LfTex)}</div><div>${fmt(Lf, 6)}</div>` : `<div>Condition</div><div>${condNote}</div>`}
          <div>Computed ${Tx('\\max_j |\\tau_j|')}</div><div>${fmtE(tauMax, 3)}</div>
          ${apriori}
          <div>Bound on ${Tx('|e_n|')}</div><div><b>${method === 'cn' && h * L > 1 ? 'does not apply (hL > 1)' : isFinite(bound) ? fmtE(bound, 3) : '∞'}</b></div>
          <div>Actual ${Tx('\\max_j|e_j|')}</div><div><b>${errs.length ? fmtE(Math.max(...errs), 3) : '—'}</b></div>
        </div>
        <p class="note">L is estimated along the ${ref.kind} solution${P.dim === 2 ? ' using the ∞-norm of the Jacobian' : ''}; τⱼ is computed by substituting the ${ref.kind} solution into the method. The bound is usually pessimistic because of the exponential factor.</p>`);
    }

    // ---------- Convergence table ----------
    const rows = [];
    const ns = [];
    for (let k = 0; k < 6 && n * 2 ** k <= 200000; k++) ns.push(n * 2 ** k);
    const yT = referenceAtT(P, T, ns[ns.length - 1]);
    let prev = null;
    for (const nk of ns) {
      let e = null;
      try {
        const s = Solvers.solveIVP(method, P.F, P.y0, T, nk);
        if (s.blowup == null && yT) e = errNorm(yT, s.Ys[s.Ys.length - 1]);
      } catch { e = null; }
      const ratio = prev != null && e ? prev / e : null;
      rows.push([nk, fmt(T / nk, 6), e != null ? fmtE(e, 3) : '—', ratio ? ratio.toFixed(2) : '', ratio ? Math.log2(ratio).toFixed(2) : '']);
      prev = e;
    }
    // shift the ratio into the row of the larger h (as e_h / e_{h/2})
    for (let i = 0; i < rows.length - 1; i++) { rows[i][3] = rows[i + 1][3]; rows[i][4] = rows[i + 1][4]; }
    if (rows.length) { rows[rows.length - 1][3] = ''; rows[rows.length - 1][4] = ''; }
    const convCard = card('Convergence as h → 0', `
      ${table(['n = T/h', 'h', '|y(T) − Y<sub>n</sub>|', 'e<sub>h</sub> / e<sub>h/2</sub>', 'observed order'], rows)}
      <p class="note">For a method of order p the ratio tends to 2<sup>p</sup> = ${2 ** M.order}. Errors near 1e-13 are dominated by rounding.</p>`);

    // ---------- Step table ----------
    const stride = Math.max(1, Math.ceil((sol.Ys.length - 1) / 400));
    const stepRows = [];
    for (let j = 0; j < sol.Ys.length; j += stride) {
      const Y = sol.Ys[j];
      stepRows.push([j, fmt(sol.ts[j], 8), ...Y.map((v) => fmt(v, 10)), ...(ref.ys[j] ? ref.ys[j].map((v) => fmt(v, 10)) : Y.map(() => '—')),
        errs[j] != null ? fmtE(errs[j], 2) : '—', j > 0 && taus[j - 1] != null ? fmtE(taus[j - 1], 2) : '']);
    }
    const yh = P.dim === 1 ? ['Y<sub>j</sub>', 'y(t<sub>j</sub>)'] : ['Y<sub>1,j</sub>', 'Y<sub>2,j</sub>', 'y<sub>1</sub>(t<sub>j</sub>)', 'y<sub>2</sub>(t<sub>j</sub>)'];
    const stepCard = card('Steps', table(['j', 't<sub>j</sub>', ...yh, '|e<sub>j</sub>|', '|τ<sub>j</sub>|'], stepRows, { scroll: true }) +
      (stride > 1 ? `<p class="note">Showing every ${stride}th step.</p>` : ''));

    UI.html($('#ivp-cards'), [result, methodCard, boundCard, convCard, stepCard].join(''));
  }

  function init() {
    const f = form();
    const run = UI.debounce(update, 350);
    f.addEventListener('input', (e) => { if (e.target.name !== 'preset') run(); });
    f.addEventListener('change', (e) => {
      if (e.target.name === 'preset' && PRESETS[e.target.value]) {
        UI.setForm(f, PRESETS[e.target.value]);
        e.target.value = '';
      }
      update();
    });
    update();
  }

  return { init, update };
})();
