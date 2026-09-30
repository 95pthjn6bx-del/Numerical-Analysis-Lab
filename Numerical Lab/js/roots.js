// Root-finding window: bisection, fixed-point iteration and Newton's method.

const RootFinding = (() => {
  const { $, fmt, fmtE, texNum, T, table, card, plot } = UI;
  const root = () => $('#tab-roots');
  const form = () => $('#roots-form');

  const CUBIC = 'x^3 + 4x^2 - 10';
  const PRESETS = {
    cubicBis: { method: 'bisection', f: CUBIC, a: '1', b: '2', tol: '1e-12', maxIter: 60 },
    cubicG3: { method: 'fixed', g: 'sqrt(10/(4+x))', fOrig: CUBIC, a: '1', b: '2', x0: '1.5', tol: '1e-12', maxIter: 50 },
    cubicG2: { method: 'fixed', g: '0.5*sqrt(10-x^3)', fOrig: CUBIC, a: '1', b: '2', x0: '1.5', tol: '1e-12', maxIter: 50 },
    cubicG1: { method: 'fixed', g: 'x - x^3 - 4x^2 + 10', fOrig: CUBIC, a: '1', b: '2', x0: '1.5', tol: '1e-12', maxIter: 50 },
    cos: { method: 'fixed', g: 'cos(x)', fOrig: 'x - cos(x)', a: '0', b: '1', x0: '1', tol: '1e-12', maxIter: 100 },
    sinNewton: { method: 'newton', f: 'sin(x)', a: '2', b: '4', x0: '3', tol: '1e-15', maxIter: 50 },
    sqrt2: { method: 'newton', f: 'x^2 - 2', a: '0', b: '2.5', x0: '1', tol: '1e-15', maxIter: 50 },
    double: { method: 'newton', f: '(x-1)^2', a: '0', b: '2.2', x0: '2', tol: '1e-12', maxIter: 60 },
    cubicNewton: { method: 'newton', f: CUBIC, a: '0.5', b: '2.5', x0: '1', tol: '1e-15', maxIter: 50 },
  };

  let state = null; // last computed run, used by the step slider

  function showError(msg) {
    const box = $('#roots-error');
    box.hidden = !msg;
    box.textContent = msg || '';
  }

  const d1 = (parsed) => (x) => { try { return Expr.derivative(parsed, 1, x); } catch { return NaN; } };

  // Estimated order of convergence from three consecutive errors.
  function orders(errs) {
    const ps = errs.map(() => null);
    for (let n = 1; n + 1 < errs.length; n++) {
      const [e0, e1, e2] = [errs[n - 1], errs[n], errs[n + 1]];
      if ([e0, e1, e2].every((e) => e != null && e > 1e-14) && e1 !== e0) ps[n + 1] = Math.log(e2 / e1) / Math.log(e1 / e0);
    }
    return ps;
  }
  function lastOrder(ps) {
    for (let i = ps.length - 1; i >= 0; i--) if (ps[i] != null && isFinite(ps[i])) return ps[i];
    return null;
  }

  function update() {
    const o = UI.readForm(form());
    UI.applyVisibility(root());
    try {
      showError('');
      const tol = Expr.number(o.tol, 'the tolerance');
      if (!(tol > 0)) throw new Error('The tolerance must be positive.');
      const maxIter = parseInt(o.maxIter, 10);
      if (!(maxIter >= 1 && maxIter <= 500)) throw new Error('Max iterations must be between 1 and 500.');
      const a = Expr.number(o.a, 'a'), b = Expr.number(o.b, 'b');
      if (!(a < b)) throw new Error('Need a < b.');
      if (o.method === 'bisection') state = runBisection(o, a, b, tol, maxIter);
      else if (o.method === 'fixed') state = runFixed(o, a, b, tol, maxIter);
      else state = runNewton(o, a, b, tol, maxIter);
      const slider = $('#roots-step');
      slider.max = state.steps;
      slider.value = state.steps;
      draw(state.steps);
      renderCards();
    } catch (e) {
      state = null;
      showError(e.message);
      console.error(e);
    }
  }

  // ---------- Bisection ----------
  function runBisection(o, a, b, tol, maxIter) {
    const fx = Expr.compile(o.f);
    const res = Solvers.bisection(fx.f, a, b, tol, maxIter);
    const xstar = res.converged ? Solvers.polish(fx.f, d1(fx), res.root) : null;
    const errs = res.steps.map((s) => (xstar == null ? null : Math.abs(s.m - xstar)));
    return { method: 'bisection', fx, a, b, res, xstar, errs, xs: res.steps.map((s) => s.m), steps: res.steps.length - 1 };
  }

  // ---------- Fixed point ----------
  function runFixed(o, a, b, tol, maxIter) {
    const gx = Expr.compile(o.g, ['x'], 'g(x)');
    const fOrig = String(o.fOrig).trim() ? Expr.compile(o.fOrig, ['x'], 'f(x)') : null;
    const x0 = Expr.number(o.x0, 'x₀');
    const res = Solvers.fixedPoint(gx.f, x0, tol, maxIter);
    const gp = d1(gx);
    const xstar = res.converged ? Solvers.polish((x) => gx.f(x) - x, (x) => gp(x) - 1, res.xs[res.xs.length - 1]) : null;
    const errs = res.xs.map((x) => (xstar == null ? null : Math.abs(x - xstar)));

    // Check the two conditions of the fixed-point theorem on [a,b].
    let gmin = Infinity, gmax = -Infinity;
    for (const x of NA.linspace(a, b, 2001)) { const v = gx.f(x); gmin = Math.min(gmin, v); gmax = Math.max(gmax, v); }
    const slack = 1e-12 * Math.max(1, Math.abs(a), Math.abs(b));
    const cond1 = isFinite(gmin) && isFinite(gmax) && gmin >= a - slack && gmax <= b + slack;
    let lam;
    try { lam = Expr.supNormDerivative(gx, 1, a, b); } catch (e) { lam = { value: NaN, finite: false, error: e.message }; }
    const cond2 = lam.finite && lam.value < 1;
    const x0in = x0 >= a && x0 <= b;
    return { method: 'fixed', gx, fOrig, a, b, x0, res, xstar, errs, xs: res.xs, steps: res.xs.length - 1, gmin, gmax, cond1, cond2, lam, x0in, gpStar: xstar != null ? gp(xstar) : null };
  }

  // ---------- Newton ----------
  function runNewton(o, a, b, tol, maxIter) {
    const fx = Expr.compile(o.f);
    const df = d1(fx);
    const x0 = Expr.number(o.x0, 'x₀');
    const res = Solvers.newton(fx.f, df, x0, tol, maxIter);
    const last = res.xs[res.xs.length - 1];
    const xstar = res.converged ? Solvers.polish(fx.f, df, last) : null;
    const errs = res.xs.map((x) => (xstar == null ? null : Math.abs(x - xstar)));
    return { method: 'newton', fx, df, a, b, x0, res, xstar, errs, xs: res.xs, steps: res.xs.length - 1, dfStar: xstar != null ? df(xstar) : null };
  }

  // ---------- Plots ----------
  function draw(k) {
    if (!state) return;
    $('#roots-step-label').textContent = `n = ${k}`;
    if (state.method === 'bisection') drawBisection(k);
    else if (state.method === 'fixed') drawCobweb(k);
    else drawNewton(k);
  }

  function curve(f, lo, hi, n = 800) {
    const xs = NA.linspace(lo, hi, n);
    return { xs, ys: xs.map((x) => { const v = f(x); return isFinite(v) ? v : null; }) };
  }

  function drawBisection(k) {
    const { fx, a, b, res } = state;
    const pad = (b - a) * 0.08;
    const c = curve(fx.f, a - pad, b + pad);
    const finite = c.ys.filter((v) => v != null);
    const span = Math.max(...finite) - Math.min(...finite) || 1;
    const s = res.steps[k];
    const traces = [
      { x: c.xs, y: c.ys, name: 'f(x)', mode: 'lines', line: { color: UI.cssVar('--c-f'), width: 2.5 } },
    ];
    // Bracket "ladder" below the axis: one bar per iteration.
    const lx = [], ly = [], levels = Math.min(k, 14);
    const base = Math.min(0, ...finite);
    for (let n = 0; n <= levels; n++) {
      const st = res.steps[n], y = base - span * 0.06 * (n + 1);
      lx.push(st.a, st.b, null); ly.push(y, y, null);
    }
    traces.push({ x: lx, y: ly, name: 'brackets [aₙ, bₙ]', mode: 'lines+markers', line: { color: UI.cssVar('--c-p'), width: 3 }, marker: { size: 5, color: UI.cssVar('--c-p') }, hoverinfo: 'x' });
    const ms = res.steps.slice(0, k + 1);
    traces.push({ x: ms.map((q) => q.m), y: ms.map((q) => q.fm), name: 'midpoints (mₙ, f(mₙ))', mode: 'markers', marker: { color: UI.cssVar('--c-node'), size: 7 }, text: ms.map((q) => `n = ${q.n}`), hovertemplate: '%{text}<br>m = %{x}<br>f(m) = %{y}<extra></extra>' });
    plot($('#roots-plot'), traces, {
      title: { text: `Bisection: step n = ${k}, bracket [${fmt(s.a, 8)}, ${fmt(s.b, 8)}]`, font: { size: 14 } },
      shapes: [{ type: 'rect', xref: 'x', yref: 'paper', x0: s.a, x1: s.b, y0: 0, y1: 1, fillcolor: UI.cssVar('--c-band'), line: { width: 0 } }],
    });
    drawConvergence();
  }

  function drawCobweb(k) {
    const { gx, a, b, xs } = state;
    const shown = xs.slice(0, k + 1).filter((x) => Math.abs(x) < 1e6);
    let lo = Math.min(a, ...shown), hi = Math.max(b, ...shown);
    const pad = (hi - lo) * 0.08 || 1;
    lo -= pad; hi += pad;
    const c = curve(gx.f, lo, hi);
    const px = [], py = [];
    for (let n = 0; n < k && n + 1 < xs.length; n++) {
      const x = xs[n], gxv = xs[n + 1];
      if (n === 0) { px.push(x); py.push(x); }
      px.push(x, gxv); py.push(gxv, gxv);
    }
    const traces = [
      { x: [lo, hi], y: [lo, hi], name: 'y = x', mode: 'lines', line: { color: UI.cssVar('--axis'), width: 1.5, dash: 'dash' } },
      { x: c.xs, y: c.ys, name: 'y = g(x)', mode: 'lines', line: { color: UI.cssVar('--c-f'), width: 2.5 } },
      { x: px, y: py, name: 'cobweb', mode: 'lines', line: { color: UI.cssVar('--c-p'), width: 1.8 } },
      { x: xs.slice(0, k + 1), y: xs.slice(0, k + 1).map(() => lo + (hi - lo) * 0.0), name: 'iterates xₙ', mode: 'markers', marker: { color: UI.cssVar('--c-node'), size: 6 }, text: xs.slice(0, k + 1).map((_, n) => `n = ${n}`), hovertemplate: '%{text}<br>x = %{x}<extra></extra>' },
    ];
    if (state.xstar != null) traces.push({ x: [state.xstar], y: [state.xstar], name: 'fixed point', mode: 'markers', marker: { color: UI.cssVar('--c-tau'), size: 11, symbol: 'star' } });
    const cy = c.ys.filter((v) => v != null);
    const ylo = Math.max(Math.min(lo, ...cy), lo - 3 * (hi - lo)), yhi = Math.min(Math.max(hi, ...cy), hi + 3 * (hi - lo));
    plot($('#roots-plot'), traces, {
      title: { text: `Fixed-point iteration (cobweb diagram), x${sub(k)} = ${fmt(xs[k], 10)}`, font: { size: 14 } },
      yaxis: { range: [ylo, yhi], title: { text: 'y' } },
      shapes: [{ type: 'rect', xref: 'x', yref: 'paper', x0: a, x1: b, y0: 0, y1: 1, fillcolor: UI.cssVar('--c-band'), line: { width: 0 } }],
    });
    drawConvergence();
  }

  function drawNewton(k) {
    const { fx, df, a, b, xs } = state;
    const shown = xs.slice(0, k + 1).filter((x) => Math.abs(x) < 1e6);
    let lo = Math.min(a, ...shown), hi = Math.max(b, ...shown);
    const pad = (hi - lo) * 0.05;
    lo -= pad; hi += pad;
    const c = curve(fx.f, lo, hi);
    const tx = [], ty = [], vx = [], vy = [];
    for (let n = 0; n < k && n + 1 < xs.length; n++) {
      const x = xs[n], y = fx.f(x), nx = xs[n + 1];
      // tangent from (x_n, f(x_n)) down to (x_{n+1}, 0), extended a little beyond
      const slope = df(x), ext = (nx - x) * 0.15;
      tx.push(x - ext, nx + ext, null); ty.push(y - slope * ext, slope * (nx + ext - x) + y, null);
      vx.push(nx, nx, null); vy.push(0, fx.f(nx), null);
    }
    const traces = [
      { x: c.xs, y: c.ys, name: 'f(x)', mode: 'lines', line: { color: UI.cssVar('--c-f'), width: 2.5 } },
      { x: tx, y: ty, name: 'tangent lines', mode: 'lines', line: { color: UI.cssVar('--c-p'), width: 1.6 } },
      { x: vx, y: vy, name: 'xₙ₊₁ → f(xₙ₊₁)', mode: 'lines', line: { color: UI.cssVar('--axis'), width: 1, dash: 'dot' }, hoverinfo: 'skip' },
      { x: xs.slice(0, k + 1), y: xs.slice(0, k + 1).map(fx.f), name: 'iterates (xₙ, f(xₙ))', mode: 'markers', marker: { color: UI.cssVar('--c-node'), size: 7 }, text: xs.slice(0, k + 1).map((_, n) => `n = ${n}`), hovertemplate: '%{text}<br>x = %{x}<br>f = %{y}<extra></extra>' },
    ];
    const cy = c.ys.filter((v) => v != null);
    plot($('#roots-plot'), traces, {
      title: { text: `Newton's method, x${sub(k)} = ${fmt(xs[k], 15)}`, font: { size: 14 } },
      yaxis: { range: [Math.min(...cy), Math.max(...cy)].map((v, i) => v + (i ? 1 : -1) * 0.05 * (Math.max(...cy) - Math.min(...cy))) },
    });
    drawConvergence();
  }

  const SUBS = '₀₁₂₃₄₅₆₇₈₉';
  const sub = (n) => String(n).split('').map((d) => SUBS[+d]).join('');

  function drawConvergence() {
    const { errs, xs } = state;
    const n = xs.map((_, i) => i);
    const pos = (v) => (v != null && v > 0 && isFinite(v) ? v : null);
    const traces = [];
    if (state.xstar != null) traces.push({ x: n, y: errs.map(pos), name: 'error eₙ = |xₙ − x*|', mode: 'lines+markers', line: { color: UI.cssVar('--c-err'), width: 2 }, marker: { size: 6 } });
    if (state.method === 'bisection') {
      const { a, b } = state;
      traces.push({ x: n, y: n.map((k) => (b - a) / 2 ** (k + 1)), name: 'bound (b − a)/2ⁿ⁺¹', mode: 'lines', line: { color: UI.cssVar('--c-w'), dash: 'dash' } });
    }
    if (state.method === 'fixed' && state.cond1 && state.cond2 && state.x0in && state.xstar != null) {
      const e0 = Math.abs(state.xstar - state.x0);
      traces.push({ x: n, y: n.map((k) => pos(state.lam.value ** k * e0)), name: 'bound λⁿ|x* − x₀|', mode: 'lines', line: { color: UI.cssVar('--c-w'), dash: 'dash' } });
    }
    const resid = state.method === 'fixed' ? xs.map((x) => Math.abs(state.gx.f(x) - x)) : xs.map((x) => Math.abs((state.fx).f(x)));
    traces.push({ x: n, y: resid.map(pos), name: state.method === 'fixed' ? '|g(xₙ) − xₙ|' : '|f(xₙ)|', mode: 'lines+markers', line: { color: UI.cssVar('--c-tau'), width: 1.5, dash: 'dot' }, marker: { size: 4 } });
    plot($('#roots-conv-plot'), traces, {
      title: { text: 'Convergence (log scale)', font: { size: 14 } },
      xaxis: { title: { text: 'iteration n' }, dtick: xs.length > 30 ? undefined : 1 },
      yaxis: { type: 'log', exponentformat: 'power' },
    });
  }

  // ---------- Cards ----------
  function statusLine(s) {
    const r = s.res;
    if (s.method === 'bisection') return r.converged ? `<span class="ok">Converged</span> after ${s.steps} bisections.` : `<span class="bad">Not converged</span> after ${s.steps + 1} bisections — increase the maximum number of iterations.`;
    if (r.converged) return `<span class="ok">Converged</span> after ${s.steps} iterations (|xₙ − xₙ₋₁| < tolerance).`;
    if (r.diverged) return `<span class="bad">Diverged</span> — the iterates grew without bound after ${s.steps} iterations.`;
    if (r.stalled) return `<span class="bad">Stopped</span> — f′(xₙ) = 0 at n = ${s.steps}, so the Newton step is undefined.`;
    return `<span class="bad">Not converged</span> after ${s.steps} iterations.`;
  }

  function renderCards() {
    const s = state;
    const ps = orders(s.errs);
    // The bisection error does not decrease monotonically, so an order estimate is meaningless there.
    const p = s.method === 'bisection' ? null : lastOrder(ps);
    const final = s.xs[s.xs.length - 1];
    const fFinal = s.method === 'fixed' ? (s.fOrig ? s.fOrig.f(final) : null) : s.fx.f(final);

    const result = card('Result', `
      <p>${statusLine(s)}</p>
      <div class="kv big">
        <div>Final approximation ${T(`x_{${s.steps}}`)}</div><div><b>${fmt(final, 15)}</b></div>
        ${fFinal != null ? `<div>${T(`f(x_{${s.steps}})`)}</div><div>${fmtE(fFinal, 3)}</div>` : ''}
        ${s.method === 'fixed' ? `<div>${T(`g(x_{${s.steps}}) - x_{${s.steps}}`)}</div><div>${fmtE(s.gx.f(final) - final, 3)}</div>` : ''}
        ${s.xstar != null ? `<div>Reference root ${T('x^*')}</div><div>${fmt(s.xstar, 16)}</div>` : ''}
        ${p != null ? `<div>Estimated order of convergence</div><div>${p.toFixed(2)}</div>` : ''}
      </div>
      ${s.xstar != null ? '<p class="note">x* is the final iterate refined to full precision, and is used to measure the error eₙ = |xₙ − x*|.</p>' : ''}`);

    let theory, rows, headers;
    if (s.method === 'bisection') {
      theory = card('The method', `
        <p>If ${T('f')} is continuous and ${T('f(a)')}, ${T('f(b)')} have opposite signs, ${T('f')} has a root in ${T('[a,b]')}. Each step halves the bracket and keeps the half where ${T('f')} changes sign:</p>
        <p>${T(`m_n = \\tfrac{1}{2}(a_n + b_n),\\qquad [a_{n+1}, b_{n+1}] = \\begin{cases}[a_n, m_n] & f(a_n)f(m_n) < 0\\\\ [m_n, b_n] & \\text{otherwise}\\end{cases}`, true)}</p>
        <p>${T(`|m_n - x^*| \\le \\frac{b-a}{2^{\\,n+1}}`, true)}</p>
        <p class="note">The error bound halves every step: linear convergence with rate ½, whatever f is. It needs about ${Math.ceil(Math.log2((s.b - s.a) / Math.max(Number(form().tol.value) || 1e-12, 1e-300)))} steps for this tolerance.</p>`);
      headers = ['n', 'a<sub>n</sub>', 'b<sub>n</sub>', 'm<sub>n</sub>', 'f(m<sub>n</sub>)', 'b<sub>n</sub> − a<sub>n</sub>', 'e<sub>n</sub> = |m<sub>n</sub> − x*|', '(b − a)/2<sup>n+1</sup>'];
      rows = s.res.steps.map((q, n) => [n, fmt(q.a, 12), fmt(q.b, 12), fmt(q.m, 14), fmtE(q.fm, 2), fmtE(q.width, 2), s.errs[n] != null ? fmtE(s.errs[n], 2) : '—', fmtE((s.b - s.a) / 2 ** (n + 1), 2)]);
    } else if (s.method === 'fixed') {
      const ok = (c) => (c ? '<span class="ok">✓ holds</span>' : '<span class="bad">✗ fails</span>');
      const allOk = s.cond1 && s.cond2;
      theory = card('Fixed-point theorem check', `
        <p>If (i) ${T('g(x)\\in[a,b]')} for all ${T('x\\in[a,b]')} and (ii) ${T("|g'(x)| \\le \\lambda < 1")} on ${T('[a,b]')}, then ${T('g')} has a unique fixed point ${T('x^*')} in ${T('[a,b]')} and for any ${T('x_0\\in[a,b]')}</p>
        <p>${T(`|x_n - x^*| \\le \\lambda^n\\,|x^* - x_0|.`, true)}</p>
        <div class="kv">
          <div>(i) range of ${T('g')} on ${T(`[${texNum(s.a, 5)}, ${texNum(s.b, 5)}]`)}</div><div>[${fmt(s.gmin, 6)}, ${fmt(s.gmax, 6)}] ${ok(s.cond1)}</div>
          <div>(ii) ${T("\\lambda = \\max_{[a,b]}|g'(x)|")}</div><div>${s.lam.error ? 'could not be computed' : s.lam.finite ? fmt(s.lam.value, 6) : '∞'} ${ok(s.cond2)}</div>
          <div>${T('x_0 \\in [a,b]')}</div><div>${ok(s.x0in)}</div>
          ${s.gpStar != null ? `<div>${T("|g'(x^*)|")}</div><div>${fmt(Math.abs(s.gpStar), 6)}</div>` : ''}
        </div>
        <p>${allOk ? (s.x0in ? 'Both conditions hold, so convergence is guaranteed.' : 'Both conditions hold, but x₀ is outside [a, b], so the theorem does not apply to this starting point.') : 'The theorem does not apply on this interval — the iteration may still converge, but nothing is guaranteed.'}</p>
        <p class="note">Near the fixed point the error shrinks by a factor of about |g′(x*)| each step (linear convergence). The range of g and λ are estimated by sampling 2001 points.</p>`);
      const bnd = allOk && s.x0in && s.xstar != null;
      headers = ['n', 'x<sub>n</sub>', '|x<sub>n</sub> − x<sub>n−1</sub>|', 'e<sub>n</sub>', 'e<sub>n</sub>/e<sub>n−1</sub>', ...(bnd ? ['λ<sup>n</sup>|x* − x<sub>0</sub>|'] : []), ...(s.fOrig ? ['f(x<sub>n</sub>)'] : [])];
      rows = s.xs.map((x, n) => [n, fmt(x, 14), n ? fmtE(Math.abs(x - s.xs[n - 1]), 2) : '', s.errs[n] != null ? fmtE(s.errs[n], 2) : '—',
        n && s.errs[n] != null && s.errs[n - 1] > 0 && s.errs[n] > 1e-15 ? fmt(s.errs[n] / s.errs[n - 1], 4) : '',
        ...(bnd ? [fmtE(s.lam.value ** n * Math.abs(s.xstar - s.x0), 2)] : []), ...(s.fOrig ? [fmtE(s.fOrig.f(x), 2)] : [])]);
    } else {
      const multiple = s.dfStar != null && Math.abs(s.dfStar) < 1e-6;
      theory = card("The method", `
        <p>Newton's method is the fixed-point iteration for ${T("g(x) = x - f(x)/f'(x)")}:</p>
        <p>${T(`x_{n+1} = x_n - \\frac{f(x_n)}{f'(x_n)}`, true)}</p>
        <p>If ${T("f(x^*) = 0")} and ${T("f'(x^*) \\ne 0")} then ${T("g'(x^*) = 0")}, and for ${T('x_0')} close enough to ${T('x^*')} the convergence is quadratic:</p>
        <p>${T(`e_{n+1} \\le K e_n^2,\\qquad e_n = |x_n - x^*|`, true)}</p>
        ${s.dfStar != null ? `<div class="kv"><div>${T("f'(x^*)")}</div><div>${fmtE(s.dfStar, 3)}</div></div>` : ''}
        ${multiple ? `<p class="warn">f′(x*) ≈ 0 — this is a multiple root, so g′(x*) ≠ 0 and convergence is only linear (see the ratio column).</p>` : ''}
        <p class="note">The number of correct digits roughly doubles each step once the iterates are close. f′ is computed exactly by automatic differentiation.</p>`);
      headers = ['n', 'x<sub>n</sub>', 'f(x<sub>n</sub>)', "f′(x<sub>n</sub>)", 'e<sub>n</sub>', 'e<sub>n</sub>/e<sub>n−1</sub>', 'e<sub>n</sub>/e<sub>n−1</sub><sup>2</sup>', 'order'];
      rows = s.xs.map((x, n) => {
        const e = s.errs[n], ep = n ? s.errs[n - 1] : null;
        const good = e != null && ep && e > 1e-15;
        return [n, fmt(x, 16), fmtE(s.fx.f(x), 2), fmt(s.df(x), 6), e != null ? fmtE(e, 2) : '—', good ? fmt(e / ep, 4) : '', good ? fmt(e / (ep * ep), 4) : '', ps[n] != null ? ps[n].toFixed(2) : ''];
      });
    }
    const iterCard = card('Iterations', `${table(headers, rows, { scroll: true })}
      <p class="note">Order of convergence r: e<sub>n+1</sub> ≈ K e<sub>n</sub><sup>r</sup>. Estimated as log(e<sub>n+1</sub>/e<sub>n</sub>) / log(e<sub>n</sub>/e<sub>n−1</sub>); r = 1 is linear, r = 2 quadratic.</p>`);
    UI.html($('#roots-cards'), [result, theory, iterCard].join(''));
  }

  function init() {
    const f = form();
    const run = UI.debounce(update, 300);
    f.addEventListener('input', (e) => { if (e.target.name !== 'preset') run(); });
    f.addEventListener('change', (e) => {
      if (e.target.name === 'preset' && PRESETS[e.target.value]) {
        UI.setForm(f, PRESETS[e.target.value]);
        e.target.value = '';
      }
      update();
    });
    $('#roots-step').addEventListener('input', (e) => draw(parseInt(e.target.value, 10)));
    update();
  }

  return { init, update };
})();
