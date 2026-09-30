// Interpolation window: linear, piecewise-linear and degree-N polynomial interpolation.

const Interpolation = (() => {
  const { $, fmt, fmtE, texNum, T, table, card, plot } = UI;
  const root = () => $('#tab-interp');
  const form = () => $('#interp-form');

  const PRESETS = {
    sqrtGood: { f: 'sqrt(x)', method: 'linear', x0: '1/4', x1: '1' },
    sqrtBad: { f: 'sqrt(x)', method: 'linear', x0: '0', x1: '3/4' },
    expPw: { f: 'exp(x^2)', method: 'piecewise', a: '0', b: '1', pwMesh: 'uniform', J: 8 },
    runge: { f: '1/(25x^2+1)', method: 'degreeN', a: '-1', b: '1', N: 10, nodes: 'uniform', form: 'lagrange' },
    rungeCheb: { f: '1/(25x^2+1)', method: 'degreeN', a: '-1', b: '1', N: 10, nodes: 'chebyshev', form: 'lagrange' },
    sinNewton: { f: 'sin(pi*x)', method: 'degreeN', a: '0', b: '2', nodes: 'custom', customNodes: '0, 0.5, 1.2, 2', form: 'newton' },
  };

  function showError(msg) {
    const box = $('#interp-error');
    box.hidden = !msg;
    box.textContent = msg || '';
  }

  // Derivative sup-norm with a readable description; never throws.
  function derivNorm(fx, k, a, b) {
    try {
      const r = Expr.supNormDerivative(fx, k, a, b);
      return r;
    } catch (e) {
      return { value: NaN, error: e.message };
    }
  }
  const normTex = (k, a, b) => `\\|f^{(${k})}\\|_{\\infty,[${texNum(a, 4)},\\,${texNum(b, 4)}]}`;
  function normText(r, k) {
    if (r.error) return `could not be computed (${r.error})`;
    if (!r.finite) return `∞ — f<sup>(${k})</sup> is unbounded near x = ${fmt(r.at, 4)}, so the bound gives no information`;
    return `${fmt(r.value, 6)} (max at x ≈ ${fmt(r.at, 4)})`;
  }

  function update() {
    const opts = UI.readForm(form());
    UI.applyVisibility(root());
    try {
      showError('');
      const fx = Expr.compile(opts.f);
      if (opts.method === 'linear') renderLinear(fx, opts);
      else if (opts.method === 'piecewise') renderPiecewise(fx, opts);
      else renderDegreeN(fx, opts);
    } catch (e) {
      showError(e.message);
      console.error(e);
    }
  }

  // ---------- common plotting ----------
  function plotMain(fx, p, lo, hi, nodes, title, extraTraces = []) {
    const xs = NA.linspace(lo, hi, 1201);
    const fy = xs.map(fx.f), py = xs.map(p);
    const traces = [
      { x: xs, y: fy, name: 'f(x)', mode: 'lines', line: { color: UI.cssVar('--c-f'), width: 2.5 } },
      { x: xs, y: py, name: 'interpolant p(x)', mode: 'lines', line: { color: UI.cssVar('--c-p'), width: 2, dash: 'dash' } },
      { x: nodes, y: nodes.map(fx.f), name: 'nodes', mode: 'markers', marker: { color: UI.cssVar('--c-node'), size: 8 } },
      ...extraTraces,
    ];
    // Keep wild polynomial oscillations from squashing the plot of f.
    const finite = fy.filter(isFinite);
    const fmin = Math.min(...finite), fmax = Math.max(...finite);
    const pmin = Math.min(...py.filter(isFinite)), pmax = Math.max(...py.filter(isFinite));
    const span = Math.max(fmax - fmin, 1e-12);
    const yr = [Math.max(pmin, fmin - 1.5 * span) , Math.min(pmax, fmax + 1.5 * span)];
    const range = [Math.min(yr[0], fmin) - 0.08 * span, Math.max(yr[1], fmax) + 0.08 * span];
    plot($('#interp-plot'), traces, { title: { text: title, font: { size: 14 } }, yaxis: { range } });
    return xs;
  }

  function plotError(fx, p, lo, hi, w) {
    const xs = NA.linspace(lo, hi, 1201);
    const traces = [{ x: xs, y: xs.map((x) => fx.f(x) - p(x)), name: 'e(x) = f(x) − p(x)', mode: 'lines', line: { color: UI.cssVar('--c-err'), width: 2 } }];
    if (w) traces.push({ x: xs, y: xs.map(w), name: 'w_{N+1}(x)', mode: 'lines', yaxis: 'y2', line: { color: UI.cssVar('--c-w'), width: 1.5, dash: 'dot' } });
    const layout = {
      title: { text: 'Error e(x) = f(x) − p(x)', font: { size: 14 } },
      margin: { l: 65, r: w ? 60 : 20, t: 30, b: 45 },
    };
    if (w) layout.yaxis2 = { overlaying: 'y', side: 'right', showgrid: false, title: { text: 'w_{N+1}(x)' } };
    plot($('#interp-err-plot'), traces, layout);
  }

  function evalCard(fx, p, xstar, lo, hi) {
    if (!String(xstar).trim()) return '';
    const x = Expr.number(xstar, 'x*');
    const fv = fx.f(x), pv = p(x);
    const warn = x < lo || x > hi ? `<p class="note">x* lies outside the interpolation interval — this is extrapolation.</p>` : '';
    return card('Evaluate at x*', `
      <div class="kv">
        <div>${T('x^*')}</div><div>${fmt(x, 10)}</div>
        <div>${T('f(x^*)')}</div><div>${fmt(fv, 12)}</div>
        <div>${T('p(x^*)')}</div><div>${fmt(pv, 12)}</div>
        <div>${T('e(x^*) = f(x^*) - p(x^*)')}</div><div>${fmtE(fv - pv, 3)}</div>
      </div>${warn}`);
  }

  function hideBasis() { $('#interp-basis-card').hidden = true; }

  // ---------- Linear ----------
  function renderLinear(fx, o) {
    hideBasis();
    const x0 = Expr.number(o.x0, 'x₀'), x1 = Expr.number(o.x1, 'x₁');
    if (x0 === x1) throw new Error('x₀ and x₁ must be different.');
    const [lo, hi] = [Math.min(x0, x1), Math.max(x0, x1)];
    const f0 = fx.f(x0), f1 = fx.f(x1);
    if (!isFinite(f0) || !isFinite(f1)) throw new Error('f is not defined at one of the nodes.');
    const lin = NA.linearInterpolant(x0, x1, f0, f1);
    const p = lin.eval;

    plotMain(fx, p, lo, hi, [x0, x1], 'Linear interpolation');
    plotError(fx, p, lo, hi);

    const maxErr = NA.maxAbsOn((x) => fx.f(x) - p(x), lo, hi);
    const d2 = derivNorm(fx, 2, lo, hi);
    const bound = ((hi - lo) ** 2 / 8) * d2.value;
    const c = f0 - lin.slope * x0;

    const cards = [
      card('Interpolant', `
        <p>${T(`p_1(x) = f(x_0) + \\left(\\frac{f(x_1)-f(x_0)}{x_1-x_0}\\right)(x-x_0)`, true)}</p>
        <p>${T(`p_1(x) = ${texNum(f0)} + ${texNum(lin.slope)}\\,(x - ${texNum(x0)}) = ${texNum(lin.slope)}\\,x ${c < 0 ? '-' : '+'} ${texNum(Math.abs(c))}`, true)}</p>`),
      card('Nodes', table(['i', 'x<sub>i</sub>', 'f(x<sub>i</sub>)'], [[0, fmt(x0, 10), fmt(f0, 10)], [1, fmt(x1, 10), fmt(f1, 10)]])),
      card('Error analysis', `
        <p>${T(`e(x) = \\frac{f''(\\xi)}{2}\\,w_2(x),\\quad w_2(x) = (x-x_0)(x-x_1)`, true)}</p>
        <p>${T(`\\|e\\|_{\\infty,[x_0,x_1]} \\le \\frac{(x_1-x_0)^2}{8}\\,\\|f''\\|_{\\infty,[x_0,x_1]}`, true)}</p>
        <div class="kv">
          <div>Actual max error ${T('\\|e\\|_\\infty')}</div><div><b>${fmtE(maxErr.value, 3)}</b> at x ≈ ${fmt(maxErr.at, 4)}</div>
          <div>${T(normTex(2, lo, hi))}</div><div>${normText(d2, 2)}</div>
          <div>Theoretical bound</div><div><b>${d2.finite ? fmtE(bound, 3) : '∞'}</b></div>
        </div>
        <p class="note">Sup norms are estimated by sampling 2001 points of [x₀, x₁]; the actual error on 4001 points.</p>`),
    ];
    UI.html($('#interp-cards'), cards.join('') + evalCard(fx, p, o.xstar, lo, hi));
  }

  // ---------- Piecewise linear ----------
  function buildMesh(o, a, b) {
    if (o.pwMesh === 'custom') {
      const inner = String(o.pwPoints).trim() ? Expr.numberList(o.pwPoints, 'mesh points') : [];
      if (inner.some((y) => y <= a || y >= b)) throw new Error('Interior mesh points must lie strictly between a and b.');
      return NA.checkDistinct([a, ...inner, b], 'Mesh points');
    }
    const J = parseInt(o.J, 10);
    if (!(J >= 1)) throw new Error('J must be a positive integer.');
    if (J > 5000) throw new Error('Please use J ≤ 5000.');
    return NA.uniformNodes(a, b, J);
  }

  function readInterval(o) {
    const a = Expr.number(o.a, 'a'), b = Expr.number(o.b, 'b');
    if (!(a < b)) throw new Error('Need a < b.');
    return [a, b];
  }

  function renderPiecewise(fx, o) {
    hideBasis();
    const [a, b] = readInterval(o);
    const ys = buildMesh(o, a, b);
    const J = ys.length - 1;
    const fs = ys.map(fx.f);
    if (fs.some((v) => !isFinite(v))) throw new Error('f is not defined at every mesh point.');
    const p = NA.piecewiseLinear(ys, fs);

    plotMain(fx, p, a, b, ys, `Piecewise-linear interpolation, J = ${J}`);
    plotError(fx, p, a, b);

    const hs = ys.slice(1).map((y, j) => y - ys[j]);
    const h = Math.max(...hs);
    const maxErr = NA.maxAbsOn((x) => fx.f(x) - p(x), a, b, Math.max(4001, 20 * J + 1));
    // e_h: maximum error at the midpoints z_j
    const eh = Math.max(...hs.map((_, j) => { const z = (ys[j] + ys[j + 1]) / 2; return Math.abs(fx.f(z) - p(z)); }));
    const d2 = derivNorm(fx, 2, a, b);
    const bound = (h * h / 8) * d2.value;

    // Per-subinterval table (only for modest J)
    let subTable = '';
    if (J <= 40) {
      const rows = hs.map((hj, j) => {
        const loc = derivNorm(fx, 2, ys[j], ys[j + 1]);
        const locMax = NA.maxAbsOn((x) => fx.f(x) - p(x), ys[j], ys[j + 1], 201).value;
        return [j + 1, `[${fmt(ys[j], 5)}, ${fmt(ys[j + 1], 5)}]`, fmt(hj, 5), fmt((fs[j + 1] - fs[j]) / hj, 6),
          fmtE(locMax, 2), loc.finite && !loc.error ? fmtE(hj * hj / 8 * loc.value, 2) : '∞'];
      });
      subTable = card('Subintervals', table(['j', '[y<sub>j−1</sub>, y<sub>j</sub>]', 'h<sub>j</sub>', 'slope', 'max |e| on sub.', 'h<sub>j</sub>²/8 ‖f″‖'], rows));
    }

    // Convergence table (uniform meshes, halving h)
    const convRows = [];
    let prev = null;
    const J0 = o.pwMesh === 'uniform' ? J : 8;
    for (let k = 0; k < 5; k++) {
      const Jk = J0 * 2 ** k;
      if (Jk > 20000) break;
      const yk = NA.uniformNodes(a, b, Jk), fk = yk.map(fx.f), pk = NA.piecewiseLinear(yk, fk);
      const hk = (b - a) / Jk;
      let e = 0;
      for (let j = 0; j < Jk; j++) { const z = (yk[j] + yk[j + 1]) / 2; e = Math.max(e, Math.abs(fx.f(z) - pk(z))); }
      convRows.push({ Jk, hk, e, bound: hk * hk / 8 * d2.value });
      if (prev) prev.ratio = prev.e / e;
      prev = convRows[convRows.length - 1];
    }
    const conv = table(['J', 'h', 'e<sub>h</sub>', 'e<sub>h</sub> / e<sub>h/2</sub>', 'bound(h) = h²/8 ‖f″‖'],
      convRows.map((r) => [r.Jk, fmt(r.hk, 5), fmtE(r.e, 2), r.ratio ? r.ratio.toFixed(2) : '', d2.finite ? fmtE(r.bound, 2) : '∞']));

    const cards = [
      card('Error analysis', `
        <p>${T(`\\|f - p_{1,J}\\|_{\\infty,[a,b]} \\le \\frac{1}{8}h^2\\,\\|f''\\|_{\\infty,[a,b]},\\qquad h := \\max_j h_j`, true)}</p>
        <div class="kv">
          <div>Mesh width ${T('h')}</div><div>${fmt(h, 8)}</div>
          <div>Actual max error ${T('\\|f-p_{1,J}\\|_\\infty')}</div><div><b>${fmtE(maxErr.value, 3)}</b> at x ≈ ${fmt(maxErr.at, 4)}</div>
          <div>${T('e_h = \\max_j |(f-p_{1,J})(z_j)|')} (midpoints)</div><div>${fmtE(eh, 3)}</div>
          <div>${T(normTex(2, a, b))}</div><div>${normText(d2, 2)}</div>
          <div>Theoretical bound</div><div><b>${d2.finite ? fmtE(bound, 3) : '∞'}</b></div>
        </div>`),
      card('Convergence as h → 0 (uniform mesh)', conv +
        `<p class="note">If e<sub>h</sub> = C h<sup>α</sup> then e<sub>h</sub>/e<sub>h/2</sub> = 2<sup>α</sup>; a ratio near 4 means O(h²).</p>`),
      subTable,
      J <= 40 ? card('Mesh points', table(['j', 'y<sub>j</sub>', 'f(y<sub>j</sub>)'], ys.map((y, j) => [j, fmt(y, 10), fmt(fs[j], 10)]))) : '',
    ];
    UI.html($('#interp-cards'), evalCard(fx, p, o.xstar, a, b) + cards.join(''));
  }

  // ---------- Degree N ----------
  function buildNodes(o, a, b) {
    if (o.nodes === 'custom') {
      const xs = Expr.numberList(o.customNodes, 'nodes');
      if (xs.length < 2) throw new Error('Enter at least two nodes.');
      return NA.checkDistinct(xs);
    }
    const N = parseInt(o.N, 10);
    if (!(N >= 1)) throw new Error('N must be a positive integer.');
    if (N > 40) throw new Error('Please use N ≤ 40.');
    return o.nodes === 'chebyshev' ? NA.chebyshevNodes(a, b, N) : NA.uniformNodes(a, b, N);
  }

  function renderDegreeN(fx, o) {
    let [a, b] = readInterval(o);
    const xs = buildNodes(o, a, b);
    const N = xs.length - 1;
    const lo = Math.min(a, xs[0]), hi = Math.max(b, xs[N]);
    const fs = xs.map(fx.f);
    if (fs.some((v) => !isFinite(v))) throw new Error('f is not defined at every node.');

    const lag = NA.lagrangeInterpolant(xs, fs);
    const newt = NA.newtonInterpolant(xs, fs);
    const p = o.form === 'newton' ? newt.eval : lag.eval;
    const w = NA.nodePolynomial(xs);
    const nodeLabel = { uniform: 'uniformly spaced', chebyshev: 'Chebyshev', custom: 'custom' }[o.nodes];

    plotMain(fx, p, lo, hi, xs, `Degree-${N} interpolation, ${N + 1} ${nodeLabel} nodes (${o.form === 'newton' ? 'Newton' : 'Lagrange'} form)`);
    plotError(fx, p, lo, hi, o.showW ? w : null);

    // Basis functions plot
    const basisCard = $('#interp-basis-card');
    if (o.form === 'lagrange' && o.showBasis) {
      basisCard.hidden = false;
      const grid = NA.linspace(lo, hi, 801);
      const traces = lag.basis.map((L, j) => ({ x: grid, y: grid.map(L), name: `L${j}`, mode: 'lines', line: { width: 1.5 } }));
      traces.push({ x: xs, y: xs.map(() => 0), mode: 'markers', name: 'nodes', marker: { color: UI.cssVar('--c-node'), size: 7 } });
      plot($('#interp-basis-plot'), traces, { title: { text: 'Lagrange basis functions L_j(x), with L_j(x_k) = δ_jk', font: { size: 14 } }, showlegend: N <= 12 });
    } else basisCard.hidden = true;

    // Error analysis
    const maxErr = NA.maxAbsOn((x) => fx.f(x) - p(x), lo, hi);
    const dN = derivNorm(fx, N + 1, lo, hi);
    const wMax = NA.maxAbsOn(w, lo, hi);
    const bound = dN.value / Expr.FACT[N + 1] * wMax.value;

    // Formula card
    let formula;
    if (o.form === 'newton') {
      const terms = newt.coeffs.map((c, k) => {
        const prod = xs.slice(0, k).map((xi) => `(x ${xi < 0 ? '+' : '-'} ${texNum(Math.abs(xi), 4)})`).join('');
        return `${texNum(c)}${prod}`;
      });
      const shown = N <= 6 ? terms.join(' + ').replace(/\+ -/g, '- ') : terms.slice(0, 3).join(' + ').replace(/\+ -/g, '- ') + ' + \\cdots';
      const dd = newt.table;
      const hdr = ['i', 'x<sub>i</sub>', 'f[x<sub>i</sub>]', ...Array.from({ length: Math.min(N, 8) }, (_, k) => `f[x<sub>i</sub>,…,x<sub>i+${k + 1}</sub>]`)];
      const rows = xs.map((x, i) => [i, fmt(x, 6), ...dd.slice(0, Math.min(N, 8) + 1).map((col) => (i < col.length ? `<span class="${i === 0 ? 'diag' : ''}">${fmt(col[i], 6)}</span>` : ''))]);
      formula = card('Interpolant — Newton divided-difference form', `
        <p>${T(`p_N(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + \\cdots`, true)}</p>
        <p>${T(`f[x_0,\\ldots,x_n] := \\frac{f[x_1,\\ldots,x_n]-f[x_0,\\ldots,x_{n-1}]}{x_n-x_0}`, true)}</p>
        <p>${T(`p_{${N}}(x) = ${shown}`, true)}</p>
        <h4>Divided-difference table</h4>
        ${table(hdr, rows)}
        <p class="note">The coefficients of p<sub>N</sub> are the top entries (row i = 0, highlighted).${N > 8 ? ' Only the first 8 difference columns are shown.' : ''}</p>`);
    } else {
      const Lterms = xs.map((xj, j) => {
        const num = xs.map((xi, i) => (i === j ? '' : `(x ${xi < 0 ? '+' : '-'} ${texNum(Math.abs(xi), 4)})`)).join('');
        const den = xs.map((xi, i) => (i === j ? '' : `(${texNum(xj, 4)} ${xi < 0 ? '+' : '-'} ${texNum(Math.abs(xi), 4)})`)).join('');
        return `L_{${j}}(x) = \\frac{${num}}{${den}}`;
      });
      const sum = fs.map((v, j) => `${texNum(v)}\\,L_{${j}}(x)`).join(' + ').replace(/\+ -/g, '- ');
      formula = card('Interpolant — Lagrange form', `
        <p>${T(`p_N(x) = \\sum_{j=0}^{N} f(x_j)\\,L_j(x),\\qquad L_j(x) = \\prod_{i\\ne j}\\frac{x-x_i}{x_j-x_i}`, true)}</p>
        <p>${T(`p_{${N}}(x) = ${N <= 8 ? sum : sum.split(' + ').slice(0, 3).join(' + ') + ' + \\cdots'}`, true)}</p>
        ${N <= 4 ? Lterms.map((t) => `<p>${T(t, true)}</p>`).join('') : '<p class="note">Individual L<sub>j</sub> are shown for N ≤ 4; tick “Plot the basis functions” to see them.</p>'}`);
    }

    // Monomial form (for reference)
    let mono = '';
    if (N <= 12) {
      const cs = NA.monomialCoefficients(xs, newt.coeffs);
      const terms = cs.map((c, k) => (Math.abs(c) < 1e-14 * Math.max(...cs.map(Math.abs)) ? null : `${texNum(c, 5)}${k === 0 ? '' : k === 1 ? 'x' : `x^{${k}}`}`)).filter(Boolean);
      mono = `<p class="note">Expanded:</p><p>${T(`p_{${N}}(x) = ${terms.join(' + ').replace(/\+ -/g, '- ') || '0'}`, true)}</p>`;
    }

    // Runge / node comparison table: max error vs N for uniform and Chebyshev nodes
    const cmpRows = [];
    for (const n of [2, 4, 6, 8, 10, 12, 16, 20]) {
      const errFor = (nodes) => {
        const pp = NA.newtonInterpolant(nodes, nodes.map(fx.f)).eval;
        return NA.maxAbsOn((x) => fx.f(x) - pp(x), a, b, 2001).value;
      };
      cmpRows.push([n, fmtE(errFor(NA.uniformNodes(a, b, n)), 2), fmtE(errFor(NA.chebyshevNodes(a, b, n)), 2)]);
    }

    const cards = [
      formula.replace('</section>', mono + '</section>'),
      card('Error analysis', `
        <p>${T(`(f-p_N)(x) = \\frac{f^{(N+1)}(\\xi)}{(N+1)!}\\,w_{N+1}(x),\\qquad w_{N+1}(x) = (x-x_0)\\cdots(x-x_N)`, true)}</p>
        <p>${T(`\\Rightarrow\\ \\|f-p_N\\|_\\infty \\le \\frac{\\|f^{(N+1)}\\|_\\infty}{(N+1)!}\\,\\|w_{N+1}\\|_\\infty`, true)}</p>
        <div class="kv">
          <div>Actual max error ${T('\\|f-p_N\\|_\\infty')} on [${fmt(lo, 4)}, ${fmt(hi, 4)}]</div><div><b>${fmtE(maxErr.value, 3)}</b> at x ≈ ${fmt(maxErr.at, 4)}</div>
          <div>${T(normTex(N + 1, lo, hi))}</div><div>${normText(dN, N + 1)}</div>
          <div>${T(`(N+1)! = ${N + 1}!`)}</div><div>${fmt(Expr.FACT[N + 1], 8)}</div>
          <div>${T('\\|w_{N+1}\\|_\\infty')}</div><div>${fmtE(wMax.value, 3)}</div>
          <div>Theoretical bound</div><div><b>${dN.finite && !dN.error ? fmtE(bound, 3) : '∞'}</b></div>
        </div>
        <p class="note">Sup norms are estimated by sampling. Derivatives are computed exactly (automatic differentiation).</p>`),
      card(`Nodes (${nodeLabel})`, table(['i', 'x<sub>i</sub>', 'f(x<sub>i</sub>)'], xs.map((x, i) => [i, fmt(x, 10), fmt(fs[i], 10)]))),
      card(`Uniform vs Chebyshev nodes on [${fmt(a, 4)}, ${fmt(b, 4)}] (Runge's phenomenon)`,
        table(['N', 'max error, uniform', 'max error, Chebyshev'], cmpRows)),
    ];
    UI.html($('#interp-cards'), evalCard(fx, p, o.xstar, lo, hi) + cards.join(''));
  }

  function init() {
    const f = form();
    const run = UI.debounce(update, 300);
    f.addEventListener('input', (e) => { if (e.target.name !== 'preset') { run(); } });
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
