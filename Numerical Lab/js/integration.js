// Numerical integration window: Newton–Cotes rules, composite rules and their error analysis.

const Integration = (() => {
  const { $, fmt, fmtE, texNum, texFrac, T, table, card, plot } = UI;
  const root = () => $('#tab-integ');
  const form = () => $('#integ-form');

  const PRESETS = {
    sqrtTrap: { f: 'sqrt(x)', interval: 'general', a: '1/4', b: '1', kind: 'single', N: '1' },
    sqrtSimp: { f: 'sqrt(x)', interval: 'general', a: '1/4', b: '1', kind: 'single', N: '2' },
    dop: { f: 'x^4', interval: 'unit', kind: 'single', N: '2' },
    comp: { f: 'exp(x^2)', interval: 'unit', kind: 'composite', N: '2', mesh: 'uniform', J: 4 },
    sqrtComp: { f: 'sqrt(x)', interval: 'unit', kind: 'composite', N: '1', mesh: 'uniform', J: 8 },
    sinComp: { f: 'sin(x)', interval: 'general', a: '0', b: 'pi', kind: 'composite', N: '1', mesh: 'uniform', J: 4 },
  };
  const RULE_NAMES = { 1: 'trapezium rule', 2: "Simpson's rule", 3: "Simpson's 3/8 rule", 4: "Boole's rule" };
  const ruleName = (N) => RULE_NAMES[N] || `${N + 1}-point Newton–Cotes rule`;

  function showError(msg) {
    const box = $('#integ-error');
    box.hidden = !msg;
    box.textContent = msg || '';
  }

  function derivNorm(fx, k, a, b, samples) {
    try { return Expr.supNormDerivative(fx, k, a, b, samples); } catch (e) { return { value: NaN, error: e.message, finite: false }; }
  }

  function update() {
    const o = UI.readForm(form());
    UI.applyVisibility(root());
    try {
      showError('');
      const fx = Expr.compile(o.f);
      let a = 0, b = 1;
      if (o.interval === 'general') {
        a = Expr.number(o.a, 'a'); b = Expr.number(o.b, 'b');
        if (!(a < b)) throw new Error('Need a < b.');
      }
      const N = parseInt(o.N, 10);
      const rule = NA.newtonCotes(N);
      let ys;
      if (o.kind === 'composite') {
        if (o.mesh === 'custom') {
          const inner = String(o.meshPoints).trim() ? Expr.numberList(o.meshPoints, 'mesh points') : [];
          if (inner.some((y) => y <= a || y >= b)) throw new Error(`Interior mesh points must lie strictly between a = ${fmt(a)} and b = ${fmt(b)}.`);
          ys = NA.checkDistinct([a, ...inner, b], 'Mesh points');
        } else {
          const J = parseInt(o.J, 10);
          if (!(J >= 1)) throw new Error('J must be a positive integer.');
          if (J > 5000) throw new Error('Please use J ≤ 5000.');
          ys = NA.uniformNodes(a, b, J);
        }
      } else ys = [a, b];
      render(fx, o, rule, a, b, ys);
    } catch (e) {
      showError(e.message);
      console.error(e);
    }
  }

  function render(fx, o, rule, a, b, ys) {
    const N = rule.N, d = rule.dop, K = rule.errorConstValue;
    const composite = o.kind === 'composite';
    const J = ys.length - 1;
    const unit = o.interval === 'unit';
    const ab = unit ? '' : '^{[a,b]}';

    // Every sample point must be finite.
    for (let j = 1; j < ys.length; j++) for (const x of rule.nodes) {
      const t = ys[j - 1] + (ys[j] - ys[j - 1]) * x;
      if (!isFinite(fx.f(t))) throw new Error(`f(${fmt(t, 6)}) is not finite — the rule cannot sample f there.`);
    }

    const Q = NA.applyComposite(rule, fx.f, ys);
    const ref = NA.referenceIntegral(fx.f, a, b);
    const E = ref.value - Q;

    drawPlot(fx, rule, ys, a, b, o.shade, composite);

    // ---------- Result card ----------
    const Qsym = composite ? `Q${ab}_{${N},J}(f)` : `Q${ab}_{${N}}(f)`;
    const result = card('Result', `
      <div class="kv big">
        <div>${T(Qsym)} (${ruleName(N)}${composite ? `, composite, J = ${J}` : ''})</div><div><b>${fmt(Q, 12)}</b></div>
        <div>${T(`I(f) = \\int_{${texNum(a, 4)}}^{${texNum(b, 4)}} f(x)\\,dx`)} (reference)</div><div>${fmt(ref.value, 14)}</div>
        <div>Error ${T(`E(f) = I(f) - Q(f)`)}</div><div><b>${fmtE(E, 3)}</b></div>
        <div>Relative error ${T('|E|/|I|')}</div><div>${ref.value !== 0 ? fmtE(Math.abs(E / ref.value), 3) : '—'}</div>
      </div>
      <p class="note">Reference value by adaptive Gauss–Kronrod quadrature (estimated accuracy ${fmtE(ref.errEstimate, 1)}).</p>`);

    // ---------- Rule card: formula, nodes and weights ----------
    const weightsTex = rule.weightsFrac.map((w, i) => `${texFrac(w)}\\,f\\!\\left(${texFrac(rule.nodesFrac[i])}\\right)`).join(' + ');
    let formula = `<p>${T(`Q_{${N}}(f) = \\sum_{i=0}^{${N}} w_i f(x_i) = ${weightsTex}`, true)}</p>
      <p class="note">Nodes ${T('x_i = i/N')} on [0, 1]; weights ${T('w_i = \\int_0^1 L_i(x)\\,dx')} (exact fractions).</p>`;
    if (!unit) formula += `<p>${T(`Q^{[a,b]}_{${N}}(f) = (b-a)\\sum_{i=0}^{${N}} w_i\\,f\\big(a + (b-a)x_i\\big)`, true)}</p>`;
    if (composite) formula += `<p>${T(`Q^{[a,b]}_{${N},J}(f) = \\sum_{j=1}^{J} Q^{[y_{j-1},y_j]}_{${N}}(f) = \\sum_{j=1}^{J} h_j \\sum_{i=0}^{${N}} w_i\\, f(y_{j-1} + h_j x_i)`, true)}</p>`;

    const hasNeg = rule.weights.some((w) => w < 0);
    let nodeTable;
    if (!composite) {
      nodeTable = table(['i', 'x<sub>i</sub> ∈ [0,1]', 'w<sub>i</sub>', ...(unit ? [] : ['t<sub>i</sub> = a + (b−a)x<sub>i</sub>', '(b−a)w<sub>i</sub>']), 'f(node)'],
        rule.nodes.map((x, i) => {
          const t = a + (b - a) * x;
          return [i, rule.nodesFrac[i].toFraction(), rule.weightsFrac[i].toFraction() + ` <span class="dim">(${fmt(rule.weights[i], 6)})</span>`,
            ...(unit ? [] : [fmt(t, 8), fmt((b - a) * rule.weights[i], 8)]), fmt(fx.f(t), 10)];
        }));
    } else {
      nodeTable = table(['i', 'x<sub>i</sub> ∈ [0,1]', 'w<sub>i</sub>'],
        rule.nodes.map((x, i) => [i, rule.nodesFrac[i].toFraction(), rule.weightsFrac[i].toFraction() + ` <span class="dim">(${fmt(rule.weights[i], 6)})</span>`]));
    }
    const ruleCard = card(`The rule: ${ruleName(N)} (N = ${N})`, formula + nodeTable +
      (hasNeg ? '<p class="note warn">Some weights are negative (N ≥ 8) — high-order Newton–Cotes rules become unstable; composite rules are preferred.</p>' : ''));

    // ---------- Error card ----------
    const errCard = composite ? compositeErrorCard(fx, rule, ys, a, b, E) : singleErrorCard(fx, rule, a, b, E, unit);

    // ---------- DoP card ----------
    const dopRows = rule.dopRows.map((r) => [r.r, r.r === 0 ? '1' : r.r === 1 ? 'x' : `x<sup>${r.r}</sup>`, r.I.toFraction(), r.Q.toFraction(),
      r.E.valueOf() === 0 ? '0' : `<b>${r.E.toFraction()}</b>`]);
    const dopCard = card(`Degree of precision`, `
      ${table(['r', 'x<sup>r</sup>', 'I(x<sup>r</sup>)', `Q<sub>${N}</sub>(x<sup>r</sup>)`, `E<sub>${N}</sub>(x<sup>r</sup>)`], dopRows, { highlight: dopRows.length - 1 })}
      <p>Hence the DoP of the ${ruleName(N)} is <b>d = ${d}</b>. ${N % 2 === 1 ? `(N odd ⇒ DoP = N.)` : `(N even ⇒ DoP = N + 1 — one better than expected.)`}</p>
      <p>${T(`E_{${N}}(f) = \\left[\\frac{E_{${N}}(x^{${d + 1}})}{${d + 1}!}\\right] f^{(${d + 1})}(\\xi) = ${texFrac(rule.errorConst)}\\, f^{(${d + 1})}(\\xi)`, true)}</p>`);

    // ---------- Convergence / comparison ----------
    const convCard = composite ? compositeConvergence(fx, rule, a, b, ref.value, o, ys) : ruleComparison(fx, a, b, ref.value);

    // ---------- Subinterval breakdown ----------
    let subCard = '';
    if (composite && J <= 30) {
      const rows = [];
      for (let j = 1; j <= J; j++) {
        const q = NA.applyRule(rule, fx.f, ys[j - 1], ys[j]);
        const I = NA.referenceIntegral(fx.f, ys[j - 1], ys[j]).value;
        rows.push([j, `[${fmt(ys[j - 1], 5)}, ${fmt(ys[j], 5)}]`, fmt(ys[j] - ys[j - 1], 5), fmt(q, 10), fmt(I, 10), fmtE(I - q, 2)]);
      }
      subCard = card('Subintervals', table(['j', '[y<sub>j−1</sub>, y<sub>j</sub>]', 'h<sub>j</sub>', `Q<sub>${N}</sub><sup>[y<sub>j−1</sub>,y<sub>j</sub>]</sup>(f)`, 'exact ∫', `E<sub>${N}</sub><sup>[y<sub>j−1</sub>,y<sub>j</sub>]</sup>(f)`], rows));
    }

    UI.html($('#integ-cards'), [result, errCard, ruleCard, subCard, dopCard, convCard].join(''));
  }

  function normRow(r, k, a, b) {
    const lbl = T(`\\|f^{(${k})}\\|_{\\infty,[${texNum(a, 4)},\\,${texNum(b, 4)}]}`);
    let val;
    if (r.error) val = `could not be computed (${r.error})`;
    else if (!r.finite) val = `∞ — f<sup>(${k})</sup> is unbounded near x = ${fmt(r.at, 4)}, so the bound gives no information`;
    else val = `${fmt(r.value, 6)} (max at x ≈ ${fmt(r.at, 4)})`;
    return `<div>${lbl}</div><div>${val}</div>`;
  }

  function singleErrorCard(fx, rule, a, b, E, unit) {
    const N = rule.N, d = rule.dop, K = rule.errorConstValue, Kabs = Math.abs(K);
    const dn = derivNorm(fx, d + 1, a, b);
    const bound = Kabs * (b - a) ** (d + 2) * dn.value;
    const Ktex = texFrac(rule.errorConst);
    const exact = unit
      ? `E_{${N}}(f) = ${Ktex}\\,f^{(${d + 1})}(\\xi),\\quad \\xi\\in[0,1]`
      : `E^{[a,b]}_{${N}}(f) = ${Ktex}\\,(b-a)^{${d + 2}}\\,f^{(${d + 1})}(\\eta),\\quad \\eta\\in[a,b]`;
    const bnd = unit
      ? `|E_{${N}}(f)| \\le ${texFrac(rule.errorConst.abs())}\\,\\|f^{(${d + 1})}\\|_{\\infty,[0,1]}`
      : `|E^{[a,b]}_{${N}}(f)| \\le ${texFrac(rule.errorConst.abs())}\\,(b-a)^{${d + 2}}\\,\\|f^{(${d + 1})}\\|_{\\infty,[a,b]}`;
    // The implied f^(d+1)(ξ) from the actual error
    const impliedDeriv = E / (K * (b - a) ** (d + 2));
    return card('Error analysis', `
      <p>${T(exact, true)}</p>
      <p>${T(bnd, true)}</p>
      <div class="kv">
        <div>Actual error ${T('E(f)')}</div><div><b>${fmtE(E, 3)}</b></div>
        ${normRow(dn, d + 1, a, b)}
        <div>Theoretical bound on ${T('|E|')}</div><div><b>${dn.finite && !dn.error ? fmtE(bound, 3) : '∞'}</b></div>
        <div>Implied ${T(`f^{(${d + 1})}(${unit ? '\\xi' : '\\eta'})`)}</div><div>${fmt(impliedDeriv, 6)}</div>
      </div>
      <p class="note">The rule is exact (E = 0) for every polynomial of degree ≤ ${d}.</p>`);
  }

  function compositeErrorCard(fx, rule, ys, a, b, E) {
    const N = rule.N, d = rule.dop, Kabs = Math.abs(rule.errorConstValue);
    const hs = ys.slice(1).map((y, j) => y - ys[j]);
    const h = Math.max(...hs), J = hs.length;
    const dn = derivNorm(fx, d + 1, a, b);
    const bound = Kabs * (b - a) * h ** (d + 1) * dn.value;
    // Sharper, subinterval-wise bound  Σ |K| h_j^{d+2} ||f^{(d+1)}||_{[y_{j-1},y_j]}
    let sharp = null;
    if (J <= 400) {
      sharp = 0;
      for (let j = 0; j < J; j++) {
        const r = derivNorm(fx, d + 1, ys[j], ys[j + 1], 101);
        if (!r.finite || r.error) { sharp = Infinity; break; }
        sharp += Kabs * hs[j] ** (d + 2) * r.value;
      }
    }
    const Ktex = texFrac(rule.errorConst);
    return card('Error analysis (composite rule)', `
      <p>${T(`E^{[a,b]}_{${N},J}(f) = \\sum_{j=1}^{J} E^{[y_{j-1},y_j]}_{${N}}(f) = ${Ktex}\\sum_{j=1}^{J} h_j^{${d + 2}}\\, f^{(${d + 1})}(\\eta_j)`, true)}</p>
      <p>${T(`\\big|E^{[a,b]}_{${N},J}(f)\\big| \\le ${texFrac(rule.errorConst.abs())}\\,(b-a)\\,\\|f^{(${d + 1})}\\|_{\\infty,[a,b]}\\,h^{${d + 1}},\\qquad h := \\max_j h_j`, true)}</p>
      <div class="kv">
        <div>Mesh width ${T('h')}, subintervals ${T('J')}</div><div>${fmt(h, 8)}, ${J}</div>
        <div>Actual error ${T('E(f)')}</div><div><b>${fmtE(E, 3)}</b></div>
        ${normRow(dn, d + 1, a, b)}
        <div>Theoretical bound</div><div><b>${dn.finite && !dn.error ? fmtE(bound, 3) : '∞'}</b></div>
        ${sharp !== null ? `<div>Sharper bound ${T(`\\sum_j ${texFrac(rule.errorConst.abs())} h_j^{${d + 2}}\\|f^{(${d + 1})}\\|_{\\infty,[y_{j-1},y_j]}`)}</div><div>${isFinite(sharp) ? fmtE(sharp, 3) : '∞'}</div>` : ''}
      </div>
      <p class="note">The error is O(h<sup>${d + 1}</sup>): halving h should divide the error by about 2<sup>${d + 1}</sup> = ${2 ** (d + 1)}.
      ${!dn.finite ? ' Here f is not smooth enough on [a,b] — try a mesh that isolates the problem point and estimate that subinterval separately.' : ''}</p>`);
  }

  function compositeConvergence(fx, rule, a, b, I, o, ys) {
    const d = rule.dop, Kabs = Math.abs(rule.errorConstValue);
    const dn = derivNorm(fx, d + 1, a, b);
    const J0 = o.mesh === 'uniform' ? ys.length - 1 : 1;
    const rows = [];
    let prev = null;
    for (let k = 0; k < 7; k++) {
      const J = J0 * 2 ** k;
      if (J > 20000) break;
      const h = (b - a) / J;
      const Q = NA.applyComposite(rule, fx.f, NA.uniformNodes(a, b, J));
      const row = { J, h, Q, E: Math.abs(I - Q), bound: Kabs * (b - a) * h ** (d + 1) * dn.value };
      if (prev) prev.ratio = prev.E / row.E;
      rows.push(row); prev = row;
    }
    const tbl = table(['J', 'h', `Q<sub>${rule.N},J</sub>(f)`, '|E|', '|E<sub>h</sub>| / |E<sub>h/2</sub>|', 'bound'],
      rows.map((r) => [r.J, fmt(r.h, 5), fmt(r.Q, 12), fmtE(r.E, 2), r.ratio && isFinite(r.ratio) ? r.ratio.toFixed(2) : '', dn.finite && !dn.error ? fmtE(r.bound, 2) : '∞']));
    return card('Convergence as h → 0 (uniform meshes)', tbl +
      `<p class="note">Expected ratio 2<sup>${d + 1}</sup> = ${2 ** (d + 1)} for smooth f. Errors near 1e-15 are at the level of rounding error.</p>`);
  }

  function ruleComparison(fx, a, b, I) {
    const rows = [];
    for (let N = 1; N <= 8; N++) {
      const r = NA.newtonCotes(N);
      const Q = NA.applyRule(r, fx.f, a, b);
      const dn = derivNorm(fx, r.dop + 1, a, b, 501);
      const bound = Math.abs(r.errorConstValue) * (b - a) ** (r.dop + 2) * dn.value;
      rows.push([N, ruleName(N), r.dop, fmt(Q, 12), fmtE(I - Q, 2), dn.finite && !dn.error ? fmtE(bound, 2) : '∞']);
    }
    return card('Comparison of Newton–Cotes rules on this integral',
      table(['N', 'rule', 'DoP', 'Q<sub>N</sub>(f)', 'E(f)', 'bound'], rows));
  }

  function drawPlot(fx, rule, ys, a, b, shade, composite) {
    const J = ys.length - 1;
    const xs = NA.linspace(a, b, 1201);
    const px = [], py = [], nx = [], ny = [];
    const perSub = Math.max(8, Math.ceil(1200 / J));
    for (let j = 1; j <= J; j++) {
      const l = ys[j - 1], r = ys[j];
      const nodes = rule.nodes.map((x) => l + (r - l) * x);
      const vals = nodes.map(fx.f);
      const p = NA.newtonInterpolant(nodes, vals).eval;
      for (const x of NA.linspace(l, r, perSub)) { px.push(x); py.push(p(x)); }
      nodes.forEach((x, i) => { if (i > 0 || j === 1) { nx.push(x); ny.push(vals[i]); } });
    }
    const pColor = UI.cssVar('--c-p');
    const traces = [];
    traces.push({ x: px, y: py, name: `p_${rule.N} on each subinterval`, mode: 'lines', line: { color: pColor, width: 1.8, dash: 'dash' },
      fill: shade ? 'tozeroy' : 'none', fillcolor: UI.cssVar('--c-fill') });
    traces.push({ x: xs, y: xs.map(fx.f), name: 'f(x)', mode: 'lines', line: { color: UI.cssVar('--c-f'), width: 2.5 } });
    if (nx.length <= 400) traces.push({ x: nx, y: ny, name: 'quadrature nodes', mode: 'markers', marker: { color: UI.cssVar('--c-node'), size: 7 } });
    if (composite && J <= 200) {
      const vx = [], vy = [];
      ys.forEach((y) => { vx.push(y, y, null); vy.push(0, fx.f(y), null); });
      traces.push({ x: vx, y: vy, name: 'mesh points y_j', mode: 'lines', line: { color: UI.cssVar('--axis'), width: 1 }, hoverinfo: 'skip' });
    }
    plot($('#integ-plot'), traces, { title: { text: `${composite ? 'Composite ' + ruleName(rule.N) : ruleName(rule.N).replace(/^./, (c) => c.toUpperCase())}${composite ? `, J = ${J}` : ''} on [${fmt(a, 4)}, ${fmt(b, 4)}]`, font: { size: 14 } } });
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
    update();
  }

  return { init, update };
})();
