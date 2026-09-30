// Small shared UI helpers: formatting, tables, maths rendering, plotting.

const UI = (() => {
  const $ = (sel, root = document) => root.querySelector(sel);
  const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];

  // Significant-figure formatting (trailing zeros trimmed).
  function fmt(v, sig = 8) {
    if (v === Infinity || v === -Infinity) return v > 0 ? '∞' : '−∞';
    if (typeof v !== 'number' || isNaN(v)) return '—';
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e-4 && a < 1e8) return String(parseFloat(v.toPrecision(sig)));
    return fmtE(v, sig - 1);
  }
  // Errors in short scientific form, e.g. 2.60e-2
  function fmtE(v, digits = 2) {
    if (v === Infinity) return '∞';
    if (typeof v !== 'number' || isNaN(v)) return '—';
    if (v === 0) return '0';
    return v.toExponential(digits).replace('e+', 'e');
  }
  // A number in LaTeX, e.g. 2.6\times10^{-2}
  function texNum(v, sig = 6) {
    if (!isFinite(v)) return v > 0 ? '\\infty' : '-\\infty';
    const a = Math.abs(v);
    if (v === 0) return '0';
    if (a >= 1e-4 && a < 1e6) return String(parseFloat(v.toPrecision(sig)));
    const [m, e] = v.toExponential(sig - 1).split('e');
    return `${parseFloat(m)}\\times 10^{${parseInt(e, 10)}}`;
  }
  const texFrac = (q) => q.toLatex();

  function tex(el, src, display = false) {
    try { katex.render(src, el, { displayMode: display, throwOnError: false }); } catch { el.textContent = src; }
  }
  // Render all $...$ / $$...$$ inside html string into an element.
  function html(el, content) {
    el.innerHTML = content;
    $$('[data-tex]', el).forEach((n) => tex(n, n.getAttribute('data-tex'), n.hasAttribute('data-display')));
  }
  const T = (src, display = false) =>
    `<span data-tex="${src.replace(/&/g, '&amp;').replace(/"/g, '&quot;')}"${display ? ' data-display' : ''}></span>`;

  function table(headers, rows, opts = {}) {
    const th = headers.map((h) => `<th>${h}</th>`).join('');
    const tr = rows.map((r, i) => `<tr${opts.highlight === i ? ' class="hl"' : ''}>${r.map((c) => `<td>${c}</td>`).join('')}</tr>`).join('');
    return `<div class="table-wrap${opts.scroll ? ' scroll' : ''}"><table class="${opts.cls || ''}"><thead><tr>${th}</tr></thead><tbody>${tr}</tbody></table></div>`;
  }

  function card(title, body, cls = '') {
    return `<section class="card ${cls}"><h3>${title}</h3>${body}</section>`;
  }

  function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }

  function plot(el, traces, layout = {}) {
    // Titles are rendered as an HTML heading above the chart so a wrapping legend never overlaps them.
    let heading = el.previousElementSibling;
    if (!heading || !heading.classList.contains('plot-title')) {
      heading = document.createElement('div');
      heading.className = 'plot-title';
      el.parentNode.insertBefore(heading, el);
    }
    heading.textContent = layout.title ? layout.title.text : '';
    heading.hidden = !layout.title;
    layout = { ...layout };
    delete layout.title;
    const base = {
      margin: { l: 55, r: 20, t: 30, b: 45 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: cssVar('--text'), size: 12 },
      xaxis: { gridcolor: cssVar('--grid'), zerolinecolor: cssVar('--axis'), title: { text: 'x' } },
      yaxis: { gridcolor: cssVar('--grid'), zerolinecolor: cssVar('--axis') },
      legend: { orientation: 'h', x: 0, y: 1.0, yanchor: 'bottom' },
      hovermode: 'closest',
    };
    const merged = { ...base, ...layout, xaxis: { ...base.xaxis, ...(layout.xaxis || {}) }, yaxis: { ...base.yaxis, ...(layout.yaxis || {}) } };
    Plotly.react(el, traces, merged, { responsive: true, displaylogo: false });
  }

  function debounce(fn, ms = 250) {
    let t;
    return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
  }

  // Show only elements whose data-show matches current select values, e.g. data-show="method=piecewise".
  function applyVisibility(root) {
    $$('[data-show]', root).forEach((el) => {
      const ok = el.getAttribute('data-show').split('&').every((cond) => {
        const [name, vals] = cond.split('=');
        const input = $(`[name="${name}"]`, root);
        return input && vals.split('|').includes(input.value);
      });
      el.hidden = !ok;
    });
  }

  function readForm(root) {
    const o = {};
    $$('input, select, textarea', root).forEach((el) => {
      if (!el.name) return;
      o[el.name] = el.type === 'checkbox' ? el.checked : el.value;
    });
    return o;
  }

  function setForm(root, values) {
    for (const [k, v] of Object.entries(values)) {
      const el = $(`[name="${k}"]`, root);
      if (!el) continue;
      if (el.type === 'checkbox') el.checked = !!v; else el.value = v;
    }
  }

  return { $, $$, fmt, fmtE, texNum, texFrac, tex, html, T, table, card, plot, debounce, applyVisibility, readForm, setForm, cssVar };
})();
