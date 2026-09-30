// Tabs, theme toggle and start-up.

(function () {
  const { $, $$ } = UI;
  const MODULES = { interp: Interpolation, integ: Integration, roots: RootFinding, ivp: IVP };
  const started = {};
  let current = null;

  // Typeset the maths in form labels (<i class="m">...</i>).
  $$('.controls .m').forEach((el) => UI.tex(el, el.textContent));

  // ---------- Theme ----------
  const root = document.documentElement;
  const media = window.matchMedia('(prefers-color-scheme: dark)');
  const effectiveTheme = () => root.dataset.theme || (media.matches ? 'dark' : 'light');

  function applyTheme(rerender) {
    const eff = effectiveTheme();
    root.dataset.effectiveTheme = eff;
    $('#theme-toggle').setAttribute('aria-label', eff === 'dark' ? 'Switch to light mode' : 'Switch to dark mode');
    if (rerender) Object.keys(started).forEach((k) => MODULES[k].update());
  }
  $('#theme-toggle').addEventListener('click', () => {
    const next = effectiveTheme() === 'dark' ? 'light' : 'dark';
    root.dataset.theme = next;
    try { localStorage.setItem('na-theme', next); } catch {}
    applyTheme(true);
  });
  media.addEventListener('change', () => { if (!root.dataset.theme) applyTheme(true); });
  applyTheme(false);

  // ---------- Tabs ----------
  function show(tab) {
    if (!MODULES[tab]) tab = 'interp';
    current = tab;
    $$('.tabs [role="tab"]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.tab === tab)));
    $$('.panel').forEach((p) => (p.hidden = p.id !== 'tab-' + tab));
    if (!started[tab]) {
      started[tab] = true;
      MODULES[tab].init();
    } else {
      // Plotly needs a resize after being un-hidden.
      $$('#tab-' + tab + ' .js-plotly-plot').forEach((el) => Plotly.Plots.resize(el));
    }
    try { localStorage.setItem('na-tab', tab); } catch {}
  }

  $$('.tabs [role="tab"]').forEach((b) => b.addEventListener('click', () => show(b.dataset.tab)));
  let initial = 'interp';
  try { initial = localStorage.getItem('na-tab') || 'interp'; } catch {}
  show(initial);
})();
