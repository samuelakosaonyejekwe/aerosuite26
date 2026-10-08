// Runs before first paint (classic script, no inline code needed under the strict CSP).
(function () {
  var t = 'auto';
  try { t = (JSON.parse(localStorage.getItem('aerosuite26.settings') || '{}').theme) || 'auto'; } catch (e) { /* storage blocked */ }
  if (t === 'auto') t = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  document.documentElement.setAttribute('data-theme', t);
})();
