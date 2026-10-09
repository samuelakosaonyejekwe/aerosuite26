// High-fidelity bridge: prepare a case for an external open-source solver, run it (own computer, free
// cloud runner through the user's GitHub repository, or a cluster), read the results back and set them
// beside the native suites.

import { h, clear, add, setKids, icon, num, btn, badge, card, toast, pickFiles, downloadBlob, downloadText } from '../dom.js';
import { state, ls, emit, on } from '../../core/store.js';
import { renderPlot } from '../plots.js';
import { dataTable, kpiGrid } from '../results.js';

let M = null;
const mods = async () => (M ||= Object.fromEntries(await Promise.all([['D', '../../core/bridge/decks.js'], ['R', '../../core/bridge/results.js'], ['Z', '../../core/bridge/zip.js'], ['G', '../../core/bridge/remote.js']].map(async ([k, p]) => [k, await import(p)]))));

// kept for the session only; written to this browser's storage only when the user ticks the box
let sessionToken = '';
const mem = { solver: 'su2', analysis: {}, settings: {}, modelId: '', route: 'local' };

const SOLVER_ICON = { su2: 'flow', openfoam: 'wave', calculix: 'truss' };
const INSTALL = {
  su2: { site: 'https://su2code.github.io/download.html', lines: (n) => ['# 1. Once: download the SU2 binaries for your system from https://su2code.github.io/download.html,', '#    unpack them and add the bin folder to your PATH.', '# 2. Unpack the case and run it (Linux, macOS, or Windows with WSL / Git Bash):', `unzip ${n}.zip -d ${n}`, `cd ${n}`, 'NP=4 bash run.sh', '', '# Plain Windows command prompt instead of step 2:', 'SU2_CFD config.cfg'] },
  calculix: { site: 'https://www.calculix.de', lines: (n) => ['# 1. Once: install CalculiX.', '#    Ubuntu / Debian / WSL:   sudo apt-get install calculix-ccx', '#    Any system with conda:   conda install -c conda-forge calculix', `unzip ${n}.zip -d ${n}`, `cd ${n}`, 'NP=4 bash run.sh', '', '# Or one job at a time (any system):', 'ccx -i static'] },
  openfoam: { site: 'https://www.openfoam.com/download', lines: (n) => ['# 1. Once: install OpenFOAM from openfoam.com (v2312 or later), e.g. on Ubuntu:', '#    curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash', '#    sudo apt-get install openfoam2412-default', '# 2. Load its environment, unpack the case and run it:', 'source /usr/lib/openfoam/openfoam2412/etc/bashrc', `unzip ${n}.zip -d ${n}`, `cd ${n}`, 'NP=4 bash run.sh'] },
};
const HPC = (n, solver) => ['#!/bin/bash', `#SBATCH --job-name=${n.slice(0, 24)}`, '#SBATCH --nodes=1', '#SBATCH --ntasks=32', '#SBATCH --time=12:00:00', '', '# load the solver the way your cluster provides it (ask its documentation for the exact name)', `module load ${solver === 'su2' ? 'su2' : solver === 'openfoam' ? 'openfoam' : 'calculix'}`, '', 'cd "$SLURM_SUBMIT_DIR"', 'NP="$SLURM_NTASKS" bash run.sh'];

export async function render(root, _p, { setCrumb }) {
  setCrumb('High-fidelity bridge');
  const host = h('div'); root.append(host);
  const { D, R, Z, G } = await mods();
  let appVersion = ''; try { const r = await fetch('version.json', { cache: 'no-cache' }); if (r.ok) appVersion = (await r.json()).version || ''; } catch { /* offline: the manifest records "unknown" */ }
  const plots = [], killPlots = () => plots.splice(0).forEach((p) => p.destroy());
  let built = null, zipBytes = null, result = null, abort = null;
  const saved = ls.get('bridgeGh', null) || {};
  if (saved.token && !sessionToken) sessionToken = saved.token;

  const model = () => state.geometry.find((g) => g.id === mem.modelId)?.model || null;
  const analysis = () => mem.analysis[mem.solver] || Object.keys(D.SOLVERS[mem.solver].analyses)[0];
  const defaults = () => D.defaultSettings(state.case, mem.solver, { analysis: analysis(), model: model(), up: state.up });
  const settings = () => (mem.settings[mem.solver] ||= defaults());
  const resetSettings = () => { mem.settings[mem.solver] = defaults(); };
  const invalidate = () => { built = null; zipBytes = null; paintBuild(); paintRun(); paintSteps(); };

  // ------------------------------------------------------------------ step strip
  const stepsHost = h('div', { class: 'steps' });
  function paintSteps() {
    const done = [!!built, !!built, !!built, !!result, !!result && !!state.up.hifi?.ts];
    setKids(stepsHost, [['Choose', 'the solver and the kind of analysis'], ['Review', 'the settings filled in from your case'], ['Download', 'the complete, ready-to-run case'], ['Run', 'on your computer, a free cloud runner or a cluster'], ['Compare', 'the results with the native suites']].map(([t, d], i) => h('a', { class: `step ${done[i] ? 'done' : ''}`, href: `#bridge-step-${i + 1}`, onclick: (ev) => { ev.preventDefault(); document.getElementById(`bridge-step-${i + 1}`)?.scrollIntoView({ behavior: 'smooth', block: 'start' }); } }, h('b', null, t), h('p', null, d))));
  }

  // ------------------------------------------------------------------ 1. choose
  const chooseHost = h('div', { class: 'stack' });
  function paintChoose() {
    const s = D.SOLVERS[mem.solver], usable = state.geometry.filter((g) => (mem.solver === 'su2' ? D.modelHasVolume(g.model) : mem.solver === 'openfoam' ? D.modelHasSurface(g.model) : D.modelHasStructure(g.model)));
    if (mem.modelId && !usable.some((g) => g.id === mem.modelId)) mem.modelId = '';
    const need = mem.solver === 'su2' ? 'volume mesh' : mem.solver === 'openfoam' ? 'closed surface' : 'structural mesh';
    setKids(chooseHost,
      h('div', { class: 'grid g3' }, Object.entries(D.SOLVERS).map(([id, v]) => h('button', { type: 'button', class: 'suite-card', style: { textAlign: 'left', cursor: 'pointer', borderColor: id === mem.solver ? 'var(--accent)' : '' }, 'aria-pressed': String(id === mem.solver), onclick: () => { mem.solver = id; paintChoose(); paintSettings(); invalidate(); } },
        h('div', { class: 'top-r' }, h('span', { class: 'num' }, icon(SOLVER_ICON[id], 20)), h('h3', null, v.label), id === mem.solver ? badge('Selected', 'accent') : null), h('p', null, v.what)))),
      h('div', { class: 'grid g2' },
        h('label', { class: 'stack' }, h('span', { class: 'small muted' }, 'Kind of analysis'), h('select', { class: 'inp', onchange: (ev) => { mem.analysis[mem.solver] = ev.target.value; resetSettings(); paintSettings(); invalidate(); } }, Object.entries(s.analyses).map(([id, label]) => h('option', { value: id, selected: id === analysis() }, label)))),
        h('label', { class: 'stack' }, h('span', { class: 'small muted' }, `Imported ${need} (optional)`), h('select', { class: 'inp', onchange: (ev) => { mem.modelId = ev.target.value; resetSettings(); paintSettings(); invalidate(); } }, h('option', { value: '' }, 'None — generate the model from the case'), usable.map((g) => h('option', { value: g.id, selected: g.id === mem.modelId }, `${g.model.name} (${g.label})`))))),
      h('p', { class: 'muted small' }, usable.length ? `An imported ${need} replaces the generated one.` : h('span', null, `No imported ${need} in this session. The model is generated from the case; to use your own, import it on the `, h('a', { href: '#/geometry' }, 'Geometry & mesh'), ' page first.')));
  }

  // ------------------------------------------------------------------ 2. settings
  const setHost = h('div', { class: 'stack' });
  const visible = (fd, s) => {
    const k = fd.key, a = analysis();
    if (mem.solver === 'su2') {
      const imp = s.mesh === 'imported';
      if (fd.group === 'Imported mesh') return imp;
      if (['dim', 'resolution', 'naca', 'chord_m', 'farfield_chords'].includes(k)) return !imp;
      if (k === 'span_m') return !imp && s.dim === '3-D extruded';
      if (k === 'yplus') return !imp && a !== 'euler';
      if (k === 'reynolds') return a !== 'euler';
      if (k === 'mach') return s.regime === 'compressible';
      if (k === 'V_ms') return s.regime === 'incompressible';
      if (k === 'ref_area_m2' || k === 'ref_length_m') return imp;
    }
    if (mem.solver === 'openfoam') {
      if (k === 'turbulence' || k === 'iterations') return a !== 'pimpleFoam';
      if (k === 'end_time_s') return a === 'pimpleFoam';
      if (k === 'T_K' || k === 'p_Pa') return a === 'rhoSimpleFoam';
    }
    if (mem.solver === 'calculix') {
      const box = s.source === 'wing box';
      if (['element', 't_spar_mm', 't_rib_mm', 'n_span', 'n_chord', 'n_height', 'rib_every'].includes(k)) return box;
      if (k === 'n_modes') return a === 'all' || a === 'modal';
      if (k === 'n_buckle') return a === 'all' || a === 'buckle';
      if (k === 'dyn_time_s') return a === 'dynamic' || a === 'explicit';
    }
    return true;
  };
  // repainting replaces the field that is being edited, which fires its blur/change again: defer and coalesce
  let pending = 0;
  const later = () => { clearTimeout(pending); pending = setTimeout(() => { paintSettings(); invalidate(); }, 0); };
  function paintSettings() {
    const s = settings(), def = defaults(), fields = D.SETTING_FIELDS[mem.solver].filter((fd) => visible(fd, s)), groups = [...new Set(fields.map((fd) => fd.group))];
    const input = (fd) => {
      const id = `bridge-${mem.solver}-${fd.key}`, edited = s[fd.key] !== def[fd.key];
      let el;
      if (fd.type === 'select') el = h('select', { class: `inp ${edited ? 'edited' : ''}`, id, onchange: (ev) => { s[fd.key] = ev.target.value; later(); } }, fd.options.map((o) => h('option', { value: o, selected: o === s[fd.key] }, o)));
      else if (fd.type === 'text') el = h('input', { class: `inp ${edited ? 'edited' : ''}`, id, type: 'text', value: s[fd.key] ?? '', autocomplete: 'off', spellcheck: false, onchange: (ev) => { s[fd.key] = ev.target.value.trim(); later(); } });
      else el = h('input', { class: `inp ${edited ? 'edited' : ''}`, id, type: 'number', step: fd.step || 'any', min: fd.min ?? null, max: fd.max ?? null, value: s[fd.key], inputMode: 'decimal', onchange: (ev) => { const v = Number(ev.target.value); if (ev.target.value === '' || !Number.isFinite(v)) { ev.target.classList.add('invalid'); return; } s[fd.key] = v; later(); } });
      return h('div', { class: 'field' }, h('label', { for: id }, fd.label, fd.unit && fd.unit !== '-' ? h('span', { class: 'u' }, fd.unit) : null), el, fd.help ? h('div', { class: 'help' }, fd.help) : null);
    };
    const nEdited = fields.filter((fd) => s[fd.key] !== def[fd.key]).length, c = state.case;
    setKids(setHost,
      h('div', { class: 'note' }, icon('link'), h('div', null, h('b', null, `Filled in from “${c.meta.name || 'your case'}”. `), `Speed ${num(c.flight.V_ms)} m/s at ${num(c.atm.alt_m)} m, angle of attack ${num(c.flight.alpha_deg)}°, section NACA ${c.wing.airfoil}, material ${c.struct.material}. Change anything below; `, h('a', { href: '#/case' }, 'edit the case'), ' to change it everywhere.')),
      h('div', { class: 'grid g2' }, groups.map((g) => h('div', null, h('h4', null, g), h('div', { class: 'form' }, fields.filter((fd) => fd.group === g).map(input))))),
      h('div', { class: 'row' }, btn('Reset to the case values', () => { resetSettings(); paintSettings(); invalidate(); }, { ic: 'refresh', kind: 'ghost sm', disabled: !nEdited }), nEdited ? h('span', { class: 'muted small' }, `${nEdited} value${nEdited > 1 ? 's' : ''} changed from the case defaults (outlined).`) : null));
  }

  // ------------------------------------------------------------------ 3. build and download
  const buildHost = h('div', { class: 'stack' });
  const kb = (n) => (n > 1e6 ? `${num(n / 1e6, 3)} MB` : `${num(n / 1e3, 3)} kB`);
  async function build() {
    setKids(buildHost, h('div', { class: 'row muted' }, h('i', { class: 'spin' }), 'Writing the decks and the mesh…'));
    await new Promise((r) => setTimeout(r, 30));
    try {
      built = D.buildCase(state.case, { solver: mem.solver, settings: { ...settings(), analysis: analysis() }, model: model(), appVersion, up: state.up });
      zipBytes = await Z.zipWrite(built.files.map((f) => ({ path: f.path, text: f.text, exec: f.exec })));
      toast('The case is ready to download.', 'ok');
    } catch (e) { built = null; zipBytes = null; setKids(buildHost, h('div', { class: 'note bad' }, icon('warn'), h('div', null, h('b', null, 'The case could not be built. '), e.message)), btn('Try again', build, { ic: 'refresh' })); paintRun(); paintSteps(); return; }
    paintBuild(); paintRun(); paintSteps();
  }
  function paintBuild() {
    if (!built) { setKids(buildHost, h('p', { class: 'muted small' }, 'One click writes every input file the solver needs — configuration, mesh, boundary conditions, a run script and a manifest that ties the files to this exact case.'), h('div', { class: 'row' }, btn('Build the case', build, { ic: 'play', kind: 'primary big' }))); return; }
    const m = built.manifest;
    setKids(buildHost,
      h('div', { class: 'note ok' }, icon('check'), h('div', null, h('b', null, `${built.title}. `), `${built.files.length} files, ${kb(zipBytes.length)} zipped.`)),
      h('div', { class: 'row' }, btn(`Download ${built.name}.zip`, () => downloadBlob(`${built.name}.zip`, new Blob([zipBytes], { type: 'application/zip' })), { ic: 'download', kind: 'primary big' }), btn('Rebuild', build, { ic: 'refresh', kind: 'ghost' })),
      built.warnings.length ? h('div', { class: 'note warn' }, icon('warn'), h('div', null, h('b', null, 'Check before relying on the result'), h('ul', { style: { margin: '4px 0 0', paddingLeft: '18px' } }, built.warnings.map((w) => h('li', null, w))))) : null,
      h('div', null, h('h4', null, 'What is in the deck'), h('ul', { class: 'small', style: { margin: '6px 0 0', paddingLeft: '18px' } }, built.notes.map((n) => h('li', null, n)))),
      h('details', null, h('summary', { class: 'small', style: { cursor: 'pointer' } }, 'Files in the archive (download any one on its own)'),
        h('div', { class: 'table-wrap', style: { marginTop: '8px' } }, h('table', { class: 'data' }, h('thead', null, h('tr', null, h('th', null, 'File'), h('th', { class: 'num' }, 'Size'), h('th', null, ''))),
          h('tbody', null, built.files.map((f) => h('tr', null, h('td', { class: 'mono' }, f.path), h('td', { class: 'num' }, kb(new Blob([f.text]).size)), h('td', null, btn('', () => downloadText(f.path.split('/').pop() || 'file', f.text), { ic: 'download', kind: 'sm ghost', title: `Download ${f.path}` }))))))),
        h('p', { class: 'muted small' }, `Traceability: case hash ${m.caseHash.slice(0, 16)}…, app version ${m.appVersion}. The manifest lists a checksum for every file.`)));
  }

  // ------------------------------------------------------------------ 4. run
  const runHost = h('div', { class: 'stack' }), ghStatus = h('div', { class: 'stack' });
  const code = (lines) => { const text = lines.join('\n'); return h('div', { class: 'stack' }, h('pre', { class: 'mono', style: { margin: 0, padding: '10px 12px', borderRadius: '10px', background: 'var(--surface-2)', overflowX: 'auto', fontSize: '.8rem', lineHeight: 1.5 } }, text), h('div', { class: 'row' }, btn('Copy', async () => { try { await navigator.clipboard.writeText(text); toast('Copied.', 'ok', 1800); } catch { toast('Copy is blocked in this browser; select the text instead.', 'bad'); } }, { ic: 'doc', kind: 'sm ghost' }))); };
  const ext = (href, text) => h('a', { href, target: '_blank', rel: 'noopener noreferrer' }, text, ' ', icon('external', 14));
  function paintRun() {
    const name = built?.name || `${mem.solver}-case`, tabs = [['local', 'On your computer'], ['github', 'Free cloud runner (GitHub)'], ['hpc', 'Cluster or HPC']];
    const body = h('div', { class: 'stack' });
    if (!built) body.append(h('p', { class: 'muted small' }, 'Build the case in step 3 first; the commands below then use its exact file name.'));
    if (mem.route === 'local') add(body, [
      h('p', { class: 'small' }, `Best for the section runs and the wing-box model: they take minutes on a laptop. ${D.SOLVERS[mem.solver].label} is free and open source. `, ext(INSTALL[mem.solver].site, 'Get it here')),
      code(INSTALL[mem.solver].lines(name)),
      h('p', { class: 'muted small' }, 'NP is the number of processor cores to use. When it finishes, come back to step 5 and drop the result files (or the whole case folder, zipped).')]);
    else if (mem.route === 'hpc') add(body, [
      h('p', { class: 'small' }, 'For fine 3-D meshes and large-eddy simulation. Copy the zip to the cluster, unpack it, and submit a job that calls the same run script. A typical SLURM script:'),
      code(HPC(name, mem.solver)),
      code([`scp ${name}.zip you@cluster.example.org:~/`, `ssh you@cluster.example.org "unzip ${name}.zip -d ${name} && cd ${name} && sbatch job.sh"`]),
      h('p', { class: 'muted small' }, 'Save the first block as job.sh inside the case folder. Module names differ between clusters. Copy the result files back and drop them in step 5.')]);
    else add(body, [githubPanel()]);
    setKids(runHost, h('div', { class: 'tabs', role: 'tablist' }, tabs.map(([id, label]) => h('button', { class: `tab ${mem.route === id ? 'on' : ''}`, role: 'tab', 'aria-selected': String(mem.route === id), type: 'button', onclick: () => { mem.route = id; paintRun(); } }, label))), body);
  }
  function githubPanel() {
    const repoIn = h('input', { class: 'inp', type: 'text', placeholder: 'owner/repository', value: saved.repo || '', autocomplete: 'off', spellcheck: false, 'aria-label': 'Your GitHub repository' });
    const tokIn = h('input', { class: 'inp', type: 'password', placeholder: 'github_pat_…', value: sessionToken, autocomplete: 'off', spellcheck: false, 'aria-label': 'Personal access token' });
    const keep = h('input', { type: 'checkbox', checked: !!saved.token }), threads = h('input', { class: 'inp', type: 'number', min: 1, max: 64, step: 1, value: 4, style: { width: '90px' }, 'aria-label': 'Processor cores' });
    const remember = () => { saved.repo = repoIn.value.trim(); sessionToken = tokIn.value.trim(); if (keep.checked) saved.token = sessionToken; else delete saved.token; ls.set('bridgeGh', { ...saved }); };
    keep.addEventListener('change', () => remember());
    const client = () => { remember(); const { owner, repo } = G.parseRepo(repoIn.value); return G.createGithub({ token: sessionToken, owner, repo }); };
    const say = (kind, ...kids) => setKids(ghStatus, h('div', { class: `note ${kind}` }, kind === 'busy' ? h('i', { class: 'spin' }) : icon(kind === 'ok' ? 'check' : kind === 'bad' || kind === 'warn' ? 'warn' : 'info'), h('div', null, kids)));
    const check = async () => { try { const gh = client(); say('busy', 'Checking…'); const info = await gh.repoInfo(), wf = await gh.workflow(); say('ok', h('b', null, `Connected to ${info.fullName}. `), `The solve workflow is ${wf.state}. Default branch: ${info.defaultBranch}.`); } catch (e) { say('bad', e.message); } };
    const go = async () => {
      if (!built || !zipBytes) { toast('Build the case in step 3 first.', 'info'); return; }
      abort?.abort(); abort = new AbortController(); const signal = abort.signal;
      try {
        const gh = client(); let runUrl = '';
        const status = (msg, run) => { if (run?.html_url) runUrl = run.html_url; say('busy', msg, runUrl ? [' ', ext(runUrl, 'Open the run on GitHub')] : null, ' ', btn('Stop watching', () => abort?.abort(), { kind: 'sm ghost' })); };
        const sub = await G.submitCase(gh, { zipBytes, name: built.name, solver: built.solver, threads: Math.max(1, Math.round(Number(threads.value) || 4)), signal, onStatus: status });
        status(sub.route === 'inline' ? 'The case was sent with the request. Waiting for the run…' : `The case was uploaded to branch ${sub.branch}. Waiting for the run…`);
        const run = await G.followRun(gh, sub, { signal, onStatus: status });
        runUrl = run.html_url || runUrl;
        let files = [];
        try { const arts = await G.fetchArtifacts(gh, run.id, { signal, onStatus: status }); for (const a of arts) for (const f of await Z.zipRead(a.bytes)) files.push({ name: f.name, bytes: f.bytes }); }
        catch (e) { if (e?.name === 'AbortError') throw e; say('warn', h('b', null, `The run finished (${run.conclusion}) but its results could not be downloaded here. `), e.message, ' ', ext(runUrl, 'Open the run'), ', download the artifact at the foot of its page and drop the zip in step 5.'); return; }
        say(run.conclusion === 'success' ? 'ok' : 'warn', h('b', null, run.conclusion === 'success' ? 'The run finished and its results are loaded below. ' : `The run ended with “${run.conclusion}”. Whatever it wrote is loaded below; the solver log explains what happened. `), ext(runUrl, 'Open the run on GitHub'));
        await ingest(files, true);
      } catch (e) { if (e?.name === 'AbortError') say('info', 'Stopped watching. The run itself continues on GitHub; download its artifact there and drop it in step 5.'); else say('bad', e.message); }
    };
    return h('div', { class: 'stack' },
      h('p', { class: 'small' }, 'GitHub runs the solver for you on its hosted machines (free for public repositories, a monthly allowance for private ones). The case goes to your own repository, the workflow there installs the solver, runs it and hands the results back. Nothing is committed to your main branch.'),
      h('details', null, h('summary', { class: 'small', style: { cursor: 'pointer' } }, 'One-time setup (about five minutes)'),
        h('ol', { class: 'small', style: { paddingLeft: '20px' } },
          h('li', null, 'Put this app in a GitHub repository of your own (fork or copy it). It already contains the workflow file ', h('span', { class: 'mono' }, '.github/workflows/solve.yml'), '; make sure it is on the default branch and that Actions are enabled.'),
          h('li', null, 'On GitHub open Settings → Developer settings → Personal access tokens → Fine-grained tokens → Generate new token. ', ext('https://github.com/settings/personal-access-tokens/new', 'Open that page')),
          h('li', null, 'Choose “Only select repositories” and pick that one repository. Under Repository permissions set:', h('ul', null, G.TOKEN_PERMISSIONS.map((p) => h('li', null, p)))),
          h('li', null, 'Give it a short expiry, generate it, and paste it below.'))),
      h('div', { class: 'grid g2' },
        h('label', { class: 'stack' }, h('span', { class: 'small muted' }, 'Your repository'), repoIn),
        h('label', { class: 'stack' }, h('span', { class: 'small muted' }, 'Personal access token'), tokIn)),
      h('div', { class: 'row small' }, h('label', { class: 'row' }, keep, 'Remember the token in this browser'), h('label', { class: 'row' }, 'Cores', threads)),
      h('div', { class: 'note' }, icon('shield'), h('div', null, 'The token is sent only to api.github.com, straight from this page. It stays in memory until you close the tab unless you tick the box, which saves it in this browser’s storage on this device only. It is never part of a saved project or report. Untick the box to erase a saved token.')),
      h('div', { class: 'row' }, btn('Send the case and run it', go, { ic: 'play', kind: 'primary', disabled: !built }), btn('Check connection', check, { ic: 'link' }), !built ? h('span', { class: 'muted small' }, 'Build the case in step 3 first.') : h('span', { class: 'muted small' }, `${kb(zipBytes.length)} — ${Z.base64Encode(zipBytes).length <= G.INLINE_LIMIT ? 'small enough to travel with the request' : `will be uploaded to branch ${G.CASE_BRANCH}`}.`)),
      ghStatus,
      h('p', { class: 'muted small' }, 'Without a token you can still use the cloud runner: commit the unpacked case as cases/<name>/ in the repository, open its Actions tab, choose “Solve”, press “Run workflow”, then download the artifact and drop it in step 5.'));
  }

  // ------------------------------------------------------------------ 5. results
  const resHost = h('div', { class: 'stack' }), busy = h('div', { class: 'row muted', hidden: true }, h('i', { class: 'spin' }), h('span', null, 'Reading…'));
  async function ingest(files, fromRun = false) {
    busy.hidden = false; await new Promise((r) => setTimeout(r, 30));
    try {
      const flat = [];
      const expand = async (name, bytes, depth) => { if (/\.zip$/i.test(name) && depth < 3) { for (const f of await Z.zipRead(bytes)) await expand(f.name, f.bytes, depth + 1); } else flat.push({ name, bytes }); };
      for (const f of files) await expand(f.name, f.bytes ?? new Uint8Array(await f.arrayBuffer()), 0);
      const r = R.readResults(flat, { manifest: built?.manifest || null });
      if (!r.kpis.length && !r.plots.length) { toast(flat.length ? 'None of these files is a result this page reads. See the list of expected files.' : 'The archive was empty.', 'bad', 7000); if (!result) paintResults(); return; }
      result = r; paintResults(); paintSteps();
      toast(`${r.sources.filter((s) => s.ok && s.kind !== 'manifest').length} result file(s) read.`, 'ok');
      if (fromRun) document.getElementById('bridge-step-5')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    } catch (e) { toast(`The results could not be read: ${e.message}`, 'bad', 8000); } finally { busy.hidden = true; }
  }
  const take = async () => ingest(await pickFiles({ multiple: true }));
  const drop = h('div', { class: 'drop', tabIndex: 0, role: 'button', onclick: take, onkeydown: (ev) => { if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); take(); } }, ondragover: (ev) => { ev.preventDefault(); drop.classList.add('over'); }, ondragleave: () => drop.classList.remove('over'), ondrop: (ev) => { ev.preventDefault(); drop.classList.remove('over'); ingest([...ev.dataTransfer.files]); } },
    icon('upload', 32), h('b', null, 'Drop the result files here, or click to choose'), h('span', { class: 'muted small' }, 'A zip of the whole case folder or the downloaded artifact works too. SU2: history.csv, surface_flow.vtu, flow.vtu · OpenFOAM: postProcessing files and log.* · CalculiX: *.dat and *.frd'), busy);
  function store() {
    const rec = R.hifiRecord(result);
    state.up.hifi = { ...(state.up.hifi || {}), ...rec };
    ls.set('up', state.up); emit('run', { suite: 'hifi', analysis: rec.solver });
    toast('High-fidelity values stored on the data bus under “hifi”.', 'ok'); paintResults(); paintSteps();
  }
  function paintResults() {
    killPlots(); clear(resHost);
    if (!result) { resHost.append(h('p', { class: 'muted small' }, built ? `Expected from this case: ${built.expected.join(', ')}.` : 'Nothing read yet.')); return; }
    const r = result, m = r.manifest, rows = R.compare(r.values, state.up, { manifest: m }), stale = m && m.caseHash !== D.caseHash(state.case);
    resHost.append(h('div', { class: 'row' }, h('h4', { class: 'grow' }, `${D.SOLVERS[r.solver]?.label || 'Solver'} results${m ? ` — ${m.title}` : ''}`), btn('Clear', () => { result = null; paintResults(); paintSteps(); }, { ic: 'close', kind: 'sm ghost' })));
    if (!m) resHost.append(h('div', { class: 'note' }, icon('info'), h('div', null, 'No manifest.json came with these files, so they cannot be tied to a particular case. Include it (it is in the case zip) for full traceability.')));
    if (stale) resHost.append(h('div', { class: 'note warn' }, icon('warn'), h('div', null, h('b', null, 'These results belong to a different version of the case. '), 'The case has been edited since the deck was built (or the files come from another case), so the comparison below is not like for like.')));
    if (r.warnings.length) resHost.append(h('div', { class: 'note warn' }, icon('warn'), h('div', null, h('b', null, r.warnings.length === 1 ? 'Note' : 'Notes'), h('ul', { style: { margin: '4px 0 0', paddingLeft: '18px' } }, r.warnings.map((w) => h('li', null, w))))));
    resHost.append(kpiGrid(r.kpis));
    if (r.plots.length) { const grid = h('div', { class: 'plots' }); resHost.append(grid); for (const p of r.plots) { try { plots.push(renderPlot(grid, p)); } catch (e) { grid.append(h('div', { class: 'note bad' }, `Chart “${p.title}” could not be drawn: ${e.message}`)); } } }
    // comparison with the native suites
    const cmp = h('div', { class: 'stack' });
    if (rows.length) {
      const have = rows.filter((x) => x.native != null);
      cmp.append(h('p', { class: 'muted small' }, have.length ? 'Each high-fidelity value beside what the built-in suites computed for the same case. The ratio is the factor that would bring the native model onto the high-fidelity answer — use it to calibrate, or treat a large difference as a prompt to check both models.' : 'The native suites have not been run for this case yet, so there is nothing to compare against. Run them (Integrated run does all at once) and come back.'),
        dataTable({ title: 'High-fidelity against native', columns: ['Quantity', 'Unit', 'High-fidelity', 'Native suite', 'Difference', 'Difference [%]', 'Ratio', 'Native from', 'Note'], rows: rows.map((x) => [x.label, x.unit, x.hifi, x.native, x.diff, x.pct, x.ratio, x.suite, x.note]) }));
      const cp = R.comparePlot(rows); if (cp) { const g = h('div', { class: 'plots' }); cmp.append(g); try { plots.push(renderPlot(g, cp)); } catch { /* table already shows it */ } }
      if (!have.length) cmp.append(h('div', { class: 'row' }, h('a', { class: 'btn', href: '#/integrated' }, icon('graph', 18), h('span', null, 'Go to Integrated run'))));
    } else cmp.append(h('p', { class: 'muted small' }, 'These files hold no quantity that a native suite also computes.'));
    const storedAt = state.up.hifi?.ts;
    resHost.append(card('Comparison with the native suites', cmp),
      card('Keep these values', h('div', { class: 'stack' }, h('p', { class: 'muted small' }, 'Storing puts the numbers on the shared data bus under the name “hifi” (for example hifi.CL, hifi.tip_deflection_m), together with the solver, the analysis and the case hash. They are saved with the project and appear in exported project files.'),
        h('div', { class: 'row' }, btn('Store on the data bus', store, { ic: 'link', kind: 'primary' }), btn('Download summary (JSON)', () => downloadText(`${m?.name || 'hifi'}-results.json`, JSON.stringify({ manifest: m, values: r.values, kpis: r.kpis, comparison: rows, sources: r.sources }, null, 1), 'application/json'), { ic: 'download' }),
          storedAt ? h('span', { class: 'muted small' }, `Last stored ${new Date(storedAt).toLocaleString()} (${state.up.hifi.solver || 'solver unknown'}).`) : null))));
    for (const t of r.tables) resHost.append(card(null, dataTable(t), { collapsible: false }));
    resHost.append(card('Files read', dataTable({ title: '', columns: ['File', 'Recognised as', 'Read', 'Detail'], rows: r.sources.map((s) => [s.name, s.kind || '—', s.ok ? 'yes' : 'no', s.note]) }), { collapsible: true, open: false }));
  }

  // ------------------------------------------------------------------ page
  const section = (n, title, body) => { const c = card(`${n}. ${title}`, body); c.id = `bridge-step-${n}`; return c; };
  host.append(
    h('div', { class: 'page-h' }, h('div', { class: 'grow' }, h('div', { class: 'eyebrow' }, 'Beyond what a browser can solve'), h('h1', null, 'High-fidelity bridge'),
      h('p', null, 'Some analyses need more than a browser tab: body-fitted 3-D RANS and large-eddy simulation, detailed shell and solid finite elements. This page prepares those cases from the same aircraft data, tells you exactly how to run them on proven open-source solvers, and reads the answers back so you can check and calibrate the built-in suites.'))),
    stepsHost, h('div', { class: 'gap' }),
    section(1, 'Choose the solver and the analysis', chooseHost), h('div', { class: 'gap' }),
    section(2, 'Review the settings', setHost), h('div', { class: 'gap' }),
    section(3, 'Build and download the case', buildHost), h('div', { class: 'gap' }),
    section(4, 'Run it', runHost), h('div', { class: 'gap' }),
    section(5, 'Bring the results back', h('div', { class: 'stack' }, drop, resHost)), h('div', { class: 'gap' }),
    h('div', { class: 'note' }, icon('shield'), h('div', null, 'The decks are starting points written from preliminary-design data. A converged run on a coarse mesh is still a coarse answer: refine the mesh until the result stops changing, and compare with test data where you have it, before using any number for a decision.')));
  paintSteps(); paintChoose(); paintSettings(); paintBuild(); paintRun(); paintResults();
  const off = on('case', () => { mem.settings = {}; paintChoose(); paintSettings(); invalidate(); }), offGeo = on('geometry', paintChoose);
  return () => { off(); offGeo(); clearTimeout(pending); abort?.abort(); killPlots(); };
}
