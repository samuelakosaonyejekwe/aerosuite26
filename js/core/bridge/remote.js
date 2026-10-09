// High-fidelity bridge — optional remote execution on GitHub Actions runners in the user's OWN repository.
// The personal access token is held by the caller and passed in; this module never stores it and attaches
// it only to requests whose URL starts with https://api.github.com/. Artifact downloads are redirected by
// GitHub to a short-lived signed storage URL; the platform fetch drops the Authorization header on that
// cross-origin redirect, so the token never leaves api.github.com.
// No DOM. `fetchImpl` can be injected for tests.

import { base64Encode } from './zip.js';

export const API = 'https://api.github.com';
export const WORKFLOW_FILE = 'solve.yml';
export const CASE_BRANCH = 'aerosuite-cases';
/** workflow_dispatch inputs are limited to 65 535 characters in total; stay well inside it. */
export const INLINE_LIMIT = 48000;
export const TOKEN_PERMISSIONS = ['Actions: Read and write (start the workflow, read its status, download results)', 'Contents: Read and write (only needed when the case is too large to send inline and is uploaded to a branch)', 'Metadata: Read-only (always included)'];

/** "owner/repo", a github.com URL or a git remote → { owner, repo }. */
export function parseRepo(input) {
  const s = String(input || '').trim().replace(/\.git$/, '').replace(/\/+$/, '');
  const m = /^(?:https?:\/\/(?:www\.)?github\.com\/|git@github\.com:)?([A-Za-z0-9](?:[A-Za-z0-9-]{0,38}))\/([A-Za-z0-9._-]{1,100})$/.exec(s);
  if (!m) throw new Error('Enter the repository as owner/name, for example jane/aerosuite.');
  return { owner: m[1], repo: m[2] };
}
export const safeCaseName = (name) => String(name || 'case').replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^[-.]+|[-.]+$/g, '').slice(0, 80) || 'case';

export class GithubError extends Error { constructor(message, status, detail) { super(message); this.name = 'GithubError'; this.status = status; this.detail = detail; } }

function explain(status, detail, what) {
  const d = detail ? ` GitHub said: “${String(detail).slice(0, 200)}”.` : '';
  if (status === 401) return `GitHub rejected the token while trying to ${what}. It may be mistyped, expired or revoked.${d}`;
  if (status === 403) return `The token is not allowed to ${what}. Give it these repository permissions: Actions (read and write)${/upload|branch|file/.test(what) ? ' and Contents (read and write)' : ''}. A rate limit also gives this answer; wait a minute and retry.${d}`;
  if (status === 404) return `GitHub could not find what was needed to ${what}. Check the repository name, that the token was created for this repository, and that .github/workflows/${WORKFLOW_FILE} is on its default branch.${d}`;
  if (status === 422) return `GitHub refused the request to ${what}.${d}`;
  return `GitHub returned an error (${status}) while trying to ${what}.${d}`;
}

/**
 * Minimal GitHub REST client for one repository.
 * @param {{token: string, owner: string, repo: string, fetchImpl?: typeof fetch}} o
 */
export function createGithub({ token, owner, repo, fetchImpl = globalThis.fetch }) {
  if (!token || !/^[A-Za-z0-9_]{20,255}$/.test(token)) throw new Error('Paste a GitHub personal access token (it starts with github_pat_ or ghp_).');
  const base = `${API}/repos/${encodeURIComponent(owner)}/${encodeURIComponent(repo)}`;
  async function call(path, { method = 'GET', body, what = 'reach GitHub', binary = false, signal, okStatus = [] } = {}) {
    const url = path.startsWith('http') ? path : base + path;
    if (!url.startsWith(`${API}/`)) throw new Error('Refusing to send the token to a host other than api.github.com.');
    let res;
    try {
      res = await fetchImpl(url, { method, signal, cache: 'no-store', redirect: 'follow', referrerPolicy: 'no-referrer', headers: { Authorization: `Bearer ${token}`, Accept: 'application/vnd.github+json', 'X-GitHub-Api-Version': '2022-11-28', ...(body ? { 'Content-Type': 'application/json' } : {}) }, body: body ? JSON.stringify(body) : undefined });
    } catch (e) {
      if (e?.name === 'AbortError') throw e;
      throw new GithubError(`The request to ${what} did not go through (${e?.message || 'network error'}). Check the connection; if this is the results download, the app's content-security policy may not yet allow the storage host.`, 0, String(e?.message || e));
    }
    if (okStatus.includes(res.status)) return { status: res.status, data: null };
    if (!res.ok) { let detail = ''; try { detail = (await res.json())?.message || ''; } catch { /* no body */ } throw new GithubError(explain(res.status, detail, what), res.status, detail); }
    if (binary) return { status: res.status, data: new Uint8Array(await res.arrayBuffer()) };
    if (res.status === 204) return { status: 204, data: null };
    let data = null; try { data = await res.json(); } catch { /* empty body */ }
    return { status: res.status, data };
  }
  const enc = (p) => p.split('/').map(encodeURIComponent).join('/');
  const api = {
    owner, repo,
    async repoInfo(signal) { const { data } = await call('', { what: 'read the repository', signal }); return { defaultBranch: data.default_branch, private: !!data.private, fullName: data.full_name, htmlUrl: data.html_url }; },
    async workflow(signal) { const { data } = await call(`/actions/workflows/${WORKFLOW_FILE}`, { what: 'find the solve workflow', signal }); return { id: data.id, state: data.state, htmlUrl: data.html_url }; },
    /** Make sure a branch exists (created from the default branch when missing). */
    async ensureBranch(branch, defaultBranch, signal) {
      const r = await call(`/git/ref/heads/${enc(branch)}`, { what: 'look up the case branch', okStatus: [404], signal });
      if (r.status !== 404) return { created: false };
      const { data } = await call(`/git/ref/heads/${enc(defaultBranch)}`, { what: 'read the default branch', signal });
      await call('/git/refs', { method: 'POST', body: { ref: `refs/heads/${branch}`, sha: data.object.sha }, what: 'create the case branch', signal });
      return { created: true };
    },
    /** Create or replace one file on a branch through the contents API. */
    async putFile(path, bytes, branch, message, signal) {
      const cur = await call(`/contents/${enc(path)}?ref=${encodeURIComponent(branch)}`, { what: 'check for an earlier upload of the case file', okStatus: [404], signal });
      const { data } = await call(`/contents/${enc(path)}`, { method: 'PUT', what: 'upload the case file to the case branch', signal, body: { message, content: base64Encode(bytes), branch, ...(cur.data?.sha ? { sha: cur.data.sha } : {}) } });
      return { path: data.content?.path || path, sha: data.content?.sha };
    },
    async dispatch(ref, inputs, signal) { await call(`/actions/workflows/${WORKFLOW_FILE}/dispatches`, { method: 'POST', body: { ref, inputs }, what: 'start the solve workflow', signal }); },
    async recentRuns(signal) { const { data } = await call(`/actions/workflows/${WORKFLOW_FILE}/runs?event=workflow_dispatch&per_page=20`, { what: 'list workflow runs', signal }); return data.workflow_runs || []; },
    async run(id, signal) { const { data } = await call(`/actions/runs/${id}`, { what: 'read the run status', signal }); return data; },
    async jobs(id, signal) { const { data } = await call(`/actions/runs/${id}/jobs?per_page=30`, { what: 'read the run steps', signal }); return data.jobs || []; },
    async artifacts(id, signal) { const { data } = await call(`/actions/runs/${id}/artifacts?per_page=50`, { what: 'list the result files', signal }); return data.artifacts || []; },
    /** Download an artifact archive (GitHub answers with a redirect to a signed storage URL valid for one minute). */
    async downloadArtifact(id, signal) { const { data } = await call(`/actions/artifacts/${id}/zip`, { what: 'download the results', binary: true, signal }); return data; },
    async cancel(id, signal) { await call(`/actions/runs/${id}/cancel`, { method: 'POST', what: 'cancel the run', okStatus: [409], signal }); },
  };
  return api;
}

const sleep = (ms, signal) => new Promise((res, rej) => { const t = setTimeout(res, ms); signal?.addEventListener('abort', () => { clearTimeout(t); rej(Object.assign(new Error('Stopped.'), { name: 'AbortError' })); }, { once: true }); });
export const newTag = () => `r${Date.now().toString(36)}${Math.random().toString(36).slice(2, 8)}`;

/**
 * Send a case archive to the repository and start the solve workflow.
 * Small archives travel inline as a workflow input; larger ones are committed as cases/<name>/case.zip on a
 * separate branch (never the default branch).
 * @returns {Promise<{tag, route: 'inline'|'branch', branch: string|null, defaultBranch, since: number}>}
 */
export async function submitCase(gh, { zipBytes, name, solver, threads = 4, signal, onStatus = () => {} }) {
  const caseName = safeCaseName(name), tag = newTag(), b64 = base64Encode(zipBytes);
  onStatus('Checking the repository and the workflow…');
  const info = await gh.repoInfo(signal), wf = await gh.workflow(signal);
  if (wf.state && wf.state !== 'active') throw new Error('The solve workflow exists but is disabled. Enable it on the repository\'s Actions tab.');
  const since = Date.now() - 15000, inputs = { solver, case: caseName, tag, threads: String(threads) };
  let route = 'inline', branch = null;
  if (b64.length <= INLINE_LIMIT) inputs.case_b64 = b64;
  else {
    route = 'branch'; branch = CASE_BRANCH;
    if (zipBytes.length > 60e6) throw new Error('The case archive is larger than 60 MB, which is too much to upload from the browser. Commit it to the repository with git and start the workflow from the Actions tab.');
    onStatus(`Uploading the case (${(zipBytes.length / 1e6).toFixed(2)} MB) to branch ${branch}…`);
    await gh.ensureBranch(branch, info.defaultBranch, signal);
    await gh.putFile(`cases/${caseName}/case.zip`, zipBytes, branch, `Add solver case ${caseName}`, signal);
    inputs.branch = branch;
  }
  onStatus('Starting the workflow…');
  await gh.dispatch(info.defaultBranch, inputs, signal);
  return { tag, route, branch, defaultBranch: info.defaultBranch, since, caseName, repoUrl: info.htmlUrl };
}

/** Find the run started by submitCase (GitHub does not return its id), then poll until it finishes. */
export async function followRun(gh, { tag, since }, { signal, onStatus = () => {}, pollMs = 8000, findTimeoutMs = 120000, sleepFn = sleep } = {}) {
  let run = null; const t0 = Date.now();
  while (!run) {
    const runs = await gh.recentRuns(signal);
    run = runs.find((r) => (r.display_title || r.name || '').includes(tag)) || null;
    if (!run && Date.now() - t0 > findTimeoutMs) throw new Error('The workflow was started but its run did not appear within two minutes. Look at the repository\'s Actions tab.');
    if (!run) { onStatus('Waiting for the run to appear in the queue…'); await sleepFn(Math.min(pollMs, 4000), signal); }
  }
  void since;
  for (;;) {
    let step = '';
    try { const jobs = await gh.jobs(run.id, signal), active = jobs.flatMap((j) => j.steps || []).find((s) => s.status === 'in_progress'); if (active) step = ` — ${active.name}`; } catch (e) { if (e?.name === 'AbortError') throw e; }
    onStatus(run.status === 'queued' ? 'Queued: waiting for a free runner…' : run.status === 'completed' ? `Finished: ${run.conclusion}` : `Running${step}…`, run);
    if (run.status === 'completed') return run;
    await sleepFn(pollMs, signal);
    run = await gh.run(run.id, signal);
  }
}

/** Download every artifact of a finished run. Returns [{ name, bytes }] (each is a ZIP archive). */
export async function fetchArtifacts(gh, runId, { signal, onStatus = () => {} } = {}) {
  const list = (await gh.artifacts(runId, signal)).filter((a) => !a.expired), out = [];
  if (!list.length) throw new Error('The run produced no result files. Open its log on GitHub to see why.');
  for (const a of list) { onStatus(`Downloading ${a.name} (${(a.size_in_bytes / 1e6).toFixed(2)} MB)…`); out.push({ name: a.name, bytes: await gh.downloadArtifact(a.id, signal) }); }
  return out;
}
