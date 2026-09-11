(() => {
  const EV_PER_HARTREE = 27.211386245988;
  const state = { models: [], xtbResults: [], selectedResults: [], selectedResultId: "", stage: null, component: null, loadToken: 0, sortDirection: "desc" };
  const $ = (selector) => document.querySelector(selector);

  function escapeHtml(value) {
    return String(value ?? "").replace(/[&<>'"]/g, (character) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "'": "&#39;", '"': "&quot;" })[character]);
  }

  async function request(url) {
    const response = await fetch(url);
    const body = await response.json().catch(() => null);
    if (!response.ok) throw new Error(body?.detail || `${response.status} ${response.statusText}`);
    return body;
  }

  function finite(value) {
    const number = Number(value);
    return Number.isFinite(number) ? number : null;
  }

  function median(values) {
    const sorted = values.filter(Number.isFinite).sort((a, b) => a - b);
    if (!sorted.length) return null;
    const middle = Math.floor(sorted.length / 2);
    return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
  }

  function range(values, digits = 3, suffix = "") {
    const numbers = values.filter(Number.isFinite);
    if (!numbers.length) return "—";
    const low = Math.min(...numbers); const high = Math.max(...numbers);
    return Math.abs(high - low) < 10 ** (-(digits + 1)) ? `${low.toFixed(digits)}${suffix}` : `${low.toFixed(digits)} – ${high.toFixed(digits)}${suffix}`;
  }

  function resolvedSource(model) {
    const source = model.source || {};
    if (source.family !== "complex") return source;
    const surface = state.models.find((item) => item.id === source.surface_model_id);
    const adsorbate = state.models.find((item) => item.id === source.adsorbate_model_id);
    const surfaceSource = surface?.source || {};
    const adsorbateSource = adsorbate?.source || {};
    const pose = source.pose || {};
    return {
      ...source,
      isComplex: true,
      polymer: surfaceSource.polymer || null,
      repeats: surfaceSource.repeats ?? null,
      pH: adsorbateSource.pH ?? surfaceSource.pH ?? null,
      microstate: adsorbateSource.microstate || null,
      adsorbateCompound: adsorbateSource.compound || null,
      poseId: source.pose_id || null,
      siteId: pose.site_id || null,
      surfaceRegion: pose.surface_region || null,
    };
  }

  function selectedModels() {
    const target = $("#analysisPolymer").value;
    const engine = $("#analysisEngine").value;
    const repeats = $("#analysisRepeats").value;
    const query = ($("#analysisSearch")?.value || "").trim().toLowerCase();
    return state.models.filter((model) => {
      const source = resolvedSource(model);
      const targetMatches = target === "all"
        || (target === "dopamine" && (source.compound === "dopamine" || source.adsorbateCompound === "dopamine"))
        || source.polymer === target;
      const searchable = [model.id, model.title, model.engine, model.formula, source.polymer, source.compound, source.adsorbateCompound, source.microstate, source.pH, source.repeats, source.poseId, source.siteId, source.surfaceRegion].join(" ").toLowerCase();
      return targetMatches
        && (engine === "all" || model.engine === engine)
        && (repeats === "all" || Number(source.repeats) === Number(repeats))
        && (!query || searchable.includes(query));
    });
  }

  function enrichResults(models) {
    const byId = new Map(models.map((model) => [model.id, model]));
    const statusFilter = $("#analysisStatus").value;
    const enriched = state.xtbResults.map((result) => {
      const model = byId.get(result.model_id);
      if (!model) return null;
      const source = resolvedSource(model);
      const energy = finite(result.energy_hartree);
      const atoms = finite(model.atom_count);
      const repeats = finite(source.repeats);
      const perAtomEh = energy !== null && atoms ? energy / atoms : null;
      return {
        ...result,
        model,
        isComplex: Boolean(source.isComplex),
        poseSetId: source.pose_set_id || null,
        poseId: source.poseId || null,
        siteId: source.siteId || null,
        surfaceRegion: source.surfaceRegion || null,
        compound: source.compound || null,
        polymer: source.polymer || "—",
        pH: finite(source.pH),
        microstate: source.microstate || null,
        repeats,
        atoms,
        energy,
        perAtomEh,
        perAtomEv: perAtomEh === null ? null : perAtomEh * EV_PER_HARTREE,
        gap: finite(result.homo_lumo_gap_ev),
        gradient: finite(result.gradient_norm_hartree),
        runtime: finite(result.wall_time_seconds),
        groupKey: source.isComplex
          ? `complex|${source.pose_set_id || "—"}|${result.optimization_mode || "full"}|${result.solvent || "none"}`
          : `${source.compound || source.polymer || "—"}|${source.polymer || "—"}|${repeats || "—"}|${source.microstate || "—"}|${source.pH ?? "—"}|${result.simulation_type || (result.optimize ? "optimization" : "singlepoint")}|${result.solvent || "none"}`,
      };
    }).filter(Boolean);
    const filtered = enriched.filter((result) => statusFilter === "all" || (statusFilter === "completed" ? result.status === "completed" : result.status !== "completed"));
    const groups = new Map();
    filtered.forEach((result) => { if (!groups.has(result.groupKey)) groups.set(result.groupKey, []); groups.get(result.groupKey).push(result); });
    filtered.forEach((result) => {
      const values = groups.get(result.groupKey).map((item) => item.perAtomEv).filter(Number.isFinite);
      const centre = median(values);
      result.deltaGroupEv = result.perAtomEv === null || centre === null ? null : result.perAtomEv - centre;
    });
    return { results: filtered, groups };
  }

  function sortResults(results) {
    const field = $("#analysisSort")?.value || "run";
    const direction = state.sortDirection === "asc" ? 1 : -1;
    const value = (result) => {
      if (field === "energy") return result.energy;
      if (field === "perAtom") return result.perAtomEv;
      if (field === "gap") return result.gap;
      if (field === "runtime") return result.runtime;
      if (field === "atoms") return result.atoms;
      // Results arrive newest-first; negate the index so the default
      // descending direction keeps the newest calculation at the top.
      return -state.xtbResults.findIndex((item) => item.model_id === result.model_id);
    };
    return results.slice().sort((left, right) => {
      const a = value(left); const b = value(right);
      if (a === null && b === null) return 0;
      if (a === null) return 1;
      if (b === null) return -1;
      return (a - b) * direction;
    });
  }

  function groupLabel(result) {
    if (result.isComplex) return `${result.polymer}–${result.microstate || "adsorbate"} · ${result.repeats || "?"}-mer · ${result.optimization_mode === "surface_fixed" ? "surface-fixed" : "full"} · ${result.solvent || "none"}`;
    if (result.compound === "dopamine") {
      const pH = result.pH === null ? "pH ?" : `pH ${result.pH}`;
      return `Dopamine · ${pH} · ${result.microstate || "microstate ?"} · ${result.simulation_type || "xTB"} · ${result.solvent || "none"}`;
    }
    return `${result.polymer} · ${result.repeats || "?"}-mer · ${result.simulation_type || "xTB"} · ${result.solvent || "none"}`;
  }

  function setupStage() {
    if (!window.NGL) return;
    try {
      state.stage = new NGL.Stage("analysisViewport", { backgroundColor: "#e9edf3", quality: "high", tooltip: true });
      state.stage.setParameters({ cameraType: "perspective", clipDist: 0, fogNear: 70, fogFar: 100 });
      window.addEventListener("resize", () => state.stage?.handleResize());
    } catch (_) { state.stage = null; }
  }

  function inspectorLabel(result) {
    if (result.isComplex) return `${result.polymer}–${result.microstate || "adsorbate"} · ${result.poseId || "pose"} · ${result.siteId || "site ?"} · ${result.model_id}`;
    if (result.compound === "dopamine") return `Dopamine · pH ${result.pH ?? "?"} · ${result.microstate || "state"} · ${result.model_id}`;
    return `${result.polymer} · ${result.repeats || "?"}-mer · ${result.model.engine} · ${result.model_id}`;
  }

  function inspectorDetails(result, coordinateFile) {
    const rows = [
      ["Model", result.model.title || result.model_id],
      ["Coordinates", coordinateFile === "xtbopt.xyz" ? "Optimized geometry" : "Input geometry"],
      ["Status", result.converged ? "Completed · converged" : (result.status || "review")],
      ["Total energy", result.energy === null ? "—" : `${result.energy.toFixed(8)} Eh`],
      ["Energy / atom", result.perAtomEv === null ? "—" : `${result.perAtomEv.toFixed(5)} eV`],
      ["HOMO–LUMO gap", result.gap === null ? "—" : `${result.gap.toFixed(5)} eV`],
      ["Gradient norm", result.gradient === null ? "—" : `${result.gradient.toExponential(3)} Eh`],
      ["Charge / UHF", `${result.charge ?? 0} / ${result.uhf ?? 0}`],
      ["Solvent", result.solvent || "none"],
      ["Runtime", result.runtime === null ? "—" : `${result.runtime.toFixed(3)} s`],
    ];
    if (result.isComplex) rows.splice(2, 0,
      ["Pose / site", `${result.poseId || "?"} · ${result.siteId || "?"} · ${result.surfaceRegion || "surface"}`],
      ["Optimization mode", result.optimization_mode === "surface_fixed" ? "Surface-fixed screen" : "Full complex refinement"],
    );
    $("#analysisResultDetails").innerHTML = rows.map(([key, value]) => `<div><dt>${escapeHtml(key)}</dt><dd title="${escapeHtml(value)}">${escapeHtml(value)}</dd></div>`).join("");
    const meaning = result.status !== "completed"
      ? "This run needs review before interpretation. Check the raw log for the failure or missing-binary reason."
      : coordinateFile === "xtbopt.xyz"
        ? "This is the lowest-energy geometry found from this starting structure. It is a local minimum, not proof of the global minimum or a binding affinity."
        : "This is the unrelaxed starting geometry. Compare it with Optimized geometry to see how xTB moved the atoms during minimization.";
    $("#analysisResultMeaning").textContent = meaning;
    const base = `/api/models/${encodeURIComponent(result.model_id)}/xtb`;
    const links = [`<a href="${base}/output?file=xtb.out" target="_blank" rel="noreferrer">Open xtb.out</a>`, `<a href="${base}/result" target="_blank" rel="noreferrer">Open JSON</a>`];
    if ((result.files || []).includes("xtbopt.xyz")) links.push(`<a href="${base}/output?file=xtbopt.xyz" target="_blank" rel="noreferrer">Download optimized XYZ</a>`);
    links.push(`<a href="${base}/output?file=input.xyz" target="_blank" rel="noreferrer">Download input XYZ</a>`);
    if ((result.files || []).includes("charges")) links.push(`<a href="${base}/output?file=charges" target="_blank" rel="noreferrer">Download charges</a>`);
    $("#analysisResultLinks").innerHTML = links.join("");
  }

  async function loadInspectorStructure(result) {
    const empty = $("#analysisViewerEmpty"); const loading = $("#analysisViewerLoading");
    if (!result || !state.stage) {
      empty.classList.remove("hidden");
      return;
    }
    let file = $("#analysisStructureType").value;
    if (file === "xtbopt.xyz" && !(result.files || []).includes("xtbopt.xyz")) file = "input.xyz";
    $("#analysisStructureType").value = file;
    const token = ++state.loadToken;
    empty.classList.add("hidden"); loading.classList.remove("hidden");
    state.stage.removeAllComponents(); state.component = null;
    try {
      // The bundled NGL build has no XYZ parser. Keep the raw XYZ downloads,
      // but use PDB coordinates for the live inspector (optimized coordinates
      // are converted server-side while preserving the model CONECT records).
      const optimized = file === "xtbopt.xyz";
      const url = optimized
        ? `/api/models/${encodeURIComponent(result.model_id)}/xtb/output?file=xtbopt.pdb&t=${Date.now()}`
        : `/api/models/${encodeURIComponent(result.model_id)}/structure?format=pdb&t=${Date.now()}`;
      state.component = await state.stage.loadFile(url, { ext: "pdb", defaultRepresentation: false });
      if (token !== state.loadToken) return;
      state.component.addRepresentation("ball+stick", { colorScheme: "element", multipleBond: "symmetric", aspectRatio: 1.7, quality: "high" });
      state.component.autoView(450);
      $("#analysisViewerLabel").textContent = `${file === "xtbopt.xyz" ? "Optimized" : "Input"} · ${result.model.title || result.model_id}`;
    } catch (_) {
      if (token === state.loadToken) { empty.classList.remove("hidden"); $("#analysisViewerLabel").textContent = "Coordinates unavailable"; }
    } finally { if (token === state.loadToken) loading.classList.add("hidden"); }
  }

  function renderInspector(results) {
    const select = $("#analysisResultSelect");
    const completed = results.filter((result) => result.status === "completed");
    if (!completed.length) {
      select.innerHTML = '<option value="">No completed result</option>'; $("#analysisResultDetails").innerHTML = "<div><dt>Result</dt><dd>No completed result</dd></div>"; $("#analysisResultMeaning").textContent = "Run xTB or choose a filter with completed results."; $("#analysisResultLinks").innerHTML = ""; $("#analysisViewerEmpty").classList.remove("hidden"); $("#analysisViewerLabel").textContent = "No structure selected"; return;
    }
    select.innerHTML = completed.map((result) => `<option value="${escapeHtml(result.model_id)}">${escapeHtml(inspectorLabel(result))}</option>`).join("");
    const selected = completed.find((result) => result.model_id === state.selectedResultId) || completed[0];
    state.selectedResultId = selected.model_id; state.selectedResults = results; select.value = selected.model_id;
    const file = $("#analysisStructureType").value;
    inspectorDetails(selected, file === "xtbopt.xyz" && (selected.files || []).includes("xtbopt.xyz") ? file : "input.xyz");
    loadInspectorStructure(selected);
  }

  function renderSummary(models, results, groups) {
    const targets = new Set(models.map((model) => { const source = resolvedSource(model); return source.isComplex ? `${source.polymer}-${source.microstate}` : source.polymer || source.compound; }).filter(Boolean));
    const engines = new Set(models.map((model) => model.engine));
    const completed = results.filter((result) => result.status === "completed");
    const review = results.length - completed.length;
    const runtimes = completed.map((result) => result.runtime).filter(Number.isFinite);
    const normalized = completed.map((result) => result.perAtomEv).filter(Number.isFinite);
    $("#analysisSummary").textContent = `${models.length} model${models.length === 1 ? "" : "s"} · ${targets.size} target${targets.size === 1 ? "" : "s"} · ${engines.size} engine${engines.size === 1 ? "" : "s"} · ${results.length} xTB result${results.length === 1 ? "" : "s"}`;
    $("#analysisResultCount").textContent = `${results.length} / ${state.xtbResults.length}`;
    $("#analysisTableMeta").textContent = `${results.length} visible · ${completed.length} completed · ${groups.size} comparison group${groups.size === 1 ? "" : "s"}`;
    const metricValues = [
      ["Jobs", results.length],
      ["Completed", completed.length],
      ["Needs review", review],
      ["Median runtime", runtimes.length ? `${median(runtimes).toFixed(2)} s` : "—"],
      ["Size groups", groups.size],
      ["E/atom spread", normalized.length > 1 ? `${(Math.max(...normalized) - Math.min(...normalized)).toFixed(3)} eV` : "—"],
    ];
    $("#xtbMetrics").innerHTML = metricValues.map(([label, value]) => `<article class="xtb-metric"><span>${label}</span><strong>${value}</strong></article>`).join("");
    const groupSpreads = [...groups.values()].map((group) => {
      const values = group.filter((result) => result.status === "completed").map((result) => result.perAtomEv).filter(Number.isFinite);
      return { label: groupLabel(group[0]), spread: values.length > 1 ? Math.max(...values) - Math.min(...values) : 0 };
    }).sort((a, b) => b.spread - a.spread);
    const broadest = groupSpreads[0];
    $("#analysisQuickRead").textContent = !completed.length
      ? "Quick read: no completed calculation is available for this filter."
      : `Quick read: ${completed.length}/${results.length} runs converged. ${broadest && broadest.spread > 0 ? `The widest within-group E/atom spread is ${broadest.spread.toFixed(3)} eV (${broadest.label}). ` : "Matched models are energetically tight. "}Total energy is a size-dependent diagnostic, not an affinity score; use matched complex–fragment calculations for adsorption. `;
  }

  function renderEngineChart(models) {
    const chart = $("#engineChart");
    if (!models.length) { chart.innerHTML = '<div class="chart-empty">No models match the current filter.</div>'; return; }
    const grouped = {};
    models.forEach((model) => { grouped[model.engine] = (grouped[model.engine] || []).concat(Number(model.wall_time_seconds || 0)); });
    const values = Object.entries(grouped).map(([label, times]) => ({ label, value: times.reduce((sum, value) => sum + value, 0) / times.length }));
    const max = Math.max(...values.map((item) => item.value), 0.01);
    chart.innerHTML = values.map((item) => `<div class="bar-row"><span title="${escapeHtml(item.label)}">${escapeHtml(item.label.replace(" gen_params + gen_coords", "").replace(" MoleculeBuilder", "").replace("RDKit ETKDGv3", "RDKit"))}</span><div class="bar-track"><i style="width:${Math.max(3, (item.value / max) * 100)}%"></i></div><b>${item.value.toFixed(2)} s</b></div>`).join("");
  }

  function renderRepeatChart(models) {
    const chart = $("#repeatChart");
    const polymerModels = models.filter((model) => model.source?.family === "polymer" && model.source?.polymer && Number.isFinite(Number(model.source.repeats)));
    if (!polymerModels.length) { chart.innerHTML = '<div class="chart-empty">Polymer models are required for this plot.</div>'; return; }
    const max = Math.max(...polymerModels.map((model) => Number(model.source.repeats)), 1);
    chart.innerHTML = polymerModels.map((model) => { const repeats = Number(model.source.repeats); const left = 8 + (repeats / max) * 84; return `<a class="dot-row" href="/?model=${encodeURIComponent(model.id)}" title="Open ${escapeHtml(model.title)}"><span>${escapeHtml(model.source.polymer)}</span><div class="dot-track"><i style="left:${left}%"></i></div><b>${repeats}-mer</b></a>`; }).join("");
  }

  function renderModelTable(models) {
    const body = $("#modelTable");
    if (!models.length) { body.innerHTML = '<tr><td colspan="8" class="table-empty">No models match the current filter.</td></tr>'; return; }
    body.innerHTML = models.map((model) => {
      const source = resolvedSource(model); const validation = model.validation || {}; const valid = validation.rdkit_sanitized;
      return `<tr><td><a class="text-link" href="/?model=${encodeURIComponent(model.id)}">${escapeHtml(model.title)}</a><small>${escapeHtml(model.id)}</small></td><td>${escapeHtml(source.polymer || "—")}</td><td>${source.repeats || "—"}</td><td>${escapeHtml(model.engine)}</td><td>${model.atom_count}</td><td class="mono">${escapeHtml(model.formula)}</td><td>${Number(model.wall_time_seconds || 0).toFixed(2)} s</td><td><span class="validation-pill ${valid ? "ok" : "pending"}">${valid ? "Basic 3D" : "Review"}</span></td></tr>`;
    }).join("");
  }

  function renderBarChart(selector, values, formatter, barClass = "") {
    const chart = $(selector);
    if (!values.length) { chart.innerHTML = '<div class="chart-empty">No completed result is available for this filter.</div>'; return; }
    const numbers = values.map((item) => item.value).filter(Number.isFinite);
    if (!numbers.length) { chart.innerHTML = '<div class="chart-empty">This descriptor is not present in the saved xTB output.</div>'; return; }
    const min = Math.min(...numbers); const max = Math.max(...numbers); const span = Math.max(max - min, 1e-12);
    chart.innerHTML = values.filter((item) => Number.isFinite(item.value)).map((item) => {
      const width = Math.max(4, ((item.value - min) / span) * 96 + 4);
      return `<div class="bar-row"><span title="${escapeHtml(item.label)}">${escapeHtml(item.label)}</span><div class="bar-track"><i class="${barClass}" style="width:${width}%"></i></div><b>${formatter(item.value)}</b></div>`;
    }).join("");
  }

  function renderXtbCharts(results) {
    const completed = results.filter((result) => result.status === "completed");
    const label = (result) => result.isComplex
      ? `${result.polymer} ${result.repeats || "?"}-mer · ${result.poseId || "pose"}`
      : result.compound === "dopamine"
      ? `DA pH${result.pH ?? "?"} · ${(result.microstate || "state").replace("DA_", "")}`
      : `${result.polymer} ${result.repeats || "?"}-mer · ${result.model.engine.split(" ")[0]}`;
    renderBarChart("#xtbEnergyChart", completed.map((result) => ({ label: label(result), value: result.energy })), (value) => value.toFixed(4));
    renderBarChart("#xtbNormalizedChart", completed.map((result) => ({ label: label(result), value: result.perAtomEv })), (value) => `${value.toFixed(3)} eV`);
    renderBarChart("#xtbRuntimeChart", results.map((result) => ({ label: label(result), value: result.runtime || 0 })), (value) => `${value.toFixed(2)} s`, "bar-ok");
    renderBarChart("#xtbGapChart", completed.map((result) => ({ label: label(result), value: result.gap })), (value) => `${value.toFixed(2)} eV`);
  }

  function renderComparison(groups) {
    const body = $("#comparisonTable");
    if (!groups.size) { body.innerHTML = '<tr><td colspan="8" class="table-empty">No xTB results for the current filter.</td></tr>'; return; }
    body.innerHTML = [...groups.values()].sort((a, b) => String(groupLabel(a[0])).localeCompare(String(groupLabel(b[0])))).map((group) => {
      const completed = group.filter((result) => result.status === "completed");
      const normalized = completed.map((result) => result.perAtomEv).filter(Number.isFinite);
      const best = normalized.length ? completed.filter((result) => Number.isFinite(result.perAtomEv)).sort((a, b) => a.perAtomEv - b.perAtomEv)[0] : null;
      const runtimes = completed.map((result) => result.runtime).filter(Number.isFinite);
      const atomCounts = group.map((result) => result.atoms).filter(Number.isFinite);
      const total = completed.map((result) => result.energy).filter(Number.isFinite);
      const spread = normalized.length > 1 ? Math.max(...normalized) - Math.min(...normalized) : 0;
      const assessment = group.some((result) => result.status !== "completed") ? "Incomplete" : spread < 0.01 ? "Tight" : spread < 0.05 ? "Moderate" : "Review geometry";
      return `<tr><td><b>${escapeHtml(groupLabel(group[0]))}</b></td><td>${group.length}</td><td>${range(atomCounts, 0)}</td><td class="mono">${range(total, 4)}</td><td class="mono">${range(normalized, 3, " eV")}</td><td>${best ? escapeHtml(best.model.engine) : "—"}</td><td>${runtimes.length ? `${median(runtimes).toFixed(2)} s` : "—"}</td><td class="assessment"><span class="validation-pill ${assessment === "Tight" ? "ok" : "pending"}">${assessment}</span></td></tr>`;
    }).join("");
  }

  function renderXtbTable(results) {
    const body = $("#xtbTable");
    if (!results.length) { body.innerHTML = '<tr><td colspan="9" class="table-empty">No xTB results for the current filter. Run a job from Simulation.</td></tr>'; return; }
    body.innerHTML = sortResults(results).map((result) => {
      const completed = result.status === "completed";
      const sizeState = result.isComplex
        ? `${result.poseId || "pose"} · ${result.siteId || "site ?"} · ${result.surfaceRegion || "surface"}`
        : result.compound === "dopamine"
        ? `pH ${result.pH ?? "?"} · ${result.microstate || "state"}`
        : (result.repeats ? `${result.repeats}-mer` : "—");
      return `<tr><td><a class="text-link" href="/?model=${encodeURIComponent(result.model_id || "")}">${escapeHtml(result.model.title || result.model_id)}</a><small>${escapeHtml(result.model_id || "")}</small></td><td><b>${escapeHtml(sizeState)}</b><small>${escapeHtml(result.model.engine || "—")}</small></td><td>${result.atoms || "—"}</td><td class="mono">${result.energy === null ? "—" : result.energy.toFixed(5)}</td><td class="mono">${result.perAtomEv === null ? "—" : result.perAtomEv.toFixed(4)}</td><td>${result.gap === null ? "—" : result.gap.toFixed(3)}</td><td>${result.runtime === null ? "—" : `${result.runtime.toFixed(2)} s`}</td><td><span class="validation-pill ${completed ? "ok" : "pending"}">${escapeHtml(result.status || "unknown")}</span></td><td><button class="table-action inspect-result" type="button" data-result-id="${escapeHtml(result.model_id || "")}">Inspect</button></td></tr>`;
    }).join("");
    body.querySelectorAll(".inspect-result").forEach((button) => button.addEventListener("click", () => {
      state.selectedResultId = button.dataset.resultId;
      renderInspector(state.selectedResults);
      $(".analysis-inspector-card")?.scrollIntoView({ behavior: "smooth", block: "nearest" });
    }));
  }

  function render() {
    const models = selectedModels();
    const { results, groups } = enrichResults(models);
    renderSummary(models, results, groups);
    renderEngineChart(models); renderRepeatChart(models); renderModelTable(models);
    renderXtbCharts(results); renderComparison(groups); renderXtbTable(results); renderInspector(results);
  }

  async function load() {
    try {
      const [health, models, xtbResults] = await Promise.all([request("/api/health"), request("/api/models?limit=500"), request("/api/xtb/results?limit=500")]);
      state.models = models; state.xtbResults = xtbResults;
      $("#apiStatus").classList.add("ready"); $("#apiStatus span").textContent = `Engine ready · v${health.version}`;
      render();
    } catch (error) {
      $("#apiStatus").classList.add("error"); $("#apiStatus span").textContent = "Service offline"; $("#analysisSummary").textContent = error.message;
    }
  }

  setupStage();
  $("#analysisPolymer").addEventListener("change", render);
  $("#analysisEngine").addEventListener("change", render);
  $("#analysisRepeats").addEventListener("change", render);
  $("#analysisStatus").addEventListener("change", render);
  $("#analysisSearch").addEventListener("input", render);
  $("#analysisSort").addEventListener("change", render);
  $("#analysisSortDirection").addEventListener("click", () => { state.sortDirection = state.sortDirection === "asc" ? "desc" : "asc"; $("#analysisSortDirection").textContent = state.sortDirection === "asc" ? "↑" : "↓"; render(); });
  $("#clearAnalysisFilters").addEventListener("click", () => {
    $("#analysisSearch").value = ""; $("#analysisPolymer").value = "all"; $("#analysisEngine").value = "all"; $("#analysisRepeats").value = "all"; $("#analysisStatus").value = "all"; $("#analysisSort").value = "run"; state.sortDirection = "desc"; $("#analysisSortDirection").textContent = "↓"; render();
  });
  document.querySelectorAll(".analysis-tab").forEach((tab) => tab.addEventListener("click", () => {
    const panel = tab.dataset.panel;
    document.querySelectorAll(".analysis-tab").forEach((item) => { const active = item === tab; item.classList.toggle("active", active); item.setAttribute("aria-selected", active ? "true" : "false"); });
    document.querySelectorAll(".analysis-panel").forEach((item) => { const active = item.dataset.panelContent === panel; item.classList.toggle("active", active); item.hidden = !active; });
  }));
  $("#analysisResultSelect").addEventListener("change", (event) => { state.selectedResultId = event.target.value; renderInspector(state.selectedResults); });
  $("#analysisStructureType").addEventListener("change", () => renderInspector(state.selectedResults));
  $("#analysisResetView").addEventListener("click", () => state.component?.autoView(400));
  $("#refreshAnalysis").addEventListener("click", load);
  load();
})();
