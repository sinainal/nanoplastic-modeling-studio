(() => {
  const state = { models: [], poseSets: [], poseSet: null, results: {}, checked: new Set(), stage: null, component: null, token: 0, activeJob: null, xtb: null, page: 0, pageSize: 8, selectedPoseIndex: 0, experimentSetIds: new Set(), xtbQueue: [] };
  const $ = (selector) => document.querySelector(selector);
  const escapeHtml = (value) => String(value ?? "").replace(/[&<>'"]/g, (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "'": "&#39;", '"': "&quot;" })[char]);

  async function request(url, options) {
    const response = await fetch(url, options);
    const body = await response.json().catch(() => null);
    if (!response.ok) throw new Error(body?.detail || `${response.status} ${response.statusText}`);
    return body;
  }

  function toast(message, error = false) {
    const node = $("#toast"); node.textContent = message; node.classList.toggle("error", error); node.classList.add("show");
    window.clearTimeout(toast.timer); toast.timer = window.setTimeout(() => node.classList.remove("show"), 3600);
  }

  function setupStage() {
    if (!window.NGL) return;
    try { state.stage = new NGL.Stage("runViewport", { backgroundColor: "#e9edf3", quality: "high", tooltip: true }); state.stage.setParameters({ cameraType: "perspective", clipDist: 0 }); window.addEventListener("resize", () => state.stage?.handleResize()); }
    catch (_) { state.stage = null; }
  }

  function representationParams(type, surface = false) {
    const ghosted = $("#runOverlayStyle")?.value !== "solid";
    const palette = $("#runColorMode")?.value || "pose";
    const color = surface ? { colorScheme: "uniform", colorValue: ghosted ? 0x8190a3 : 0x6f7b8c } : { colorScheme: palette === "pose" ? "chainid" : "element" };
    const transparency = surface && ghosted ? { opacity: 0.28, depthWrite: false } : {};
    if (type === "spacefill") return { ...color, ...transparency, radiusScale: surface ? 0.48 : 0.68, quality: "high" };
    if (type === "line") return { ...color, ...transparency, multipleBond: true, linewidth: surface ? 1 : 2 };
    if (type === "licorice") return { ...color, ...transparency, multipleBond: "symmetric", radiusScale: surface ? 0.34 : 0.58, quality: "high" };
    return { ...color, ...transparency, multipleBond: "symmetric", aspectRatio: surface ? 2.8 : 2.1, radiusScale: surface ? 0.42 : 0.78, quality: "high" };
  }

  function applyRepresentation() {
    if (!state.component) return;
    const type = $("#runRepresentation").value;
    state.component.removeAllRepresentations();
    state.component.addRepresentation(type, { ...representationParams(type, true), sele: ":A" });
    state.component.addRepresentation(type, { ...representationParams(type, false), sele: "not :A" });
  }

  function effectivePageSize() { const total = availableCount() || 1; return state.pageSize === "all" ? total : Math.max(1, Number(state.pageSize) || 8); }
  function pagePoses() { const start = state.page * effectivePageSize(); return (state.poseSet?.poses || []).slice(start, start + effectivePageSize()); }

  function availableCount() { return state.poseSet?.poses?.length || 0; }
  function requestedCount() { return Math.max(1, Math.min(availableCount() || 1, Number.parseInt($("#runExperimentPoseCount")?.value || "8", 10) || 8)); }

  function selectedModels() { return (state.poseSet?.poses || []).filter((pose) => state.checked.has(pose.model_id)); }

  function renderStatus() {
    const ready = state.xtb?.status === "ready"; const badge = $("#runXtbBadge");
    badge.textContent = `xTB · ${ready ? "ready" : "unavailable"}`; badge.classList.toggle("ready", ready); badge.classList.toggle("error", !ready);
    $("#runSelectionBadge").textContent = `${state.checked.size} selected`;
    $("#runXtbButton").disabled = !ready || !state.checked.size || Boolean(state.activeJob);
  }

  function renderPreview() {
    const selected = state.checked.size; const total = availableCount(); const mode = $("#runOptimizationMode").value === "surface_fixed" ? "Surface-fixed screen" : "Full complex refinement";
    $("#runBatchPreview").innerHTML = `<div><span>Batch</span><b>${selected}/${total} selected</b></div><div><span>Mode</span><b>${mode}</b></div><div><span>Output</span><b>Optimized local minima; Eads needs matched references</b></div>`;
  }

  function renderPoseList() {
    const poses = state.poseSet?.poses || []; const holder = $("#runPoseList");
    if (!holder) return;
    const count = $("#runPoseCount"); if (count) count.textContent = `${poses.length} poses`;
    if (!poses.length) { holder.innerHTML = '<div class="table-empty">Open a saved pose set to inspect its members.</div>'; return; }
    holder.innerHTML = poses.map((pose, index) => {
      const result = state.results[pose.model_id]; const parameters = pose.parameters || {}; const status = result?.status || "not run";
      const energy = Number(result?.energy_hartree); const site = parameters.site_id || "site"; const orientation = parameters.label || parameters.orientation || "orientation";
      return `<label class="run-pose-row ${state.checked.has(pose.model_id) ? "selected" : ""}"><input type="checkbox" data-model-id="${escapeHtml(pose.model_id)}" ${state.checked.has(pose.model_id) ? "checked" : ""}><span><b>${escapeHtml(pose.pose_id || `P${index + 1}`)}</b><small>${escapeHtml(site)} · ${escapeHtml(orientation)}</small></span><span class="run-pose-energy">${Number.isFinite(energy) ? energy.toFixed(5) + " Eh" : "—"}<small>${escapeHtml(status)}</small></span></label>`;
    }).join("");
    holder.querySelectorAll("input").forEach((input) => input.addEventListener("change", () => { if (input.checked) state.checked.add(input.dataset.modelId); else state.checked.delete(input.dataset.modelId); render(); loadViewer(); }));
  }

  function renderSet() {
    const set = state.poseSet;
    if (!set) { $("#runSetSummary").innerHTML = "<span>No saved pose set selected.</span>"; $("#runViewerTitle").textContent = "Choose a pose set"; return; }
    const surface = state.models.find((model) => model.id === set.surface_model_id); const ligand = state.models.find((model) => model.id === set.adsorbate_model_id);
    const polymer = surface?.source?.polymer || "Polymer"; const repeats = surface?.source?.repeats || "?"; const adsorbate = ligand?.source?.microstate || "adsorbate";
    $("#runSetSummary").innerHTML = `<span><b>${escapeHtml(polymer)} ${repeats}-mer</b></span><i>+</i><span><b>${escapeHtml(adsorbate)}</b></span><i>·</i><span>${set.site_count || "?"} sites × ${set.orientations_per_site || "?"} orientations · charge <b>${Number(set.formal_charge || 0) >= 0 ? "+" : ""}${set.formal_charge || 0}</b></span>`;
    $("#runViewerTitle").textContent = `${polymer} ${repeats}-mer · ${adsorbate}`;
    $("#runViewerMeta").textContent = `${availableCount()} pose members · ${state.checked.size} selected`;
  }

  function renderViewerPagination() {
    const total = availableCount(); const pageSize = effectivePageSize(); const pages = Math.max(1, Math.ceil(total / pageSize)); state.page = Math.min(state.page, pages - 1);
    const start = total ? state.page * pageSize + 1 : 0; const end = Math.min((state.page + 1) * pageSize, total);
    $("#runPosePageLabel").textContent = total ? `${start}–${end} of ${total} poses` : "No poses";
    $("#runPreviousPosePage").disabled = !total || state.page === 0; $("#runNextPosePage").disabled = !total || state.page >= pages - 1;
  }

  function renderViewerPoseCards() {
    const grid = $("#runPoseCardGrid"); const poses = pagePoses();
    if (!poses.length) { grid.innerHTML = '<div class="pose-grid-empty">No poses loaded.</div>'; return; }
    const start = state.page * effectivePageSize();
    grid.innerHTML = poses.map((pose, localIndex) => {
      const index = start + localIndex; const result = state.results[pose.model_id]; const parameters = pose.parameters || {}; const site = parameters.site_id ? `${parameters.site_id} · ${parameters.surface_region || "surface"}` : pose.label || "surface pose"; const orientation = parameters.label || parameters.orientation || pose.method || "orientation";
      return `<article class="pose-card ${index === state.selectedPoseIndex ? "active" : ""}" data-pose-index="${index}"><label><input class="run-viewer-pose-check" type="checkbox" data-model-id="${escapeHtml(pose.model_id)}" ${state.checked.has(pose.model_id) ? "checked" : ""}><span>${escapeHtml(pose.pose_id || `P${index + 1}`)}</span></label><b>${escapeHtml(site)}</b><small>${escapeHtml(orientation)}</small><i class="pose-status-dot ${result?.status === "completed" ? "complete" : result?.status === "failed" ? "failed" : ""}"></i></article>`;
    }).join("");
    grid.querySelectorAll(".run-viewer-pose-check").forEach((input) => input.addEventListener("change", () => { if (input.checked) state.checked.add(input.dataset.modelId); else state.checked.delete(input.dataset.modelId); render(); loadViewer(); }));
    grid.querySelectorAll(".pose-card").forEach((card) => card.addEventListener("click", (event) => { if (event.target.closest("label")) return; state.selectedPoseIndex = Number(card.dataset.poseIndex); renderViewerPoseCards(); }));
  }

  function renderExperimentDesigner() {
    const list = $("#runExperimentSetList"); if (!list) return;
    if (!state.experimentSetIds.size && state.poseSet?.id) state.experimentSetIds.add(state.poseSet.id);
    list.innerHTML = state.poseSets.length ? state.poseSets.map((set) => `<label><input type="checkbox" value="${escapeHtml(set.id)}" ${state.experimentSetIds.has(set.id) ? "checked" : ""}>${escapeHtml(set.title || set.id)} · ${set.poses?.length || 0}</label>`).join("") : "<em>No saved pose sets.</em>";
    list.querySelectorAll("input").forEach((input) => input.addEventListener("change", () => { if (input.checked) state.experimentSetIds.add(input.value); else state.experimentSetIds.delete(input.value); renderExperimentDesigner(); }));
    const representatives = Math.max(1, Number($("#runExperimentPoseCount").value || 8)); const total = state.experimentSetIds.size * representatives;
    $("#runExperimentSummary").textContent = `${state.experimentSetIds.size} sets · ${total} jobs`;
    $("#runExperimentDetail").textContent = state.experimentSetIds.size ? `${state.experimentSetIds.size} saved set(s) × up to ${representatives} representatives. The manifest records current solvent and optimization mode; it does not start xTB.` : "Select one or more saved pose sets to make a simulation plan.";
    $("#saveRunExperiment").disabled = !state.experimentSetIds.size;
  }

  function runExperimentDesign() { return { pose_set_ids: [...state.experimentSetIds], representatives_per_set: Math.max(1, Number($("#runExperimentPoseCount").value || 8)), simulation_type: "optimization", optimization_mode: $("#runOptimizationMode").value, solvent: $("#runSolvent").value, timeout_seconds: Number($("#runTimeout").value || 300), selection_rule: "select first N saved pose members; replace with energy-ranked representatives only after a completed preliminary screen" }; }

  async function saveRunExperiment() {
    if (!state.experimentSetIds.size) return toast("Select at least one saved pose set.", true);
    const button = $("#saveRunExperiment"); button.disabled = true;
    try { const record = await request("/api/experiments", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ title: $("#runExperimentName").value || "GFN2-xTB experiment", stage: "simulation", design: runExperimentDesign(), notes: "Batch xTB simulation plan. No calculations were started when this manifest was saved." }) }); toast(`Simulation plan saved: ${record.id}`); }
    catch (error) { toast(`Could not save simulation plan: ${error.message}`, true); }
    finally { renderExperimentDesigner(); }
  }

  function render() { renderSet(); renderPoseList(); renderStatus(); renderPreview(); renderViewerPagination(); renderViewerPoseCards(); renderExperimentDesigner(); }

  async function loadViewer() {
    const empty = $("#runViewerEmpty"), loading = $("#runViewerLoading"); const set = state.poseSet; const chosen = selectedModels();
    if (!state.stage || !set || !chosen.length) { empty.classList.remove("hidden"); $("#runViewerLabel").textContent = "No pose selected"; return; }
    const pageIds = new Set(pagePoses().map((pose) => pose.model_id)); const selected = chosen.filter((pose) => pageIds.has(pose.model_id)); const ids = selected.map((pose) => pose.model_id);
    if (!selected.length) { empty.classList.remove("hidden"); $("#runViewerLabel").textContent = "No checked poses on this page"; return; }
    const token = ++state.token; empty.classList.add("hidden"); loading.classList.remove("hidden"); state.stage.removeAllComponents();
    try {
      const url = `/api/pose-sets/${encodeURIComponent(set.id)}/overlay?pose_ids=${encodeURIComponent(ids.join(","))}&t=${Date.now()}`;
      state.component = await state.stage.loadFile(url, { ext: "pdb", defaultRepresentation: false }); if (token !== state.token) return;
      applyRepresentation(); state.component.autoView(420);
      $("#runViewerLabel").textContent = chosen.length > selected.length ? `This page · ${selected.length} shown / ${chosen.length} selected` : `Selected poses · ${selected.length}`;
    } catch (_) { if (token === state.token) { empty.classList.remove("hidden"); $("#runViewerLabel").textContent = "Coordinates unavailable"; } }
    finally { if (token === state.token) loading.classList.add("hidden"); }
  }

  async function openSet(id, preserveSelection = false) {
    if (!id) return;
    state.poseSet = await request(`/api/pose-sets/${encodeURIComponent(id)}`);
    if (!preserveSelection) { state.checked = new Set((state.poseSet.poses || []).map((pose) => pose.model_id)); state.page = 0; state.selectedPoseIndex = 0; }
    state.experimentSetIds.add(state.poseSet.id);
    $("#runPoseSetSelect").value = state.poseSet.id; render(); await loadViewer();
  }

  function chooseBatch() {
    const poses = state.poseSet?.poses || []; const limit = requestedCount();
    const ordered = poses.slice().sort((a, b) => {
      const ea = Number(state.results[a.model_id]?.energy_hartree), eb = Number(state.results[b.model_id]?.energy_hartree);
      return Number.isFinite(ea) && Number.isFinite(eb) ? ea - eb : 0;
    });
    state.checked = new Set(ordered.slice(0, limit).map((pose) => pose.model_id)); state.page = 0; render(); loadViewer();
  }

  function showProgress(job) {
    const node = $("#runJobProgress"); if (!job) { node.classList.add("hidden"); return; }
    const percent = job.total ? Math.round((job.completed || 0) / job.total * 100) : 0; node.classList.remove("hidden"); node.querySelector("span").textContent = `GFN2-xTB · ${job.completed || 0}/${job.total || 0}`; node.querySelector("i").style.width = `${percent}%`; node.querySelector("b").textContent = job.status;
  }

  function currentXtbTask() { return state.poseSet && state.checked.size ? { poseSetId: state.poseSet.id, modelIds: [...state.checked], charge: state.poseSet.formal_charge, mode: $("#runOptimizationMode").value, solvent: $("#runSolvent").value, timeout: Number($("#runTimeout").value || 300) } : null; }
  function renderXtbQueue() { const box = $("#xtbTaskQueue"), run = $("#runXtbQueue"); run.disabled = !state.xtbQueue.length; box.innerHTML = state.xtbQueue.length ? state.xtbQueue.map((task, index) => `<button type="button" data-xtb-queue="${index}"><b>${escapeHtml(`${task.poseSetId.slice(-8)} · ${task.modelIds.length} poses · ${task.mode}`)}</b><i>×</i></button>`).join("") : "<span>No queued xTB batches.</span>"; document.querySelectorAll("[data-xtb-queue]").forEach(button => button.addEventListener("click", () => { state.xtbQueue.splice(Number(button.dataset.xtbQueue), 1); renderXtbQueue(); })); }
  async function runXtbQueue() { $("#runXtbQueue").disabled = true; while (state.xtbQueue.length) { const task = state.xtbQueue.shift(); renderXtbQueue(); await runBatch(task); } toast("xTB queue completed."); }
  async function runBatch(task = null) {
    task = task || currentXtbTask(); if (!task) return;
    state.activeJob = "xtb"; render();
    try {
      const job = await request("/api/xtb/jobs", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ model_ids: task.modelIds, simulation_type: "optimization", optimization_mode: task.mode, solvent: task.solvent, charge: task.charge, uhf: 0, timeout_seconds: task.timeout }) });
      let current = job; while (current.status !== "completed") { showProgress(current); await new Promise((resolve) => setTimeout(resolve, 850)); current = await request(`/api/xtb/jobs/${current.job_id}`); }
      showProgress(current); const failed = (current.results || []).filter((item) => item.status !== "completed").length; toast(failed ? `Batch completed with ${failed} failed run(s).` : "GFN2-xTB batch completed.", Boolean(failed)); await load(true);
    } catch (error) { toast(`xTB batch failed: ${error.message}`, true); }
    finally { state.activeJob = null; showProgress(null); render(); }
  }

  async function load(preserveSelection = false) {
    try {
      const [health, models, sets, results] = await Promise.all([request("/api/health"), request("/api/models?limit=500"), request("/api/pose-sets?limit=100"), request("/api/xtb/results?limit=500")]);
      state.models = models; state.poseSets = sets; state.xtb = health.xtb; state.results = Object.fromEntries(results.filter((item) => item.model_id).map((item) => [item.model_id, item]));
      $("#apiStatus").classList.add("ready"); $("#apiStatus span").textContent = `Engines ready · v${health.version}`;
      const select = $("#runPoseSetSelect"); select.innerHTML = sets.map((set) => `<option value="${escapeHtml(set.id)}">${escapeHtml(set.title || set.id)} · ${set.poses?.length || 0} poses</option>`).join("") || '<option value="">No saved pose sets</option>';
      const active = preserveSelection && state.poseSet?.id ? state.poseSet.id : (state.poseSet?.id || sets[0]?.id); if (active) await openSet(active, preserveSelection); else render();
    } catch (error) { $("#apiStatus").classList.add("error"); $("#apiStatus span").textContent = "Service offline"; toast(error.message, true); }
  }

  setupStage();
  $("#runPoseSetSelect").addEventListener("change", (event) => openSet(event.target.value));
  $("#runXtbButton").addEventListener("click", () => runBatch());
  $("#addXtbTask").addEventListener("click", () => { const task = currentXtbTask(); if (!task) return toast("Select at least one pose first.", true); state.xtbQueue.push(task); renderXtbQueue(); toast("xTB batch added to queue."); });
  $("#runXtbQueue").addEventListener("click", runXtbQueue);
  $("#runOptimizationMode").addEventListener("change", () => { renderPreview(); renderExperimentDesigner(); }); $("#runSolvent").addEventListener("change", renderExperimentDesigner); $("#runTimeout").addEventListener("input", renderExperimentDesigner);
  $("#runOverlayCount").addEventListener("change", () => { const value = $("#runOverlayCount").value; state.pageSize = value === "all" ? "all" : Number(value); state.page = 0; render(); loadViewer(); });
  ["#runRepresentation", "#runColorMode", "#runOverlayStyle"].forEach((id) => $(id).addEventListener("change", applyRepresentation));
  $("#runPreviousPosePage").addEventListener("click", () => { if (state.page > 0) { state.page -= 1; render(); loadViewer(); } }); $("#runNextPosePage").addEventListener("click", () => { const max = Math.ceil(availableCount() / effectivePageSize()) - 1; if (state.page < max) { state.page += 1; render(); loadViewer(); } });
  $("#runExperimentPoseCount")?.addEventListener("input", renderExperimentDesigner); $("#saveRunExperiment")?.addEventListener("click", saveRunExperiment);
  $("#openRunExperiment")?.addEventListener("click", () => { renderExperimentDesigner(); $("#runExperimentDialog")?.showModal(); });
  $("#runResetView").addEventListener("click", () => state.component?.autoView(420));
  renderXtbQueue(); load();
})();
