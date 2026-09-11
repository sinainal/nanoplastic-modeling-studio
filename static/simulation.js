(() => {
  const state = {
    models: [], poseSets: [], poseSet: null, xtbResults: {}, xtbStatus: null, crestStatus: null,
    stage: null, component: null, loadToken: 0, page: 0, pageSize: 8, selectedPoseIndex: 0,
    checked: new Set(), activeJob: null, expandedPoseSetId: null, poseQueue: [],
  };
  const $ = (selector) => document.querySelector(selector);

  const escapeHtml = (value) => String(value ?? "").replace(/[&<>'"]/g, (character) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "'": "&#39;", '"': "&quot;" })[character]);
  function toast(message, isError = false) { const el = $("#toast"); el.textContent = message; el.classList.toggle("error", isError); el.classList.add("show"); clearTimeout(toast.timer); toast.timer = setTimeout(() => el.classList.remove("show"), 3600); }
  async function request(url, options = {}) { const response = await fetch(url, options); const body = await response.json().catch(() => null); if (!response.ok) throw new Error(body?.detail || `${response.status} ${response.statusText}`); return body; }

  function setupStage() {
    if (!window.NGL) return;
    try {
      state.stage = new NGL.Stage("simulationViewport", { backgroundColor: "#e9edf3", quality: "high", tooltip: true });
      state.stage.setParameters({ backgroundColor: "#e9edf3", cameraType: "perspective", clipDist: 0, fogNear: 70, fogFar: 100 });
      window.addEventListener("resize", () => state.stage?.handleResize());
    } catch (_) { state.stage = null; }
  }

  function representationParams(type, surface = false) {
    const style = $("#overlayStyle")?.value || "ghost";
    const ghosted = style !== "solid";
    const palette = $("#poseColorMode")?.value || "pose";
    // In pose mode NGL uses the PDB chain assigned to each adsorbate copy,
    // giving every selected pose a stable colour. Element mode restores the
    // conventional C/O/N/H palette for chemical inspection.
    const color = surface
      ? { colorScheme: "uniform", colorValue: ghosted ? 0x8190a3 : 0x6f7b8c }
      : { colorScheme: palette === "pose" ? "chainid" : "element" };
    const transparency = surface && ghosted ? { opacity: 0.28, depthWrite: false } : {};
    if (type === "spacefill") return { ...color, ...transparency, radiusScale: surface ? 0.48 : 0.68, quality: "high" };
    if (type === "line") return { ...color, ...transparency, multipleBond: true, linewidth: surface ? 1 : 2 };
    if (type === "licorice") return { ...color, ...transparency, multipleBond: "symmetric", radiusScale: surface ? 0.34 : 0.58, quality: "high" };
    return { ...color, ...transparency, multipleBond: "symmetric", aspectRatio: surface ? 2.8 : 2.1, radiusScale: surface ? 0.42 : 0.78, quality: "high" };
  }

  function applyRepresentation() {
    if (!state.component) return;
    const type = $("#simulationRepresentation").value;
    state.component.removeAllRepresentations();
    state.component.addRepresentation(type, { ...representationParams(type, true), sele: ":A" });
    state.component.addRepresentation(type, { ...representationParams(type, false), sele: "not :A" });
  }

  function pagePoses() {
    if (!state.poseSet) return [];
    const pageSize = effectivePageSize();
    const start = state.page * pageSize;
    return state.poseSet.poses.slice(start, start + pageSize);
  }

  function effectivePageSize() {
    const total = state.poseSet?.poses?.length || 1;
    return state.pageSize === "all" ? total : Math.max(1, Number(state.pageSize) || 8);
  }

  async function loadViewer() {
    if (!state.stage || !state.poseSet?.poses?.length) return;
    const loading = $("#simulationViewerLoading"); const empty = $("#simulationViewerEmpty");
    loading.classList.remove("hidden"); empty.classList.add("hidden");
    const token = ++state.loadToken; state.stage.removeAllComponents(); state.component = null;
    try {
      const selected = state.poseSet.poses.filter((pose) => state.checked.has(pose.model_id));
      const pageIds = new Set(pagePoses().map((pose) => pose.model_id));
      const visible = selected.filter((pose) => pageIds.has(pose.model_id));
      if (!visible.length) {
        state.stage.removeAllComponents(); state.component = null;
        empty.classList.remove("hidden");
        $("#simulationViewerLabel").textContent = selected.length ? "No checked poses on this page" : "Select poses to preview";
        return;
      }
      const ids = visible.map((pose) => pose.model_id).map((id) => encodeURIComponent(id)).join(",");
      const url = `/api/pose-sets/${encodeURIComponent(state.poseSet.id)}/overlay?pose_ids=${ids}&t=${Date.now()}`;
      const label = selected.length > visible.length
        ? `This page · ${visible.length} shown / ${selected.length} selected`
        : `Selected poses · ${visible.length}`;
      state.component = await state.stage.loadFile(url, { ext: "pdb", defaultRepresentation: false });
      if (token !== state.loadToken) return;
      applyRepresentation(); state.component.autoView(420); $("#simulationViewerLabel").textContent = label;
    } catch (error) {
      if (token === state.loadToken) { empty.classList.remove("hidden"); $("#simulationViewerLabel").textContent = "Coordinates unavailable"; toast(`3D preview failed: ${error.message}`, true); }
    } finally { if (token === state.loadToken) loading.classList.add("hidden"); }
  }

  function preferredSurface(models) {
    const candidates = models.filter((m) => m.source?.family === "polymer");
    return candidates.find((m) => m.source?.polymer === "PE" && Number(m.source?.repeats) === 9 && m.engine === "RDKit ETKDGv3") || candidates.find((m) => m.source?.polymer === "PE") || candidates[0];
  }
  function preferredAdsorbate(models) {
    const candidates = models.filter((m) => m.source?.compound === "dopamine");
    return candidates.find((m) => m.source?.microstate === "DAH_plus" && Number(m.source?.pH) === 7) || candidates.find((m) => m.source?.microstate === "DAH_plus") || candidates[0];
  }

  function populateModelSelectors() {
    const surfaces = state.models.filter((model) => model.source?.family === "polymer");
    const adsorbates = state.models.filter((model) => model.source?.compound === "dopamine");
    const surface = preferredSurface(state.models); const adsorbate = preferredAdsorbate(state.models);
    $("#surfaceModel").innerHTML = surfaces.map((model) => `<option value="${escapeHtml(model.id)}">${escapeHtml(model.source?.polymer || "Polymer")} · ${model.source?.repeats || "?"}-mer · ${escapeHtml(model.engine)}</option>`).join("") || '<option value="">No polymer model</option>';
    $("#adsorbateModel").innerHTML = adsorbates.map((model) => `<option value="${escapeHtml(model.id)}">${escapeHtml(model.source?.microstate || model.title)} · pH ${escapeHtml(model.source?.pH ?? "—")}</option>`).join("") || '<option value="">No dopamine model</option>';
    if (surface) $("#surfaceModel").value = surface.id; if (adsorbate) $("#adsorbateModel").value = adsorbate.id;
    renderSystemSummary();
  }

  function renderPoseSetTable() {
    const body = $("#poseSetResults");
    if (!body) return;
    const active = state.poseSet?.id || "";
    const expanded = state.expandedPoseSetId;
    $("#poseSetLibraryMeta").textContent = `${state.poseSets.length} saved · ${active ? "active set selected" : "select a set"}`;
    if (!state.poseSets.length) {
      body.innerHTML = '<tr><td colspan="7" class="table-empty">No saved pose sets.</td></tr>';
      return;
    }
    const strategyLabel = (strategy) => strategy === "surface_scan" ? "SASA surface scan" : strategy === "crest" ? "CREST local refinement" : "Legacy systematic";
    const createdLabel = (value) => value ? new Date(value).toLocaleString(undefined, { year: "numeric", month: "short", day: "numeric", hour: "2-digit", minute: "2-digit" }) : "—";
    body.innerHTML = state.poseSets.map((item) => {
      const isActive = item.id === active;
      const isExpanded = item.id === expanded;
      const charge = Number(item.formal_charge || 0);
      const poseCount = item.poses?.length || 0;
      const detail = isExpanded ? `<tr class="pose-set-detail-row"><td colspan="7"><div class="pose-set-detail">
        <div class="pose-set-detail-head">
          <div><span class="workspace-kicker">POSE MEMBERS</span><h3>${escapeHtml(item.title || "Pose set")}</h3><p id="poseTableMeta">${poseCount} poses · ${escapeHtml(item.strategy || "systematic")} · charge ${charge >= 0 ? "+" : ""}${charge}</p></div>
          <div class="pose-member-actions"><button id="selectPageButton" class="secondary-button compact-button" type="button">Select best N</button><button id="runCrestButton" class="secondary-button compact-button" type="button">Run CREST on selected</button><a class="text-button" href="/analysis">Open analysis</a></div>
        </div>
        <div class="table-scroll pose-set-detail-scroll"><table class="pose-set-detail-table"><thead><tr><th></th><th>Pose</th><th>Method</th><th>Geometry</th><th>xTB energy</th><th>H–L gap</th><th>Status</th><th></th></tr></thead><tbody id="poseResults"><tr><td colspan="8" class="table-empty">Loading poses…</td></tr></tbody></table></div>
      </div></td></tr>` : "";
      return `<tr class="pose-set-row ${isActive ? "active" : ""} ${isExpanded ? "expanded" : ""}" data-pose-set-id="${escapeHtml(item.id)}">
        <td class="pose-set-toggle-cell"><button class="pose-set-toggle" type="button" data-pose-set-id="${escapeHtml(item.id)}" aria-expanded="${isExpanded ? "true" : "false"}" aria-label="${isExpanded ? "Collapse" : "Expand"} ${escapeHtml(item.title || "pose set")}">${isExpanded ? "⌄" : "›"}</button></td>
        <td><b>${escapeHtml(item.title || "Untitled pose set")}</b><small class="table-subtext">${isActive ? "Active workspace" : "Saved pose set"}</small></td>
        <td><span class="method-pill">${escapeHtml(strategyLabel(item.strategy))}</span></td>
        <td class="mono">${poseCount}</td><td class="mono">${charge >= 0 ? "+" : ""}${charge}</td><td>${escapeHtml(createdLabel(item.created_at))}</td>
        <td><button class="table-action pose-set-open" type="button" data-pose-set-id="${escapeHtml(item.id)}">${isActive ? "Open" : "Load"}</button></td>
      </tr>${detail}`;
    }).join("");
    body.querySelectorAll(".pose-set-toggle").forEach((button) => button.addEventListener("click", () => {
      const id = button.dataset.poseSetId;
      if (state.expandedPoseSetId === id) { state.expandedPoseSetId = null; renderAll(); }
      else { openPoseSet(id); }
    }));
    body.querySelectorAll(".pose-set-open").forEach((button) => button.addEventListener("click", () => openPoseSet(button.dataset.poseSetId)));
    $("#selectPageButton")?.addEventListener("click", selectBestPoses);
    $("#runCrestButton")?.addEventListener("click", () => runJob("crest"));
  }

  function renderSystemSummary() {
    const surface = state.models.find((m) => m.id === $("#surfaceModel").value); const adsorbate = state.models.find((m) => m.id === $("#adsorbateModel").value);
    if (!surface || !adsorbate) { $("#systemSummary").innerHTML = "<span>Select valid models.</span>"; return; }
    const charge = Number(surface.formal_charge || 0) + Number(adsorbate.formal_charge || 0);
    $("#systemSummary").innerHTML = `<span><b>${escapeHtml(surface.formula)}</b> surface</span><i>+</i><span><b>${escapeHtml(adsorbate.formula)}</b> adsorbate</span><i>·</i><span>complex charge <b>${charge >= 0 ? "+" : ""}${charge}</b></span>`;
    renderCli();
  }

  function syncScanUi() {
    const sites = Math.max(1, Math.min(25, Number($("#surfaceSiteCount").value || 5)));
    const orientations = Math.max(1, Math.min(8, Number($("#orientationsPerSite").value || 4)));
    const total = sites * orientations;
    $("#generatePoseSet").textContent = `Generate ${total} surface poses`;
    renderRunPlan();
    renderPoseExperiment();
    renderCli();
  }

  function checkedExperimentValues(containerId) {
    return Array.from(document.querySelectorAll(`${containerId} input:checked`)).map((input) => input.value);
  }

  function poseExperimentDesign() {
    const polymers = checkedExperimentValues("#poseExperimentPolymers");
    const repeats = checkedExperimentValues("#poseExperimentRepeats").map(Number);
    const conformers = Math.max(1, Number($("#poseExperimentConformers")?.value || 1));
    const seeds = Math.max(1, Number($("#poseExperimentSeeds")?.value || 1));
    const sites = Math.max(1, Number($("#surfaceSiteCount")?.value || 5));
    const orientations = Math.max(1, Number($("#orientationsPerSite")?.value || 4));
    const adsorbate = state.models.find((model) => model.id === $("#adsorbateModel")?.value);
    return {
      polymers, repeats, conformers_per_composition: conformers, independent_seeds: seeds,
      poses_per_conformer: sites * orientations,
      pose_generation: { strategy: "surface_scan", site_count: sites, orientations_per_site: orientations, distance_A: Number($("#poseDistance")?.value || 3.2), probe_radius_A: Number($("#surfaceProbeRadius")?.value || 1.4) },
      adsorbate: { model_id: adsorbate?.id || null, title: adsorbate?.title || $("#adsorbateModel option:checked")?.textContent || "selected adsorbate", formal_charge: adsorbate?.formal_charge ?? null },
    };
  }

  function renderPoseExperiment() {
    const summary = $("#poseExperimentSummary");
    if (!summary) return;
    const design = poseExperimentDesign();
    const compositions = design.polymers.length * design.repeats.length;
    const modelCount = compositions * design.conformers_per_composition * design.independent_seeds;
    const poses = modelCount * design.poses_per_conformer;
    summary.textContent = `${modelCount || 0} models · ${poses || 0} poses`;
    const detail = $("#poseExperimentDetail");
    detail.textContent = compositions
      ? `${compositions} polymer–size compositions × ${design.conformers_per_composition} conformer(s) × ${design.independent_seeds} seed(s) × ${design.poses_per_conformer} poses. Saving creates a planning manifest only.`
      : "Select at least one polymer and oligomer size to define a matrix.";
    $("#savePoseExperiment").disabled = !compositions;
  }

  async function savePoseExperiment() {
    const design = poseExperimentDesign();
    if (!design.polymers.length || !design.repeats.length) return toast("Select at least one polymer and oligomer size.", true);
    const button = $("#savePoseExperiment"); button.disabled = true;
    try {
      const record = await request("/api/experiments", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ title: $("#poseExperimentName").value || "Pose generation experiment", stage: "pose_generation", design, notes: "Batch pose-generation matrix. No structures were started by this plan." }) });
      toast(`Pose plan saved: ${record.id}`);
    } catch (error) { toast(`Could not save pose plan: ${error.message}`, true); }
    finally { button.disabled = false; renderPoseExperiment(); }
  }

  function syncOptimizationUi() {
    const constrained = $("#optimizationMode").value === "surface_fixed";
    $("#optimizationModeNote").textContent = constrained
      ? "Surface-fixed screening keeps every polymer atom fixed so surface sites remain comparable. Use full refinement only for selected representatives."
      : "Full refinement releases the complete PE–dopamine complex. Run it on a small, diverse set after the surface-fixed screen.";
    $("#runXtbButton").textContent = constrained ? "Run constrained xTB" : "Run full xTB refinement";
    renderRunPlan();
    renderCli();
  }

  function renderEngineStatus() {
    const setBadge = (id, name, status) => { const el = $(id); const ready = status?.status === "ready"; el.textContent = `${name} · ${ready ? "ready" : "unavailable"}`; el.classList.toggle("ready", ready); el.classList.toggle("error", !ready); };
    setBadge("#crestBadge", "CREST", state.crestStatus); setBadge("#xtbBadge", "xTB", state.xtbStatus);
  }

  function requestedSimulationCount() {
    const total = state.poseSet?.poses?.length || 0;
    const requested = Number.parseInt($("#simulationPoseCount")?.value || "8", 10);
    return total ? Math.min(total, Math.max(1, Number.isFinite(requested) ? requested : 8)) : Math.max(1, Number.isFinite(requested) ? requested : 8);
  }

  function renderRunPlan() {
    const target = $("#simulationRunPlan");
    if (!target) return;
    const sites = Math.max(1, Math.min(25, Number($("#surfaceSiteCount")?.value || 5)));
    const orientations = Math.max(1, Math.min(8, Number($("#orientationsPerSite")?.value || 4)));
    const planned = sites * orientations;
    const poseCount = state.poseSet?.poses?.length || 0;
    const selected = state.checked?.size || 0;
    const mode = $("#optimizationMode")?.value === "surface_fixed" ? "surface-fixed" : "full refinement";
    const coverage = poseCount ? `${state.poseSet.site_count} sites × ${state.poseSet.orientations_per_site} orientations = ${poseCount} poses` : `${sites} sites × ${orientations} orientations = ${planned} planned poses`;
    const batch = poseCount ? `${selected}/${poseCount} selected · ${mode}` : "Generate a pose set before launching xTB";
    target.innerHTML = `<div><span>Pose coverage</span><b>${coverage}</b></div><div><span>xTB batch</span><b>${batch}</b></div><div><span>Scientific output</span><b>Optimized local minima; Eads requires matched references</b></div>`;
  }

  function selectAllPoseMembers() {
    state.checked = new Set((state.poseSet?.poses || []).map((pose) => pose.model_id));
  }

  function renderSelection() {
    const count = state.checked.size; $("#selectionBadge").textContent = `${count} selected`;
    $("#runXtbButton").disabled = !state.poseSet || !count || state.xtbStatus?.status !== "ready" || Boolean(state.activeJob);
    const requested = requestedSimulationCount();
    const selectButton = $("#selectBestPosesButton");
    if (selectButton) { selectButton.textContent = `Select best ${requested}`; selectButton.disabled = !state.poseSet || !(state.poseSet.poses || []).length || Boolean(state.activeJob); }
    const crestButton = $("#runCrestButton");
    if (crestButton) crestButton.disabled = !state.poseSet || !count || state.crestStatus?.status !== "ready" || Boolean(state.activeJob);
    const pageButton = $("#selectPageButton");
    if (pageButton) { pageButton.textContent = `Select best ${requested}`; pageButton.disabled = !state.poseSet || !(state.poseSet.poses || []).length || Boolean(state.activeJob); }
    renderRunPlan();
    renderCli();
  }

  function selectBestPoses() {
    const poses = state.poseSet?.poses || [];
    if (!poses.length) return;
    const requested = requestedSimulationCount();
    const ranked = poses.map((pose, index) => ({ pose, index, result: state.xtbResults[pose.model_id] })).sort((left, right) => {
      const leftEnergy = Number(left.result?.energy_hartree); const rightEnergy = Number(right.result?.energy_hartree);
      const leftHasEnergy = left.result?.status === "completed" && Number.isFinite(leftEnergy);
      const rightHasEnergy = right.result?.status === "completed" && Number.isFinite(rightEnergy);
      if (leftHasEnergy && rightHasEnergy && leftEnergy !== rightEnergy) return leftEnergy - rightEnergy;
      if (leftHasEnergy !== rightHasEnergy) return leftHasEnergy ? -1 : 1;
      const leftCrest = String(left.pose.method || "").startsWith("crest");
      const rightCrest = String(right.pose.method || "").startsWith("crest");
      if (leftCrest !== rightCrest) return leftCrest ? -1 : 1;
      return left.index - right.index;
    });
    state.checked = new Set(ranked.slice(0, requested).map(({ pose }) => pose.model_id));
    renderAll();
    toast(`Selected ${Math.min(requested, poses.length)} best-available pose(s); completed xTB energies are ranked first.`);
  }

  function renderPoseCards() {
    const grid = $("#poseCardGrid"); const poses = pagePoses();
    if (!poses.length) { grid.innerHTML = '<div class="pose-grid-empty">No poses on this page.</div>'; return; }
    const start = state.page * effectivePageSize();
    grid.innerHTML = poses.map((pose, localIndex) => {
      const globalIndex = start + localIndex; const result = state.xtbResults[pose.model_id]; const active = globalIndex === state.selectedPoseIndex;
      const parameters = pose.parameters || {}; const site = parameters.site_id ? `${parameters.site_id} · ${parameters.surface_region}` : pose.label;
      const orientation = parameters.label || parameters.name || pose.method || "surface pose";
      return `<article class="pose-card ${active ? "active" : ""}" data-pose-index="${globalIndex}"><label><input class="pose-check" type="checkbox" data-model-id="${escapeHtml(pose.model_id)}" ${state.checked.has(pose.model_id) ? "checked" : ""}><span>${escapeHtml(pose.pose_id)}</span></label><b>${escapeHtml(site)}</b><small>${escapeHtml(orientation)}</small><i class="pose-status-dot ${result?.status === "completed" ? "complete" : result?.status === "failed" ? "failed" : ""}"></i></article>`;
    }).join("");
    grid.querySelectorAll(".pose-check").forEach((input) => input.addEventListener("change", (event) => { event.stopPropagation(); if (input.checked) state.checked.add(input.dataset.modelId); else state.checked.delete(input.dataset.modelId); renderSelection(); renderPoseTable(); loadViewer(); }));
    grid.querySelectorAll(".pose-card").forEach((card) => card.addEventListener("click", (event) => { if (event.target.closest("label")) return; state.selectedPoseIndex = Number(card.dataset.poseIndex); renderPoseCards(); loadViewer(); }));
  }

  function renderPagination() {
    const total = state.poseSet?.poses?.length || 0; const pageSize = effectivePageSize(); const pages = Math.max(1, Math.ceil(total / pageSize)); state.page = Math.min(state.page, pages - 1);
    const start = total ? state.page * pageSize + 1 : 0; const end = Math.min((state.page + 1) * pageSize, total);
    $("#posePageLabel").textContent = total ? `${start}–${end} of ${total} poses` : "No poses";
    $("#previousPosePage").disabled = !total || state.page === 0; $("#nextPosePage").disabled = !total || state.page >= pages - 1;
  }

  function renderPoseTable() {
    const body = $("#poseResults"); const poses = state.poseSet?.poses || [];
    if (!body) return;
    if (!poses.length) { body.innerHTML = '<tr><td colspan="8" class="table-empty">No active pose set.</td></tr>'; return; }
    $("#poseTableMeta").textContent = `${poses.length} poses · ${state.poseSet.strategy} · charge ${state.poseSet.formal_charge >= 0 ? "+" : ""}${state.poseSet.formal_charge}`;
    body.innerHTML = poses.map((pose, index) => {
      const result = state.xtbResults[pose.model_id]; const parameters = pose.parameters || {}; const status = result?.status || "not run";
      const orientation = parameters.tilt == null ? parameters.kind || "—" : `${Number(parameters.tilt).toFixed(0)}° tilt · ${Number(parameters.azimuth || 0).toFixed(0)}° az`;
      const site = parameters.site_id ? `${parameters.site_id} · ${parameters.surface_region || "surface"}` : "Legacy centre";
      const mode = result?.optimization_mode === "surface_fixed" ? "fixed" : result?.optimization_mode === "full" ? "full" : "";
      return `<tr class="${index === state.selectedPoseIndex ? "selected-row" : ""}"><td><input class="table-pose-check" type="checkbox" data-model-id="${escapeHtml(pose.model_id)}" ${state.checked.has(pose.model_id) ? "checked" : ""}></td><td><b>${escapeHtml(pose.pose_id)}</b><small>${escapeHtml(site)}</small></td><td>${escapeHtml(pose.method)}</td><td><b>${escapeHtml(parameters.label || parameters.name || "orientation")}</b><small>${escapeHtml(orientation)}</small></td><td class="mono">${result?.energy_hartree == null ? "—" : Number(result.energy_hartree).toFixed(7)}</td><td>${result?.homo_lumo_gap_ev == null ? "—" : Number(result.homo_lumo_gap_ev).toFixed(3)}</td><td><span class="validation-pill ${status === "completed" ? "ok" : "pending"}">${escapeHtml(status)}${mode ? ` · ${mode}` : ""}</span></td><td><button class="table-action pose-view-action" type="button" data-pose-index="${index}">Focus</button></td></tr>`;
    }).join("");
    body.querySelectorAll(".table-pose-check").forEach((input) => input.addEventListener("change", () => { if (input.checked) state.checked.add(input.dataset.modelId); else state.checked.delete(input.dataset.modelId); renderSelection(); renderPoseCards(); loadViewer(); }));
    body.querySelectorAll(".pose-view-action").forEach((button) => button.addEventListener("click", () => { state.selectedPoseIndex = Number(button.dataset.poseIndex); state.page = Math.floor(state.selectedPoseIndex / effectivePageSize()); renderAll(); loadViewer(); }));
  }

  function renderAll() {
    const active = Boolean(state.poseSet);
    const surface = state.models.find((model) => model.id === state.poseSet?.surface_model_id);
    const adsorbate = state.models.find((model) => model.id === state.poseSet?.adsorbate_model_id);
    const compactTitle = state.poseSet
      ? `${surface?.source?.polymer || "Polymer"} · ${surface?.source?.repeats || "?"}-mer + ${adsorbate?.source?.microstate || "adsorbate"} · ${state.poseSet.site_count || "?"} sites / ${state.poseSet.poses.length} poses`
      : "Create or open a pose set";
    $("#simulationModelTitle").textContent = compactTitle;
    renderPagination(); renderPoseCards(); renderPoseSetTable(); renderPoseTable(); renderSelection();
  }

  function cliCommand() {
    const selected = Array.from(state.checked); const poseIds = state.poseSet?.poses?.filter((p) => state.checked.has(p.model_id)).map((p) => p.pose_id) || [];
    if (state.poseSet && selected.length) return `python -m modeling.app.cli pose-xtb --pose-set-id ${state.poseSet.id} --poses ${poseIds.join(" ")} --solvent ${$("#simulationSolvent").value} --optimization-mode ${$("#optimizationMode").value}`;
    const surface = $("#surfaceModel").value || "<surface-model>"; const adsorbate = $("#adsorbateModel").value || "<adsorbate-model>";
    return `python -m modeling.app.cli pose-create --surface-model ${surface} --adsorbate-model ${adsorbate} --sites ${Number($("#surfaceSiteCount").value || 5)} --orientations ${Number($("#orientationsPerSite").value || 4)} --distance ${Number($("#poseDistance").value || 3.2)} --probe-radius ${Number($("#surfaceProbeRadius").value || 1.4)} --seed ${Number($("#poseSeed").value || 20260827)} --strategy surface_scan`;
  }
  function renderCli() { $("#cliPreview").textContent = cliCommand(); }

  async function loadAll(keepPoseSet = true) {
    try {
      const hadPoseSet = Boolean(state.poseSet);
      const [health, models, poseSets, results] = await Promise.all([request("/api/health"), request("/api/models?limit=500"), request("/api/pose-sets?limit=100"), request("/api/xtb/results?limit=500")]);
      state.models = models; state.poseSets = poseSets; state.xtbStatus = health.xtb; state.crestStatus = health.crest; state.xtbResults = {}; results.forEach((item) => { if (item.model_id) state.xtbResults[item.model_id] = item; });
      $("#apiStatus").classList.add("ready"); $("#apiStatus span").textContent = `Engines ready · v${health.version}`;
      populateModelSelectors(); renderEngineStatus();
      renderPoseExperiment();
      if (keepPoseSet && state.poseSet) state.poseSet = await request(`/api/pose-sets/${state.poseSet.id}`);
      else if (!state.poseSet && poseSets.length) state.poseSet = poseSets[0];
      if (!hadPoseSet && state.poseSet) selectAllPoseMembers();
      state.expandedPoseSetId = state.poseSet?.id || null;
      renderAll(); if (state.poseSet) loadViewer();
    } catch (error) { $("#apiStatus").classList.add("error"); $("#apiStatus span").textContent = "Service offline"; toast(error.message, true); }
  }

  function posePayloadFromForm() {
    return { surface_model_id: $("#surfaceModel").value, adsorbate_model_id: $("#adsorbateModel").value, site_count: Number($("#surfaceSiteCount").value || 5), orientations_per_site: Number($("#orientationsPerSite").value || 4), distance_A: Number($("#poseDistance").value || 3.2), probe_radius_A: Number($("#surfaceProbeRadius").value || 1.4), seed: Number($("#poseSeed").value || 20260827), strategy: "surface_scan" };
  }
  function renderPoseQueue() {
    const box = $("#poseTaskQueue"), run = $("#runPoseQueue"); run.disabled = !state.poseQueue.length;
    box.innerHTML = state.poseQueue.length ? state.poseQueue.map((p, i) => `<button type="button" data-pose-queue="${i}"><b>${escapeHtml(`${p.surface_model_id.slice(-8)} + ${p.adsorbate_model_id.slice(-8)} · ${p.site_count}×${p.orientations_per_site}`)}</b><i>×</i></button>`).join("") : "<span>No queued pose jobs.</span>";
    document.querySelectorAll("[data-pose-queue]").forEach(button => button.addEventListener("click", () => { state.poseQueue.splice(Number(button.dataset.poseQueue), 1); renderPoseQueue(); }));
  }
  async function runPoseQueue() {
    $("#runPoseQueue").disabled = true;
    while (state.poseQueue.length) { const payload = state.poseQueue.shift(); renderPoseQueue(); await createPoseSet(payload); }
    toast("Pose-generation queue completed.");
  }

  async function createPoseSet(queuedPayload = null) {
    const button = $("#generatePoseSet"); button.disabled = true;
    try {
      const payload = queuedPayload || posePayloadFromForm();
      state.poseSet = await request("/api/pose-sets", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) }); state.poseSets = [state.poseSet, ...state.poseSets.filter((item) => item.id !== state.poseSet.id)]; state.expandedPoseSetId = state.poseSet.id; state.page = 0; state.selectedPoseIndex = 0; selectAllPoseMembers(); renderAll(); await loadViewer();
      toast(`${state.poseSet.site_count} surface sites × ${state.poseSet.orientations_per_site} orientations generated (${state.poseSet.poses.length} poses).`);
    } catch (error) { toast(`Pose generation failed: ${error.message}`, true); }
    finally { state.activeJob = null; setProgress(null); button.disabled = false; renderSelection(); }
  }

  async function openPoseSet(id) {
    if (!id) return;
    try { state.poseSet = await request(`/api/pose-sets/${encodeURIComponent(id)}`); state.expandedPoseSetId = id; state.page = 0; state.selectedPoseIndex = 0; selectAllPoseMembers(); renderAll(); await loadViewer(); }
    catch (error) { toast(error.message, true); }
  }

  function setProgress(job, label) {
    const element = $("#jobProgress"); if (!job) { element.classList.add("hidden"); return; }
    const percent = job.total ? Math.round((job.completed || 0) / job.total * 100) : 0; element.classList.remove("hidden"); element.querySelector("span").textContent = `${label} · ${job.completed || 0}/${job.total || 0}`; element.querySelector("i").style.width = `${percent}%`; element.querySelector("b").textContent = job.status;
  }

  async function pollJob(kind, jobId) {
    while (true) {
      const job = await request(`/api/${kind}/jobs/${jobId}`); setProgress(job, kind === "xtb" ? "GFN2-xTB" : "CREST");
      if (job.status === "completed") return job;
      await new Promise((resolve) => setTimeout(resolve, 900));
    }
  }

  async function runJob(kind) {
    const modelIds = Array.from(state.checked); if (!modelIds.length) return;
    state.activeJob = kind; renderSelection();
    try {
      const payload = kind === "xtb" ? { model_ids: modelIds, simulation_type: "optimization", optimization_mode: $("#optimizationMode").value, solvent: $("#simulationSolvent").value, charge: state.poseSet.formal_charge, uhf: 0, timeout_seconds: Number($("#simulationTimeout").value || 300) } : { model_ids: modelIds, top_k: Number($("#crestTopK").value || 3), solvent: $("#crestSolvent").value, timeout_seconds: 600 };
      const job = await request(`/api/${kind}/jobs`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) }); const finished = await pollJob(kind, job.job_id); const failed = (finished.results || []).filter((item) => item.status !== "completed").length;
      if (kind === "crest") { state.poseSet = await request(`/api/pose-sets/${state.poseSet.id}`); state.poseSets = [state.poseSet, ...state.poseSets.filter((item) => item.id !== state.poseSet.id)]; state.expandedPoseSetId = state.poseSet.id; selectAllPoseMembers(); }
      await loadAll(true); toast(failed ? `${kind.toUpperCase()} completed with ${failed} failed run(s).` : `${kind.toUpperCase()} job completed.`, Boolean(failed));
    } catch (error) { toast(`${kind.toUpperCase()} job failed: ${error.message}`, true); }
    finally { state.activeJob = null; setProgress(null); renderSelection(); }
  }

  function togglePageSelection() { selectBestPoses(); }

  setupStage();
  const refreshButton = $("#refreshSimulation");
  if (refreshButton) refreshButton.addEventListener("click", () => loadAll(true));
  $("#surfaceModel").addEventListener("change", () => { renderSystemSummary(); renderPoseExperiment(); }); $("#adsorbateModel").addEventListener("change", () => { renderSystemSummary(); renderPoseExperiment(); });
  ["#surfaceSiteCount", "#orientationsPerSite", "#poseDistance", "#surfaceProbeRadius", "#poseSeed"].forEach((id) => $(id).addEventListener("change", syncScanUi));
  ["#crestSolvent", "#crestTopK", "#simulationSolvent"].forEach((id) => $(id).addEventListener("change", renderCli));
  $("#optimizationMode").addEventListener("change", syncOptimizationUi);
  $("#generatePoseSet").addEventListener("click", () => createPoseSet());
  $("#addPoseTask").addEventListener("click", () => { state.poseQueue.push(posePayloadFromForm()); renderPoseQueue(); toast("Pose setup added to queue."); });
  $("#runPoseQueue").addEventListener("click", runPoseQueue);
  $("#selectBestPosesButton").addEventListener("click", selectBestPoses);
  $("#simulationPoseCount").addEventListener("input", renderSelection);
  $("#simulationPoseCount").addEventListener("change", () => { const value = requestedSimulationCount(); $("#simulationPoseCount").value = String(value); renderSelection(); });
  $("#simulationRepresentation").addEventListener("change", applyRepresentation); $("#poseColorMode").addEventListener("change", applyRepresentation); $("#simulationResetView").addEventListener("click", () => state.component?.autoView(350));
  $("#overlayStyle").addEventListener("change", applyRepresentation);
  $("#overlayPageSize").addEventListener("change", () => { const value = $("#overlayPageSize").value; state.pageSize = value === "all" ? "all" : Number(value); state.page = 0; renderAll(); loadViewer(); });
  $("#previousPosePage").addEventListener("click", () => { if (state.page > 0) { state.page -= 1; renderAll(); loadViewer(); } });
  $("#nextPosePage").addEventListener("click", () => { const max = Math.ceil((state.poseSet?.poses?.length || 0) / effectivePageSize()) - 1; if (state.page < max) { state.page += 1; renderAll(); loadViewer(); } });
  $("#runXtbButton").addEventListener("click", () => runJob("xtb"));
  $("#copyCliButton").addEventListener("click", async () => { const command = cliCommand(); try { await navigator.clipboard.writeText(command); toast("CLI command copied."); } catch (_) { toast(command); } });
  ["#poseExperimentConformers", "#poseExperimentSeeds"].forEach((id) => $(id)?.addEventListener("input", renderPoseExperiment));
  document.querySelectorAll("#poseExperimentPolymers input, #poseExperimentRepeats input").forEach((input) => input.addEventListener("change", renderPoseExperiment));
  $("#savePoseExperiment")?.addEventListener("click", savePoseExperiment);
  $("#openPoseExperiment")?.addEventListener("click", () => { renderPoseExperiment(); $("#poseExperimentDialog")?.showModal(); });
  syncScanUi(); syncOptimizationUi(); renderPoseQueue(); loadAll(false);
})();
