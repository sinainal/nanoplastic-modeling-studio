(() => {
  const state = {
    catalog: null,
    stage: null,
    component: null,
    activeModel: null,
    xtbStatus: null,
    xtbResults: {},
    spinning: false,
    history: [],
    modelQueue: [],
  };

  const $ = (selector) => document.querySelector(selector);
  const $$ = (selector) => Array.from(document.querySelectorAll(selector));

  function toast(message, isError = false) {
    const element = $("#toast");
    element.textContent = message;
    element.classList.toggle("error", isError);
    element.classList.add("show");
    clearTimeout(toast.timer);
    toast.timer = setTimeout(() => element.classList.remove("show"), 3200);
  }

  async function request(url, options = {}) {
    const response = await fetch(url, options);
    let body = null;
    try { body = await response.json(); } catch (_) { body = null; }
    if (!response.ok) throw new Error(body?.detail || `${response.status} ${response.statusText}`);
    return body;
  }

  function setupStage() {
    if (!window.NGL) {
      toast("The NGL viewer could not be loaded.", true);
      return;
    }
    try {
      state.stage = new NGL.Stage("viewport", {
        backgroundColor: "#e9edf3",
        quality: "high",
        tooltip: true,
      });
      state.stage.setParameters({ cameraType: "perspective", clipDist: 0, fogNear: 70, fogFar: 100 });
      window.addEventListener("resize", () => state.stage?.handleResize());
    } catch (error) {
      state.stage = null;
      toast("WebGL is unavailable; model metadata remains accessible.", true);
    }
  }

  function representationParams(type) {
    if (type === "spacefill") return { colorScheme: "element", radiusScale: 0.82, quality: "high" };
    if (type === "line") return { colorScheme: "element", multipleBond: true, linewidth: 2 };
    if (type === "licorice") return { colorScheme: "element", multipleBond: "symmetric", radiusScale: 0.75, quality: "high" };
    return { colorScheme: "element", multipleBond: "symmetric", aspectRatio: 1.8, quality: "high" };
  }

  function applyRepresentation() {
    if (!state.component) return;
    const type = $("#representation").value;
    state.component.removeAllRepresentations();
    state.component.addRepresentation(type, representationParams(type));
  }

  async function loadStructure(model) {
    if (!state.stage) return;
    $("#viewerLoading").classList.remove("hidden");
    $("#viewerEmpty").classList.add("hidden");
    state.stage.removeAllComponents();
    state.component = null;
    try {
      const component = await state.stage.loadFile(`/api/models/${model.id}/structure?format=pdb&t=${Date.now()}`, {
        ext: "pdb",
        defaultRepresentation: false,
      });
      state.component = component;
      applyRepresentation();
      component.autoView(450);
      state.stage.setSpin(state.spinning);
    } catch (error) {
      $("#viewerEmpty").classList.remove("hidden");
      toast(`Structure could not be displayed: ${error.message}`, true);
    } finally {
      $("#viewerLoading").classList.add("hidden");
    }
  }

  function showModel(model) {
    state.activeModel = model;
    $("#modelTitle").textContent = model.title;
    $("#metricFormula").textContent = model.formula;
    $("#metricAtoms").textContent = String(model.atom_count);
    $("#metricCharge").textContent = model.formal_charge > 0 ? `+${model.formal_charge}` : String(model.formal_charge);
    $("#metricEngine").textContent = model.engine;
    $("#metricTime").textContent = `${model.wall_time_seconds.toFixed(2)} s`;

    const source = model.source || {};
    const details = [
      ["Model ID", model.id],
      ["Model level", model.model_kind === "surface_patch" ? "2×2 surface proxy" : "molecule / oligomer"],
      ["Optimization", model.optimization],
      ["Conformers", `${model.conformers_generated}/${model.conformers_requested}`],
      ["Seed", model.seed],
    ];
    if (source.microstate) details.splice(2, 0, ["Microstate", `${source.microstate} · pH ${source.pH}`]);
    if (source.polymer) {
      details.splice(2, 0, ["Polymer", `${source.polymer} · ${source.repeats} repeat`]);
      if (source.pH !== undefined) details.splice(3, 0, ["Environment", `pH ${source.pH} · ${Number(source.temperature_K || 310.15).toFixed(2)} K`]);
    }
    $("#modelDetails").innerHTML = details.map(([key, value]) => `<div><dt>${escapeHtml(key)}</dt><dd>${escapeHtml(String(value))}</dd></div>`).join("");

    const warnings = model.warnings?.length ? model.warnings : ["This model passed basic 3D geometry validation."];
    $("#warnings").innerHTML = warnings.map((warning) => `<p>${escapeHtml(warning)}</p>`).join("");
    if ($("#runXtbButton")) $("#runXtbButton").disabled = false;
    const savedXtb = state.xtbResults[model.id];
    if (savedXtb) renderXtbResult(savedXtb);
    else renderXtbStatus(state.xtbStatus);
    renderHistory();
    loadStructure(model);
  }

  function renderXtbStatus(status) {
    const element = $("#xtbStatus");
    const output = $("#xtbOutput");
    if (!element || !output) return;
    const available = status?.status === "ready";
    element.classList.toggle("ready", available);
    element.classList.toggle("error", !available);
    element.querySelector("span").textContent = available
      ? `Ready · ${status.method || "GFN2-xTB"}`
      : (status?.reason || "xTB is not available on this machine.");
    if (!available) output.textContent = "Install xtb 6.7+ and set XTB_BIN if the executable is outside PATH.";
  }

  function renderXtbResult(result) {
    const element = $("#xtbStatus");
    const output = $("#xtbOutput");
    if (!element || !output) return;
    const complete = result?.status === "completed";
    element.classList.toggle("ready", complete);
    element.classList.toggle("error", !complete);
    element.querySelector("span").textContent = complete
      ? `Completed · ${Number(result.wall_time_seconds || 0).toFixed(2)} s`
      : (result?.reason || "xTB did not run");
    output.textContent = JSON.stringify(result, null, 2);
  }

  function escapeHtml(value) {
    return value.replace(/[&<>'"]/g, (character) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "'": "&#39;", '"': "&quot;" })[character]);
  }

  function renderHistory() {
    const container = $("#modelHistory");
    if (!state.history.length) {
      container.innerHTML = '<div class="history-empty">Build your first model to add it here.</div>';
      return;
    }
    container.innerHTML = state.history.map((model) => {
      const active = model.id === state.activeModel?.id ? " active" : "";
      const shortName = model.source?.polymer || model.source?.microstate || "MODEL";
      return `<button class="history-item${active}" data-model-id="${model.id}">
        <div class="history-top"><b>${escapeHtml(model.title)}</b><span>${escapeHtml(shortName)}</span></div>
        <div class="history-meta"><span>${model.atom_count} atoms</span><i></i><span>${escapeHtml(model.formula)}</span><i></i><span>${model.wall_time_seconds.toFixed(2)} s</span></div>
      </button>`;
    }).join("");
    $$(".history-item").forEach((button) => button.addEventListener("click", () => {
      const model = state.history.find((item) => item.id === button.dataset.modelId);
      if (model) showModel(model);
    }));
  }

  async function loadHistory(autoLoad = false) {
    try {
      state.history = await request("/api/models?limit=40");
      renderHistory();
      if (autoLoad && state.history.length) {
        const requested = new URLSearchParams(window.location.search).get("model");
        const selected = state.history.find((item) => item.id === requested) || state.history[0];
        showModel(selected);
      }
    } catch (error) {
      toast(`Model history could not be loaded: ${error.message}`, true);
    }
  }

  function updateFamily(family) {
    $("#family").value = family;
    $$(".segment").forEach((button) => button.classList.toggle("active", button.dataset.family === family));
    $("#smallMoleculeFields").classList.toggle("hidden", family !== "small_molecule");
    $("#polymerFields").classList.toggle("hidden", family !== "polymer");
    updateEngineAvailability();
  }

  function updatePolymerLimits() {
    if (!state.catalog) return;
    const definition = state.catalog.polymers.find((item) => item.id === $("#polymer").value);
    if (!definition) return;
    const slider = $("#repeats");
    slider.max = definition.max_repeats;
    if (Number(slider.value) > definition.max_repeats) slider.value = definition.default_repeats;
    $("#repeatsValue").textContent = slider.value;
    updateEngineAvailability();
  }

  function updateEngineAvailability() {
    const select = $("#engine");
    if (!select || !state.catalog?.engines) return;
    const family = $("#family").value;
    const polymer = $("#polymer").value;
    const modelKind = $("#modelKind").value;
    const statuses = Object.fromEntries(state.catalog.engines.map((item) => [item.id, item]));
    Array.from(select.options).forEach((option) => {
      if (option.value === "rdkit_etkdg") {
        option.disabled = false;
        return;
      }
      const status = statuses[option.value];
      const unsupported = family !== "polymer"
        || status?.status !== "ready"
        || (option.value === "polyply" && (polymer === "PET" || modelKind === "surface_patch"))
        || (option.value === "psp" && polymer === "PET");
      option.disabled = unsupported;
      if (status?.status !== "ready") option.title = status?.reason || "Engine unavailable";
      else if (option.value === "polyply" && polymer === "PET") option.title = "PET force-field block is not validated yet";
      else if (option.value === "polyply" && modelKind === "surface_patch") option.title = "Polyply produces a condensed-phase box, not a surface patch";
      else option.title = "";
    });
    if (select.options[select.selectedIndex]?.disabled) select.value = "rdkit_etkdg";
  }

  function payloadFromForm() {
    const family = $("#family").value;
    return {
      family,
      model_kind: family === "polymer" ? $("#modelKind").value : "molecule",
      pH: Number(family === "polymer" ? $("#polymerPH").value : $("#pH").value),
      temperature_K: Number($("#temperatureK").value || 310.15),
      microstate: $("#microstate").value,
      engine: $("#engine").value,
      polymer: $("#polymer").value,
      repeats: Number($("#repeats").value),
      conformers: Number($("#conformers").value),
      seed: Number($("#seed").value),
      smiles: $("#customSmiles")?.value.trim() || null,
      chembl_id: $("#customChembl")?.value.trim() || null,
      name: $("#customName")?.value.trim() || null,
    };
  }

  async function build(payload) {
    const button = $("#buildButton");
    button.disabled = true;
    $("#viewerLoading").classList.remove("hidden");
    try {
      const model = await request("/api/models", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      state.history = [model, ...state.history.filter((item) => item.id !== model.id)];
      showModel(model);
      toast(`${model.title} was built successfully.`);
    } catch (error) {
      $("#viewerLoading").classList.add("hidden");
      toast(`Model generation failed: ${error.message}`, true);
    } finally {
      button.disabled = false;
    }
  }

  function modelTaskLabel(payload) {
    return payload.family === "polymer"
      ? `${payload.polymer} · ${payload.repeats}-mer · ${payload.engine}`
      : `Dopamine · pH ${payload.pH} · ${payload.microstate}`;
  }

  function renderModelQueue() {
    const holder = $("#modelTaskQueue");
    const run = $("#runModelQueue");
    run.disabled = !state.modelQueue.length;
    holder.innerHTML = state.modelQueue.length
      ? state.modelQueue.map((payload, index) => `<button type="button" data-queue-index="${index}"><b>${escapeHtml(modelTaskLabel(payload))}</b><i title="Remove">×</i></button>`).join("")
      : "<span>No queued model builds.</span>";
    $$("[data-queue-index]").forEach((button) => button.addEventListener("click", () => {
      state.modelQueue.splice(Number(button.dataset.queueIndex), 1);
      renderModelQueue();
    }));
  }

  async function runModelQueue() {
    const run = $("#runModelQueue");
    run.disabled = true;
    while (state.modelQueue.length) {
      const payload = state.modelQueue.shift();
      renderModelQueue();
      await build(payload);
    }
    toast("Model build queue completed.");
  }

  async function runXtbOnActiveModel() {
    if (!state.activeModel) {
      toast("Build or select a model first.", true);
      return;
    }
    const button = $("#runXtbButton");
    button.disabled = true;
    $("#xtbOutput").textContent = "Submitting xTB job…";
    try {
      const result = await request(`/api/models/${state.activeModel.id}/xtb`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ solvent: $("#xtbSolvent").value, optimize: true, timeout_seconds: 300 }),
      });
      state.xtbResults[state.activeModel.id] = result;
      renderXtbResult(result);
      toast(result.status === "completed" ? "GFN2-xTB finished." : result.reason, result.status !== "completed");
    } catch (error) {
      renderXtbResult({ status: "failed", reason: error.message });
      toast(`xTB failed: ${error.message}`, true);
    } finally {
      button.disabled = false;
    }
  }

  function applyPreset(name) {
    const presets = {
      dopamine4: { family: "small_molecule", pH: 4, microstate: "auto", model_kind: "molecule", conformers: 6, seed: 20260826 },
      dopamine9: { family: "small_molecule", pH: 9, microstate: "auto", model_kind: "molecule", conformers: 6, seed: 20260827 },
      pe9: { family: "polymer", polymer: "PE", repeats: 9, model_kind: "molecule", conformers: 5, seed: 20260828 },
      pe15: { family: "polymer", polymer: "PE", repeats: 15, model_kind: "molecule", conformers: 5, seed: 20260829 },
      pet: { family: "polymer", polymer: "PET", repeats: 2, model_kind: "molecule", conformers: 5, seed: 20260829 },
    };
    const preset = presets[name];
    if (!preset) return;
    updateFamily(preset.family);
    if (preset.pH) $("#pH").value = String(preset.pH);
    if (preset.pH) $("#polymerPH").value = String(preset.pH);
    if (preset.microstate) $("#microstate").value = preset.microstate;
    if (preset.polymer) $("#polymer").value = preset.polymer;
    if (preset.repeats) $("#repeats").value = String(preset.repeats);
    $("#conformers").value = String(preset.conformers);
    $("#seed").value = String(preset.seed);
    updatePolymerLimits();
    build({ ...payloadFromForm(), ...preset });
  }

  function bindEvents() {
    $$(".segment").forEach((button) => button.addEventListener("click", () => updateFamily(button.dataset.family)));
    $("#repeats").addEventListener("input", () => $("#repeatsValue").textContent = $("#repeats").value);
    $("#polymer").addEventListener("change", updatePolymerLimits);
    $("#modelKind").addEventListener("change", updateEngineAvailability);
    $("#engine").addEventListener("change", updateEngineAvailability);
    $("#modelForm").addEventListener("submit", (event) => { event.preventDefault(); build(payloadFromForm()); });
    $("#addModelTask").addEventListener("click", () => { state.modelQueue.push(payloadFromForm()); renderModelQueue(); toast("Model setup added to queue."); });
    $("#runModelQueue").addEventListener("click", runModelQueue);
    $$(".preset").forEach((button) => button.addEventListener("click", () => applyPreset(button.dataset.preset)));
    $("#representation").addEventListener("change", applyRepresentation);
    $("#resetViewButton").addEventListener("click", () => state.component?.autoView(400));
    $("#spinButton").addEventListener("click", () => {
      state.spinning = !state.spinning;
      state.stage?.setSpin(state.spinning);
      $("#spinButton").classList.toggle("active", state.spinning);
    });
    $("#refreshHistory").addEventListener("click", () => loadHistory(false));
    const xtbButton = $("#runXtbButton");
    if (xtbButton) xtbButton.addEventListener("click", runXtbOnActiveModel);
  }

  async function initialize() {
    setupStage();
    bindEvents();
    renderModelQueue();
    try {
      const [health, catalog] = await Promise.all([request("/api/health"), request("/api/catalog")]);
      state.catalog = catalog;
      state.xtbStatus = health.xtb || await request("/api/xtb/status");
      $("#apiStatus").classList.add("ready");
      $("#apiStatus span").textContent = `Engine ready · v${health.version}`;
      updateEngineAvailability();
      updatePolymerLimits();
      renderXtbStatus(state.xtbStatus);
      await loadHistory(true);
    } catch (error) {
      $("#apiStatus").classList.add("error");
      $("#apiStatus span").textContent = "Service offline";
      toast(error.message, true);
    }
  }

  initialize();
})();
