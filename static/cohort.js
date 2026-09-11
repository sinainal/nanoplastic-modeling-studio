(() => {
  const $ = (selector) => document.querySelector(selector);
  const polymers = ["PE", "PP", "PS", "PET"];
  const colours = { PE: "#477998", PP: "#658777", PS: "#ad624e", PET: "#76638a" };

  async function request(url) {
    const response = await fetch(url);
    const body = await response.json().catch(() => null);
    if (!response.ok) throw new Error(body?.detail || `${response.status} ${response.statusText}`);
    return body;
  }

  function n(value, digits = 2) {
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric.toFixed(digits) : "—";
  }

  function dataFor(summary, polymer, repeats) {
    return summary.find((item) => item.polymer === polymer && Number(item.repeats) === Number(repeats));
  }

  function renderMetrics(data) {
    const qc = data.quality_control;
    const values = [
      [qc.completed, "completed complexes"],
      [qc.polymers.length, "polymer chemistries"],
      [qc.sizes.length, "oligomer sizes"],
      ["24", "poses / composition"],
      ["Matched", "reference cycle"],
    ];
    $("#cohortMetrics").innerHTML = values.map(([value, label]) => `<article><strong>${value}</strong><span>${label}</span></article>`).join("");
  }

  function renderEnergy(data) {
    const sizes = data.quality_control.sizes;
    const cells = [
      '<div class="energy-header empty"></div>',
      ...sizes.map((size) => `<div class="energy-header">${size}-mer</div>`),
    ];
    polymers.forEach((polymer) => {
      cells.push(`<div class="energy-polymer"><i style="background:${colours[polymer]}"></i>${polymer}</div>`);
      sizes.forEach((size) => {
        const item = dataFor(data.summary, polymer, size);
        if (!item) { cells.push('<div class="energy-cell">—</div>'); return; }
        const intensity = Math.min(1, Math.max(0.10, Math.abs(item.median_kcal_mol) / 5.0));
        cells.push(`<button class="energy-cell" type="button" title="${polymer} ${size}-mer: median ${n(item.median_kcal_mol, 3)}, IQR ${n(item.iqr_kcal_mol, 3)} kcal mol⁻¹" style="--polymer:${colours[polymer]}; --intensity:${intensity}"><b>${n(item.median_kcal_mol, 2)}</b><small>IQR ${n(item.iqr_kcal_mol, 2)}</small></button>`);
      });
    });
    $("#cohortEnergyGrid").innerHTML = cells.join("");
  }

  function renderRanks(data) {
    const rows = data.quality_control.sizes.map((size) => {
      const order = data.ranks[String(size)] || [];
      const confidence = data.rank_confidence[String(size)];
      const labels = order.map((polymer, index) => `<span class="rank-pill rank-${index + 1}" style="--polymer:${colours[polymer]}">${index + 1} ${polymer}</span>`).join("");
      const note = confidence == null
        ? "Separate 9-mer cohort — bootstrap rank pending"
        : `Full-order bootstrap support: ${(confidence * 100).toFixed(1)}%`;
      return `<div class="rank-row"><div><b>${size}-mer</b><small>${note}</small></div><div class="rank-order">${labels}</div></div>`;
    });
    $("#cohortRanks").innerHTML = rows.join("");
  }

  function renderQC(data) {
    const qc = data.quality_control;
    const checks = [
      [qc.completed === qc.expected, "Completed calculation matrix", `${qc.completed}/${qc.expected}`],
      [qc.reference_cycle.includes("matched"), "Matched fragment reference cycle", "Fixed polymer geometry"],
      [true, "Comparable energy metric", "Static Eads only — not ΔGbind"],
      [false, "9-mer rank uncertainty", "Bootstrap analysis pending"],
    ];
    $("#cohortQC").innerHTML = checks.map(([passed, title, detail]) => `<div class="qc-row"><i class="${passed ? "ok" : "review"}">${passed ? "✓" : "!"}</i><div><b>${title}</b><small>${detail}</small></div></div>`).join("");
    const badge = $("#cohortQcBadge");
    badge.className = "validation-pill pending";
    badge.textContent = "Usable with caution";
  }

  function render(data) {
    $("#cohortTitle").textContent = data.title;
    $("#cohortProtocol").textContent = `${data.protocol.adsorbate} · ${data.protocol.medium} · ${data.protocol.method}`;
    renderMetrics(data); renderEnergy(data); renderRanks(data); renderQC(data);
  }

  async function loadCohort(id) {
    try { render(await request(`/api/analysis/cohorts/${encodeURIComponent(id)}`)); }
    catch (error) {
      $("#cohortProtocol").textContent = error.message;
      $("#cohortEnergyGrid").innerHTML = '<div class="chart-empty">Cohort data unavailable.</div>';
    }
  }

  async function init() {
    try {
      const cohorts = await request("/api/analysis/cohorts");
      const select = $("#cohortSelect");
      if (!cohorts.length) throw new Error("No completed analysis cohort is available.");
      select.innerHTML = cohorts.map((item) => `<option value="${item.id}">${item.title} · ${item.completed} completed</option>`).join("");
      select.addEventListener("change", () => loadCohort(select.value));
      await loadCohort(select.value);
    } catch (error) { $("#cohortProtocol").textContent = error.message; }
  }

  init();
})();
