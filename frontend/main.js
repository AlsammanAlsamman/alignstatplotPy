const state = { fastaFile: null, annoFile: null, lastResult: null, activeTab: null };

function setupDropzone(zoneId, inputId, filenameId, onFile) {
  const zone = document.getElementById(zoneId);
  const input = document.getElementById(inputId);
  const filenameEl = document.getElementById(filenameId);

  const setFile = (file) => {
    if (!file) return;
    onFile(file);
    filenameEl.textContent = file.name;
  };

  zone.addEventListener("click", () => input.click());
  input.addEventListener("change", (e) => setFile(e.target.files[0]));

  ["dragenter", "dragover"].forEach((evt) =>
    zone.addEventListener(evt, (e) => {
      e.preventDefault();
      zone.classList.add("dragover");
    })
  );
  ["dragleave", "drop"].forEach((evt) =>
    zone.addEventListener(evt, (e) => {
      e.preventDefault();
      zone.classList.remove("dragover");
    })
  );
  zone.addEventListener("drop", (e) => setFile(e.dataTransfer.files[0]));
}

setupDropzone("dropzone", "fasta-input", "fasta-filename", (f) => (state.fastaFile = f));
setupDropzone("anno-dropzone", "anno-input", "anno-filename", (f) => (state.annoFile = f));

const statusEl = document.getElementById("status");
const runBtn = document.getElementById("run-btn");
const resultsSection = document.getElementById("results");
const tabsEl = document.getElementById("tabs");
const tabContentEl = document.getElementById("tab-content");

function setStatus(message, isError = false) {
  statusEl.textContent = message;
  statusEl.classList.toggle("error", isError);
}

function renderTable(records) {
  if (!records || records.length === 0) return "<p>No data.</p>";
  const columns = Object.keys(records[0]);
  const header = columns.map((c) => `<th>${c}</th>`).join("");
  const rows = records
    .map((row) => `<tr>${columns.map((c) => `<td>${formatCell(row[c])}</td>`).join("")}</tr>`)
    .join("");
  return `<div class="table-wrap"><table><thead><tr>${header}</tr></thead><tbody>${rows}</tbody></table></div>`;
}

function formatCell(value) {
  if (typeof value === "number" && !Number.isInteger(value)) return value.toFixed(4);
  return value;
}

function humanize(name) {
  return name.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase());
}

function showTab(name) {
  state.activeTab = name;
  [...tabsEl.children].forEach((btn) => btn.classList.toggle("active", btn.dataset.tab === name));

  const result = state.lastResult;
  if (name === "stats") {
    tabContentEl.innerHTML = renderTable(result.stats);
  } else if (name === "diversity") {
    tabContentEl.innerHTML = renderTable(result.diversity);
  } else if (name === "distance") {
    tabContentEl.innerHTML = renderTable(result.distance);
  } else if (name === "impact") {
    tabContentEl.innerHTML = renderTable(result.impact);
  } else {
    const b64 = result.plots[name];
    tabContentEl.innerHTML = `<img src="data:image/png;base64,${b64}" alt="${name}" />`;
  }
}

function renderResults(result) {
  state.lastResult = result;
  resultsSection.hidden = false;

  const tabs = [
    ["stats", "Sequence stats"],
    ["diversity", "Diversity"],
    ["distance", "Distance matrix"],
  ];
  if (result.impact && result.impact.length) tabs.push(["impact", "SNP impact"]);
  Object.keys(result.plots).forEach((name) => tabs.push([name, humanize(name)]));

  tabsEl.innerHTML = "";
  tabs.forEach(([id, label], i) => {
    const btn = document.createElement("button");
    btn.className = "tab-btn";
    btn.dataset.tab = id;
    btn.textContent = label;
    btn.addEventListener("click", () => showTab(id));
    tabsEl.appendChild(btn);
    if (i === 0) showTab(id);
  });
}

runBtn.addEventListener("click", async () => {
  if (!state.fastaFile) {
    setStatus("Please choose a FASTA file first.", true);
    return;
  }

  const form = new FormData();
  form.append("fasta", state.fastaFile);
  if (state.annoFile) form.append("annotation", state.annoFile);
  form.append("max_miss_per", document.getElementById("max-miss").value);
  form.append("n_clusters", document.getElementById("n-clusters").value);
  form.append("min_cluster_length", document.getElementById("min-cluster-length").value);
  form.append("variant_window", document.getElementById("variant-window").value);

  runBtn.disabled = true;
  setStatus("Running analysis...");

  try {
    const response = await fetch("/api/run", { method: "POST", body: form });
    const data = await response.json();
    if (!response.ok) throw new Error(data.error || "Analysis failed.");
    setStatus("Done.");
    renderResults(data);
  } catch (err) {
    setStatus(err.message, true);
  } finally {
    runBtn.disabled = false;
  }
});
