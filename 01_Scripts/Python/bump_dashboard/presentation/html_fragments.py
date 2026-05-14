"""
presentation.html_fragments
============================

Static HTML/CSS/JS fragments assembled by ``DashboardHtmlRenderer``.

Each fragment is a plain string constant.  No data embedding happens here —
that is the renderer's responsibility.

Structure mirrors the rendered page:
  ``HEAD``    → ``<head>`` tag  (CSS design-tokens, layout, component styles)
  ``BODY``    → ``<body>`` tag  (sidebar controls + main chart area)
  ``SCRIPT``  → ``<script>`` tag (Plotly initialisation + all JS logic)

Token substitution uses unambiguous sentinel strings:
  ``%%PLOTLY_JS_URL%%``   – CDN URL for Plotly
  ``%%PATHWAYS_JSON%%``   – JSON array of pathway objects
  ``%%METADATA_JSON%%``   – JSON metadata object

These sentinels are chosen to be 100% safe inside an HTML raw-string context
(double-percent prefix is never valid JS) and are easily greppable.

Audit fixes applied (2026-05)
------------------------------
Bug 1  – Tooltip clipping at top of viewport: replaced ``bottom:125%`` with
         smart ``_positionTooltip()`` JS that flips to ``top:`` when the
         tooltip would overflow above the viewport.
Bug 2  – "Display Mode" had no ? help icon: added one with explanatory text.
Bug 3  – "Visual Style" tooltip did not explain Bézier amplitude origin:
         tooltip now explains that curvature magnitude = |TrajDev NES| and
         direction (up/down) = sign of TrajDev NES, and only sig. pathways
         are curved.
Bug 4/5 – Status-bar "Showing N pathways" did not react to significance
         toggles: root cause was that ``filteredData`` was set by the global
         scope filter before per-mutation visibility predicates were applied.
         The status-bar now counts the union of pathways visible in *either*
         mutation after all filters are applied.
Bug 6  – Single ``<input type="number">`` replaced by dual-handle range
         sliders for: Min |NES|, |Early NES|, |Late NES|, |TrajDev NES|.
         Each slider shows live "X – Y" readout and filters inclusively.
Bug 7  – Highlight search: (a) now uses the active Bézier path when curved
         mode is on; (b) draws highlighted traces on top with high z-order,
         bold colour, and reduced opacity on all non-matching lines so
         matches are clearly visible even in dense charts.

Architecture notes
------------------
* JS is organised into named IIFE-free sections mirroring the Python layer:
  CONFIG → STATE → INIT → CONTROLS → FILTERS → GEOMETRY → COLOURS →
  TOOLTIPS → TRACES → RENDER → UPDATE.
* ``_readSettings()`` is the single source of truth for UI state per cycle.
* All filter predicates are pure functions; ``updateCharts()`` is the sole
  orchestrator that mutates global STATE.
"""

# ---------------------------------------------------------------------------
# HEAD — styles
# ---------------------------------------------------------------------------

HEAD: str = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Interactive Bump Chart Dashboard</title>
  <script src="%%PLOTLY_JS_URL%%"></script>
  <style>
    /* ===== Design tokens ===== */
    :root {
      --primary:        #2c3e50;
      --accent:         #3498db;
      --accent-hover:   #2980b9;
      --bg:             #f8f9fa;
      --sidebar-width:  360px;
      --border:         #ddd;
      --text-muted:     #666;
      --radius:         4px;
      --shadow-sm:      0 2px 4px rgba(0,0,0,.05);
    }

    /* ===== Reset / base ===== */
    *, *::before, *::after { box-sizing: border-box; }
    body {
      margin: 0; padding: 0;
      font-family: -apple-system, system-ui, sans-serif;
      font-size: 14px;
      background: var(--bg);
      display: flex;
      height: 100vh;
      overflow: hidden;
    }

    /* ===== Sidebar ===== */
    .sidebar {
      width: var(--sidebar-width);
      flex-shrink: 0;
      background: white;
      border-right: 1px solid var(--border);
      padding: 16px;
      overflow-y: auto;
      display: flex;
      flex-direction: column;
      gap: 14px;
      box-shadow: 2px 0 5px rgba(0,0,0,.05);
    }

    /* ===== Control groups ===== */
    .ctrl-group {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    .ctrl-group label.group-label {
      font-weight: 600;
      color: var(--primary);
      display: flex;
      align-items: center;
      justify-content: space-between;
    }

    /* ===== Radio / checkbox rows ===== */
    .check-item {
      display: flex;
      align-items: center;
      gap: 6px;
      padding: 2px 4px;
      border-radius: var(--radius);
      cursor: pointer;
    }
    .check-item:hover { background: #f0f8ff; }

    /* ===== Scrollable checkbox list ===== */
    .check-list {
      display: flex;
      flex-direction: column;
      gap: 3px;
      max-height: 150px;
      overflow-y: auto;
      border: 1px solid #eee;
      padding: 5px;
      border-radius: var(--radius);
    }

    /* ===== Inputs ===== */
    input[type="text"],
    select {
      padding: 7px 8px;
      border: 1px solid var(--border);
      border-radius: var(--radius);
      width: 100%;
    }

    /* ===== Dual-handle range slider ===== */
    .range-wrap {
      display: flex;
      flex-direction: column;
      gap: 3px;
    }
    .range-track {
      position: relative;
      height: 28px;
      display: flex;
      align-items: center;
    }
    /* Both range inputs stacked; the lower one is inverted */
    .range-track input[type="range"] {
      position: absolute;
      width: 100%;
      height: 4px;
      appearance: none;
      -webkit-appearance: none;
      background: transparent;
      pointer-events: none;
    }
    .range-track input[type="range"]::-webkit-slider-thumb {
      -webkit-appearance: none;
      pointer-events: all;
      width: 16px; height: 16px;
      border-radius: 50%;
      background: var(--accent);
      border: 2px solid white;
      box-shadow: 0 1px 3px rgba(0,0,0,.25);
      cursor: pointer;
    }
    .range-track input[type="range"]::-moz-range-thumb {
      pointer-events: all;
      width: 14px; height: 14px;
      border-radius: 50%;
      background: var(--accent);
      border: 2px solid white;
      cursor: pointer;
    }
    .range-fill {
      position: absolute;
      height: 4px;
      background: var(--accent);
      border-radius: 2px;
      pointer-events: none;
    }
    .range-track-bg {
      position: absolute;
      width: 100%;
      height: 4px;
      background: #ddd;
      border-radius: 2px;
    }
    .range-readout {
      font-size: 0.8em;
      color: var(--text-muted);
      text-align: right;
    }

    /* ===== Button ===== */
    .btn {
      padding: 7px 12px;
      background: var(--accent);
      color: white;
      border: none;
      border-radius: var(--radius);
      cursor: pointer;
    }
    .btn:hover { background: var(--accent-hover); }
    .btn-sm { padding: 3px 8px; font-size: 0.82em; }

    /* ===== Tooltip (Bug 1 fix: positioning is handled by JS) ===== */
    .info-wrap {
      position: relative;
      display: inline-block;
    }
    .info-icon {
      cursor: help;
      color: #999;
      border: 1px solid #ccc;
      border-radius: 50%;
      width: 16px;
      height: 16px;
      font-size: 11px;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      margin-left: 4px;
      background: white;
      user-select: none;
    }
    .info-icon:hover { color: var(--accent); border-color: var(--accent); }
    /*
     * Tooltip box:  default position is ABOVE the icon (bottom:130%).
     * When JS detects viewport overflow the class .tip-below is added,
     * which flips it to BELOW (top:130%, bottom:auto).
     * This fixes Bug 1: tooltip no longer clips at the top of the window.
     */
    .tooltip-text {
      display: none;              /* hidden by default */
      width: 240px;
      background: #222;
      color: #f0f0f0;
      font-size: 0.83em;
      font-weight: normal;
      line-height: 1.5;
      text-align: left;
      border-radius: 6px;
      padding: 9px 10px;
      position: fixed;            /* fixed so it escapes any overflow:hidden parent */
      z-index: 9999;
      box-shadow: 0 3px 10px rgba(0,0,0,.35);
      pointer-events: none;
    }

    /* ===== Divider ===== */
    hr.divider {
      border: 0;
      border-top: 1px solid #eee;
      margin: 2px 0;
      width: 100%;
    }

    /* ===== Filter label row ===== */
    .filter-label-row {
      display: flex;
      align-items: center;
      gap: 4px;
      font-size: 0.9em;
      font-weight: 500;
      color: #444;
      margin-bottom: 2px;
    }

    /* ===== Main content ===== */
    .main {
      flex: 1;
      display: flex;
      flex-direction: column;
      padding: 16px;
      gap: 10px;
      min-width: 0;
    }

    .page-header {
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding-bottom: 10px;
      border-bottom: 1px solid var(--border);
    }
    .page-title { font-size: 1.15em; font-weight: bold; color: var(--primary); }
    .status-bar { font-size: 0.88em; color: var(--text-muted); }

    /* ===== Charts ===== */
    .charts-row {
      flex: 1;
      display: flex;
      gap: 10px;
      min-height: 0;
    }
    .chart-card {
      flex: 1;
      background: white;
      border-radius: 8px;
      box-shadow: var(--shadow-sm);
      display: flex;
      flex-direction: column;
      min-width: 0;
    }
    .chart-card-header {
      padding: 10px;
      font-weight: bold;
      text-align: center;
      border-bottom: 1px solid #eee;
    }
    .chart-plotly {
      flex: 1;
      min-height: 0;
      width: 100%;
    }

    /* ===== GSVA replicate-level click-through panel ===== */
    .gsva-modal-overlay {
      display: none;
      position: fixed;
      inset: 0;
      background: rgba(0,0,0,.32);
      z-index: 10000;
      align-items: center;
      justify-content: center;
    }
    .gsva-modal-overlay.is-open { display: flex; }
    .gsva-modal {
      background: white;
      border-radius: 8px;
      box-shadow: 0 12px 32px rgba(0,0,0,.25);
      width: min(880px, 92vw);
      max-height: 88vh;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .gsva-modal-header {
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 12px;
      padding: 14px 18px;
      border-bottom: 1px solid #eee;
    }
    .gsva-modal-title {
      font-size: 1.02em;
      font-weight: 600;
      color: var(--primary);
      line-height: 1.3;
    }
    .gsva-modal-subtitle {
      font-size: 0.83em;
      color: var(--text-muted);
      margin-top: 3px;
    }
    .gsva-modal-pattern-row {
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
      margin-top: 6px;
      font-size: 0.84em;
    }
    .gsva-pattern-chip {
      padding: 2px 8px;
      border-radius: 10px;
      background: #f0f0f0;
      color: #222;
      font-weight: 500;
    }
    .gsva-modal-close {
      background: transparent;
      border: none;
      cursor: pointer;
      font-size: 1.4em;
      color: #888;
      padding: 0 6px;
      line-height: 1;
    }
    .gsva-modal-close:hover { color: #222; }
    .gsva-modal-body {
      padding: 12px 18px 18px 18px;
      overflow-y: auto;
    }
    .gsva-plot {
      width: 100%;
      height: 360px;
    }
    .gsva-caption {
      font-size: 0.78em;
      color: var(--text-muted);
      margin-top: 8px;
      line-height: 1.45;
    }
    /* GSVA driver verdict block — appears between the plot and the caption */
    .gsva-driver-verdict {
      font-size: 0.82em;
      margin: 10px 0 6px 0;
      line-height: 1.7;
      border-top: 1px solid var(--border);
      padding-top: 8px;
    }
    .gsva-driver-verdict .verdict-deltas { color: var(--text-muted); margin-bottom: 3px; }
    .gsva-driver-verdict .verdict-line   { font-weight: 600; }
    .gsva-empty {
      padding: 28px 14px;
      text-align: center;
      color: var(--text-muted);
      font-size: 0.92em;
    }
  </style>
</head>
"""

# ---------------------------------------------------------------------------
# BODY — sidebar + main layout
# ---------------------------------------------------------------------------

BODY: str = """\
<body>

<!-- Floating tooltip element (singleton, repositioned by JS) -->
<div id="global-tooltip" class="tooltip-text" role="tooltip" aria-hidden="true"></div>

<!-- ═══════════════════════════════════ SIDEBAR ═══════════════════════════ -->
<nav class="sidebar" aria-label="Dashboard controls">

  <!-- ── Display Mode ─────────────────────────────────────── -->
  <!--
    Bug 2 fix: Display Mode now has a ? info icon.
  -->
  <div class="ctrl-group">
    <label class="group-label">
      Display Mode
      <span class="info-wrap"
            data-tip="<b>Uniform:</b> every pathway line has the same thickness and opacity — easier to spot individual paths in dense charts.<br><br><b>Weighted:</b> line thickness and opacity reflect how <em>rare</em> the pattern is across the dataset. Dominant patterns (most frequent) are thin and faint; rare patterns are thick and fully opaque. This lets uncommon biology stand out against the background crowd.">
        <span class="info-icon" aria-label="info about Display Mode">?</span>
      </span>
    </label>
    <div class="check-item">
      <input type="radio" name="mode" value="uniform" onchange="updateCharts()">
      Uniform (all lines same width)
    </div>
    <div class="check-item">
      <input type="radio" name="mode" value="weighted" checked onchange="updateCharts()">
      Weighted (rarer pattern = thicker / more opaque)
    </div>
  </div>

  <!-- ── Y-Axis Metric ─────────────────────────────────────── -->
  <div class="ctrl-group">
    <label class="group-label">
      Y-Axis Metric
      <span class="info-wrap"
            data-tip="<b>NES (Normalised Enrichment Score):</b> raw effect size and direction. Positive = gene-set up-regulated vs control; negative = down-regulated. Scale is shared across both charts so amplitudes are comparable.<br><br><b>Rank:</b> replaces the raw NES with the pathway's ordinal position (1 = highest NES). Useful for comparing <em>relative ordering</em> across time-points without being distracted by absolute magnitudes.">
        <span class="info-icon" aria-label="info about Y-Axis Metric">?</span>
      </span>
    </label>
    <div class="check-item">
      <input type="radio" name="ytype" value="nes" checked onchange="updateCharts()">
      NES (effect size)
    </div>
    <div class="check-item">
      <input type="radio" name="ytype" value="rank" onchange="updateCharts()">
      Rank (relative position)
    </div>
  </div>

  <!-- ── Visual Style ──────────────────────────────────────── -->
  <!--
    Bug 3 fix: tooltip now fully explains Bézier amplitude source.
  -->
  <div class="ctrl-group">
    <label class="group-label">
      Visual Style
      <span class="info-wrap"
            data-tip="<b>Straight lines</b> connect each pathway's Early (D35) and Late (D65) NES values directly.<br><br><b>Curved lines</b> add a Bézier arc whose properties encode the Trajectory Deviation (TrajDev) contrast:<br>• <em>Curvature magnitude</em> = |TrajDev NES| — a larger absolute NES produces a more pronounced bow.<br>• <em>Curvature direction</em> = sign of TrajDev NES — positive TrajDev bows the line upward; negative bows it downward.<br>• Curves are drawn <em>only</em> for pathways with significant TrajDev (p.adj &lt; 0.05); non-significant pathways remain straight even when this option is on.">
        <span class="info-icon" aria-label="info about Visual Style">?</span>
      </span>
    </label>
    <div class="check-item">
      <input type="checkbox" id="chk-curved" onchange="updateCharts()">
      Curved lines (encode trajectory deviation)
    </div>
  </div>

  <!-- ── Color By ──────────────────────────────────────────── -->
  <div class="ctrl-group">
    <label class="group-label">
      Color By
      <span class="info-wrap"
            data-tip="<b>Pattern:</b> each trajectory-classification pattern (Compensation, Progressive, etc.) has its own project-canonical colour.<br><br><b>NES options:</b> apply a Blue–White–Orange diverging scale to the selected NES value. Blue = strongly down-regulated; white = near zero; orange = strongly up-regulated. Useful for spotting directional gradients across the chart.">
        <span class="info-icon" aria-label="info about Color By">?</span>
      </span>
    </label>
    <select id="color-by" onchange="updateCharts()">
      <option value="pattern" selected>Pattern (default)</option>
      <option value="nes_early">Early NES</option>
      <option value="nes_late">Late NES</option>
      <option value="nes_trajdev">TrajDev NES</option>
    </select>
  </div>

  <hr class="divider">

  <!-- ── Significance quick-filters ───────────────────────── -->
  <!--
    Bugs 4 & 5 fix: the status-bar count now reflects per-mutation visibility
    (see JS updateCharts); these checkboxes correctly drive that count.
  -->
  <div class="ctrl-group">
    <label class="group-label">
      Significance Filters
      <span class="info-wrap"
            data-tip="<b>TrajDev only:</b> keep pathways where the Maturation-specific (TrajDev) contrast is significant. These are the pathways whose trajectory differs between mutant and control during maturation.">
        <span class="info-icon" aria-label="info about significance filters">?</span>
      </span>
    </label>
    <div class="check-item">
      <input type="checkbox" id="filter-sig-traj" onchange="updateCharts()">
      Significant TrajDev only
    </div>
  </div>

  <!-- ── NES range sliders ─────────────────────────────────── -->
  <!--
    Bug 6 fix: single-number inputs replaced by dual-handle range sliders.
    Each slider pair defines a min–max window; pathways must fall inside
    the window to be displayed (for that mutation).
  -->
  <div class="ctrl-group">
    <label class="group-label">
      NES Range Filters
      <span class="info-wrap"
            data-tip="Each slider defines an <em>inclusive</em> min–max window. A pathway is shown only when its |NES| for that contrast falls inside the window.<br><br>Drag the <b>left knob</b> to set the minimum; drag the <b>right knob</b> to set the maximum. This lets you isolate pathways in a specific effect-size band — for example, moderate effects (|NES| 1.0–2.0) or very strong effects (|NES| ≥ 2.5).">
        <span class="info-icon" aria-label="info about NES range filters">?</span>
      </span>
    </label>

    <!-- Min |NES| (any contrast) -->
    <div class="filter-label-row">Min |NES| (any contrast)</div>
    <div class="range-wrap" id="rw-min-nes">
      <div class="range-track">
        <div class="range-track-bg"></div>
        <div class="range-fill" id="rf-min-nes"></div>
        <input type="range" id="rs-min-nes-lo" min="0" max="4" step="0.05" value="0"
               oninput="syncSlider('min-nes')">
        <input type="range" id="rs-min-nes-hi" min="0" max="4" step="0.05" value="4"
               oninput="syncSlider('min-nes')">
      </div>
      <div class="range-readout" id="rr-min-nes">0 – 4</div>
    </div>

    <details style="margin-top:6px;">
      <summary style="cursor:pointer;font-weight:500;color:#555;font-size:0.92em;">
        Per-contrast NES windows
      </summary>
      <div style="padding:8px 0 0 0;display:flex;flex-direction:column;gap:10px;">

        <!-- |Early NES| -->
        <div>
          <div class="filter-label-row">|Early NES| range</div>
          <div class="range-wrap" id="rw-early-nes">
            <div class="range-track">
              <div class="range-track-bg"></div>
              <div class="range-fill" id="rf-early-nes"></div>
              <input type="range" id="rs-early-nes-lo" min="0" max="4" step="0.05" value="0"
                     oninput="syncSlider('early-nes')">
              <input type="range" id="rs-early-nes-hi" min="0" max="4" step="0.05" value="4"
                     oninput="syncSlider('early-nes')">
            </div>
            <div class="range-readout" id="rr-early-nes">0 – 4</div>
          </div>
        </div>

        <!-- |Late NES| -->
        <div>
          <div class="filter-label-row">|Late NES| range</div>
          <div class="range-wrap" id="rw-late-nes">
            <div class="range-track">
              <div class="range-track-bg"></div>
              <div class="range-fill" id="rf-late-nes"></div>
              <input type="range" id="rs-late-nes-lo" min="0" max="4" step="0.05" value="0"
                     oninput="syncSlider('late-nes')">
              <input type="range" id="rs-late-nes-hi" min="0" max="4" step="0.05" value="4"
                     oninput="syncSlider('late-nes')">
            </div>
            <div class="range-readout" id="rr-late-nes">0 – 4</div>
          </div>
        </div>

        <!-- |TrajDev NES| -->
        <div>
          <div class="filter-label-row">|TrajDev NES| range</div>
          <div class="range-wrap" id="rw-trajdev-nes">
            <div class="range-track">
              <div class="range-track-bg"></div>
              <div class="range-fill" id="rf-trajdev-nes"></div>
              <input type="range" id="rs-trajdev-nes-lo" min="0" max="4" step="0.05" value="0"
                     oninput="syncSlider('trajdev-nes')">
              <input type="range" id="rs-trajdev-nes-hi" min="0" max="4" step="0.05" value="4"
                     oninput="syncSlider('trajdev-nes')">
            </div>
            <div class="range-readout" id="rr-trajdev-nes">0 – 4</div>
          </div>
        </div>

        <!-- Stage-specific significance -->
        <div class="check-item">
          <input type="checkbox" id="filter-sig-early" onchange="updateCharts()">
          Sig. Early (p.adj &lt; 0.05)
        </div>
        <div class="check-item">
          <input type="checkbox" id="filter-sig-late" onchange="updateCharts()">
          Sig. Late (p.adj &lt; 0.05)
        </div>

      </div>
    </details>
  </div>

  <hr class="divider">

  <!-- ── Pattern filter ────────────────────────────────────── -->
  <div class="ctrl-group">
    <label class="group-label">Filter: Patterns</label>
    <div class="check-list" id="pattern-list" role="group"></div>
    <button class="btn btn-sm" style="margin-top:4px;" onclick="toggleAll('pattern-list')">Toggle all</button>
  </div>

  <!-- ── GSVA driver filter ────────────────────────────────── -->
  <div class="ctrl-group">
    <label class="group-label">
      Filter: GSVA driver
      <span class="info-wrap"
            data-tip="The GSEA pattern label tells you the two trajectories differ; the GSVA driver tells you <b>which arm is moving</b>. Computed per (pathway, mutation) from per-sample GSVA medians (D65−D35) using ε = 0.10.<br><br><b>mutant_driven:</b> only the mutant arm changes substantially.<br><b>ctrl_driven:</b> only the Ctrl arm changes — the GSEA contrast is real, but it's the Ctrl developmental trajectory carrying it (apparent rescue / closing gap).<br><b>both_moving:</b> both arms move (often in opposite directions — true crossover).<br><b>neither_moving:</b> both arms stay flat (rare; GSEA contrast may be noise or pattern-shape effect).">
        <span class="info-icon" aria-label="info about GSVA driver filter">?</span>
      </span>
    </label>
    <div class="check-list" id="driver-list" role="group"></div>
    <button class="btn btn-sm" style="margin-top:4px;" onclick="toggleAll('driver-list')">Toggle all</button>
  </div>

  <!-- ── Database filter ───────────────────────────────────── -->
  <div class="ctrl-group">
    <label class="group-label">Filter: Databases</label>
    <div class="check-list" id="db-list" role="group"></div>
    <button class="btn btn-sm" style="margin-top:4px;" onclick="toggleAll('db-list')">Toggle all</button>
  </div>

  <!-- ── Highlight search ──────────────────────────────────── -->
  <!--
    Bug 7 fix: highlight now uses the active Bézier geometry; non-matching
    lines are dimmed so the highlight stands out in dense charts.
  -->
  <div class="ctrl-group">
    <label class="group-label">
      Highlight Pathway
      <span class="info-wrap"
            data-tip="Type any fragment of a pathway description. Matching pathways are drawn on top in bold black (using the same straight or curved geometry as the rest of the chart). All other pathways are simultaneously dimmed so matches are easily spotted even in dense views. Clear the box to restore normal colours.">
        <span class="info-icon" aria-label="info about highlight search">?</span>
      </span>
    </label>
    <input type="text" id="search-box" placeholder="Search description…" oninput="updateHighlights()">
  </div>

  <!-- ── Reset ─────────────────────────────────────────────── -->
  <div class="ctrl-group" style="margin-top:auto;">
    <button class="btn" onclick="resetView()">Reset view</button>
  </div>

</nav><!-- end sidebar -->


<!-- ═══════════════════════════════════ MAIN ══════════════════════════════ -->
<main class="main">
  <header class="page-header">
    <div class="page-title">Pathway Trajectory Dashboard</div>
    <div class="status-bar" id="status-bar" aria-live="polite">Loading…</div>
  </header>

  <div class="charts-row">
    <div class="chart-card">
      <div class="chart-card-header" style="color:#0072B2;">G32A Mutation</div>
      <div id="chart-g32a" class="chart-plotly"></div>
    </div>
    <div class="chart-card">
      <div class="chart-card-header" style="color:#D55E00;">R403C Mutation</div>
      <div id="chart-r403c" class="chart-plotly"></div>
    </div>
  </div>
</main>

<!-- ═══════════════════════ GSVA REPLICATE-LEVEL PANEL ═══════════════════════ -->
<!--
  Phase 2b: Click any pathway line to see its per-sample GSVA enrichment
  scores grouped by genotype × day. Jittered dots are individual replicates;
  bars are group means ± approximate 95% CI (mean ± 1.96 × SE, small n).
-->
<div class="gsva-modal-overlay" id="gsva-modal" role="dialog" aria-modal="true"
     aria-labelledby="gsva-modal-title" onclick="if(event.target===this) closeGsvaModal()">
  <div class="gsva-modal">
    <div class="gsva-modal-header">
      <div>
        <div id="gsva-modal-title" class="gsva-modal-title">Pathway</div>
        <div id="gsva-modal-subtitle" class="gsva-modal-subtitle"></div>
        <div id="gsva-modal-patterns" class="gsva-modal-pattern-row"></div>
      </div>
      <button class="gsva-modal-close" type="button" aria-label="Close panel"
              onclick="closeGsvaModal()">×</button>
    </div>
    <div class="gsva-modal-body">
      <div id="gsva-plot" class="gsva-plot"></div>
      <div id="gsva-driver-verdict" class="gsva-driver-verdict"></div>
      <div id="gsva-caption" class="gsva-caption"></div>
    </div>
  </div>
</div>
"""

# ---------------------------------------------------------------------------
# SCRIPT — full JS application
# ---------------------------------------------------------------------------
# Internal section order (each section is clearly delimited):
#   §CONFIG   – injected constants + style tables
#   §STATE    – mutable runtime state (minimal surface)
#   §INIT     – bootstrap and DOM wiring
#   §SLIDER   – dual-handle range-slider mechanics  (Bug 6)
#   §TOOLTIP  – viewport-aware tooltip positioning  (Bug 1)
#   §CONTROLS – UI event handlers (reset, toggleAll)
#   §SETTINGS – _readSettings() single-source-of-truth reader
#   §FILTERS  – pure predicate functions (no side-effects)
#   §GEOMETRY – Bézier + text-wrap helpers
#   §COLOUR   – NES diverging colour scale helpers
#   §TRACES   – pattern / NES-colour / highlight trace builders (Bug 7)
#   §RENDER   – Plotly chart assembly
#   §UPDATE   – updateCharts() orchestrator (Bugs 4 & 5 fix)

SCRIPT: str = r"""
<script>
"use strict";

// ╔══════════════════════════════════════════════════════════════╗
// ║  §CONFIG — constants injected by Python renderer            ║
// ╚══════════════════════════════════════════════════════════════╝

const RAW_DATA = %%PATHWAYS_JSON%%;
const METADATA  = %%METADATA_JSON%%;

const PATTERN_COLORS = METADATA.pattern_colors;
const WEIGHT_CATS    = METADATA.weight_categories;
const PATTERN_DEFS   = METADATA.pattern_definitions;

/**
 * Line style per weight category.
 *
 * "dominant" patterns are the most common and rendered thin + translucent
 * so they form a background crowd. "rare" patterns are thick + opaque
 * so uncommon biology always stands out on top.
 */
const WEIGHT_STYLES = {
  dominant:  { width: 0.8, opacity: 0.22 },
  common:    { width: 1.5, opacity: 0.45 },
  uncommon:  { width: 2.5, opacity: 0.75 },
  rare:      { width: 3.5, opacity: 1.00 },
};
const UNIFORM_STYLE = { width: 1.2, opacity: 0.45 };

/** Dimmed style applied to non-matching lines when a highlight search is active. */
const DIMMED_STYLE = { width: 0.6, opacity: 0.10 };

/** Diverging colour palette — mirrors Python domain.rules */
const CLR_NEG  = '#2166AC';   // blue  (negative NES)
const CLR_MID  = '#F7F7F7';   // white (zero NES)
const CLR_POS  = '#B35806';   // orange (positive NES)
const NES_CMAX = 3.5;

/** Bézier control-point offset factor in NES mode (fraction of |TrajDev NES|). */
const BEZIER_NES_SCALE  = 0.55;
/** Bézier offset factor in Rank mode (fraction of visible rank range). */
const BEZIER_RANK_SCALE = 0.12;

// ╔══════════════════════════════════════════════════════════════╗
// ║  §STATE — mutable runtime (minimal surface area)           ║
// ╚══════════════════════════════════════════════════════════════╝

/** Pathways passing the global (database + pattern) scope filter. */
let scopeData = [];
/** Pathways visible in G32A after all per-mutation filters. */
let visibleG32A = [];
/** Pathways visible in R403C after all per-mutation filters. */
let visibleR403C = [];

let highlightSearch = '';
let globalYMin = 0;
let globalYMax = 0;

// ╔══════════════════════════════════════════════════════════════╗
// ║  §INIT — bootstrap                                          ║
// ╚══════════════════════════════════════════════════════════════╝

// Fixed display order for GSVA driver labels (matches the decision-table in the spec).
const GSVA_DRIVER_ORDER = ['mutant_driven', 'ctrl_driven', 'both_moving', 'neither_moving'];

// Human-readable expansions used in the modal verdict and filter list.
const GSVA_DRIVER_LABELS = {
  mutant_driven:   'mutant-driven (mutant carries the change)',
  ctrl_driven:     'ctrl-driven (Ctrl developmental trajectory closes the gap)',
  both_moving:     'both moving (both arms change substantially — often a true crossover)',
  neither_moving:  'neither moving (both arms flat — contrast may reflect shape/noise)',
};

function init() {
  _populateCheckList('pattern-list', METADATA.patterns, true);
  _populateCheckList('db-list',      METADATA.databases, false);
  _populateDriverList();
  _initAllSliders();
  _initTooltips();
  updateCharts();
}

/**
 * Populate #driver-list with checkboxes in the canonical GSVA_DRIVER_ORDER.
 * Only labels that actually appear in the data (across both mutations) are shown;
 * always start with all boxes checked.
 */
function _populateDriverList() {
  const present = new Set();
  RAW_DATA.forEach(d => {
    if (d.Driver_G32A  != null) present.add(d.Driver_G32A);
    if (d.Driver_R403C != null) present.add(d.Driver_R403C);
  });
  const el = document.getElementById('driver-list');
  el.innerHTML = '';
  GSVA_DRIVER_ORDER.forEach(key => {
    if (!present.has(key)) return;        // skip labels absent from this build
    const div = document.createElement('div');
    div.className = 'check-item';
    const shortLabel = key.replace(/_/g, '\u202f');   // thin-space for readability
    div.innerHTML =
      `<input type="checkbox" value="${_esc(key)}" checked onchange="updateCharts()"> ${_esc(shortLabel)}`;
    el.appendChild(div);
  });
  // Always add a "(no data)" bucket so pathways without GSVA can be included.
  const noData = document.createElement('div');
  noData.className = 'check-item';
  noData.innerHTML =
    `<input type="checkbox" value="__null__" checked onchange="updateCharts()"> (no data)`;
  el.appendChild(noData);
}

function _populateCheckList(id, items, showColor) {
  const el = document.getElementById(id);
  el.innerHTML = '';
  items.forEach(item => {
    const div = document.createElement('div');
    div.className = 'check-item';
    let inner = `<input type="checkbox" value="${_esc(item)}" checked onchange="updateCharts()"> ${_esc(item)}`;
    if (showColor && PATTERN_COLORS[item]) {
      inner += `<span style="width:10px;height:10px;background:${_esc(PATTERN_COLORS[item])};
                             border-radius:50%;display:inline-block;margin-left:4px;flex-shrink:0;"></span>`;
    }
    div.innerHTML = inner;
    el.appendChild(div);
  });
}

function _esc(s) {
  // Minimal HTML-entity escape for dynamic content inserted into innerHTML
  return String(s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')
    .replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §SLIDER — dual-handle range slider logic  (Bug 6)          ║
// ╚══════════════════════════════════════════════════════════════╝

/**
 * Each logical slider is identified by a key (e.g. "min-nes").
 * Elements expected:
 *   #rs-{key}-lo   lower range input
 *   #rs-{key}-hi   upper range input
 *   #rf-{key}      fill bar div
 *   #rr-{key}      readout span
 */
const SLIDER_KEYS = ['min-nes', 'early-nes', 'late-nes', 'trajdev-nes'];

function _initAllSliders() {
  SLIDER_KEYS.forEach(key => syncSlider(key, /*initial*/ true));
}

/**
 * Called oninput on either thumb.  Enforces lo ≤ hi, updates the fill bar
 * and the readout label, then triggers a chart update (unless initialising).
 */
function syncSlider(key, isInit) {
  const lo = document.getElementById(`rs-${key}-lo`);
  const hi = document.getElementById(`rs-${key}-hi`);
  if (!lo || !hi) return;

  let vLo = parseFloat(lo.value);
  let vHi = parseFloat(hi.value);

  // Enforce ordering
  if (vLo > vHi) {
    // Figure out which thumb was moved by comparing to last known state
    // Simple heuristic: clamp the one that would cross the other
    if (document.activeElement === lo) { lo.value = vHi; vLo = vHi; }
    else                               { hi.value = vLo; vHi = vLo; }
  }

  // Update fill bar
  const min    = parseFloat(lo.min);
  const max    = parseFloat(lo.max);
  const range  = max - min;
  const pLo    = range > 0 ? (vLo - min) / range * 100 : 0;
  const pHi    = range > 0 ? (vHi - min) / range * 100 : 100;
  const fill   = document.getElementById(`rf-${key}`);
  if (fill) {
    fill.style.left  = pLo + '%';
    fill.style.width = (pHi - pLo) + '%';
  }

  // Update readout
  const rr = document.getElementById(`rr-${key}`);
  if (rr) rr.textContent = `${vLo.toFixed(2)} – ${vHi.toFixed(2)}`;

  if (!isInit) updateCharts();
}

/** Read a slider pair and return { lo, hi }. */
function _sliderValues(key) {
  return {
    lo: parseFloat(document.getElementById(`rs-${key}-lo`).value) || 0,
    hi: parseFloat(document.getElementById(`rs-${key}-hi`).value) || 4,
  };
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §TOOLTIP — viewport-aware positioning  (Bug 1)             ║
// ╚══════════════════════════════════════════════════════════════╝

/**
 * We use a single floating <div id="global-tooltip"> that is moved to the
 * hovered icon position and shown/hidden via JS.  The tooltip is
 * position:fixed so it escapes any overflow:hidden ancestor (sidebar scroll).
 *
 * Direction logic:
 *   - Preferred: appear ABOVE the icon (tooltip bottom edge = icon top - 8px).
 *   - If that would clip above viewport top, flip BELOW instead.
 */
function _initTooltips() {
  const tip = document.getElementById('global-tooltip');
  document.querySelectorAll('.info-wrap').forEach(wrap => {
    const icon = wrap.querySelector('.info-icon');
    const html = wrap.getAttribute('data-tip') || '';

    icon.addEventListener('mouseenter', () => {
      tip.innerHTML = html;
      tip.style.display = 'block';

      const iconRect = icon.getBoundingClientRect();
      const tipH    = tip.offsetHeight;
      const tipW    = tip.offsetWidth;
      const margin  = 8;

      // Horizontal: right-align with icon, clamp to viewport left edge
      let left = iconRect.right - tipW;
      if (left < 4) left = 4;

      // Vertical: try above first
      let top = iconRect.top - tipH - margin;
      if (top < 4) {
        // Flip below
        top = iconRect.bottom + margin;
      }

      tip.style.left = left + 'px';
      tip.style.top  = top  + 'px';
    });

    icon.addEventListener('mouseleave', () => {
      tip.style.display = 'none';
    });
  });
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §CONTROLS — UI event handlers                              ║
// ╚══════════════════════════════════════════════════════════════╝

function toggleAll(id) {
  const inputs = document.querySelectorAll(`#${id} input[type="checkbox"]`);
  const allOn  = Array.from(inputs).every(i => i.checked);
  inputs.forEach(i => { i.checked = !allOn; });
  updateCharts();
}

function resetView() {
  // Checkboxes in lists → all on
  document.querySelectorAll('.check-list input[type="checkbox"]').forEach(i => i.checked = true);

  // Standalone filter checkboxes → off
  ['filter-sig-traj','filter-sig-early','filter-sig-late',
   'chk-curved'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.checked = false;
  });

  // Radio defaults
  document.querySelector('input[name="mode"][value="weighted"]').checked = true;
  document.querySelector('input[name="ytype"][value="nes"]').checked     = true;

  // Select default
  document.getElementById('color-by').value = 'pattern';

  // Search
  document.getElementById('search-box').value = '';
  highlightSearch = '';

  // Sliders → full range
  SLIDER_KEYS.forEach(key => {
    const lo = document.getElementById(`rs-${key}-lo`);
    const hi = document.getElementById(`rs-${key}-hi`);
    if (lo) lo.value = lo.min;
    if (hi) hi.value = hi.max;
    syncSlider(key, /*isInit*/ true);
  });

  updateCharts();
}

function updateHighlights() {
  highlightSearch = document.getElementById('search-box').value.toLowerCase().trim();
  updateCharts();
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §SETTINGS — single-source-of-truth UI state reader         ║
// ╚══════════════════════════════════════════════════════════════╝

function _readSettings() {
  const nesRange     = _sliderValues('min-nes');
  const earlyRange   = _sliderValues('early-nes');
  const lateRange    = _sliderValues('late-nes');
  const trajdevRange = _sliderValues('trajdev-nes');

  return {
    mode:          document.querySelector('input[name="mode"]:checked').value,
    yType:         document.querySelector('input[name="ytype"]:checked').value,
    curved:        document.getElementById('chk-curved').checked,
    colorBy:       document.getElementById('color-by').value,
    patterns:      _checkedValues('#pattern-list input:checked'),
    dbs:           _checkedValues('#db-list input:checked'),
    drivers:       _checkedValues('#driver-list input:checked'),

    // NES range sliders  {lo, hi}
    nesRange,
    earlyRange,
    lateRange,
    trajdevRange,

    // Significance toggles
    sigTraj:   document.getElementById('filter-sig-traj').checked,
    sigEarly:  document.getElementById('filter-sig-early').checked,
    sigLate:   document.getElementById('filter-sig-late').checked,

    // Highlight search (already lowercased in updateHighlights)
    highlight: highlightSearch,
  };
}

function _checkedValues(selector) {
  return Array.from(document.querySelectorAll(selector)).map(i => i.value);
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §FILTERS — pure predicate functions (no side-effects)      ║
// ╚══════════════════════════════════════════════════════════════╝

/**
 * Global scope filter: retains rows whose database AND at least one mutation
 * pattern are in the selected sets.  Applied once per cycle; result is stored
 * in scopeData.
 */
function _applyGlobalFilter(s) {
  return RAW_DATA.filter(d =>
    s.dbs.includes(d.database) && (
      s.patterns.includes(d['Pattern_G32A']) ||
      s.patterns.includes(d['Pattern_R403C'])
    )
  );
}

/**
 * Per-mutation visibility predicate.
 *
 * Returns true iff *d* should be drawn for mutation *mut* under settings *s*.
 *
 * Logic order (fastest checks first):
 *   1. Pattern in selected set
 *   2. NES range windows (any-contrast + per-contrast)
 *   3. Significance gates
 *
 * Bugs 4 & 5 fix: this function is now used to compute the status-bar count,
 * not just to filter drawing.  The count reflects the union of pathways
 * visible in *either* mutation.
 */
function _isVisible(d, mut, s) {
  const pat = d[`Pattern_${mut}`];
  if (!pat || !s.patterns.includes(pat)) return false;

  const absEarly   = Math.abs(d[`NES_Early_${mut}`]   ?? 0);
  const absLate    = Math.abs(d[`NES_Late_${mut}`]    ?? 0);
  const absTrajdev = Math.abs(d[`NES_TrajDev_${mut}`] ?? 0);

  // Global |NES| range: passes if ANY contrast is inside the [lo, hi] window
  const nesLo = s.nesRange.lo, nesHi = s.nesRange.hi;
  const anyInNesWindow = (
    (absEarly   >= nesLo && absEarly   <= nesHi) ||
    (absLate    >= nesLo && absLate    <= nesHi) ||
    (absTrajdev >= nesLo && absTrajdev <= nesHi)
  );
  if (!anyInNesWindow) return false;

  // Per-contrast NES windows (both bounds must be met, independently)
  const {lo: eLo, hi: eHi} = s.earlyRange;
  if (absEarly   < eLo || absEarly   > eHi)   return false;
  const {lo: lLo, hi: lHi} = s.lateRange;
  if (absLate    < lLo || absLate    > lHi)    return false;
  const {lo: tLo, hi: tHi} = s.trajdevRange;
  if (absTrajdev < tLo || absTrajdev > tHi) return false;

  // GSVA driver filter — per-mutation gate.
  // null driver → treated as '__null__' bucket.
  const drv = d[`Driver_${mut}`];
  const drvKey = (drv == null) ? '__null__' : drv;
  if (!s.drivers.includes(drvKey)) return false;

  // Significance gates
  if (s.sigTraj  && !d[`Sig_TrajDev_${mut}`])                          return false;
  if (s.sigEarly && !(d[`padj_Early_${mut}`] < 0.05))                  return false;
  if (s.sigLate  && !(d[`padj_Late_${mut}`]  < 0.05))                  return false;

  return true;
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §GEOMETRY — Bézier + text helpers                          ║
// ╚══════════════════════════════════════════════════════════════╝

/**
 * Sample n+1 points along a quadratic Bézier curve through three control
 * points (x0,y0)→(x1,y1)→(x2,y2).
 */
function _bezier(x0, y0, x1, y1, x2, y2, n) {
  const xs = [], ys = [];
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    const a = (1-t)*(1-t), b = 2*(1-t)*t, c = t*t;
    xs.push(a*x0 + b*x1 + c*x2);
    ys.push(a*y0 + b*y1 + c*y2);
  }
  return { x: xs, y: ys };
}

function _wrapText(str, w) {
  if (!str) return '';
  w = w || 35;
  const words = str.split(' ');
  let line = '', result = '';
  words.forEach(word => {
    if ((line + word).length > w) { result += line + '<br>'; line = word + ' '; }
    else line += word + ' ';
  });
  return (result + line).trim();
}

/**
 * Compute the Bézier control-point vertical offset from the segment midpoint.
 *
 * Bug 3 (tooltip explanation) documents this:
 *   - Offset magnitude = |trajdevNES| × scale factor.
 *   - Offset direction = sign of trajdevNES (positive → upward bow in NES
 *     space; in rank mode the y-axis is reversed so sign is also flipped).
 */
function _curveOffset(trajdevNES, yType) {
  if (yType === 'rank') {
    const rangeSpan = globalYMax - globalYMin || 1;
    return -1 * trajdevNES * rangeSpan * BEZIER_RANK_SCALE;
  }
  return trajdevNES * BEZIER_NES_SCALE;
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §COLOUR — NES diverging palette helpers                    ║
// ╚══════════════════════════════════════════════════════════════╝

function _hex2rgb(h) {
  h = h.replace('#','');
  return [parseInt(h.slice(0,2),16), parseInt(h.slice(2,4),16), parseInt(h.slice(4,6),16)];
}
function _lerp(c1, c2, t) {
  const [r1,g1,b1] = _hex2rgb(c1), [r2,g2,b2] = _hex2rgb(c2);
  const r = Math.round(r1+(r2-r1)*t), g = Math.round(g1+(g2-g1)*t), b = Math.round(b1+(b2-b1)*t);
  return '#' + [r,g,b].map(v => v.toString(16).padStart(2,'0')).join('');
}
function _nesColor(nes) {
  if (nes == null || isNaN(nes)) return '#999';
  const t = Math.max(-1, Math.min(1, nes / NES_CMAX));
  return t <= 0 ? _lerp(CLR_NEG, CLR_MID, 1+t) : _lerp(CLR_MID, CLR_POS, t);
}
function _lineStyle(settings, pattern) {
  if (settings.mode === 'uniform') return UNIFORM_STYLE;
  const cat = WEIGHT_CATS[pattern] || 'rare';
  return WEIGHT_STYLES[cat] || UNIFORM_STYLE;
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §TRACES — chart trace builders                             ║
// ╚══════════════════════════════════════════════════════════════╝

// ── Tooltip text ────────────────────────────────────────────────

function _buildTooltip(d, mut) {
  const fmt  = v => (v != null && !isNaN(v)) ? (+v).toFixed(2) : 'N/A';
  const fmtP = v => (v != null && !isNaN(v)) ? (+v).toExponential(2) : 'N/A';
  const sig  = v => (v != null && v < 0.05)  ? ' *' : '';
  const pat  = d[`Pattern_${mut}`] || '—';
  const desc = _wrapText(d.description || d.pathway_id, 38);
  const def  = _wrapText(PATTERN_DEFS[pat] || '', 48);

  return (
    `<b>${desc}</b><br>` +
    `<span style="color:#aaa;">DB: ${d.database} | Pattern: ${pat}</span><br><br>` +
    `<b>NES  [p.adj]:</b><br>` +
    `Early:   ${fmt(d[`NES_Early_${mut}`])}  [${fmtP(d[`padj_Early_${mut}`])}]${sig(d[`padj_Early_${mut}`])}<br>` +
    `TrajDev: ${fmt(d[`NES_TrajDev_${mut}`])} [${fmtP(d[`padj_TrajDev_${mut}`])}]${sig(d[`padj_TrajDev_${mut}`])}<br>` +
    `Late:    ${fmt(d[`NES_Late_${mut}`])}   [${fmtP(d[`padj_Late_${mut}`])}]${sig(d[`padj_Late_${mut}`])}<br>` +
    (def ? `<br><i style="color:#ccc;">${def}</i>` : '')
  );
}

// ── Segment geometry ─────────────────────────────────────────────

/**
 * Build the (x, y, hover) arrays for one pathway segment.
 *
 * Bug 7 fix: curved segments now reuse _bezier() via _buildCurvedSegment,
 * so highlighted pathways follow the same geometry as the regular traces.
 */
function _buildSegment(d, mut, settings, forceStyle) {
  const yKey = settings.yType === 'rank' ? 'Rank' : 'NES';
  const y1   = d[`${yKey}_Early_${mut}`];
  const y2   = d[`${yKey}_Late_${mut}`];
  if (y1 == null || y2 == null) return null;

  const hover = _buildTooltip(d, mut);
  const pid   = d.pathway_id;

  if (settings.curved && d[`Sig_TrajDev_${mut}`]) {
    const traj = d[`NES_TrajDev_${mut}`] || 0;
    const yMid = (y1 + y2) / 2;
    const off  = _curveOffset(traj, settings.yType);
    const pts  = _bezier(0, y1, 0.5, yMid + off, 1, y2, 24);
    return {
      xs:    [...pts.x, null],
      ys:    [...pts.y, null],
      hover: [...pts.x.map(() => hover), null],
      cd:    [...pts.x.map(() => pid), null],
    };
  }
  // Straight line — customdata carries pathway_id at each real (non-null) point
  return {
    xs: [0, 1, null],
    ys: [y1, y2, null],
    hover: [hover, hover, null],
    cd: [pid, pid, null],
  };
}

// ── Pattern-colored traces ────────────────────────────────────────

/**
 * Returns one Plotly trace per pattern group.
 *
 * Bug 7 fix: when a highlight search is active, all traces in this function
 * are dimmed to DIMMED_STYLE.  A separate highlight overlay is added by
 * _buildHighlightTraces().
 */
function _buildPatternTraces(rows, mut, settings) {
  const byPattern = {};
  rows.forEach(d => {
    const p = d[`Pattern_${mut}`];
    if (!byPattern[p]) byPattern[p] = [];
    byPattern[p].push(d);
  });

  // Render dominant patterns first (background) → rare patterns last (foreground)
  const catRank = { dominant:0, common:1, uncommon:2, rare:3 };
  const sorted  = Object.keys(byPattern).sort((a, b) => {
    const da = catRank[WEIGHT_CATS[a]||'rare'] - catRank[WEIGHT_CATS[b]||'rare'];
    return da !== 0 ? da : byPattern[b].length - byPattern[a].length;
  });

  const dimmed = !!settings.highlight;

  return sorted.map(pat => {
    const color            = PATTERN_COLORS[pat] || '#999';
    const { width, opacity } = dimmed ? DIMMED_STYLE : _lineStyle(settings, pat);
    const xs=[], ys=[], texts=[], cds=[];
    byPattern[pat].forEach(d => {
      const seg = _buildSegment(d, mut, settings);
      if (!seg) return;
      xs.push(...seg.xs); ys.push(...seg.ys); texts.push(...seg.hover); cds.push(...seg.cd);
    });
    if (!xs.length) return null;
    return {
      x: xs, y: ys, mode: 'lines', type: 'scattergl',
      line: { color, width }, opacity,
      name: pat, text: texts, hoverinfo: 'text',
      customdata: cds,
    };
  }).filter(Boolean);
}

// ── NES-diverging colored traces ──────────────────────────────────

function _buildNesColorTraces(rows, mut, settings) {
  const dimmed = !!settings.highlight;
  return rows.map(d => {
    const pat = d[`Pattern_${mut}`];
    let nesVal = 0;
    if (settings.colorBy === 'nes_early')   nesVal = d[`NES_Early_${mut}`]   || 0;
    if (settings.colorBy === 'nes_late')    nesVal = d[`NES_Late_${mut}`]    || 0;
    if (settings.colorBy === 'nes_trajdev') nesVal = d[`NES_TrajDev_${mut}`] || 0;
    const color            = _nesColor(nesVal);
    const { width, opacity } = dimmed ? DIMMED_STYLE : _lineStyle(settings, pat);
    const seg = _buildSegment(d, mut, settings);
    if (!seg) return null;
    return {
      x: seg.xs, y: seg.ys, mode: 'lines', type: 'scattergl',
      line: { color, width }, opacity,
      name: (d.description || '').slice(0, 30),
      text: seg.hover, hoverinfo: 'text', showlegend: false,
      customdata: seg.cd,
    };
  }).filter(Boolean);
}

// ── Highlight overlay traces  (Bug 7) ────────────────────────────

/**
 * Build highlight overlay traces for pathways whose description matches the
 * search term.
 *
 * Fixes for Bug 7:
 *   a) Uses _buildSegment() so curved geometry is preserved when curved
 *      mode is active.
 *   b) Drawn as a separate high-z trace; non-matching lines are dimmed (see
 *      DIMMED_STYLE applied in the pattern/NES trace builders above).
 *   c) Uses bold black lines with markers so matches are unmistakable.
 */
function _buildHighlightTraces(rows, mut, settings) {
  if (!settings.highlight) return [];

  // Colour-by-pattern uses 'pattern' mode; NES modes show single colour
  const matchRows = rows.filter(d =>
    (d.description || '').toLowerCase().includes(settings.highlight) ||
    (d.pathway_id  || '').toLowerCase().includes(settings.highlight)
  );
  if (!matchRows.length) return [];

  // Group by pattern so each match can use its pattern colour overlaid with a
  // thick black stroke for maximum contrast
  const byPattern = {};
  matchRows.forEach(d => {
    const p = d[`Pattern_${mut}`] || '_unknown';
    if (!byPattern[p]) byPattern[p] = [];
    byPattern[p].push(d);
  });

  const traces = [];
  Object.entries(byPattern).forEach(([pat, group]) => {
    const baseColor = (settings.colorBy === 'pattern')
      ? (PATTERN_COLORS[pat] || '#E74C3C')
      : '#E74C3C';
    const xs=[], ys=[], texts=[], cds=[];
    group.forEach(d => {
      const seg = _buildSegment(d, mut, settings);
      if (!seg) return;
      xs.push(...seg.xs); ys.push(...seg.ys); texts.push(...seg.hover); cds.push(...seg.cd);
    });
    if (!xs.length) return;

    // Thick black outline pass
    traces.push({
      x: xs, y: ys, mode: 'lines', type: 'scatter',
      line: { color: 'black', width: 5 }, opacity: 1.0,
      showlegend: false, hoverinfo: 'skip',
      customdata: cds,
    });
    // Coloured inner pass
    traces.push({
      x: xs, y: ys, mode: 'lines+markers', type: 'scatter',
      line: { color: baseColor, width: 3 },
      marker: { size: 7, color: baseColor, symbol: 'circle',
                line: { color: 'black', width: 1.5 } },
      opacity: 1.0, name: `Match: ${pat}`,
      text: texts, hoverinfo: 'text',
      customdata: cds,
    });
  });

  return traces;
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §RENDER — Plotly chart assembly                            ║
// ╚══════════════════════════════════════════════════════════════╝

function renderChart(divId, mut, settings, visibleRows) {
  // Base traces (dimmed when highlight active)
  const traces = settings.colorBy === 'pattern'
    ? _buildPatternTraces(visibleRows, mut, settings)
    : _buildNesColorTraces(visibleRows, mut, settings);

  // Highlight overlay (Bug 7: drawn on top)
  const hlTraces = _buildHighlightTraces(visibleRows, mut, settings);
  traces.push(...hlTraces);

  const isRank = settings.yType === 'rank';

  const layout = {
    margin:     { t: 30, b: 30, l: 55, r: 20 },
    showlegend: false,
    xaxis: {
      tickvals:   [0, 1],
      ticktext:   ['Early (D35)', 'Late (D65)'],
      range:      [-0.15, 1.15],
      zeroline:   false,
      fixedrange: true,
    },
    yaxis: {
      title:     { text: isRank ? 'Rank (1 = highest NES)' : 'NES', standoff: 8 },
      range:     [globalYMin, globalYMax],
      autorange: isRank ? 'reversed' : false,
      zeroline:  !isRank,
      zerolinecolor: '#ccc',
    },
    hovermode:  'closest',
    hoverlabel: {
      bgcolor:    'white',
      bordercolor: '#555',
      font:       { size: 11, color: '#222', family: 'system-ui,sans-serif' },
      namelength: -1,
      align:      'left',
    },
  };

  Plotly.react(divId, traces, layout, { displayModeBar: false, responsive: true })
    .then(() => {
      // §GSVA_PANEL: register click handler once per chart div (idempotent).
      const div = document.getElementById(divId);
      if (div && !div._gsvaClickWired) {
        div.on('plotly_click', evt => {
          const pt = evt && evt.points && evt.points[0];
          const pid = pt && pt.customdata;
          if (pid) openGsvaModal(pid);
        });
        div._gsvaClickWired = true;
      }
    });
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §UPDATE — orchestrator  (Bugs 4 & 5 fix)                   ║
// ╚══════════════════════════════════════════════════════════════╝

/**
 * updateCharts() is the single orchestration function called by every UI
 * control.  It:
 *   1. Reads all settings.
 *   2. Applies the global scope filter (database + pattern).
 *   3. Applies per-mutation visibility filters to get visibleG32A / visibleR403C.
 *   4. Computes the shared Y axis from the union of visible rows.
 *   5. Updates the status-bar with the count of UNIQUE pathways visible in
 *      either mutation (Bugs 4 & 5 fix: count is now correct for all filters).
 *   6. Renders both charts.
 */
function updateCharts() {
  const s = _readSettings();

  // Step 1: global scope filter
  scopeData = _applyGlobalFilter(s);

  // Step 2: per-mutation visibility
  visibleG32A  = scopeData.filter(d => _isVisible(d, 'G32A',  s));
  visibleR403C = scopeData.filter(d => _isVisible(d, 'R403C', s));

  // Step 3: union count for status bar
  const unionIds = new Set([
    ...visibleG32A.map(d  => d.pathway_id),
    ...visibleR403C.map(d => d.pathway_id),
  ]);

  // Step 4: shared Y axis
  if (s.yType === 'nes') {
    let mx = 0;
    [...visibleG32A, ...visibleR403C].forEach(d => {
      ['G32A','R403C'].forEach(mut => {
        ['NES_Early_','NES_Late_','NES_TrajDev_'].forEach(pf => {
          const v = d[pf + mut];
          if (v != null && !isNaN(v)) mx = Math.max(mx, Math.abs(v));
        });
      });
    });
    mx = Math.max(mx * 1.08, 1.0);
    globalYMin = -mx; globalYMax = mx;
  } else {
    // Rank: 1 at top (reversed axis), total = largest rank visible
    const maxRank = Math.max(
      ...visibleG32A.map(d  => d['Rank_Late_G32A']  || 0),
      ...visibleR403C.map(d => d['Rank_Late_R403C'] || 0),
      1,
    );
    globalYMin = 1; globalYMax = maxRank;
  }

  // Step 5: status bar (Bugs 4 & 5 fix)
  const hlCount = unionIds.size > 0 && s.highlight
    ? Array.from(unionIds).filter(id => {
        const d = RAW_DATA.find(r => r.pathway_id === id);
        return d && (
          (d.description || '').toLowerCase().includes(s.highlight) ||
          (d.pathway_id  || '').toLowerCase().includes(s.highlight)
        );
      }).length
    : null;

  let statusText = `Showing ${unionIds.size.toLocaleString()} pathway${unionIds.size !== 1 ? 's' : ''}`;
  if (hlCount !== null) statusText += ` · ${hlCount} highlighted`;
  document.getElementById('status-bar').textContent = statusText;

  // Step 6: render
  renderChart('chart-g32a',  'G32A',  s, visibleG32A);
  renderChart('chart-r403c', 'R403C', s, visibleR403C);
}

// ╔══════════════════════════════════════════════════════════════╗
// ║  §GSVA_PANEL — replicate-level click-through modal           ║
// ╚══════════════════════════════════════════════════════════════╝
//
// Clicking a pathway line opens a modal showing its per-sample GSVA
// enrichment scores. The shared sample index in METADATA.gsva_sample_index
// gives canonical column order; each pathway's d.gsva_scores is positionally
// aligned. Group means use approximate 95% CI = mean ± 1.96·SE (n is small;
// surfaced for qualitative inspection, not formal inference).

const GSVA_GROUP_COLORS = {
  Ctrl:  '#666666',
  G32A:  '#0072B2',
  R403C: '#D55E00',
};
const GSVA_DAY_ORDER      = ['D35', 'D65'];
const GSVA_GENOTYPE_ORDER = ['Ctrl', 'G32A', 'R403C'];

/** Build a Set of pathway_ids that have any GSVA data, for empty-state copy. */
const _GSVA_AVAILABLE = (() => {
  const s = new Set();
  RAW_DATA.forEach(d => {
    if (Array.isArray(d.gsva_scores) && d.gsva_scores.some(v => v != null)) {
      s.add(d.pathway_id);
    }
  });
  return s;
})();

function _gsvaGroupKey(g, day) { return `${g}_${day}`; }

function _gsvaGroupOrder() {
  // Ctrl-D35, G32A-D35, R403C-D35, Ctrl-D65, G32A-D65, R403C-D65
  const out = [];
  GSVA_DAY_ORDER.forEach(day => {
    GSVA_GENOTYPE_ORDER.forEach(g => out.push({ genotype: g, day, key: _gsvaGroupKey(g, day) }));
  });
  return out;
}

function _gsvaIndexByGroup() {
  const index = METADATA.gsva_sample_index || [];
  const byGroup = {};
  index.forEach((s, i) => {
    const k = _gsvaGroupKey(s.genotype, s.day);
    if (!byGroup[k]) byGroup[k] = [];
    byGroup[k].push({ idx: i, sample_id: s.sample_id });
  });
  return byGroup;
}

function _mean(arr) {
  const v = arr.filter(x => x != null && !isNaN(x));
  if (!v.length) return null;
  return v.reduce((a,b) => a+b, 0) / v.length;
}

// Geometry helpers ported from 02_Analysis/revision/supplements/Supp10.replicate_level_gsva.py
// so the modal mirrors the Supp10 figure (Early/Late x-stops with within-day
// genotype offsets, deterministic per-sample spread, median+IQR aggregation).
const _GSVA_GENO_OFFSET = { Ctrl: -0.28, G32A: 0.0, R403C: +0.28 };
const _GSVA_DAY_X       = { D35: 1.0, D65: 2.0 };

function _xPosition(day, genotype) {
  return (_GSVA_DAY_X[day] ?? 0) + (_GSVA_GENO_OFFSET[genotype] ?? 0);
}

function _deterministicSpread(n, halfWidth = 0.07) {
  if (n <= 1) return [0];
  const step = (2 * halfWidth) / (n - 1);
  const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = -halfWidth + i * step;
  return out;
}

function _median(arr) {
  const v = arr.filter(x => x != null && !isNaN(x)).slice().sort((a,b) => a-b);
  if (!v.length) return null;
  const mid = v.length / 2;
  return v.length % 2 ? v[Math.floor(mid)] : 0.5 * (v[mid - 1] + v[mid]);
}

function _quantile(arr, q) {
  const v = arr.filter(x => x != null && !isNaN(x)).slice().sort((a,b) => a-b);
  if (!v.length) return null;
  if (v.length === 1) return v[0];
  const pos = (v.length - 1) * q;
  const lo  = Math.floor(pos), hi = Math.ceil(pos);
  return v[lo] + (v[hi] - v[lo]) * (pos - lo);
}

function _hexToRgba(hex, alpha) {
  const h = hex.replace('#', '');
  const r = parseInt(h.substr(0, 2), 16);
  const g = parseInt(h.substr(2, 2), 16);
  const b = parseInt(h.substr(4, 2), 16);
  return `rgba(${r},${g},${b},${alpha})`;
}

function openGsvaModal(pid) {
  const pathway = RAW_DATA.find(d => d.pathway_id === pid);
  if (!pathway) return;
  const modal = document.getElementById('gsva-modal');

  // Header content
  document.getElementById('gsva-modal-title').textContent =
    pathway.description || pathway.pathway_id;
  document.getElementById('gsva-modal-subtitle').textContent =
    `${pathway.database} · ${pathway.pathway_id}`;

  const patEl = document.getElementById('gsva-modal-patterns');
  patEl.innerHTML = '';
  const pG = pathway.Pattern_G32A, pR = pathway.Pattern_R403C;
  function _chip(label, pattern) {
    if (!pattern) return null;
    const span = document.createElement('span');
    span.className = 'gsva-pattern-chip';
    const c = PATTERN_COLORS[pattern] || '#ccc';
    span.style.background = c + '22';
    span.style.borderLeft = `3px solid ${c}`;
    span.textContent = `${label}: ${pattern}`;
    return span;
  }
  [_chip('G32A GSEA pattern', pG), _chip('R403C GSEA pattern', pR)]
    .filter(Boolean).forEach(c => patEl.appendChild(c));

  const captionEl = document.getElementById('gsva-caption');
  const plotEl    = document.getElementById('gsva-plot');

  modal.classList.add('is-open');

  const scores = pathway.gsva_scores;
  const sampleIndex = METADATA.gsva_sample_index || [];
  if (!Array.isArray(scores) || !scores.some(v => v != null) || !sampleIndex.length) {
    plotEl.innerHTML =
      '<div class="gsva-empty">No per-sample GSVA data is available for this pathway ' +
      '(typically because the gene set fell outside the size filter of the GSVA pipeline).</div>';
    captionEl.textContent = '';
    return;
  }

  // Supp10-style layout: two x-stops (Early D35, Late D65); within each day,
  // genotypes are offset ±0.28 around their day-x so the Ctrl/G32A/R403C
  // clusters separate cleanly. Per-genotype trajectory connects the Early
  // and Late medians; the IQR band fills between Q1 and Q3 along that path.
  const byGroup = _gsvaIndexByGroup();
  const dotTraces   = [];
  const bandTraces  = [];
  const lineTraces  = [];
  const ringTraces  = [];

  // Track per-genotype day-level summaries for the trajectory + band layer.
  // genoStats[geno] = { D35: {median,q1,q3,n}, D65: {...} }
  const genoStats = {};
  GSVA_GENOTYPE_ORDER.forEach(g => { genoStats[g] = {}; });

  GSVA_DAY_ORDER.forEach(day => {
    GSVA_GENOTYPE_ORDER.forEach(geno => {
      const members = byGroup[_gsvaGroupKey(geno, day)] || [];
      const pairs = members
        .map(m => ({ v: scores[m.idx], sid: m.sample_id }))
        .filter(p => p.v != null && !isNaN(p.v));
      if (!pairs.length) return;

      const vals    = pairs.map(p => p.v);
      const xCentre = _xPosition(day, geno);
      const offsets = _deterministicSpread(vals.length, 0.07);
      const color   = GSVA_GROUP_COLORS[geno] || '#666';

      dotTraces.push({
        x: pairs.map((_, i) => xCentre + offsets[i]),
        y: vals,
        text: pairs.map(p => p.sid),
        mode: 'markers',
        type: 'scatter',
        showlegend: false,
        marker: {
          size: 9,
          color,
          opacity: 0.85,
          line: { color: 'white', width: 1 },
        },
        hovertemplate: `<b>%{text}</b><br>${geno} · ${day}<br>GSVA = %{y:.3f}<extra></extra>`,
      });

      genoStats[geno][day] = {
        median: _median(vals),
        q1: _quantile(vals, 0.25),
        q3: _quantile(vals, 0.75),
        n: vals.length,
      };
    });
  });

  // Build the band, trajectory line, and median rings per genotype.
  GSVA_GENOTYPE_ORDER.forEach(geno => {
    const s35 = genoStats[geno]['D35'];
    const s65 = genoStats[geno]['D65'];
    if (!s35 || !s65) return;

    const color    = GSVA_GROUP_COLORS[geno] || '#666';
    const xPair    = [_GSVA_DAY_X['D35'], _GSVA_DAY_X['D65']];
    const fillRgba = _hexToRgba(color, 0.15);

    // IQR band: lower (Q1) trace then upper (Q3) trace with fill:'tonexty'.
    bandTraces.push({
      x: xPair, y: [s35.q1, s65.q1],
      mode: 'lines', type: 'scatter',
      line: { color: 'rgba(0,0,0,0)', width: 0 },
      hoverinfo: 'skip', showlegend: false,
    });
    bandTraces.push({
      x: xPair, y: [s35.q3, s65.q3],
      mode: 'lines', type: 'scatter',
      line: { color: 'rgba(0,0,0,0)', width: 0 },
      fill: 'tonexty', fillcolor: fillRgba,
      hoverinfo: 'skip', showlegend: false,
    });

    // Trajectory line connecting Early -> Late medians.
    lineTraces.push({
      x: xPair, y: [s35.median, s65.median],
      mode: 'lines', type: 'scatter',
      line: { color, width: 2.8 },
      hoverinfo: 'skip', showlegend: false,
    });

    // White-ring median markers — one per (genotype, day).
    ringTraces.push({
      x: xPair, y: [s35.median, s65.median],
      mode: 'markers', type: 'scatter',
      marker: { size: 14, color: 'white', symbol: 'circle',
                line: { color, width: 2.6 } },
      text: [
        `${geno} · D35 — median ${s35.median?.toFixed(3)} (n=${s35.n}, IQR ${s35.q1?.toFixed(3)}–${s35.q3?.toFixed(3)})`,
        `${geno} · D65 — median ${s65.median?.toFixed(3)} (n=${s65.n}, IQR ${s65.q1?.toFixed(3)}–${s65.q3?.toFixed(3)})`,
      ],
      hovertemplate: '%{text}<extra></extra>',
      showlegend: false,
    });
  });

  const layout = {
    margin: { t: 10, b: 50, l: 60, r: 12 },
    xaxis: {
      tickmode: 'array',
      tickvals: [_GSVA_DAY_X['D35'], _GSVA_DAY_X['D65']],
      ticktext: ['Early (D35)', 'Late (D65)'],
      range: [0.4, 2.6],
      zeroline: false,
      fixedrange: true,
    },
    yaxis: {
      title: { text: 'GSVA enrichment score', standoff: 8 },
      zeroline: false,
      fixedrange: true,
    },
    showlegend: false,
    hovermode: 'closest',
    shapes: [
      // Dashed y = 0 reference line (Supp10 convention).
      { type: 'line', x0: 0.4, x1: 2.6, y0: 0, y1: 0, xref: 'x', yref: 'y',
        line: { color: 'gray', width: 1, dash: 'dash' } },
    ],
  };

  // Draw order (back -> front): bands, trajectory lines, sample dots, rings.
  const traces = [...bandTraces, ...lineTraces, ...dotTraces, ...ringTraces];
  Plotly.react(plotEl, traces, layout,
    { displayModeBar: false, responsive: true });

  // ── GSVA driver verdict ─────────────────────────────────────────────────
  // Compute per-genotype Δ (D65 − D35) from the genoStats we already built.
  const deltas = {};
  GSVA_GENOTYPE_ORDER.forEach(g => {
    const a = genoStats[g]?.D35, b = genoStats[g]?.D65;
    deltas[g] = (a && b && a.median != null && b.median != null)
      ? b.median - a.median
      : null;
  });

  const verdictEl = document.getElementById('gsva-driver-verdict');
  const _fmt = v => (v == null) ? 'n/a' : (v >= 0 ? '+' : '') + v.toFixed(2);
  const _color = geno => GSVA_GROUP_COLORS[geno] || '#666';

  // Delta rows (Ctrl grey, mutants coloured).
  let deltaHtml = '<div class="verdict-deltas">';
  GSVA_GENOTYPE_ORDER.forEach(g => {
    const col = (g === 'Ctrl') ? '#888' : _color(g);
    deltaHtml +=
      `<span style="color:${col};margin-right:14px;">Δ<sub>${_esc(g)}</sub>&nbsp;=&nbsp;${_fmt(deltas[g])}</span>`;
  });
  deltaHtml += '</div>';

  // Verdict lines — one per mutation.
  let verdictHtml = '';
  ['G32A', 'R403C'].forEach(mut => {
    const rawDriver = pathway[`Driver_${mut}`];
    const col = _color(mut);
    if (rawDriver == null) {
      verdictHtml +=
        `<div class="verdict-line" style="color:${col};">${_esc(mut)}: ` +
        `<span style="font-weight:400;color:#888;">(no GSVA driver — pathway has no per-sample data)</span></div>`;
    } else {
      const humanLabel = GSVA_DRIVER_LABELS[rawDriver] || rawDriver;
      verdictHtml +=
        `<div class="verdict-line" style="color:${col};">${_esc(mut)}: ${_esc(humanLabel)}</div>`;
    }
  });

  verdictEl.innerHTML = deltaHtml + verdictHtml;

  captionEl.innerHTML =
    'Each dot is one biological replicate (deterministic horizontal spread around its genotype cluster). ' +
    'White rings mark group medians; shaded bands span the IQR (Q1–Q3); thin lines connect Early → Late medians per genotype. ' +
    '<b>Reading the verdict</b>: Δ = median(D65) − median(D35) per group (ε = 0.10). ' +
    'Sign-reversal patterns appear as a Ctrl ↔ mutant crossover; ' +
    'Compensation can reflect either a mutant rebound (<em>mutant-driven</em>) or a closing gap driven by Ctrl developmental dynamics ' +
    '(<em>ctrl-driven</em>) — the Δ values and verdict line resolve which.';
}

function closeGsvaModal() {
  const modal = document.getElementById('gsva-modal');
  modal.classList.remove('is-open');
}

// Close on Escape
document.addEventListener('keydown', e => {
  if (e.key === 'Escape') closeGsvaModal();
});

// ════════════════════════════════════════
//  Bootstrap
// ════════════════════════════════════════
document.addEventListener('DOMContentLoaded', init);
// Fallback for synchronous script execution
if (document.readyState !== 'loading') init();
</script>
</body>
</html>
"""
