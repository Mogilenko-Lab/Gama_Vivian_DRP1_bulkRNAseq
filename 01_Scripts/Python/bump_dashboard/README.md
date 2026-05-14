# `bump_dashboard` — Interactive Pathway-Trajectory Bump-Chart Dashboard

Self-contained Python package that generates a **single, self-contained HTML
dashboard** visualising how gene-set enrichment scores (NES) evolve across two
differentiation time-points (Early D35 → Late D65) for two DRP1 mutations
(G32A, R403C).

Output: `03_Results/02_Analysis/Plots/Trajectory_Flow/interactive_bump_dashboard.html`

---

## What it does

1. **Loads** the wide-format GSEA master table produced upstream by
   `viz_bump_charts.load_data()`.
2. **Enriches** the table with per-contrast adjusted p-values (pivoted from
   long format) for Early, Late, and TrajDev contrasts of both mutations.
3. **Filters** to the configured scope union (`"focused"` + `"significant"` by
   default) and removes pathways with no trajectory pattern classification.
4. **Computes** frequency-based weight categories per pattern (`dominant /
   common / uncommon / rare`) and NES-descending ranks per mutation per stage.
5. **Validates** every row against Pydantic `PathwayRecord` / `MutationStats`
   schemas; invalid rows are skipped with a warning.
6. **Serialises** domain objects to JSON-safe flat dicts and injects them into
   a self-contained HTML template (Plotly CDN, no server required).
7. **Renders** two side-by-side Plotly bump charts (one per mutation) with a
   shared interactive sidebar and writes the result to disk.

The dashboard supports:
- **Y-axis toggle**: NES value vs. rank.
- **Colour-by mode**: trajectory pattern, Early NES, Late NES, or TrajDev NES.
- **Visual style**: straight lines or quadratic Bézier curves (curvature
  magnitude = |TrajDev NES|; direction = sign of TrajDev NES).
- **Dual-handle range sliders** for |NES|, Early NES, Late NES, TrajDev NES.
- **Database and pattern checkboxes**, significance toggles, and a keyword
  highlight search that draws matching traces on top.
- **Smart tooltips** that flip above/below the viewport edge.
- **Live status bar** that counts visible pathways after all active filters.

---

## Architecture

The package uses a strict four-layer dependency hierarchy with no circular
imports:

```
domain  ←  application  ←  infrastructure
                        ←  presentation
```

| Layer | Package | Key class | Responsibility |
|-------|---------|-----------|----------------|
| Domain | `domain/` | `PathwayRecord`, `MutationStats`, `DashboardConfig` | Pydantic value objects, domain rules, geometry helpers |
| Application | `application/` | `DashboardDataService`, `DashboardSerializer` | Use-case orchestration — load, enrich, validate, serialise |
| Infrastructure | `infrastructure/` | `DashboardOutputWriter` | Sole file-write concern |
| Presentation | `presentation/` | `DashboardHtmlRenderer`, `html_fragments` | HTML/CSS/JS template assembly via sentinel token substitution |
| Entry point | `application/` | `DashboardPipeline` | Wires all layers; single `run()` call |

### Key files

| File | Purpose |
|------|---------|
| `domain/schema.py` | Pydantic contracts (`PathwayRecord`, `MutationStats`, `DashboardConfig`, `DashboardMetadata`); schema version |
| `domain/rules.py` | Pure business rules: significance threshold, NES→colour mapping, rank computation, weight-category thresholds |
| `domain/geometry.py` | Stateless geometric helpers: quadratic Bézier sampling, control-point offset, text wrapping |
| `application/data_service.py` | Loads raw data, pivots p-values, filters by scope, computes ranks/weight categories, assembles domain objects |
| `application/serializer.py` | Flattens domain objects to JSON-safe dicts mirroring JS field names |
| `application/pipeline.py` | Orchestrates the full build → serialise → render → write workflow |
| `presentation/html_fragments.py` | Static HTML/CSS/JS strings (`HEAD`, `BODY`, `SCRIPT`) with `%%TOKEN%%` sentinels |
| `presentation/html_renderer.py` | Combines fragments, validates tokens, injects JSON payloads |
| `infrastructure/output_writer.py` | Resolves output path via `config.get_project_root()`, creates dirs, writes UTF-8 HTML |

---

## Quick start

```python
from Python.bump_dashboard import DashboardPipeline

pipeline = DashboardPipeline()
output_path = pipeline.run()
print(f"Dashboard written to: {output_path}")
```

Custom configuration:

```python
from Python.bump_dashboard import DashboardPipeline
from Python.bump_dashboard.domain.schema import DashboardConfig

cfg = DashboardConfig(
    scopes=["focused"],                        # only focused pathways
    output_filename="bump_focused_only.html",  # custom filename
)
pipeline = DashboardPipeline(config=cfg)
pipeline.run()
```

---

## Design notes

- **Domain objects are frozen** (`model_config = {"frozen": True}`) — all
  mutation happens in the application layer before assembly.
- **All dependency injection** — every collaborator can be replaced with a
  stub, so the pipeline is fully unit-testable without touching the filesystem.
- **Colours are never hardcoded** — NES diverging colours (Blue `#2166AC` ↔
  White ↔ Orange `#B35806`) and pattern colours are sourced from the
  project-canonical `color_config.py` / `pattern_definitions.py`.
- **JS is purely data-driven** — the rendered HTML contains no hard-coded
  pathway data; all state is driven from the two JSON payloads injected at
  build time.
