# AGENTS.md - Universal AI Agent Instructions for Gama_Vivian_DRP1_bulkRNAseq

**Version:** 2.2.0
**Generated:** 2026-05-11
**Purpose:** Single source of truth for AI agent behavior in the project

---

## 1. Rigorous AI Validation & Cross-Checking Rules

**CRITICAL: NEVER guess or assume technical implementation details. Always verify empirically.** 


1. **Verify R API and Default Parameters:**
   - **Do not guess** what default parameters are
   - Before confirming or denying a claim about parameters, check the actual function signature

2. **Verify Experimental Design and Metadata:**
   - You MUST read and parse the actual `metadata.csv` (e.g., using `head` and R scripts) to determine the true group sizes. 
   - Example: Run `Rscript -e 'meta <- read.csv("03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv", sep=";"); table(meta$days, meta$genotype)'` to find the smallest group size

3. **Verify Code Invocation:**
   - Check if the code actually passes the claimed arguments, or if it relies on defaults

4. **Verify Final Data States (Checkpoints):**
   - Use R to load the intermediate RDS checkpoints and check the dimensions: 
     `Rscript -e 'obj <- readRDS("03_Results/02_Analysis/checkpoints/model_objects.rds"); print(nrow(obj$DGE))'`

## 2. Workspace & Environment Clutter Rules (CRITICAL FOR AGENTS)

**NO CLUTTER:**
- Do NOT leave stale test artifacts, temporary scripts (e.g., `test.py`, `debug.R`), or diagnostic output files scattered around the root or analysis directories.
- If you create a temporary file to debug an environment issue (like PDF rendering or library imports), you MUST delete it after extracting the necessary insight.
- The repository must remain clean and strictly follow the documented directory architecture.

**Environment Navigation:**
- Assume Python scripts are run using `python3` with standard data science libraries (`pandas`, `numpy`, `matplotlib`, `seaborn`) available.
- R scripts are executed via `Rscript` or sourced within an R session.
- If dependencies are missing, do not attempt to globally install them without asking, or use a temporary virtual environment that you clean up afterward.
- Be careful with `matplotlib` outputs to PDF; prefer rasterizing heatmaps (`rasterized=True` in seaborn) and embedding fonts (`mpl.rcParams['pdf.fonttype'] = 42`) to avoid viewer rendering bugs.

## 3. Critical Coding Rules
1. **Annotate genes BEFORE filtering** - Never lose gene IDs by filtering first.
2. **Use `filterByExpr()`** - Never use manual count thresholds like `rowSums >= 10`.
3. **Cache anything >1 minute** - Always use the `load_or_compute()` pattern for GSEA, DE fitting, etc.
4. **Never hardcode colors** - Load project colors via `source("02_analysis/config/color_config.R")` (or `01_Scripts/R_scripts/color_config.R`).
5. **Separate concerns** - Processing scripts never plot; viz scripts never compute.

## 3a. Python Style — No Inlining Inside `main()`

**Rule: Never define helper functions *inside* `main()` (or any other function). Every function must be a top-level named module-level block.**

When writing or refactoring Python scripts in `02_Analysis/`:

- **Extract** lambdas, nested `def`s, and ad-hoc inline logic into their own top-level `def func_name(...)`. Each function should have its own ``"""docstring"""`` and be placed in a clearly labeled section of the file (e.g., `# == Low-level helpers`, `# ==== High-level orchestration`).
- **`main()` should be a thin orchestrator** - it reads config/ data, calls helpers, and writes outputs. No domain logic belongs inline in the `if __name__ == "__main__":` guard or `main()`.  
- **`_-prefixed` names are "private to the module" and **should appear _after_ `main()`** in the source file to make the public API surface immediately visible to the reader.  

Anti-pattern:
```python
def main():
    def inline_helper(x):  # ← NEVER
        pass
    pass
```

Pattern:
```python
def _inline_helper(x):  # ← ALWAYS top-level
    """Compute something useful."""
    pass

def main():  # thin orchestrator
    _inline_helper(data)
```

**Key Locations:**
- `00_Data/`: Input data (read-only)
- `01_Scripts/`: Shared code (RNAseq-toolkit, R_scripts helper functions)
- `02_Analysis/`: Project-specific scripts and execution pipelines
- `03_Results/`: All outputs (checkpoints, tables, plots)

## 4. Configuration System

All hardcoded values/thresholds go in configuration files, usually referenced centrally:
**Importing Config in R:**
```r
config <- yaml::read_yaml("02_analysis/config/analysis_config.yaml") # if used
DE_FDR_CUTOFF <- config$thresholds$de_fdr
```

## 6. Current Status & Task Execution

When asked to perform analysis:
1. Review this file (`AGENTS.md`) for methodology and verification rules.
2. Load relevant checkpoints (`03_Results/02_Analysis/checkpoints/`) instead of re-running raw data processing.
3. Keep track of progress in `tasks.md`.
