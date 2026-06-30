# Provenance & Reproducibility Lock

Single source of truth for the exact code, tools, and environment behind the
published results. Every numerical result and figure in the manuscript was
produced from the pinned commits below.

## Paper release

| Item | Value |
|------|-------|
| Repository tag | **`v2.1.1`** (commit `b1fa558`) — post-revision final |
| GitHub | `MogilenkoLabVUMC/Gama_Vivian_DRP1_bulkRNAseq` |
| Zenodo DOI | [10.5281/zenodo.20213748](https://doi.org/10.5281/zenodo.20213748) |

## Pinned submodules (the lock)

Both submodules are pinned by exact commit (the gitlink recorded in this repo's
tree). These are the commits that generated the published figures and tables —
**not** later development tips.

| Submodule | Commit (pinned) | Tag | Affects results? |
|-----------|-----------------|-----|------------------|
| `01_Scripts/RNAseq-toolkit` | `532982df1fd04ce9cc0b1717c75721dcdfe846c7` | `drp1-paper-v2.1.1` (to create) / `v0.1.0-beta-42-g532982d` | **Yes** — DE, GSEA, volcano, IO helpers |
| `01_Scripts/SciAgent-toolkit` | `fb6012a731fa0aaaef417932960f8d923f7caa43` | `drp1-paper-v2.1.1` (to create) | **No** — AI-agent tooling (skills/prompts); not on the results path, pinned for full provenance only |

> **Do not** run `git submodule update --remote` or "bump to latest" on these.
> Doing so silently moves the analysis onto untested toolkit code. The pin is
> the gitlink commit, not a branch.

### Verifying the lock

```bash
git submodule status
#  532982df... 01_Scripts/RNAseq-toolkit  (v0.1.0-beta-42-g532982d)
#  fb6012a7... 01_Scripts/SciAgent-toolkit
```

## Computational environment

| Layer | Pin |
|-------|-----|
| Dev container | `scdock-r-dev:v0.5.1` ([scbio-docker v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)) |
| R packages | `R_session_info.txt` (full `sessionInfo`), `R_packages.txt` |
| Python packages | `python_requirements_freeze.txt` (exact freeze), `requirements.txt` (core) |
| Upstream alignment/quant | [bulkRNAseq_pipeline_scripts](https://github.com/tony-zhelonkin/bulkRNAseq_pipeline_scripts) on [scbio-docker v0.2.0](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.2.0) |
| Reference genome | GRCh38.p14 (GCF_000001405.40) |

## Reproduce from a clean clone

```bash
git clone --recurse-submodules \
  https://github.com/MogilenkoLabVUMC/Gama_Vivian_DRP1_bulkRNAseq.git
cd Gama_Vivian_DRP1_bulkRNAseq
git checkout v2.1.1
git submodule update --init --recursive   # lands on the pinned commits above

# Open in VS Code → "Dev Containers: Reopen in Container" (scdock-r-dev:v0.5.1)
pip install -r python_requirements_freeze.txt   # inside container

# Core pipeline (see README "Reproduce" section for the full ordered list)
Rscript 02_Analysis/1.1.main_pipeline.R
```

## Provenance history (why the lock matters)

Local `main` briefly carried four un-pushed commits that bumped the submodules
to later dev tips (`RNAseq f6f8aff`, `SciAgent 9149092`) — those tips were never
used for the paper. The `repro-lock-v2.1.1` branch reverts the gitlinks back to
the `v2.1.1` published state recorded above. `origin/main` always pointed at the
correct `v2.1.1` pins.
