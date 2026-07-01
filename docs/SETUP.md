# Setup

This file used to carry standalone setup steps that have since drifted from the
code. Setup and reproducibility now live in two authoritative places:

- **Reproduce the published analysis** → [`../PROVENANCE.md`](../PROVENANCE.md)
  — pinned toolkit commits, container image, environment freezes, and the exact
  clean-clone reproduce steps for tag `v2.1.1`.
- **Install / build the environment** → [`INSTALL.md`](INSTALL.md)
  — prerequisites, Docker image build, submodule init, Dev Container launch,
  runtime package installs, and troubleshooting.

## Quick start

```bash
# 1. Clone with submodules (pins land on the published toolkit commits)
git clone --recurse-submodules \
  https://github.com/Mogilenko-Lab/Gama_Vivian_DRP1_bulkRNAseq.git
cd Gama_Vivian_DRP1_bulkRNAseq

# 2. Open in VS Code → "Dev Containers: Reopen in Container"
#    (image: scdock-r-dev:v0.5.1 — see INSTALL.md to build it)

# 3. Run the main pipeline inside the container
Rscript 02_Analysis/1.1.main_pipeline.R
```

See [`CLAUDE.md`](CLAUDE.md) for the full command reference and architecture.
