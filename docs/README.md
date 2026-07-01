# Documentation Index

This directory contains comprehensive documentation for the DRP1 bulk RNA-seq analysis project.

**Last Updated:** 2026-06-30

---

## Quick Reference

| Need | Documentation |
|------|---------------|
| Project overview | [../README.md](../README.md) |
| Setup instructions | [SETUP.md](SETUP.md) or [INSTALL.md](INSTALL.md) |
| Analysis architecture | [CLAUDE.md](CLAUDE.md) |
| Script reference | [../02_Analysis/SCRIPTS.md](../02_Analysis/SCRIPTS.md) |
| **Pattern classification** | [PATTERN_CLASSIFICATION.md](PATTERN_CLASSIFICATION.md) |
| Scientific background | [bio_notes.md](bio_notes.md) |
| Reproducibility lock | [../PROVENANCE.md](../PROVENANCE.md) |
| Data lineage | [DATA_LINEAGE.md](DATA_LINEAGE.md) |
| Migration notes | [MIGRATION.md](MIGRATION.md) |

---

## Core Documentation

### Pattern Classification System

- **[PATTERN_CLASSIFICATION.md](PATTERN_CLASSIFICATION.md)** - Canonical reference for pattern classification:
  - 8-pattern taxonomy (Compensation, Sign_reversal, Progressive, etc.)
  - Super-category system for simplified reporting
  - Threshold definitions and biological interpretation
  - Algorithm pseudocode
  - Python/R implementation synchronization

**Implementation files:**
- `01_Scripts/Python/pattern_definitions.py` - Python canonical source
- `02_Analysis/1.7.create_master_gsva_table.R` - R implementation

### Scientific Context

- **[bio_notes.md](bio_notes.md)** - Curated biological context:
  - DRP1 mutations and clinical phenotype
  - Mitochondrial dynamics and bioenergetics
  - Synaptic ribosomes and local protein synthesis
  - Energy-translation crisis model
  - Calcium dysregulation, epilepsy mechanisms, and critical period plasticity

### Setup and Installation

- **[SETUP.md](SETUP.md)** - Quick setup guide:
  - Container setup via VS Code Dev Containers
  - Environment verification
  - Quick start commands

- **[INSTALL.md](INSTALL.md)** - Detailed installation instructions:
  - Prerequisites and system requirements
  - Step-by-step container setup
  - Troubleshooting common issues

### Analysis Architecture

- **[CLAUDE.md](CLAUDE.md)** - Claude Code instructions and architecture:
  - Repository structure and pipeline overview
  - Essential commands and workflows
  - Pattern classification system synchronization
  - Common tasks and troubleshooting

### Project Documentation

- **[../PROVENANCE.md](../PROVENANCE.md)** - Reproducibility lock: pinned submodule commits, container, environment, and clean-clone reproduce steps
- **[MIGRATION.md](MIGRATION.md)** - Technical migration notes from original workstation setup
- **[DATA_LINEAGE.md](DATA_LINEAGE.md)** - Data provenance and processing history

> Version history lives in the git tags and commit log (`git tag`, `git log`); the
> per-release Zenodo archive is linked from [../PROVENANCE.md](../PROVENANCE.md).

---

## Getting Started

### For New Users

1. Start with [../README.md](../README.md) - Project overview and key findings
2. Then [SETUP.md](SETUP.md) - Quick setup guide
3. Optionally [INSTALL.md](INSTALL.md) - Detailed installation instructions

### For Developers

1. Start with [CLAUDE.md](CLAUDE.md) - Claude Code instructions and architecture
2. Reference [../02_Analysis/SCRIPTS.md](../02_Analysis/SCRIPTS.md) - Complete script inventory
3. Review [PATTERN_CLASSIFICATION.md](PATTERN_CLASSIFICATION.md) - Pattern system documentation
4. Check [../01_Scripts/Python/README.md](../01_Scripts/Python/README.md) - Python module architecture

### For Manuscript Reviewers

1. Start with [../README.md](../README.md) - Key findings and overview
2. Explore [../03_Results/02_Analysis/Plots/README.md](../03_Results/02_Analysis/Plots/README.md) - Visualization guide
3. Reference [bio_notes.md](bio_notes.md) - Detailed mechanistic context
4. See [../PROVENANCE.md](../PROVENANCE.md) - Exact code, versions, and environment behind the published figures

---

## Analysis Documentation

### Results Overview

- [../03_Results/02_Analysis/README.md](../03_Results/02_Analysis/README.md) - Results overview and master tables:
  - Master GSEA table (109K pathway enrichments with patterns)
  - Master GSVA tables (focused 7 modules + comprehensive all pathways)
  - Pattern storage and querying

### Plot Documentation

Individual plot folders contain detailed README files explaining:
- Generating scripts and dependencies
- Methods and parameters
- Figure descriptions
- Regeneration commands

**Key visualization folders:**
- [Publication_Figures/](../03_Results/02_Analysis/Plots/Publication_Figures/) - Manuscript-ready figures
- [Trajectory_Flow/](../03_Results/02_Analysis/Plots/Trajectory_Flow/) - Bump charts and alluvial diagrams
- [Critical_period_trajectories/](../03_Results/02_Analysis/Plots/Critical_period_trajectories/) - GSVA trajectory analysis

---

## Version History

- **`v1.0`** (commit `d6ec164`) — original manuscript submission (frozen).
- **`v2.1.1`** (commit `b1fa558`) — post-revision final; the published state. See
  [../PROVENANCE.md](../PROVENANCE.md) for the reproducibility lock (pinned toolkit
  commits, container, environment).
- Deprecated plots from earlier revisions: `../03_Results/02_Analysis/Plots/.deprecated/`

Full history is in the git tags and commit log (`git tag`, `git log`).

---

## Documentation Map

```
docs/
├── README.md                      # This index
├── CLAUDE.md                      # Claude Code instructions and architecture
├── SETUP.md                       # Setup pointer → PROVENANCE.md + INSTALL.md
├── INSTALL.md                     # Detailed installation instructions
├── PATTERN_CLASSIFICATION.md      # Pattern system (canonical)
├── bio_notes.md                   # Curated biological notes
├── MIGRATION.md                   # Migration notes
├── DATA_LINEAGE.md                # Data provenance
└── REPRODUCE_AGENT_PROMPT.md      # Clean-clone reproduce hand-off

Root documentation:
├── README.md                      # Project overview
└── PROVENANCE.md                  # Reproducibility lock (pins, container, env)
```

---

## Related Documentation

- **[../01_Scripts/Python/README.md](../01_Scripts/Python/README.md)** - Python module architecture
- **[../02_Analysis/SCRIPTS.md](../02_Analysis/SCRIPTS.md)** - Script inventory and usage
- **[../03_Results/02_Analysis/Plots/README.md](../03_Results/02_Analysis/Plots/README.md)** - Plot folder overview
