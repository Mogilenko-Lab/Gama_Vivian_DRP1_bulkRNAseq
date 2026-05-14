# DEPRECATION NOTICE

## Status: DEPRECATED (2026-01-13)

This directory (`Explorer/`) and its contents have been **deprecated** and replaced.

### Deprecated File
- `DRP1_Pathway_Explorer.html` - Old Streamlit-based pathway explorer

### Replacement
**New interactive tool:**
- **Location:** `03_Results/02_Analysis/Plots/Trajectory_Flow/interactive_bump_dashboard.html`
- **Type:** Self-contained HTML dashboard (no server required)
- **Size:** Optimized for performance
- **Features:** Enhanced filtering, pattern exploration, trajectory visualization

### Why Deprecated?

1. **Performance:** Old explorer required Streamlit server, new tool is standalone HTML
2. **Usability:** New dashboard has improved UI/UX with better filtering
3. **Maintenance:** Consolidated into single interactive visualization system
4. **Documentation:** New tool better integrated with paper figures

### Migration Guide

**Old reference:**
```
See interactive explorer at 03_Results/02_Analysis/Explorer/DRP1_Pathway_Explorer.html
```

**New reference:**
```
See interactive dashboard at 03_Results/02_Analysis/Plots/Trajectory_Flow/interactive_bump_dashboard.html
```

### Scripts Affected

- `02_Analysis/Supp5.prepare_explorer_data.py` - Data preparation (deprecated)
- `02_Analysis/Supp6.app_bump_chart_explorer.py` - Streamlit app (deprecated)

**Replacement scripts:**
- `02_Analysis/3.8.viz_interactive_bump_dashboard.py` - Generates new dashboard
- `01_Scripts/Python/viz_bump_charts.py` - Visualization module

---

**For questions:** See `03_Results/02_Analysis/Plots/Trajectory_Flow/README.md`
