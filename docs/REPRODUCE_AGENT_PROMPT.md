# Hand-off prompt — finalize the reproducibility lock (run OUTSIDE the dev container)

Paste the block below to an agent (Claude Code / Codex / Gemini) running on a
host that **has network access, SSH keys, and `gh` authenticated** for both
`github.com/tony-zhelonkin` (toolkits) and
`github.com/Mogilenko-Lab` (paper repo). The in-container session prepared the
`repro-lock-v2.1.1` branch but cannot reach the network to push tags or test a
clean clone.

---

You are finalizing the reproducibility lock for the `Gama_Vivian_DRP1_bulkRNAseq`
paper repo. The branch `repro-lock-v2.1.1` already pins both toolkit submodules to
the exact commits behind the published `v2.1.1` results. Your job: make those pins
durably fetchable by anyone, then prove a clean clone reproduces.

**Pinned commits (do not change these):**
- `RNAseq-toolkit`  → `532982df1fd04ce9cc0b1717c75721dcdfe846c7`
- `SciAgent-toolkit` → `fb6012a731fa0aaaef417932960f8d923f7caa43`

**Tasks, in order:**

1. **Push the paper branch.** From the repo, `git push -u origin repro-lock-v2.1.1`.
   Review the diff vs `origin/main` first — it should only revert the submodule
   gitlinks, switch `.gitmodules` to HTTPS, and add `PROVENANCE.md` /
   `scripts/push_toolkit_tags.sh` / README edits. Open a PR into `main` and, after
   review, merge so `main` no longer carries the un-pushed "bump to latest dev"
   drift commits.

2. **Tag the toolkit commits.** Run `scripts/push_toolkit_tags.sh` (creates and
   pushes annotated tag `drp1-paper-v2.1.1` at each pinned commit). Confirm with
   `git -C 01_Scripts/RNAseq-toolkit ls-remote --tags origin | grep drp1-paper`.

3. **Check visibility — reproducibility blocker.** A stranger must be able to fetch
   the toolkits. For each toolkit run `gh repo view tony-zhelonkin/<repo> --json visibility`.
   - If either is **private**, the clone test below WILL fail for outsiders. Either
     make it public (`gh repo edit tony-zhelonkin/<repo> --visibility public`,
     confirm with the owner first) **or** plan a Zenodo source-archive of the
     toolkits and record that DOI in `PROVENANCE.md`. Report which path you took.

4. **Clean-clone reproducibility test (the real proof).** In a fresh directory with
   NO existing credentials cached for the toolkits:
   ```bash
   git clone --recurse-submodules \
     https://github.com/Mogilenko-Lab/Gama_Vivian_DRP1_bulkRNAseq.git /tmp/repro-test
   cd /tmp/repro-test && git checkout v2.1.1 && git submodule update --init --recursive
   git submodule status   # MUST show 532982d... and fb6012a...
   ```
   Assert the two SHAs match exactly. If `git submodule update` errors with
   "reference is not a tree" / cannot fetch, the tags from step 2 are missing or the
   repo is private — go back and fix.

5. **Run the pipeline in the pinned container.** Build/open `scdock-r-dev:v0.5.1`
   (scbio-docker v0.5.1), then inside it:
   ```bash
   pip install -r python_requirements_freeze.txt
   pytest -q                                   # conftest.py present
   Rscript 02_Analysis/1.1.main_pipeline.R     # core DE+GSEA pipeline
   ```
   Compare regenerated outputs against the committed checkpoints in
   `03_Results/02_Analysis/checkpoints/` and the master tables
   (`master_gsea_table.csv`, `master_gsva_*.csv`). Spot-check a few DE/GSEA values.

6. **Report.** Summarize: tags pushed (✓/✗), toolkit visibility, clone-test SHA
   match, pytest result, and whether key numbers reproduced. Flag any drift.

Do not bump submodules to newer tips, and do not `git submodule update --remote`.
The pin is the point.
