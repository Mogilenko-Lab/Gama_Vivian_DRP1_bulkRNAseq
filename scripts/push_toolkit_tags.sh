#!/usr/bin/env bash
# One-time maintainer script — run OUTSIDE this dev container, where you have
# SSH/gh access to the toolkit repos. Creates and pushes annotated tags at the
# exact commits behind the published v2.1.1 results, so a fresh clone can always
# fetch them (a bare SHA not reachable from any ref may fail to fetch).
#
# Safe to re-run: tag creation is skipped if the tag already exists.
set -euo pipefail

TAG="drp1-paper-v2.1.1"
MSG="Code state behind Gama/Vivian DRP1 bulkRNAseq paper (repo tag v2.1.1)"

declare -A PIN=(
  ["01_Scripts/RNAseq-toolkit"]="532982df1fd04ce9cc0b1717c75721dcdfe846c7"
  ["01_Scripts/SciAgent-toolkit"]="fb6012a731fa0aaaef417932960f8d923f7caa43"
)

for path in "${!PIN[@]}"; do
  sha="${PIN[$path]}"
  echo "== $path @ $sha =="
  git -C "$path" fetch origin
  if ! git -C "$path" cat-file -e "${sha}^{commit}" 2>/dev/null; then
    echo "FATAL: commit $sha not found in $path — is the local clone complete?" >&2
    exit 1
  fi
  if git -C "$path" rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
    echo "  tag $TAG already exists locally — skipping create"
  else
    git -C "$path" tag -a "$TAG" "$sha" -m "$MSG"
    echo "  created $TAG -> $sha"
  fi
  git -C "$path" push origin "refs/tags/$TAG"
  echo "  pushed $TAG to origin"
done

echo
echo "Done. Verify a third party can fetch the pins from scratch:"
echo "  git clone --recurse-submodules https://github.com/MogilenkoLabVUMC/Gama_Vivian_DRP1_bulkRNAseq.git /tmp/repro-test"
