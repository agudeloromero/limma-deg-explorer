#!/usr/bin/env bash
# Checks that the version is the same everywhere it is recorded,
# and (optionally) that it matches a release tag such as v1.0.0.
#   bash .github/scripts/check_version.sh            # consistency only
#   bash .github/scripts/check_version.sh v1.0.0     # also compare with the tag
set -euo pipefail
cd "$(dirname "$0")/../.."

r_ver=$(sed -n 's/^APP_VERSION <- "\(.*\)"$/\1/p' Shiny_app/R/helpers.R)
py_ver=$(sed -n 's/^APP_VERSION = "\(.*\)"$/\1/p' Shiny_app/python/report.py)
cff_ver=$(sed -n 's/^version: *//p' CITATION.cff | tr -d '"')
log_ver=$(sed -n 's/^## \[\([0-9][^]]*\)\].*/\1/p' CHANGELOG.md | head -n 1)

printf '%-30s %s\n' "Shiny_app/R/helpers.R" "$r_ver" \
  "Shiny_app/python/report.py" "$py_ver" "CITATION.cff" "$cff_ver" \
  "CHANGELOG.md (latest release)" "$log_ver"

status=0
for v in "$py_ver" "$cff_ver" "$log_ver"; do
  [[ "$v" == "$r_ver" ]] || status=1
done
[[ -n "$r_ver" ]] || status=1
if [[ $# -gt 0 ]]; then
  printf '%-30s %s\n' "Release tag" "$1"
  [[ "${1#v}" == "$r_ver" ]] || status=1
fi

if [[ $status -ne 0 ]]; then
  echo "::error::Version numbers do not match. Update them all to the same value."
  exit 1
fi
echo "Version $r_ver is consistent."
