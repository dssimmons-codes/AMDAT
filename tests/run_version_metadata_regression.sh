#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

header="src/generated/version.h"
expected_version="$(tr -d '[:space:]' < VERSION)"
expected_commit="$(git rev-parse --short HEAD)"
expected_describe="$(git describe --tags --dirty --always)"

verify_header()
{
  grep -Fq "SEMVER=\"$expected_version\";" "$header"
  grep -Fq "GIT_COMMIT=\"$expected_commit\";" "$header"
  grep -Fq "GIT_DESCRIBE=\"$expected_describe\";" "$header"
}

make "$header"
verify_header

# Simulate a header left behind from an older Git/version state. Removing the
# current metadata stamp is equivalent to Make observing a new state key.
printf '#pragma once\n// deliberately stale metadata\n' > "$header"
rm -f src/generated/.version-*
make "$header"
verify_header

make clean
test ! -e "$header"
make "$header"
verify_header
