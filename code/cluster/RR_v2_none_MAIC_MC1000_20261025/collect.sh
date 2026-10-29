#!/bin/bash
set -euo pipefail
source "$RRNMW_HOME/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RRNMW_HOME"
sha256sum --quiet -c SHA256SUMS
Rscript --vanilla collect.R "$RRNMW_HOME" "$RRNMW_OUT"
