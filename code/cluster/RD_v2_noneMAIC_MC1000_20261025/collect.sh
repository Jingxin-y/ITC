#!/bin/bash
set -euo pipefail
source "$RDNM_HOME/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RDNM_HOME"
sha256sum --quiet -c SHA256SUMS
Rscript --vanilla collect.R "$RDNM_HOME" "$RDNM_OUT"
