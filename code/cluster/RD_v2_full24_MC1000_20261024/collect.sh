#!/bin/bash
set -euo pipefail
source "$RDFULL_HOME/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RDFULL_HOME"
sha256sum --quiet -c SHA256SUMS
Rscript --vanilla collect.R "$RDFULL_HOME" "$RDFULL_OUT"
