#!/bin/bash
set -euo pipefail
[[ "$#" == 0 ]] || { echo 'Usage: bash submit_auto.sh (all eight none+MAIC cases)'; exit 1; }
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$here/submit.sh"
bash "$here/auto/submit_auto_none.sh"
