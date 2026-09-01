#!/usr/bin/env bash
# Run by Snakemake right after the peddy conda env is created.
#
# peddy 0.4.8 (2021) is the last release and reads its bundled 1000G PCA
# background with np.fromstring(<bytes>, dtype=...). NumPy 2.0 removed the
# binary mode of fromstring, so on a current solve peddy dies with
#
#   ValueError: The binary mode of fromstring is removed, use frombuffer instead
#
# in het_check -> pca(). The failure is nasty because ped_check.csv is written
# first and looks fine -- only sex_check.csv, het_check.csv,
# background_pca.json and the HTML report go missing.
#
# np.frombuffer is the documented drop-in replacement; the result is immediately
# .astype(np.int32)'d into a writable copy, so read-only buffer semantics are
# not a problem here.
set -euo pipefail

pca_py="$(python -c 'import os.path, peddy; print(os.path.join(os.path.dirname(peddy.__file__), "pca.py"))')"

if grep -q 'np\.fromstring(gzip\.open(f' "$pca_py"; then
    sed -i "s/np\.fromstring(gzip\.open(f, 'rb')\.read()/np.frombuffer(gzip.open(f, 'rb').read()/" "$pca_py"
    echo "peddy.post-deploy.sh: patched np.fromstring -> np.frombuffer in $pca_py"
fi

# Fail env creation loudly rather than at job time if the patch did not take.
python -c "
import inspect, peddy.pca
assert 'np.fromstring' not in inspect.getsource(peddy.pca), 'peddy pca.py still uses np.fromstring'
"
