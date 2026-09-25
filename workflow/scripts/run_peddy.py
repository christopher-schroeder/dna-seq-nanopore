"""Run peddy's CLI with a NumPy 2 compatibility shim, then hand over to click.

peddy 0.4.8 (2021) is the last release and reads its bundled 1000G PCA
background with np.fromstring(<bytes>, dtype=...).  NumPy 2.0 removed the
binary mode of fromstring, so het_check -> pca() dies with

    ValueError: The binary mode of fromstring is removed, use frombuffer instead

The failure is nasty because ped_check runs first and writes a plausible
ped_check.csv -- only sex_check.csv, het_check.csv, background_pca.json and the
HTML report go missing.

Two fixes that do *not* work here:

  * Pinning numpy <2 in envs/peddy.yaml: every available cyvcf2 build is
    compiled against the numpy 2 ABI.
  * Patching site-packages from envs/peddy.post-deploy.sh: `snakemake
    --containerize` emits nothing but `conda env create` per environment, so
    post-deploy scripts never run when the envs are baked into container.sif.
    (Do not "clean up" that script -- its bytes feed the conda env hash that
    names /conda-envs/<hash> inside the image, so touching it invalidates the
    baked environment and forces a container rebuild.)

Shimming at call time keeps the fix with the rule that needs it and works
identically for a containerised run and a plain --use-conda run.  np.frombuffer
is the documented drop-in replacement; peddy immediately .astype(np.int32)'s the
result into a writable copy, so the read-only buffer semantics do not matter.
"""

import sys

import numpy as np

_np_fromstring = np.fromstring


def _fromstring(string, *args, **kwargs):
    """np.fromstring, but routed to frombuffer for the removed binary mode."""
    if isinstance(string, (bytes, bytearray, memoryview)) and "sep" not in kwargs:
        return np.frombuffer(string, *args, **kwargs)
    return _np_fromstring(string, *args, **kwargs)


np.fromstring = _fromstring

from peddy.__main__ import cli  # noqa: E402  -- must be imported after the shim

if __name__ == "__main__":
    sys.exit(cli())
