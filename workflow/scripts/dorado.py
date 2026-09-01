"""Run dorado on a GPU claimed with a lock file.

Only needed for local (non-SLURM) basecalling; under SLURM the scheduler hands
out the device and rules/basecalling.smk calls dorado directly.
"""

import os
import random
import subprocess
import sys
import time
from pathlib import Path

import GPUtil

LOCK_DIR = Path("/local/tmp")
DORADO = "/projects/humgen/pipelines/dna-seq-nanopore/workflow/tools/dorado-2.0.0-linux-x64/bin/dorado"


def claim_gpu(poll_seconds=30):
    """Claim a free GPU by creating its lock file, waiting until one frees up."""
    while True:
        available = GPUtil.getAvailable(order="memory", includeNan=False, limit=6)
        for gid in available:
            lock = LOCK_DIR / f"LCK_gpu_{gid}.lock"
            try:
                # O_EXCL makes the check and the claim a single atomic step, so
                # two jobs starting together cannot both take the same GPU
                fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError:
                continue
            os.close(fd)
            return gid, lock
        print("no free GPU, retrying", file=sys.stderr)
        time.sleep(poll_seconds)


# stagger concurrent jobs so they do not all inspect the GPUs at once
time.sleep(random.randint(2, 10))

gpu_id, lock = claim_gpu()
print("allocated gpu", gpu_id, file=sys.stderr)

try:
    command = " ".join([
        DORADO, "basecaller",
        str(snakemake.params.model),
        str(snakemake.input.fast5_dir),
        str(snakemake.params.remora_args),
        str(snakemake.params.basecaller_args),
        "--recursive",
        f"--device cuda:{gpu_id}",
    ])
    print(command, file=sys.stderr)
    with open(snakemake.output.ubam, "wb") as out:
        # check=True: a silent failure here would otherwise leave an empty ubam
        # that the workflow accepts as a finished basecall
        subprocess.run(command, shell=True, stdout=out, check=True)
finally:
    lock.unlink(missing_ok=True)
