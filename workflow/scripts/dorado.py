from contextlib import contextmanager
from itertools import cycle
import GPUtil
import os
from pathlib import Path
import time
import random

allocated_gid = None
# avail_gpus = GPUtil.getAvailable(order='memory', maxLoad=.1, maxMemory=.1, includeNan=False, limit=6)
avail_gpus = GPUtil.getAvailable(order='memory', includeNan=False, limit=6)


amount = random.randint(2,10)
print(f"waiting for {amount} seconds")
time.sleep(amount)


while not allocated_gid:
    for gid in cycle(avail_gpus):
        # then we've successfully created the lockfile
        if os.path.isfile(f"/local/tmp/LCK_gpu_{gid}.lock"):
            continue
        Path(f"/local/tmp/LCK_gpu_{gid}.lock").touch()
        allocated_gid = gid
        print("allocated_gid", allocated_gid)
        break

try:
    command = """
        /projects/humgen/pipelines/dna-seq-nanopore/workflow/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller \
        {params.model} \
        {input.fast5_dir} \
        {params.remora_args} \
        {params.basecaller_args} \
        --recursive \
        --device cuda:{gpu_id} \
        > {output.ubam}
    """.format(params=snakemake.params, input=snakemake.input, output=snakemake.output, gpu_id=allocated_gid)
    print(command)
    os.system(command)
finally:
    os.system(f"rm /local/tmp/LCK_gpu_{allocated_gid}.lock")