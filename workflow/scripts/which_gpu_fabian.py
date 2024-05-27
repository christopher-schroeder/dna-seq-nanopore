#!/usr/bin/env python

from io import StringIO
import os
import pandas as pd
import subprocess
import sys
import time


try:
    n = int(sys.argv[1])
except:
    n = 1


try:
    gpu_txt = sys.argv[2]
except:
    gpu_txt = None


## print("Requested GPUs:", n)



## print("{:#^60}".format(" Running jobs "))

jobs = []

if gpu_txt:
    ## print(f"gpu_txt: '{gpu_txt}'")

    if os.path.isfile(gpu_txt):
        ##print(f"Reading running jobs from: '{gpu_txt}'")

        with open(gpu_txt, "r") as f:
            #jobs = [line.strip() for line in f.readlines()]
            jobs = []
            for line in f:
                line = line.strip()
                if line:
                    job = jobs.append(line)
    else:
        ## print(f"Saving to NEW file:, '{gpu_txt}'")
        pass


if jobs:
    ## print("jobs:", jobs, len(jobs))
    pass



in_use = []
for job in jobs:
    in_use.extend(job.split(":")[0].split(","))
in_use = [ int(x.strip()) for x in in_use ]
# in_use = sorted(list(set(in_use)))
##print("GPUs in use:", in_use, len(in_use))


if len(jobs) > 0:
    time.sleep(5)


##print("{:#^60}".format(" GPUs "))

##cmdl = "nvidia-smi --format=csv --query-gpu=index,utilization.gpu,memory.used"
cmdl = "nvidia-smi --format=csv --query-gpu=index,utilization.gpu,utilization.memory"
output = subprocess.check_output(cmdl.split()).decode()
## print(output)

df = pd.read_csv(StringIO(output), sep=",")
df.columns = [column.strip() for column in df.columns] # strip whitespaces from column names
##df = df.rename(columns={'index': 'gpu'})
df['utilization.gpu [%]'] = pd.to_numeric(df['utilization.gpu [%]'].str.replace(' %',''))
df['utilization.memory [%]'] = pd.to_numeric(df['utilization.memory [%]'].str.replace(' %',''))
df['utilization.memory_available [%]'] = 100 - df['utilization.memory [%]']
##print(df.to_string(index=False))
##print(df.dtypes)

##print("{:#^60}".format(" GPUs by usage "))
df = df.sort_values(by=['utilization.gpu [%]','utilization.memory [%]'], ascending=False)
##print(df.to_string(index=False))

##print("{:#^60}".format(" GPUs by usage (not used) "))
df = df[~df['index'].isin(in_use)]
##print(df.to_string(index=False))

##print("{:#^60}".format(f" GPUs by usage (not used, n={n}) "))
##print(df.tail(n))

memory_min = df.tail(n)["utilization.memory_available [%]"].min()
##print("memory_min:", memory_min)


if memory_min > 50:
    gpu_str = ",".join( sorted([str(x) for x in df.tail(n).index.to_list()]) ) #+ ":" + str(memory_min) + "%"
    print(gpu_str)

    if gpu_str and gpu_txt:
        with open(gpu_txt, "a") as f:
            f.write(gpu_str + "\n")

