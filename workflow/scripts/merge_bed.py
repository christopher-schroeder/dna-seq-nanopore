filename = "/projects/humgen/science/depienne/project818/dna-seq-nanopore/results/methylation/818-13_1.bed"
last_split = None
for i, curr in enumerate(open(filename)):
    curr_split = curr.strip().split("\t")
    if i % 2 == 0:
        last_split = curr_split
    else:
        assert(curr_split[0] == curr_split[0])
        assert(curr_split[1] == curr_split[1])
        assert(curr_split[2] == curr_split[2])
        last_values = last_split[-1].split(" ")
        curr_values = curr_split[-1].split(" ")
        print("last", last_values)
        print("curr", curr_values)