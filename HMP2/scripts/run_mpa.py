#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil

workflow = Workflow(version="0.1", description="MPA 4 workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=45000)
workflow.add_argument(name="input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument("bowtie", desc="Metaphlan4 bowtie")
args = workflow.parse_args()
input_extension = args.input_extension

# output
output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

os.makedirs(output + "/merged", exist_ok = True)

# scratch directory
scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/"

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
max_time = args.time

# list the input fastq files
in_dir = args.input

paths = glob.glob(os.path.abspath(in_dir.rstrip("/")) + "/" + '*.' + input_extension)

files = []
for path in paths:
    files.append(path)

names = set(file.split("." + input_extension)[0] for file in files)


#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step):
    time = 1
    n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
    if step == "metaphlan4":
        time = 20 * n_gigabytes
    if time > max_time:
        time = max_time
    return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step):
    return [str(name+"."+input_extension)]

############################
# function to list targets #
############################

def list_targets(name, step):
	return [str(output + name.split("/")[-1] + "_taxonomic_profile.tsv")]

##############################
# function to run metaphlan4 #
##############################

def metaphlan4(name):
    command = '''{a}'''.format(
        a = "metaphlan " + name + "." + input_extension + " --input_type fastq --output_file " + \
        scratch + name.split("/")[-1] + "_taxonomic_profile.tsv --nproc " + str(cores) + \
        " --unclassified_estimation --bowtie2db " + args.bowtie + " --tmp_dir " + scratch
        )
    return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="metaphlan4")[0]):
		workflow.add_task_gridable(actions=metaphlan4(name),
			depends=list_depends(name=name, step="metaphlan4"),
			targets=list_targets(name=name, step="metaphlan4"),
			time=calculate_time(name=name, step="metaphlan4"),
			mem=memory,
			cores=cores,
			partition=partition
			)

merge_depends = []
for name in names:
	merge_depends.append(list_targets(name=name, step="metaphlan4")[0])

final_output = [output + "/merged/metaphlan4_taxonomic_profiles.tsv"]

command = '''{a}'''.format(
	a = "merge_metaphlan_tables.py " + output + "*.tsv > [targets[0]]"
	)

workflow.add_task_gridable(actions=command,
	depends=merge_depends,
	targets=final_output,
	time=5,
	mem=8000,
	cores=1,
	partition=partition
	)

####################
# run the workflow #
####################

workflow.go()

#
