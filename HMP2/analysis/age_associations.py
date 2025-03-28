from anadama2 import Workflow
import os
import itertools
import copy

workflow = Workflow(version="0.1", description="Age associations")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The memory in megabytes allocated to run the command", type=int, default=10000)
workflow.add_argument("time", desc="The time in minutes allocated to run the command", type=int, default=120)
workflow.add_argument("workingDirectory", desc="The directory with Maaslin3", type=str)
args = workflow.parse_args()
this_directory = str(os.path.dirname(os.path.realpath(__file__))).rstrip('/') + '/'

# output, set this to the GitHub repository
output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

# scratch directory
scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/"

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
time = args.time

def compute_outputs(dataset, version, mversion):
    if dataset == 'under16':
        return output + 'results/' + 'v' + str(version) + '_under16_Maaslin' + str(mversion) + '.tsv'
    else:
        return output + 'results/' + 'v' + str(version) + '_atleast16_Maaslin' + str(mversion) + '.tsv'

for mversion in [2, 3]:
    for dataset in ['under16', 'atleast16']:
        for version in [3,4]:
            if mversion == 2:
                script = 'age_associations_Maaslin2.R'
            else:
                script = 'age_associations_Maaslin3.R'

            new_command = 'Rscript run_scripts/' + script + ' --nCores ' + str(cores) + ' --workingDirectory ' + args.workingDirectory + ' --analysisDirectory ' + output + ' --dataset ' + dataset + ' --version ' + str(version)

            if not os.path.exists(compute_outputs(dataset, version, mversion)):
                workflow.add_task_gridable(actions=new_command,
                depends=[output + 'data/hmp2_metadata_2018-08-20.csv', output + 'data/metaphlan4_taxonomic_profiles.tsv', output + 'data/metaphlan3_taxonomic_profiles.tsv'],
                targets=compute_outputs(dataset, version, mversion),
                time=time,
                mem=memory,
                cores=cores,
                partition=partition
                )

workflow.go()
