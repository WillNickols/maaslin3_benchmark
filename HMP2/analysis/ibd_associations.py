from anadama2 import Workflow
import os
import itertools
import copy

workflow = Workflow(version="0.1", description="Get IBD associations")
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

scripts = [file for file in os.listdir(output + 'run_scripts/') if file.endswith('.R') and file.startswith('ibd')]

def compute_outputs(script, dataset, version):
    if dataset == 'taxa':
        return output + 'results/' + 'v' + str(version) + '_' + script.rstrip('.R') + '.tsv'
    else:
        return output + 'results/mbx_' + 'v' + str(version) + '_' + script.rstrip('.R') + '.tsv'

for dataset in ['taxa']: # 'mbx'
    for version in [3,4]:
        for script in scripts:
            new_command = 'Rscript run_scripts/' + script + ' --nCores ' + str(cores) + ' --workingDirectory ' + args.workingDirectory + ' --analysisDirectory ' + output + ' --dataset ' + dataset + ' --version ' + str(version)

            if not os.path.exists(compute_outputs(script, dataset, version)):
                workflow.add_task_gridable(actions=new_command,
                depends=[output + 'data/hmp2_metadata_2018-08-20.csv', output + 'data/metaphlan3_taxonomic_profiles.tsv', output + 'data/metaphlan4_taxonomic_profiles.tsv', output + 'data/annotations_hmp2.csv', output + 'data/intensities_hmp2.csv'],
                targets=compute_outputs(script, dataset, version),
                time=time,
                mem=memory,
                cores=cores,
                partition=partition
                )

workflow.go()
