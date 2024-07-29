from anadama2 import Workflow
import os
import itertools
import copy

workflow = Workflow(version="0.1", description="Run randomization procedure")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The memory in megabytes allocated to run the command", type=int, default=10000)
workflow.add_argument("time", desc="The time in minutes allocated to run the command", type=int, default=120)
workflow.add_argument('tmp', desc="Whether to use reduced parameters", action="store_true")
args = workflow.parse_args()
this_directory = str(os.path.dirname(os.path.realpath(__file__))).rstrip('/') + '/'

# output, set this to the GitHub repository
output = os.path.abspath(args.output.rstrip("randomization_test/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

# scratch directory
scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/"

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
time = args.time

if args.tmp:
    nIterations = 5
else:
    nIterations = 100

#############
# Run tools #
#############

run_tools_directory = this_directory + 'run_tools/'
working_directory = output

def compute_running_outputs(iteration, dataset):
    generation_outputs = [output + 'randomization_test/associations/' + str(dataset) + '_' + str(iteration) + '.tsv']
    return generation_outputs

tools = [file for file in os.listdir(run_tools_directory) if file.startswith('run_null_')]
tool = tools[0]

datasets = [file.split('_metadata.tsv')[0].split('_ASVs_table.tsv')[0] for file in os.listdir(this_directory + 'raw_data/')]
datasets = set(datasets)

for dataset in datasets:
    for iteration in range(nIterations):
        iteration = iteration + 1
        new_command = 'Rscript ' + run_tools_directory + tool + \
        ' --workingDirectory ' + working_directory + \
        ' --dataset ' + str(dataset) + \
        ' --iter ' + str(iteration)

        if not all(os.path.exists(file_path) for file_path in compute_running_outputs(iteration, dataset)):
            workflow.add_task_gridable(actions=new_command,
            targets=compute_running_outputs(iteration, dataset),
            time=time,
            mem=memory,
            cores=cores,
            partition=partition
            )

tools = [file for file in os.listdir(run_tools_directory) if file.startswith('run_nonnull_')]
tool = tools[0]

for dataset in datasets:

    new_command = 'Rscript ' + run_tools_directory + tool + \
    ' --workingDirectory ' + working_directory + \
    ' --dataset ' + str(dataset)

    if not all(os.path.exists(file_path) for file_path in [output + 'randomization_test/associations/' + dataset + '_non_null.tsv']):
        workflow.add_task_gridable(actions=new_command,
        targets=output + 'randomization_test/associations/' + dataset + '_non_null.tsv',
        time=time,
        mem=memory,
        cores=cores,
        partition=partition
        )

workflow.go()
