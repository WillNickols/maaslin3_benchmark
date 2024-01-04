from anadama2 import Workflow
import os

workflow = Workflow(version="0.1", description="MPA 4 workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The memory in megabytes allocated to run the command", type=int, default=10000)
workflow.add_argument("time", desc="The time in minutes allocated to run the command", type=int, default=120)
workflow.add_argument('parameters', type=str, desc="Specify the parameters")
args = workflow.parse_args()
this_directory = os.path.dirname(os.path.realpath(__file__))

# This workflow should:
# Create everything on the benchmark options spreadsheet
# Run the three methods on each
# Run the analysis scripts on those results
# 
# Can run multiple simulation parameters on multiple cores (yes) and can run multiple iterations of a single simulation on multiple cores (no)
# Can run multiple simulation types or evaluations with this script

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

#import parameters from
with open(str(args.parameters), 'r') as file:
    lines = file.readlines()

param_dict = {}
for line in lines:
    parts = line.strip().split(':')
    if len(parts) == 2:
        key = parts[0].strip()
        value = parts[1].strip().split()
        param_dict[key] = value
    else:
        raise ValueError("Wrong input format")

def compute_outputs(generator, step):
    if step == 'generate':
        return output + 'Input/general_evaluations/' + generator + '/done'
    if step == 'run_tools':
        return output + 'Output/general_evaluations/' + generator + '/done'

for generator in param_dict['generators']:
    generate_command = '''{a} && {b}'''.format(
        a = 'python3 general_evaluations/data_generation/' + generator + '_workflow.py --parameters general_evaluations/data_generation/' + generator + '.txt --working-directory ' + output + ' --cores ' + str(cores),
        b = 'touch [targets[0]]'
        )
    
    workflow.add_task_gridable(actions=generate_command,
        targets=compute_outputs(generator, 'generate'),
        time=time,
        mem=memory,
        cores=cores,
        partition=partition
        )

run_tools_command = '''{a} && {b}'''.format(
    a = 'python3 general_evaluations/run_tools/workflow.py --generators ' + ','.join(param_dict['generators']) + ' --working-directory ' + output + ' --cores ' + str(cores),
    b = 'touch [targets[0]]'
    )

workflow.add_task_gridable(actions=run_tools_command,
    depends=compute_outputs(generator, 'generate'),
    targets=compute_outputs(generator, 'run_tools'),
    time=time,
    mem=memory,
    cores=cores,
    partition=partition
    )

workflow.go()
