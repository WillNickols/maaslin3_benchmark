from anadama2 import Workflow
import os
import itertools
import copy

workflow = Workflow(version="0.1", description="Diet associations")
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

foods = ['soft_drinks', 'diet_soft_drinks', 'fruit_juice', 'water', 'alcohol', 'yogurt', 'dairy', 'probiotic', 'fruits', 
        'vegetables', 'beans', 'whole_grains', 'starch', 'eggs', 'processed_meat', 'red_meat', 'white_meat', 'shellfish',
        'fish', 'sweets', 'tea']

def compute_outputs(food, ordered, dataset):
    if dataset == 'taxa':
        return output + 'results/' + ordered + "_food_associations_" + food + ".tsv"
    else:
        return output + 'results/mbx_' + ordered + "_food_associations_" + food + ".tsv"

dataset = 'taxa'
for food in foods:
    for ordered in ['ordered', 'group']:
        new_command = 'Rscript run_scripts/' + 'diet_associations_Maaslin3.R' + ' --nCores ' + str(cores) + ' --workingDirectory ' + args.workingDirectory + ' --analysisDirectory ' + output + ' --food ' + food + ' --ordered ' + ordered  + ' --dataset ' + dataset

        if not os.path.exists(compute_outputs(food, ordered, dataset)):
            workflow.add_task_gridable(actions=new_command,
            depends=[output + 'data/hmp2_metadata_2018-08-20.csv', output + 'data/metaphlan4_taxonomic_profiles.tsv', output + 'data/annotations_hmp2.csv', output + 'data/intensities_hmp2.csv'],
            targets=compute_outputs(food, ordered, dataset),
            time=time,
            mem=memory,
            cores=cores,
            partition=partition
            )

workflow.go()
