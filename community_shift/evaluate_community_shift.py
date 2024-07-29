from anadama2 import Workflow
import os
import itertools
import copy

workflow = Workflow(version="0.1", description="Evaluate community_shift")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The memory in megabytes allocated to run the command", type=int, default=10000)
workflow.add_argument("time", desc="The time in minutes allocated to run the command", type=int, default=120)
workflow.add_argument('parameters', type=str, desc="Specify the parameters")
workflow.add_argument('tmp', desc="Whether to use reduced parameters", action="store_true")
args = workflow.parse_args()
this_directory = str(os.path.dirname(os.path.realpath(__file__))).rstrip('/') + '/'

# output, set this to the GitHub repository
output = os.path.abspath(args.output.rstrip("community_shift/")) + "/"
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

def compute_generation_outputs(generator):
    if args.tmp:
        parameters_file = 'community_shift/data_generation/' + generator + '_tmp.txt'
        nIterations = 5
    else:
        parameters_file = 'community_shift/data_generation/' + generator + '.txt'
        nIterations = 100
    
    with open(str(parameters_file), 'r') as file:
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

    metadata_types = param_dict.pop('metadataType', None)
    working_directory = output

    param_list_generation = [dict([(key, value[0]) for key, value in param_dict.items()]) for _ in range(len([item for sublist in param_dict.values() for item in sublist]))]
    counter = 0
    for key, values in param_dict.items():
        for value in values:
            param_list_generation[counter][key] = value
            counter += 1

    param_list_final = []
    for param_dict in param_list_generation:
        for metadata_type in metadata_types:
            new_param_dict = copy.deepcopy(param_dict) # Create a deep copy
            new_param_dict['metadataType'] = metadata_type
            if metadata_type == 'binary':
                new_param_dict['nMetadata'] = str(1)
                new_param_dict['nPerSubject'] = str(1)
            param_list_final.append(new_param_dict)

    param_list_final = set([frozenset(param_single_final.items()) for param_single_final in param_list_final])

    generation_outputs = []

    for param in param_list_final:
        param = dict(param)
        if int(param['nPerSubject']) > 1:
            inputSubString = 'RandomEffect'
        else:
            inputSubString = 'noRandomEffect'

        new_depends_folder = output + 'Input/community_shift/' + generator + '/' + '_'.join([inputSubString, 
        param['metadataType'], param['nSubjects'], param['nPerSubject'], param['nMicrobes'], 
        param['spikeMicrobes'], param['nMetadata'], param['effectSize'], param['effectPos'], param['readDepth']]) + '/'

        generation_outputs.extend([new_depends_folder + file_type + '_' + str(file_number) + '.tsv' for file_type, file_number in list(itertools.product(['metadata', 'abundance', 'truth', 'unscaled'], [i for i in range(1, nIterations + 1)]))])

    return generation_outputs

for generator in param_dict['generators']:
    if args.tmp:
        generate_command = 'python3 community_shift/data_generation/' + generator + '_workflow.py --parameters community_shift/data_generation/' + generator + '_tmp.txt --working-directory ' + output + ' --cores ' + str(cores) + ' --tmp'
    else:
        generate_command = 'python3 community_shift/data_generation/' + generator + '_workflow.py --parameters community_shift/data_generation/' + generator + '.txt --working-directory ' + output + ' --cores ' + str(cores)

    if not all(os.path.exists(file_path) for file_path in compute_generation_outputs(generator)):
        workflow.add_task_gridable(actions=generate_command,
            targets=compute_generation_outputs(generator),
            time=time,
            mem=memory,
            cores=cores,
            partition=partition
            )

#############
# Run tools #
#############

run_tools_directory = this_directory + 'run_tools/'
working_directory = output

def compute_running_outputs(generator, tool, param):
    tool = tool.rstrip('.R').lstrip('run_')
    if args.tmp:
        nIterations = 5
    else:
        nIterations = 100

    if int(param['nPerSubject']) > 1:
        inputSubString = 'RandomEffect'
    else:
        inputSubString = 'noRandomEffect'

    new_depends_folder = output + 'Output/community_shift/' + generator + '/' + '_'.join([inputSubString, 
    param['metadataType'], param['nSubjects'], param['nPerSubject'], param['nMicrobes'], 
    param['spikeMicrobes'], param['nMetadata'], param['effectSize'], param['effectPos'], param['readDepth'], tool]) + '/'

    generation_outputs = [new_depends_folder + 'associations_' + str(file_number) + '.tsv' for file_number in [i for i in range(1, nIterations + 1)]]

    return generation_outputs

for generator in param_dict['generators']:
    if args.tmp:
        with open(str(working_directory + 'community_shift/data_generation/' + generator + '_tmp.txt'), 'r') as file:
            lines = file.readlines()
    else:
        with open(str(working_directory + 'community_shift/data_generation/' + generator + '.txt'), 'r') as file:
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

    metadata_types = param_dict.pop('metadataType', None)

    param_list_generation = [dict([(key, value[0]) for key, value in param_dict.items()]) for _ in range(len([item for sublist in param_dict.values() for item in sublist]))]
    counter = 0
    for key, values in param_dict.items():
        for value in values:
            param_list_generation[counter][key] = value
            counter += 1

    param_list_final = []
    for param_dict in param_list_generation:
        for metadata_type in metadata_types:
            new_param_dict = copy.deepcopy(param_dict) # Create a deep copy
            new_param_dict['metadataType'] = metadata_type
            if metadata_type == 'binary':
                new_param_dict['nMetadata'] = str(1)
                new_param_dict['nPerSubject'] = str(1)
            param_list_final.append(new_param_dict)

    param_list_final = set([frozenset(param_single_final.items()) for param_single_final in param_list_final])

    tools = [file for file in os.listdir(run_tools_directory) if file.startswith('run_')]

    for param in param_list_final:
        for tool in tools:
            param = dict(param)
            if int(param['nPerSubject']) > 1:
                random_effect_string = ' --RandomEffect'
            else:
                random_effect_string = ''
            new_command = 'Rscript ' + run_tools_directory + tool + \
            random_effect_string + \
            ' --metadataType ' + param['metadataType'] + \
            ' --nSubjects ' + param['nSubjects'] + \
            ' --nPerSubject ' + param['nPerSubject'] + \
            ' --nMicrobes ' + param['nMicrobes'] + \
            ' --spikeMicrobes ' + param['spikeMicrobes'] + \
            ' --nMetadata ' + param['nMetadata'] + \
            ' --effectSize ' + param['effectSize'] + \
            ' --effectPos ' + param['effectPos'] + \
            ' --readDepth ' + param['readDepth'] + \
            ' --nCores ' + str(cores) + \
            ' --workingDirectory ' + working_directory + \
            ' --generator ' + generator
            if args.tmp:
                new_command = new_command + ' --nIterations 5'

            if not all(os.path.exists(file_path) for file_path in compute_running_outputs(generator, tool, param)):
                workflow.add_task_gridable(actions=new_command,
                depends=compute_generation_outputs(generator),
                targets=compute_running_outputs(generator, tool, param),
                time=time,
                mem=memory,
                cores=cores,
                partition=partition
                )

workflow.go()
