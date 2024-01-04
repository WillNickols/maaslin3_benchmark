import argparse
import subprocess
import concurrent.futures
import os
import copy

parser = argparse.ArgumentParser(description="")
parser.add_argument('--working-directory', type=str, help="Working directory")
parser.add_argument('--generators', type=str, help="Comma separated generators")
parser.add_argument('--cores', type=int, help="Number of cores to use")
args = parser.parse_args()

this_directory = os.getcwd()

generators = str(args.generators).split(',')
commands = ['pwd']

for generator in generators:
    with open(str(this_directory + generator + '.txt'), 'r') as file:
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
    working_directory = args.working_directory

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

    tools = [file for file in os.listdir() if file.startswith('run_')]

    for param in param_list_final:
        for tool in tools:
            param = dict(param)
            if int(param['nPerSubject']) > 1:
                random_effect_string = ' --RandomEffect'
            else:
                random_effect_string = ''
            new_command = 'Rscript ' + this_directory + tool + \
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
            ' --nCores 1' + \
            ' --workingDirectory ' + working_directory + \
            ' --generator ' + generator
            commands.append(new_command)

max_concurrent_processes = int(args.cores)

def run_command(command):
    try:
        process = subprocess.Popen(command, shell=True, universal_newlines=True)
        stdout, stderr = process.communicate()
        return_code = process.returncode
        return (command, stdout, stderr, return_code)
    except Exception as e:
        return (command, str(e), '', 1)

print("Running " + str(len(commands)) + " variations")

with concurrent.futures.ThreadPoolExecutor(max_concurrent_processes) as executor:
    results = list(executor.map(run_command, commands))