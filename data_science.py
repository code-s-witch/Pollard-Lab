import seaborn as sns
import numpy as np
import pandas as pd

# FUNCTIONS
#reads in all the sample names used as keys in the dictionary
def getFileNames(file, exp):
    samples_list = []
    inFile = open(file, 'r')
    for line in inFile:
        splitLine=line.strip().split('\n')
        name=splitLine[0]
        if exp in name:
            samples_list.append(name)
#         elif exp not in name:
#             print(name)
    inFile.close()
    return samples_list
#stores counts and length of enhancer as tuples into a dictionary using samples as key
def counts_with_length(rel_path, file, exp_dict):
    counts = []
    inFile = open(rel_path + file, 'r')
    for line in inFile:
        splitLine=line.strip().split('\t')
        start_pos = int(splitLine[1])
        end_pos = int(splitLine[2])
        length = (end_pos - start_pos)+1
        val=int(splitLine[3])
        counts.append((np.log(val+1),length))
    inFile.close()
    exp_dict[file] = counts

#normalizes for sequencing depth CPM counts per million and for enhancer length KB (kilobase)
def normalize_CPKM(exp_dict):
    for name, counts in exp_dict.items():
        count_per_ehancer = [x[0] for x in counts]
        total_counts_per_sample = sum(count_per_ehancer)
#         print(count_per_ehancer)
#         enhancer_lengths = [x[1] for x in counts]
#         print(enhancer_lengths)
        #update key with new values normalized by sequencing depth and by enhancer length
        exp_dict[name] = [(tup[0] / (tup[1] * total_counts_per_sample)) * 10**9 for tup in counts]
#         print(exp_dict[name])

#stores counts into dictionary
def countsIntoDict(rel_path, file, exp_dict):
    counts = []
    inFile = open(rel_path + file, 'r')
    for line in inFile:
        splitLine=line.strip().split('\t')
        val=int(splitLine[3])
        counts.append(np.log(val+1))
    inFile.close()
    exp_dict[file] = counts

#normalizes using counts per million/ total sample reads as scaling factor
#normalization for sequencing depth
def normalize_CPM(exp_dict):
    for name, counts in exp_dict.items():
        exp_dict[name] = [(x / sum(counts)) * 10**6 for x in counts] #update key with new normalized list

#takes a boolean that tell is to either call normalize_CPM() or normalize_CPKM()
#if CPM is true will use that method, else uses CPKM method
def normalize(CPM_bool, rel_path, samples, sampleDict):
    if CPM_bool:
        #store log transformed counts in dict
        for sample in samples:
            countsIntoDict(rel_path, sample, sampleDict)
        normalize_CPM(sampleDict)
    elif not CPM_bool:
        for sample in samples:
            counts_with_length(rel_path, sample, sampleDict)
        normalize_CPKM(sampleDict)
