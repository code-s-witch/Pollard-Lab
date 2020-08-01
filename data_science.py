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
        counts.append((val,length))
    inFile.close()
    exp_dict[file] = counts

#normalizes for sequencing depth CPM counts per million and for enhancer length KB (kilobase)
def normalize_CPKM(exp_dict):
    for name, counts in exp_dict.items():
        count_per_ehancer = [x[0] for x in counts]
        total_counts_per_sample = sum(count_per_ehancer)
        #print(count_per_ehancer)
        #print(enhancer_lengths)
        #update key with new values normalized by sequencing depth and by enhancer length
        cpkm_vals = [(tup[0] / (tup[1] * total_counts_per_sample)) * 10**9 for tup in counts]
        log_cpkm_vals = [np.log(x+1) for x in cpkm_vals]
        exp_dict[name] = log_cpkm_vals


#stores counts into dictionary
def countsIntoDict(rel_path, file, exp_dict):
    counts = []
    inFile = open(rel_path + file, 'r')
    for line in inFile:
        splitLine=line.strip().split('\t')
        val=int(splitLine[3])
        counts.append(val)
    inFile.close()
    exp_dict[file] = counts

#normalizes using counts per million/ total sample reads as scaling factor
#normalization for sequencing depth
def normalize_CPM(exp_dict):
    for name, counts in exp_dict.items():
        cpm_vals = [(x / sum(counts)) * 10**6 for x in counts]
        log_cpm_vals = [np.log(x+1) for x in cpm_vals]
        exp_dict[name] = log_cpm_vals

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

def run_ds_on_exp(samples_dir, rel_path, exp, cpm_bool, repl_bool):
    print('exp: ', exp)
    sampleDict = {}
    samples = getFileNames(samples_dir, exp)
    if cpm_bool:
        for sample in samples:
            countsIntoDict(rel_path, sample, sampleDict)
    else:
        for sample in samples:
            counts_with_length(rel_path, sample, sampleDict)

    #rename dictionary keys to specify conditions and replicates         
    rename_keys(sampleDict, samples)

    #normalize the data using CPM or CPKM
    #if CPM is True then we normalize by CPM, else we normalize by CPKM
    #take a dictionary with condition: counts
    #takes a bool for either cpm or cpkm normalization
    #takes a bool for either combined replicates or separate replicates
    normDict = norm(sampleDict, repl_bool, cpm_bool, exp)
    
    return normDict

    
#TODO: use parserArg to take arguments for normalization, replicates, exp, etc
def combine_replicates(sampleDict, cpm_bool, exp):
    replicateDict = {}
    if cpm_bool:
        print('booean is ',  cpm_bool , ' :cpm normalization')
        #combined sham replicates into a list
        sham_combined = []
        for s1,s2 in zip(sampleDict['Sham_r1'],sampleDict['Sham_r2']):
            sham_combined.append(s1+s2)    
        replicateDict[exp+'_Sham'] = sham_combined

        #combine tac replicates into a list
        tac_combined = []
        for t1, t2 in zip(sampleDict['TAC_r1'],sampleDict['TAC_r2']):
            tac_combined.append(t1+t2)
        replicateDict[exp+'_TAC'] = tac_combined

        #combine jq1-withdrawn replicates
        jq1w_combined = []
        for jw1, jw2 in zip(sampleDict['JQ1_w_r1'],sampleDict['JQ1_w_r2']):
            jq1w_combined.append(jw1+jw2)
        replicateDict[exp+'_JQ1_w'] = jq1w_combined

        #combine jq1 replicates
        jq1_combined = []
        for j1, j2 in zip(sampleDict['JQ1_r1'],sampleDict['JQ1_r2'] ):
            jq1_combined.append(j1+j2)
        replicateDict[exp+'_JQ1'] = jq1_combined    
    else:
        print('cpm_bool is ' ,cpm_bool, ' :cpkm normalization')
        #if CPM_bool is false then the dict is a list of tuples holding repl count, enhancer length
        sham_combined = []
        for s1,s2 in zip(sampleDict['Sham_r1'],sampleDict['Sham_r2']):
            sham_combined.append((s1[0]+s2[0], s1[1]))
        replicateDict[exp+'_Sham'] = sham_combined
        
        #combine tac replicates into a list
        tac_combined = []
        for t1, t2 in zip(sampleDict['TAC_r1'],sampleDict['TAC_r2']):
            tac_combined.append((t1[0]+t2[0], t1[1]))
        replicateDict[exp+'_TAC'] = tac_combined

        #combine jq1-withdrawn replicates
        jq1w_combined = []
        for jw1, jw2 in zip(sampleDict['JQ1_w_r1'],sampleDict['JQ1_w_r2']):
            jq1w_combined.append((jw1[0]+jw2[0], jw1[1]))
        replicateDict[exp+'_JQ1_w'] = jq1w_combined

        #combine jq1 replicates
        jq1_combined = []
        for j1, j2 in zip(sampleDict['JQ1_r1'],sampleDict['JQ1_r2'] ):
#             print(sampleDict['JQ1_r1'])
#             print(sampleDict['JQ1_r2'])
            jq1_combined.append((j1[0]+j2[0], j1[1]))
#         print(jq1_combined)
        replicateDict[exp+'_JQ1'] = jq1_combined    
        
    return replicateDict

def norm(sampleDict, repl_bool, cpm_bool, exp):
    replCombinedDict = {}
    if repl_bool:
        print('replicates combined')
        replCombinedDict = combine_replicates(sampleDict, cpm_bool, exp)
        if cpm_bool:
            normalize_CPM(replCombinedDict)
        else:
            normalize_CPKM(replCombinedDict)
        return replCombinedDict    
    else:
        print('replicates separate')
        if cpm_bool:
            normalize_CPM(sampleDict)
        else:
            normalize_CPKM(sampleDict)
        return sampleDict
    
    
def rename_keys(sampleDict, samples):
        #dictionary[new_key] = dictionary.pop(old_key)
        #https://stackoverflow.com/questions/4406501/change-the-name-of-a-key-in-dictionary
        #TODO: read in a list of ordered conditions with replicates and make this a for loop 
        sampleDict['Sham_r1'] = sampleDict.pop(samples[0])
        sampleDict['Sham_r2'] = sampleDict.pop(samples[1])
        sampleDict['TAC_r1'] = sampleDict.pop(samples[2])
        sampleDict['TAC_r2'] = sampleDict.pop(samples[3])
        sampleDict['JQ1_w_r1'] = sampleDict.pop(samples[4])
        sampleDict['JQ1_w_r2'] = sampleDict.pop(samples[5])
        sampleDict['JQ1_r1'] = sampleDict.pop(samples[6])
        sampleDict['JQ1_r2'] = sampleDict.pop(samples[7])
        
