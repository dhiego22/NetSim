#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 10:52:12 2018

@author: dhiego
"""

import pandas as pd
import networkx as nx
import numpy as np
import os
import itertools
from scipy.stats.stats import pearsonr 
import matplotlib.pyplot as plt
#from IPython import get_ipython
#get_ipython().run_line_magic('matplotlib', 'qt')


""" MY FUNCTIONS """
from myfunctions import createGraph, drawGraph, updateGraph, Kauffman_function, linkedGenes, gene_states_list
from myfunctions import calculateCorrelations, runGenerations, atractorsCounts

""" DIRECTORY """
os.chdir('/home/dhiego/Documents/NetSim')


""" READ DATA """
data = pd.read_csv('graphNeurotransmissor.csv')
#pega os genes
list1 = data['gene1'].tolist()
list2 = data['gene2'].tolist()
list3 = list1 + list2
genes = set(list3)


""" CREATE GRAPH """
G = createGraph(genes, data)


""" DRAW GRAPH """
pos = nx.spring_layout(G, k = 0.7)
#drawGraph(G,pos)


""" ADD MUTATIONS """      
mutations = []
#np.random.seed(seed =  1)
#for gene in genes:
#    if np.random.choice([True, False, False, False, False, False, False, False, False, False]) == True:
#        mutations.append(gene)
        
#mutations.append('SLC32A1')
#mutations.append('PLCB2')
#mutations.append('GRIN2B')
mutations.append('CACNA1C')
#mutations.append('GRIA1')
#mutations.append('GRIA4')
#mutations.append('ATF4')
#mutations.append('GAD1')
#mutations.append('MAOA')
#mutations.append('GRM3')
#mutations.append('DRD2')
#mutations.append('COMT')
#mutations.append('PRKACA')
#mutations.append('ARAF')
#mutated.append('FYN')
#print(mutations)



#for gene in genes:
    
#mutation = []
##generations        
#np.random.seed(seed =  1)    
#states = runGenerations(100, G, mutation)
##calculate correlations between each epoch
#correlation = calculateCorrelations(states)
#plt.plot(range(len(correlation)), correlation, 'r')
#plt.xlabel('Generations')
#plt.ylabel('Correlation')
#plt.title('GRIN2B + SLC32A1 mutation')
   


#np.random.seed(seed =  1)
#for generation in range(100):
#    #atualiza o grafo
#    G = updateGraph(G, mutation)
#    #desenha o grafo   
#    #plt.title('Generation: ' + str(generation))      
#    #drawGraph(G,pos,genes)
#    #salva o grafo
#    #plt.savefig("Graph" + str(generation) + ".png", format="PNG")
#    state = gene_states_list(G)
#    states.append(state)
#
#
#plt.matshow(states)







"""INITIAL STATE"""
state1 = gene_states_list(G)

"""SET SEED"""
#np.random.seed(seed =  1)

"""ATRACTORS LIST"""
atractors = []

"""SET GENERATIONS """
generations = 50

#divide mutations in many subsets
for subsetNumber in range(len(mutations)):
    #combinations = itertools.combinations(mutations, subsetNumber+1)
    #get each subset of the combinations
    #for subset in combinations:
    for s in range(30):
        #np.random.seed(seed =  s)
        #print(list(subset))
        #mutated = list(subset)
        mutated = []
        states = [] 
        states.append(state1)
        for i in range(generations):
            np.random.seed(seed =  s)
#            if i == 25:
#                mutated.append('MAOA')
#                mutated.append('GRIA4')
#                mutated.append('ATF4')
#                mutated.append('GAD1')
            #atualiza o grafo
            G = updateGraph(G, mutated)
       
            #lista de estados
            state = gene_states_list(G)
            states.append(state)
        #plt.title(list(subset))
        atractors.append(states[-1])
        #plt.matshow(states)  
plt.matshow(atractors)  
#print(atractorsCounts(atractors))
#generations        
#np.random.seed(seed =  2)    
#states = runGenerations(10, G, mutations)
#plt.matshow(states)    

#
# #calculate correlation    
#        anterior = state1
#        for j in states[1:]:
#            correlation.append(pearsonr(anterior,j)[0])
#            anterior = j























""" WRITE DICTIONARY """
#with open('graphPositions.txt', 'w') as file:
#     for k, v in pos.items():
#         file.write(str(k) + ' : ' + str(v) + '\n')
     #file.write(pos)


#pos1 = readDict('graphPositions.txt', ':')
#def replace_value_with_definition(current_dict, key_to_find, definition):
#    for key in current_dict.keys():
#        if key == key_to_find:
#            current_dict[key] = definition
#            
#for (k,v), (k2,v2) in zip(pos.items(), pos1.items()):
#    replace_value_with_definition(pos,k,v2.split(',')[0])

#def readDict(filename, sep):
#    with open(filename, "r") as f:
#        dict = {}
#        for line in f:
#            values = line.split(sep)
#            dict[values[0]] = np.asarray(values[1].rstrip("\n"))
#      
#        return(dict)




#permut = list(itertools.permutations(mutations))