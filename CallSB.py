#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:19:06 2018

@author: dhiego
"""


""" IMPORTS """
import os
import sys
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr 
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt')


""" FUNCTIONS """
# função do Kauffman
def Kauffman_function(networkState):
    state = networkState[0]
    for i in networkState[1:]:
        if np.random.choice([True,False]) == True:
            state = state and i
        else:
            state = state or i

    return state

# função do Kauffman (deterministic)
def Kauffman_function2(networkState):
    state = networkState[0]
    for i in networkState[1:]:
        state = state and i

    return state

# recebe a rede de genes e retorna uma lista com seus estados de expressao
def gene_states_list(G):
    X_list = []
    for gene in G:
        if G.nodes[gene]['expression'] == 1:
            X_list.append(1)
        else:
            X_list.append(0)
    return X_list

# retorna os genes que tem aresta em comum ao gene do input
def linkedGenes(G, gene):
    linked_genes = []
    for i in list(G[gene]):
        if  G.nodes[i]['expression'] == 1:
            linked_genes.append(1)
        else:
            linked_genes.append(0)
            
    return linked_genes 

#calcula a similaridade entre vetores
def vectorSimilarity(vec1, vec2):
    cont = 0
    for i, j in zip(vec1, vec2):
        if i == j:
            cont = cont + 1
    return float(cont) / len(vec1)

#Update graph
def updateGraph(G, mutated):
    #atualiza a rede de genes    
    for gene in G:
        if gene in mutated:
            G.nodes[gene]['expression'] = 0
        else:
            if Kauffman_function(linkedGenes(G,gene)):
                G.nodes[gene]['expression'] = 1
            else:
                G.nodes[gene]['expression'] = 0
    return G
    



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
## GRAFO G
G = nx.Graph() # inicia o grafo

# cria os nos do grafo
half = len(genes) / 2
i = 0
for gene in genes:
    if i < half:
        G.add_node(gene, expression = 0)
    else:
        G.add_node(gene, expression = 1)
    i = i + 1

    
#gera as arestas do grafo
for i, j in zip(data['gene1'], data['gene2']) :
    #if (i != 'CACNA1C') and (j != 'CACNA1C'):
        G.add_edge(i,j)


       
""" ADD MUTATIONS """      
mutated = []
#for cont in range(len(sys.argv)):
#    mutated.append(sys.argv[cont + 1])
#mutated.append('SLC32A1')
#mutated.append('PLCB2')
#mutated.append('GRIN2B')
#mutated.append('CACNA1C')
#mutated.append('GRIA1')
#mutated.append('GRIA4')
#mutated.append('ATF4')
#mutated.append('GAD1')
#mutated.append('MAOA')
#mutated.append('GRM3')mutated.append('DRD2')
#mutated.append('COMT')
#mutated.append('PRKACA')
#mutated.append('ARAF')
#mutated.append('FYN')

   

# INITIAL STATE
state1 = gene_states_list(G)

# VECTOR OF COLORS
colors = ['r', 'g', 'b', 'c', 'm', 'y']

# VECTOR OF LABELS
labels = ['seed = 1', 'seed = 2', 'seed = 3', 'seed = 4', 'seed = 5', 'seed = 6']


#f = PdfPages('foo.pdf')
correlations = []
#allStates = []

# LOOP FOR SEED CONTROL

states = []  
states.append(state1)
#G = G1

# main loop

generation = 20
for i in range(generation):
    np.random.seed(seed =  1)

    #atualiza o grafo
    G = updateGraph(G, mutated)
       
    #lista de estados
    state = gene_states_list(G)
    states.append(state)
        
#     
#   
#        
#    correlation = []
#    anterior = state1
#    for j in states[1:]:
#        correlation.append(pearsonr(anterior,j)[0])
#        anterior = j
#    correlations.append(correlation)
         
      
  
    #plt.matshow(states)
    #allStates.append(states)
    #plt.pause(0.6)
    #plt.close()
        #plt.title(mutated)
        
#plt.matshow(states)   
np.random.seed(seed =  6)
#fig = plt.figure()
corr0 = state
correlations.append(1)
mutated.append('SLC32A1')
for i in range(generation):
   

    #atualiza o grafo
    G = updateGraph(G, mutated)
       
    #lista de estados
    state = gene_states_list(G)
    states.append(state)
corr1 = state


correlations.append(pearsonr(corr0,corr1)[0])

np.random.seed(seed =  1)
mutated.append('PLCB2')
for i in range(generation):
    #np.random.seed(seed =  1)

    #atualiza o grafo
    G = updateGraph(G, mutated)
       
    #lista de estados
    state = gene_states_list(G)
    states.append(state)
corr2 = state


correlations.append(pearsonr(corr0,corr2)[0])
 

mutated.append('GRIN2B')
for i in range(generation):
    #np.random.seed(seed =  1)

    #atualiza o grafo
    G = updateGraph(G, mutated)
       
    #lista de estados
    state = gene_states_list(G)
    states.append(state)
corr3 = state


correlations.append(pearsonr(corr0,corr3)[0])

mutated.append('CACNA1C')
for i in range(generation):
    #np.random.seed(seed =  1)

    #atualiza o grafo
    G = updateGraph(G, mutated)
       
    #lista de estados
    state = gene_states_list(G)
    states.append(state)
corr4 = state


correlations.append(pearsonr(corr0,corr4)[0])


mutated.append('GRIA1')
for i in range(generation):
    #np.random.seed(seed =  1)

    #atualiza o grafo
    G = updateGraph(G, mutated)
       
    #lista de estados
    state = gene_states_list(G)
    states.append(state)
corr5 = state


correlations.append(pearsonr(corr0,corr5)[0])

mutated.append('GRIA2')
for i in range(generation):
    #np.random.seed(seed =  1)

    #atualiza o grafo
    G = updateGraph(G, mutated)
       
    #lista de estados
    state = gene_states_list(G)
    states.append(state)
corr6 = state


correlations.append(pearsonr(corr0,corr6)[0])


#print(correlations)
#plt.subplot(6, 1, 1)
plt.plot( ['No Mut','1 Mut','2 Mut','3 Mut','4 Mut', '5 Mut','6 Mut'],correlations, colors[0])
plt.xlabel('Number of mutations')
plt.ylabel('Correlation')
#fig.savefig("plot.pdf")
#plt.subplot(6, 1, 2)
#plt.plot(range(len(correlations[1])), correlations[1], colors[1])
#plt.subplot(6, 1, 3)
#plt.plot(range(len(correlations[2])), correlations[2], colors[2])
#plt.subplot(6, 1, 4)
#plt.plot(range(len(correlations[3])), correlations[3], colors[3])
#plt.subplot(6, 1, 5)
#plt.plot(range(len(correlations[4])), correlations[4], colors[4])
#plt.subplot(6, 1, 6)
#plt.plot(range(len(correlations[5])), correlations[5], colors[5])
#plt.subplot(6, 1, 1)
#plt.plot(allStates[0])
#plt.subplot(6, 1, 2)
#plt.plot(allStates[1])
#plt.subplot(6, 1, 3)
#plt.plot(allStates[2])
#plt.subplot(6, 1, 4)
#plt.plot(allStates[3])
#plt.subplot(6, 1, 5)
#plt.plot(allStates[4])
#plt.subplot(6, 1, 6)
#plt.plot(allStates[5])


#fig = plt.figure()
## lista de nos on e off         
#node_list_on = []
#node_list_off = []      
#
#for gene in genes:
#    if G.nodes[gene]['expression'] == 0:
#        node_list_off.append(gene)
#    else:
#        node_list_on.append(gene)
#        
#        
##Desenha o grafo
#pos = nx.spring_layout(G, fixed=None)
#nx.draw_networkx_nodes(G, pos,node_list_on,node_color='g',node_size=20, alpha=0.8)
#nx.draw_networkx_nodes(G, pos,node_list_off,node_color='r',node_size=20, alpha=0.8)
#nx.draw_networkx_edges(G, pos, width=0.2, alpha=0.5)
#
##labels = {}
##for gene in genes:
##    labels[gene] = gene
##nx.draw_networkx_labels(G, pos, labels, font_size=16)
#
#plt.axis('off')
#plt.show()
#fig.savefig("plot.pdf")

#from collections import Counter
#print(Counter(state).keys()) # equals to list(set(words))
#print(Counter(state).values()) # counts the elements' frequency

   
              




#plt.hist(correlations)
#dic = {}
#for i, j in zip(state, genes):
#    dic[j] = i 
#
#print(dic)




#print(correlation)

