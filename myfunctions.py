#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 15:22:49 2018

@author: dhiego
"""

import pandas as pd
import networkx as nx
import numpy as np
import os
import itertools
from scipy.stats.stats import pearsonr 
import matplotlib.pyplot as plt

""" FUNCTIONS """
def createGraph(genes, data):
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
    return G

def drawGraph(G,pos,genes):
    # lista de nos on e off         
    node_list_on = []
    node_list_off = []      

    for gene in genes:
        if G.nodes[gene]['expression'] == 0:
            node_list_off.append(gene)
        else:
            node_list_on.append(gene)

    nx.draw_networkx_nodes(G, pos,node_list_on,node_color='c',node_size=20, alpha=0.8)
    nx.draw_networkx_nodes(G, pos,node_list_off,node_color='m',node_size=20, alpha=0.8)
    nx.draw_networkx_edges(G, pos, width=0.05, alpha=0.8)

#Update graph
def updateGraph(G, mutated):
    new_G = G
    #atualiza a rede de genes    
    for gene in G:
        if gene in mutated:
            new_G.nodes[gene]['expression'] = 0
        else:
            if Kauffman_function(linkedGenes(G,gene)):
                new_G.nodes[gene]['expression'] = 1
            else:
                new_G.nodes[gene]['expression'] = 0
    return new_G
    
# função do Kauffman
def Kauffman_function(networkState):
    state = networkState[0]
    for i in networkState[1:]:
        if np.random.choice([True,False]) == True:
            state = state and i
        else:
            state = state or i

    return state

# retorna os genes que tem aresta em comum ao gene do input
def linkedGenes(G, gene):
    linked_genes = []
    for i in list(G[gene]):
        if  G.nodes[i]['expression'] == 1:
            linked_genes.append(1)
        else:
            linked_genes.append(0)
            
    return linked_genes     

# recebe a rede de genes e retorna uma lista com seus estados de expressao
def gene_states_list(G):
    X_list = []
    for gene in G:
        if G.nodes[gene]['expression'] == 1:
            X_list.append(1)
        else:
            X_list.append(0)
    return X_list

#run n generation on graph G with m mutations
def runGenerations(n, G, m):
    # INITIAL STATE
    state1 = gene_states_list(G)
    states = []  
    states.append(state1)
    for generation in range(n):
        #atualiza o grafo
        G = updateGraph(G, m)
       
        #lista de estados
        state = gene_states_list(G)
        states.append(state)  
    
    return states

#return the correlations between each state transition
def calculateCorrelations(states):
    correlation = []
    anterior = states[0]
    for j in states[1:]:
            correlation.append(pearsonr(anterior,j)[0])
            anterior = j
    
    return correlation
    
#receive a list of atractors and return a dictionary of their counts
def atractorsCounts(atractors):
     at_list = []
     at_dict = {}
     for at in atractors:
         if at in at_list:
            at_dict[list2string(at)] += 1
         else:
             at_dict[list2string(at)] = 1
             at_list.append(at)
   
     return at_dict
    
#transform a list of int into a string
def list2string(lista):
    string = ''
    for l in lista:
        string = string + str(l)
    return string
        
#return a dictionary of the counts of the states transitions       
def TransitionDict(states):
    mm = {}
    state1 = states[0]
    for state2 in states[1:]:
        if list2string(state1)+ '->' +list2string(state2) in mm:
            mm[list2string(state1)+ '->' +list2string(state2)] += 1
        else:
            mm[list2string(state1)+ '->' +list2string(state2)] = 1
            
        state1 = state2
    
    return mm
    
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
