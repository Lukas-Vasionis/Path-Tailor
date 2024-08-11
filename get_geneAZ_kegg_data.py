"""
Script returns three objects and stores them into pickle files
These objects hold
    Data about the:
        query sent to the KEGG database
        collected data:
            Pathways holding GeneA and GeneZ (PathwaysA and PathwaysZ respectively)
            The netowrk data (nodes and edges) of these pathways
        processed data:
            Table of Pathway genes marking which genes are common between PathwaysA and PathwaysZ
    Methods (functions) for these data processes
"""

import pickle

# Import custom classes
from utils import GeneSheets, SheetsSeam

input_GeneA = 'hmgcr'
input_GeneZ = 'ace'
# Getting two sets of pathways: ones with gene A and ones with gene Z
GeneA = GeneSheets(gene=input_GeneA, org='dme').get_gene_pathways()
GeneZ = GeneSheets(gene=input_GeneZ, org='dme').get_gene_pathways()

# Getting nodes and edges of these pathways (networks)
GeneA = GeneA.get_gene_pathways_network()
GeneZ = GeneZ.get_gene_pathways_network()

# Finding common nodes (gene products, compounds) within PathwaysA (holding geneA)  and PathwaysZ (hodling geneZ)
seamAZ = SheetsSeam(GeneA, GeneZ).get_common_nodes()

# In GeneSheets.gene_pathways_meta['nodes']:
# Marking which pathways have nodes from GeneA and GeneZ sets of pathways
# by adding column 'belongs_to_pathways_with_AZ' to the table (stored in gene_pathways_meta['nodes'])

GeneA = GeneA.flag_nodes_of_AZ_pathways(seamAZ)
GeneZ = GeneZ.flag_nodes_of_AZ_pathways(seamAZ)

# Open a file and use dump()
with open('outputs/geneA.pkl', 'wb') as file:
    # A new file will be created
    pickle.dump(GeneA, file)

# Open a file and use dump()
with open('outputs/geneZ.pkl', 'wb') as file:
    # A new file will be created
    pickle.dump(GeneZ, file)

# Open a file and use dump()
with open('outputs/seamAZ.pkl', 'wb') as file:
    # A new file will be created
    pickle.dump(seamAZ, file)
