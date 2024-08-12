# Path-Tailor

## The vision
Path Tailor is a computational biology tool that looks for gene interactions between 
GeneA and GeneZ across known pathways of KEGG databases. The tool returns a set of pathways 
that share common nodes (genes or metabolites) and link GeneA to GeneZ across this set. Additionally,
it visualises these interactions in the network graph **[to be developed]**. The tool is planned to have user-friendly app interface
of Streamlit. **[to be developed]**

## Overview of data gathering and processing
### Data gathering 
Path Tailor uses KEGG pathway databases as a source data for pathways, genes and graph elements of KEGG pathways
For data extraction and processing it uses libraries of `bioservices`, `Bio` that provides an interface to access 
KEGG through its API.

### Data processing
The tool looks for the pathways contaminating GeneA (PathwaysA) or GeneZ(PathwaysZ). 
This leads to one of the three scenarios:
- Tool detects pathways that hold both geneA and geneZ (pathwaysAZ)

<img src='docs/img/CASE_AZ_pathway.png' width=500vw>

- If there are no pathways with genes A and Z, then it looks for other common nodes (gene products or metabolites) 
that would link each pair of pathways (Sew them like Sheets over the Seams - hence the names of Classes within the code)

<img src='docs/img/CASE_intermediate_genes.png' width=500vw>

- If pathwaysA and pathwaysZ have no common nodes, the algorythm looks for intermediate pathways that 
share other genes with pathwaysA or pathwaysZ **[to be developed]**.
  - (Alternatively, it should be possible to introduce additional genes from omics experiments. However, 
  this idea is yet to be explored.)**[to be developed]**
  
  
<img src='docs/img/CASE_intermediate_pathways.png' width=500vw>

### Data visualisation **[to be developed]**
For data visualisation it will use combination of `networkx`, `pyviz` to create interactive network graphs with 
kegg pathway graph aesthetics. All the preference options, inputs, and graphs will be wrapped into a 
quick app of `streamlit`.




