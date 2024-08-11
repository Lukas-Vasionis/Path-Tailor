import warnings
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from bioservices.kegg import KEGG

# Import Pandas, so we can use dataframes
import pandas as pd


class GeneSheets: 
    """
    Object to store data of GeneA and GeneZ. This includes query variables and geneA/geneZ related pathway data
    
    """
    def __init__(self, gene, org_id, geneA_aliases=None, geneZ_aliases=None):
        self.gene = gene # queried name
        self.org_id = org_id #queried organism
        self.gene_aliases = self.get_gene_aliases() # found gene aliases
        self.gene_pathways = None # pathway ids where gene was found
        self.gene_pathways_meta = None #network data of found pathways (nodes and edges)

    
    def get_gene_aliases(self):
       
        print(f"Getting aliases of {self.gene}")
              
        request_aliases = KEGG().find(self.org_id, self.gene)
        
        aliases = request_aliases.strip('\n').split('\t')[1]
        list_aliases = aliases.split(',')
        list_aliases = [x.strip(' ') for x in list_aliases]

        print(list_aliases)
        return list_aliases

    def get_gene_pathways(self):
        def get_tbl_pathways(gene_aliases, org_id):
            
            tbls=[]
            for g_alias in gene_aliases:
           
                pathways = KEGG().get_pathway_by_gene(g_alias, org_id)
                
                if pathways: #if any pathways found
                    df = pd.DataFrame(list(pathways.items()), columns=['path_id', 'path_name'])
                    tbls.append(df)
            
            return pd.concat(tbls)
        
        self.gene_pathways = get_tbl_pathways(self.gene_aliases, self.org_id)
        
        return self
        
    # Decorator to suppress specific warning in GeneSheets.get_gene_pathways() method
    def suppress_the_warning(func):
        def inner(*args, **kwargs):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                
                result = func(*args, **kwargs)
                
                # Filter out the specific warning
                w = [warning for warning in w if 'XML document using an HTML parser' not in str(warning.message)]
                # If there are other warnings - print them.
                if w:
                    for warning in w:
                        warnings.showwarning(warning.message, warning.category, warning.filename, warning.lineno)
            return result
        return inner
        
    @suppress_the_warning    
    def get_gene_pathways_network(self):
        """
        From list of pathways related to gene (A or Z), it gets metadata of pathway's edges and nodes
        
        returns {"nodes":pd.DataFrame(), "edges":pd.DataFrame()} 
        """
        print("Parsing pathway ids for netowrks data")
    
        paths = self.gene_pathways["path_id"].to_list()
        paths = list(set(paths))
        tbls_pathways_nodes=[]
        tbls_pathways_edges=[]
        
        for path in paths:
            print(f'''\tPath: {path}''')
            pathway_meta = get_pathway_network_REST(path)
            # pathway_meta = get_pathway_network_bioservices(path)
            
            pathway_meta_nodes=pd.DataFrame(pathway_meta['entries'])
            pathway_meta_nodes['path_id'] = path
            tbls_pathways_nodes.append(pathway_meta_nodes)
            
            pathway_meta_edges = pd.DataFrame(pathway_meta['relations'])
            pathway_meta_edges['path_id'] = path
            tbls_pathways_edges.append(pathway_meta_edges)
    
        pathways_network_meta={
            "nodes" : pd.concat(tbls_pathways_nodes), 
            "edges" : pd.concat(tbls_pathways_edges)
        }
        
        self.gene_pathways_meta=pathways_network_meta
        return self

    def flag_nodes_of_AZ_pathways(self, Seam):
        
        print("Marking common nodes between GeneA and GeneZ sets of pathways")
        
        gene_nodes = self.gene_pathways_meta['nodes']
        
        """
        Adding this code in case the method is re-executed. 
        It prevents from making columns additional columns of belongs_to_pathways_with_AZ with suffix _x,_y, etc.
        """
        cols_to_drop=[x for x in gene_nodes.columns if 'belongs_to_pathways_with_AZ' in x]
        if cols_to_drop:
            gene_nodes = gene_nodes.drop(cols_to_drop, axis=1)

        gene_nodes = gene_nodes.merge(Seam.common_nodes, on='name', how='left')
        self.gene_pathways_meta['nodes']=gene_nodes
            
        return self

    
class SheetsSeam:
    """
    Seam - a junction formed by sewing together two pieces of material
    In this case, pathways are like sheets which we try to sew together over the common compounds or genes.
    This class stores data about common points of pathways with Gene A and Gene Z
    """
    def __init__(self, geneA_data, geneZ_data):
        self.geneA_data=geneA_data
        self.geneZ_data=geneZ_data
        self.common_nodes=None
        
    def get_common_nodes(self):
    # self.common_nodes stores common metabolates/genes/pathways betwen A and Z pathways

        geneA_nodes = self.geneA_data.gene_pathways_meta['nodes']
        geneZ_nodes = self.geneZ_data.gene_pathways_meta['nodes']
        try:
            common_nodes=geneA_nodes.merge(geneZ_nodes, on='name',how='outer', indicator=True)
        except KeyError:
            geneA_nodes.head()
            # geneZ_nodes.head()
            
            
        common_nodes['_merge']=common_nodes['_merge'].map({"left_only":"A", 'right_only':'Z', 'both':'A,Z'})
        common_nodes=common_nodes.rename(columns={'_merge':'belongs_to_pathways_with_AZ'})
        
        self.common_nodes=common_nodes[['name', 'belongs_to_pathways_with_AZ']].drop_duplicates().sort_values('name')
        
        return self


def get_pathway_network_REST(kegg_path_id):
    """
    kegg_path_id - id that looks like this dme00010 (dme - keg organism id for d melanogaster)
    
    Gets pathway's nodes (entries) and edges (relations) using Bio.KEGG module. 
    This way of getting such data:
        seems to be more reliable as sometimes bioservices dont return relations (try example of dme00010) :
        returns graphical data of node position in kegg graph. This may be very useful in further steps of graphing.
        seems to be faster (this is purely perceptional - no benchmarking was made)
    
    returns:: dict of enries and relations:
    {'entries': [{entry_id: '...',
                '...':'...'},
                {entry_id: '...',
                '...':'...'},...],
     'relations':[{'entry1': '112',
                   'entry2': '119',
                   'link': 'ECrel',
                   'value': '95',
                   'name': 'compound'},
                  {'entry1': '112',
                   'entry2': '119',
                   'link': 'ECrel',
                   'value': '95',
                   'name': 'compound'},...]
    }
    """

    def process_entries(rest_kgml_pathway):
        
        def process_entry(entry_value):
            entry_value = entry_value.__dict__
            
            # Unpacking graphics key
            """
            Each entry is a dict. The entry dict has a graphics key which is another dict.
            Here I move graphics dict to the same level as other keys in the entry dict
                FROM:
                    entry_value={
                        'id':'...',
                        'graphics':{'x' : '...', 'y':'...',...}
                    }
                TO:
                    entry_value={
                        'id':'...',
                        'graphics_x':'...',
                        'graphics_y':'...',
                        'graphics_...':'...'
                    }
            """
            
            graphics_data = entry_value['graphics'][0].__dict__
            graphics_data = {f"graphics_{k}":v for k,v in graphics_data.items()} # Adding graphics_ prefix to each key of graphics dict 
            graphics_data ['graphics_data_len'] = len(entry_value['graphics']) # adding this in case the entry has multiple graphics elements. Its strange that graphics entry is a list - why did developers do that?
            entry_value.update(graphics_data) # Adding dict of graphics key
     
            # Renaming some keys (old keys will be removed at the end of this f-tion.)
            """
            I do this in case I will need to join this 
            data with data gathered using bioservices.kegg.KEGG module
            """
            entry_value['id']=entry_value['_id']
            entry_value['name']=entry_value['_names'][0]
            entry_value['gene_names']=" ".join(entry_value['_names'])
            
            # Selecting keys of the entry dict
            entry_graphics_keys=[x for x in entry_value.keys() if x.startswith('graphics') and x !='graphics']
            keys_to_select=['id', 'name', 'type', 'link', 'gene_names']+entry_graphics_keys
            entry_value={k:v for k,v in entry_value.items() if k in keys_to_select}
            
            return entry_value
    
        entries = {k:process_entry(v) for k,v in rest_kgml_pathway.entries.items()}
        entries = entries.values()
        
        return entries
    
    def process_relations(rest_kgml_pathway):
        def process_relation(relation):
            
            #renaming in case will need to join with data retrieved with bioservises module
            relation['entry1'] = relation['_entry1'] 
            relation['entry2'] = relation['_entry2']
            relation['link'] = relation['type']
        
            relation['value'] = relation['subtypes'][0][1]
            relation['name'] = relation['subtypes'][0][0]

            # selecting relevant keys
            keys_to_keep=['entry1','entry2','link','value', 'name']
            relation={k: v for k, v in relation.items() if k in keys_to_keep}
        
            return relation
            
        relations=rest_kgml_pathway._relations
        relations=[x.__dict__ for x in relations]
        relations = [process_relation(x) for x in relations]
        
        return relations
            
    # Retrieving raw data of kegg pathway
    rest_kgml_pathway = REST.kegg_get(kegg_path_id, "kgml").read()
    rest_kgml_pathway = KGML_parser.read(rest_kgml_pathway)

    # Processing the raw data of nodes (entries) and edges (relations)
    entries = process_entries(rest_kgml_pathway)
    relations = process_relations(rest_kgml_pathway)
    
    return {'entries':entries,
           'relations':relations}

def get_pathway_network_bioservices(kegg_path_id):
    """
    Seems like a slower alternative providing less data than the get_pathway_network_REST()
    Get pathway's nodes (entries) and edges (relations) using bioservices.kegg module:
    returns:
    {'entries': [{entry_id: '...',
                '...':'...'},
                {entry_id: '...',
                '...':'...'},...],
     'relations':[{'entry1': '112',
                   'entry2': '119',
                   'link': 'ECrel',
                   'value': '95',
                   'name': 'compound'},
                  {'entry1': '112',
                   'entry2': '119',
                   'link': 'ECrel',
                   'value': '95',
                   'name': 'compound'},...]
    }
    """
    network_data=KEGG().parse_kgml_pathway(kegg_path_id)
    return network_data

