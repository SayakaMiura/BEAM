
from parsers.ancestorParser import ancestorParser
from parsers.nodeMapParser import nodeMapParser
from alignments.MegaAlignment import MegaAlignment
import os
import tempfile

"""
    A wrapper around MEGA-CC for Maximum Parsimony tree construction
    
    Given a MEGA sequence alignment, an instance of this class can generate
    an MP tree (or multiple equally parsimonious trees) using MEGA-CC and 
    construct ancestral sequences for all trees generated. The ancestral sequences
    are stored in the ancestral_states_list property, where each element in the list
    also stores the newick tree for that set of ancestral states
"""

class MegaAncestor(object):
    
    def __init__(self):
        self.mao_file = 'ancestral_seqs_ML_nucleotideNuc.mao' #this is ML
        self._alignment_file = ''
        self._input_tree_file = ''        
        self._summary = ''
        self._ancestor_file = ''
  
        self._mega_id = ''
        self._ancestor_states_dic = {}
        self.anc2dec = {}		
  
        self._temp_dir = tempfile.gettempdir() + os.sep
        	
        
    def __del__(self):
        self._cleanup_temp_files()
        
    def do_mega_ancestor(self):
        
        print('infer ancestral sequences')
        result = False
        alignment_builder = self._alignment_file
        tree_builder = self._input_tree_file		
        self._update_file_names('Ancestor')
	
        Align = MegaAlignment()
        Align.save_mega_alignment_to_file(self._alignment_file, alignment_builder) ### 
        self.save_str_to_file(tree_builder, self._input_tree_file)		

        cl = self._command_line_string()
        os.system(cl)
        if os.path.isfile(self._ancestor_file) == True:
            result = True
    
        return result
        
    def _update_file_names(self, mega_id):        
        
        print('executing megacc ancestral sequence construction in ' + self._temp_dir)
        self._mega_id = mega_id
        self._alignment_file = self._temp_dir + mega_id + '.meg'
        self._input_tree_file = self._temp_dir + mega_id + '.nwk'		
        self._ancestor_file = self._temp_dir + mega_id + '_ances.csv'#'.nwk'
        self._nodeMap_file = self._temp_dir + mega_id + '_ances_nodeMap.txt'		
        
    def _command_line_string(self):
         
        return 'megacc -a ' + self._mao_file + ' -d ' + self._alignment_file + ' -t ' + self._input_tree_file + ' -o ' + self._ancestor_file
 
    def retrieve_ancestor_states(self): 
            self.do_mega_ancestor()
            nodeMap_file = self._nodeMap_file
            nodeMap_parser = nodeMapParser()
            nodeMap_parser.input_file_name = nodeMap_file
            if not nodeMap_parser.parse() == True:
                IOError('failed to parse nodeMap file')
            self.offspring2ancestor, self.code2cell, self.cell2code = nodeMap_parser.get_nodeMap_states()
            self.A2D, self.A2lin = nodeMap_parser.get_anc2dec_lin()			

            file = self._ancestor_file	

            parser = ancestorParser()
            parser.input_file_name = file
            if not parser.parse() == True:
                IOError('failed to parse ancestral states file')
            self._ancestor_states_list = parser.get_ancestor_states()
       	   
            return 	self._ancestor_states_list, self.offspring2ancestor, self.cell2code, self.code2cell	
    def report_anc2dec_lin(self):			
        return self.A2D, self.A2lin        
    def _cleanup_temp_files(self):    
        os.remove(self._alignment_file)
        os.remove(self._input_tree_file)		
        os.remove(self._ancestor_file)
        os.remove(self._nodeMap_file)		
        summary_file = self._temp_dir + self._mega_id + '_ances_summary.txt'
        os.remove(summary_file)
        os.remove(self._temp_dir + self._mega_id + '_ances.nwk')		
       
        
    @property
    def mao_file(self):
        return self._mao_file
    
    @mao_file.setter
    def mao_file(self, value):
        self._mao_file = value
        
    @property 
    def alignment_file(self):
        return self._alignment_file
    
    @alignment_file.setter
    def alignment_file(self, value):
        self._alignment_file = value

    @property 
    def input_tree_file(self):
        return self._input_tree_file
    
    @input_tree_file.setter
    def input_tree_file(self, value):
        self._input_tree_file = value

        

    
    
    def _get_ancestral_states_files(self):
        result = []
        filename = self._temp_dir + self._mega_id + '_PP.csv'
        if os.path.isfile(filename) == True:
            result.append(filename)
  
        return result
    
    @property
    def num_trees(self):
        return len(self._pp_states_list)
    
    @property 
    def ancestral_states_list(self):
        return self._pp_states_list
        

    
    def alignment_least_back_parallel_muts(self, remove_duplicates = True):
        print('finding alignment with least parallel and back mutations...')
        files = self._get_ancestral_states_files()
        seq_maker = MakeAncSeqMPMin()
	
        result = seq_maker.get_best_alignment(files, self._mega_id, remove_duplicates, self.newick_trees)           
        return result

    def save_str_to_file(self, String, Out_file_name):
         OutF=open(Out_file_name, 'w')
         OutF.write(String)
         OutF.close()		 