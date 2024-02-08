#from alignments.FreqToMegaSeq import FreqToMegaSeq
from parsers.PPParser import PPParser
#from MakeAncSeqMPMin import MakeAncSeqMPMin
from alignments.MegaAlignment import MegaAlignment
import os
import tempfile
import shutil
import glob
"""
    A wrapper around MEGA-CC for Maximum Parsimony tree construction
    
    Given a MEGA sequence alignment, an instance of this class can generate
    an MP tree (or multiple equally parsimonious trees) using MEGA-CC and 
    construct ancestral sequences for all trees generated. The ancestral sequences
    are stored in the ancestral_states_list property, where each element in the list
    also stores the newick tree for that set of ancestral states
"""

class MegaPP(object):
    
    def __init__(self):
        self.mao_file = 'Living_seq_Nuc.mao' #this is ML
        self._alignment_file = ''
        self._input_tree_file = ''        
        self._summary = ''
        self._pp_file = ''
        self._num_trees = 0
        self._mega_id = ''
        self._pp_states_list = []
    
        self._temp_dir = tempfile.gettempdir() + os.sep
        	
        
  #  def __del__(self):
  #      print('deletion')	
  #      self._cleanup_temp_files()
        
    def do_mega_pp(self, alignment_builder, tree_builder, mega_id):
        
        print('computing PP')
        result = False
        self._update_file_names(mega_id)
    
        Align = MegaAlignment()
       # Align.save_mega_alignment_to_file(self._alignment_file, alignment_builder) ### 
       # self.save_str_to_file(tree_builder, self._input_tree_file)		
        self._input_tree_file=tree_builder
        cl = 'megacc -a Living_seq_Nuc.mao -d '+alignment_builder+' -t '+self._input_tree_file+' -o '+self._input_tree_file[:-4]+'_PP.csv --all-seqs'#self._command_line_string()
        self.PPfileLs=glob.glob(self._input_tree_file[:-4]+'_PPseq-*.csv')	
        if len(self.PPfileLs)==0:        
            os.system(cl)
            print (cl)
            self.PPfileLs=glob.glob(self._input_tree_file[:-4]+'_PPseq-*.csv')	
       # open('a','r').readlines()        
	
 		
     
        if self.PPfileLs!=[]:
            result = True
        else:  result = False   
      #      for PPfile in PPfileLs:			
                
      #           shutil.copyfile(PPfile,PPfile.split('\\')[-1])					 
    #    print result			
    
        return result
        
    def _update_file_names(self, mega_id):        
        
        print('executing megacc parsimony tree construction in ' + self._temp_dir)
        self._mega_id = mega_id
        self._alignment_file = self._temp_dir + mega_id + '.meg'
        self._input_tree_file = self._temp_dir + mega_id + '.nwk'		
        self._pp_file = self._temp_dir + mega_id + '_PP.csv'#'.nwk'
        
    def _command_line_string(self):
    
        print(self._pp_file)	   
        return 'megacc -a ' + self._mao_file + ' -d ' + self._alignment_file + ' -t ' + self._input_tree_file + ' -o ' + self._pp_file+' --all-seqs'
  
    def retrieve_pp_states(self):  
  
       # PPfileLs=glob.glob(self._pp_file[:-4]+'seq-*.csv')
        Cell2states_list={}	
        Cell2Tstates_list={} 
        Cell2Astates_list={}        
        print('parsing pp file...',len(self.PPfileLs))
        for file in self.PPfileLs:
            parser = PPParser()
            parser.input_file_name = file
			
            if not parser.parse() == True:
                IOError('failed to parse ancestral states file')
            Cell=parser.get_cellName()				
           # states_list = parser.get_pp_states()
            Tpp_list=parser.get_Tpp_states()
           # Cell2states_list[Cell]=states_list	
            Cell2Tstates_list[Cell]= Tpp_list 
            Cell2Astates_list[Cell]=parser.get_App_states()            
   
      #  self.__del__()		   		
      #  return 	Cell2states_list,Cell2Tstates_list	
        return 	Cell2Astates_list,Cell2Tstates_list	        
    def _cleanup_temp_files(self):
      
        PPfileLs=glob.glob(self._pp_file[:-4]+'seq-*.csv')
        for PPfile in PPfileLs:
           	os.remove(PPfile)	
        FileLs=glob.glob(self._temp_dir+'*.nwk')
        for File in FileLs:
           	os.remove(File)				
        FileLs=glob.glob(self._temp_dir+'*.meg')
        for File in FileLs:
           	os.remove(File)				

        summary_file = self._temp_dir + self._mega_id + '_PP_summary.txt'
        os.remove(summary_file)

        
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