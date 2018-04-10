import os
import shutil
from shutil import copy
from alignments.MegaAlignment import MegaAlignment
from ML.MegaML import MegaML
from ML.MegaPP import MegaPP
from ML.MegaAncestor import MegaAncestor
from ML.TreeAnalizer import TreeAnalizer

class PredictCellGenotype():
    def __init__(self, id0, seqs, PPcut): #seqs is list format of mega alignment
      self.PPcut=PPcut	  
      Align = MegaAlignment()
      tree_builder = MegaML() 
      tree_analyzer = TreeAnalizer()  
      self.InMeg = Align.AddNormal(seqs)
      status = tree_builder.do_mega_ml(self.InMeg, id0)
      if status == True:
               tree = tree_builder.newick_trees
      else:
                print 'failed to run megaML'	  
      self.Tree_rooted = tree_analyzer.RootTree(tree, 'Normal')
   
    def Compute_PP(self):
        Align = MegaAlignment()	
        tree_analyzer = TreeAnalizer()  
        PP_builder = MegaPP() 		
        Input=self.Tree_rooted
        print 'input for inferring missing',Input
        Meg=self.InMeg
        NameOrder, self.Cell2megSeq=Align.name2seq_with_normal(Meg)
        id = 'All_alignment' 
        status = PP_builder.do_mega_pp(self.InMeg, self.Tree_rooted, id)
        if status == True:
                    self.Cell2PP = PP_builder.retrieve_pp_states() 		

        else:
                    print 'failed to run megaPP'						
              			
    def get_PP_for_selected_nuc_corr(self):
           return self.Cell2PPsel 	

    def Correct_error5(self):
        print 'correct errors'	  	
        self.Compute_PP()
        New_seq=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']
        self.Cell2PPsel={}		
        for Cell in self.Cell2megSeq:
          Cell=Cell[1:]
          if Cell!='Normal':		  
           self.Cell2PPsel[Cell]=[]		
           original_seq=self.Cell2megSeq['#'+Cell]
           if self.Cell2PP.has_key(Cell)!=True: seq=original_seq
           else:			 		 
             self.normal_seq = self.Cell2megSeq['#Normal']			 	 
             Nuc2PPlist=self.Cell2PP[Cell]
             seq=''
             c=0
             Len=len(original_seq)		 
             while c<Len:			 
                NucPP=	Nuc2PPlist[c]					 
                Nuc=NucPP.split('\t')[0]
                PP=float(NucPP.split('\t')[1])
                self.Cell2PPsel[Cell].append(PP)				
                if PP>self.PPcut: seq+=Nuc
                else: seq+=original_seq[c]
                c+=1				 
           New_seq+=['#'+Cell,seq]
        return New_seq		
  		