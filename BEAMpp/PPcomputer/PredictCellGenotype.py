import os
import shutil
from shutil import copy
from alignments.MegaAlignment import MegaAlignment
from ML.MegaML import MegaML
from ML.MegaPP import MegaPP
from ML.MegaAncestor import MegaAncestor
from ML.TreeAnalizer import TreeAnalizer

class PredictCellGenotype():
    def __init__(self, id0, seqs, PPcut, Rtree, Normal): #seqs is list format of mega alignment
      self.PPcut=PPcut	  
      Align = MegaAlignment()
      tree_builder = MegaML() 
      tree_analyzer = TreeAnalizer()  
      self.InMeg = seqs#Align.AddNormal(seqs)
   #   status = tree_builder.do_mega_ml(self.InMeg, id0)
   #   if status == True:
   #            tree = tree_builder.newick_trees
   #   else:
   #             print('failed to run megaML')	  
      self.Tree_rooted = Rtree# tree_analyzer.RootTree(tree, 'Normal')
      self.Normal=Normal
    def Compute_PP(self):
        Align = MegaAlignment()	
        tree_analyzer = TreeAnalizer()  
        PP_builder = MegaPP() 		
        Input=self.Tree_rooted
        print('input for inferring missing',Input)
        Meg=self.InMeg
        NameOrder, self.Cell2megSeq=Align.name2seq_with_normal(open(Meg,'r').readlines())
        id = 'All_alignment' 
        status = PP_builder.do_mega_pp(self.InMeg, self.Tree_rooted, id)
        if status == True:
                   # self.Cell2PP,self.Cell2Tpp = PP_builder.retrieve_pp_states() 
                    self.Cell2App,self.Cell2Tpp = PP_builder.retrieve_pp_states()                     

        else:
                    print('failed to run megaPP')						
              			
    def get_PP_for_selected_nuc_corr(self):
           return self.Cell2PPsel 	
    def OutTpp(self, OutFname):
       # print('compute PP for observing mutant base...')	  	
        self.Compute_PP() 
        print('output mutant PP...')	        
        out=[]
        for Cell in self.Cell2Tpp:
            out.append(Cell+'\t'+'\t'.join(map(str,self.Cell2Tpp[Cell]))+'\n')
        OutF=open(OutFname+'_mut.txt','w')
        OutF.write(''.join(out))
        OutF.close()
        print('output non-mutant PP...')	        
        out=[]
        for Cell in self.Cell2App:
            out.append(Cell+'\t'+'\t'.join(map(str,self.Cell2App[Cell]))+'\n')
        OutF=open(OutFname+'_nonmut.txt','w')
        OutF.write(''.join(out))
        OutF.close()        
        
    def Correct_error5(self):
        print('correct errors')	  	
        self.Compute_PP()
        New_seq=[]
        self.Cell2PPsel={}		
        for Cell in self.Cell2megSeq:
          Cell=Cell[1:]
          if Cell!=self.Normal:		  
           self.Cell2PPsel[Cell]=[]		
           original_seq=self.Cell2megSeq['#'+Cell]
           if (Cell in self.Cell2PP)!=True: seq=original_seq
           else:			 		 
            # self.normal_seq = self.Cell2megSeq['>'+self.Normal]			 	 
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
           New_seq+=['>'+Cell,seq]
        New_seq+=['>'+self.Normal,('A'*Len)]   
        return New_seq		
  		