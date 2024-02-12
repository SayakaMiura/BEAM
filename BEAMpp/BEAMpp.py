from PPcomputer.PredictCellGenotype import PredictCellGenotype
from alignments.MegaAlignment import MegaAlignment
from ML.TreeAnalizer import TreeAnalizer
import sys
import os
import shutil
from shutil import copy
import glob

Align=MegaAlignment()
tree_analyzer = TreeAnalizer()  
InATGC=sys.argv[1]
Type=4
Normal=sys.argv[2]#'normal'#'Normal'
Cut2=0.7 #PP cut-off
dir = os.getcwd()
CSV=glob.glob(InATGC[:-6]+'_rooted_PP*') 
for i in CSV:
    os.remove(i)
if Type==4:
    Tree=InATGC[:-6]+'.nwk'
    Go='y'
    if os.path.exists(Tree)!=True and os.path.exists('FastTree.exe')==True:
             os.system('FastTree.exe -nt '+InATGC+ ' > '+Tree)
    elif os.path.exists(Tree)!=True:
            Go='n'
            print ('please provide tree file',Tree)
    if Go=='y':            
        tree_analyzer.root_tree(Tree,Normal)
        Rtree=Tree[:-4]+'_rooted.nwk'
        SNVc,ExtraLs=Align.Nuc2BEAMinWnNoFil(InATGC,InATGC[:-6]+'.meg',Normal)
        tree_analyzer.prune_tree(Rtree,ExtraLs)
        Rtree=Rtree[:-4]+'_prune.nwk'
        In=InATGC[:-6]+'.meg' 
        OutMegFile=In[:-4]+'_BEAM.meg'
        InFile=In
        In = open(In,'r').readlines()
        In.append('#'+Normal+'\n'+('A'*SNVc))
        InMeg=InFile[:-4]+'1.meg'
        Align.save_mega_alignment_to_file(InMeg, In)
        print('correct FPs and FNs')
        PP2 = PredictCellGenotype('Correct', InMeg, Cut2, Rtree, Normal)
        PP2.OutTpp(OutMegFile[:-4])
       
      #  MEGAseqs_Corrected = PP2.Correct_error5()
      #  Align.save_mega_alignment_to_file(OutMegFile[:-4]+'Correct1.fasta', MEGAseqs_Corrected)
      #  Tree=OutMegFile[:-4]+'Correct1.nwk'
      #  if os.path.exists('FastTree.exe')==True:
      #      os.system('FastTree.exe -nt '+OutMegFile[:-4]+'Correct1.fasta'+ ' > '+Tree)   
      #  tree_analyzer.root_tree(Tree,Normal)             

CSV=glob.glob(InATGC[:-6]+'_rooted_prune_PP*') 
for i in CSV:
    os.remove(i)
    

