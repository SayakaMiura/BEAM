from PPcomputer.PredictCellGenotype import PredictCellGenotype
from alignments.MegaAlignment import MegaAlignment
from clone_annotation.ClusterQuality3 import ClusterQuality3
import sys
import os
import shutil
from shutil import copy
import glob

Align=MegaAlignment()
In=sys.argv[1]
OutMegFile=In[:-4]+'_BEAM.meg'

Cut2=0.7 #PP cut-off

dir = os.getcwd()
InFile=In
In = open(In,'r').readlines()

print('correct FPs and FNs')
PP2 = PredictCellGenotype('Correct', In, Cut2)
MEGAseqs_Corrected = PP2.Correct_error5()
Cell2PPselected = PP2.get_PP_for_selected_nuc_corr()
Align.save_mega_alignment_to_file(OutMegFile[:-4]+'Correct1.meg', MEGAseqs_Corrected)

print('correct FPs and FNs 2')
PP2 = PredictCellGenotype('Correct', MEGAseqs_Corrected, Cut2)
MEGAseqs_Corrected_1 = PP2.Correct_error5()
Cell2PPselected = PP2.get_PP_for_selected_nuc_corr()
Align.save_mega_alignment_to_file(OutMegFile[:-4]+'Correct2.meg', MEGAseqs_Corrected_1)

print('Compute final PP')
PP2 = PredictCellGenotype('Correct', MEGAseqs_Corrected_1, Cut2)
MEGAseqs_Corrected_2 = PP2.Correct_error5()
Cell2PPselected = PP2.get_PP_for_selected_nuc_corr()

print('clone annotation')
In0=InFile
OutMegFile=In0[:-4]+'_BEAM.meg'
dir = os.getcwd()
MEGAseqs_Corrected=OutMegFile[:-4]+'Correct2.meg'
MEGAseqs_Corrected_1 = open(MEGAseqs_Corrected,'r').readlines()
In=open(In0,'r').readlines()
CellLs, Cell2Seq = Align.name2seq(In)
Cell2PPls={}
dir=os.getcwd()
CellC=1
for Cell in CellLs:
     PPoutF=dir+'\\All_alignment_PPseq-'+str(CellC)+'.csv'	 
     PPout=open(PPoutF,'r').readlines()
     CellN=PPout[0].split('\"')[1]	 
     PPout=PPout[3:]	 
     PPls=[]	 
     for i in PPout:
         i=i.strip().split(',')
         PP={'A':float(i[1]),'T':float(i[4])}
         PPls.append(PP)
     Cell2PPls[CellN]=PPls
     os.remove(PPoutF)	
     CellC+=1

print('adjust cell sequences')	 
Cut3=0.02*len(Cell2Seq[CellLs[0]])
Clones = ClusterQuality3(MEGAseqs_Corrected_1,Cut3,Cell2PPls,In,OutMegFile)
Clones.adjust_cell_genotype1()
Clones.fill_cloneseq()

os.remove(OutMegFile[:-4]+'Correct1.meg')
os.remove(OutMegFile[:-4]+'Correct2.meg')
os.remove(OutMegFile[:-4]+'_withCloneID.meg')
os.remove('newtree.nwk')
os.remove('test.nwk')
os.remove('test1.nwk')
CSV=glob.glob('*.csv')
for i in CSV:
    os.remove(i)