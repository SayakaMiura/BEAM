from PPcomputer.PredictCellGenotype import PredictCellGenotype
from alignments.MegaAlignment import MegaAlignment
#from clone_annotation.ClusterQuality3 import ClusterQuality3
import sys
import os
import shutil
from shutil import copy
import glob
import Functions3

Align=MegaAlignment()
In=sys.argv[1]
if In[-6:]=='.fasta':
    Normal='normal'
    SNVc,ExtraLs=Align.Nuc2BEAMinWnNoFil(In,In[:-6]+'.meg',Normal)
    In=In[:-6]+'.meg'
#TreeInf=In[:-4]+'.nwk'
TreeInf=In[:-4]+'_rooted.nwk'
if os.path.exists(TreeInf)!=True:
    Functions3.FastTree(In,'normal')
OutMegFile=In[:-4]+'_BEAM.meg'

Cut2=0.7 #PP cut-off

dir = os.getcwd()
CSV=glob.glob('*.csv')
for i in CSV:
    os.remove(i)

InFile=In
In = open(In,'r').readlines()

print('correct FPs and FNs')
PP2 = PredictCellGenotype('Correct', In, Cut2, TreeInf)
MEGAseqs_Corrected = PP2.Correct_error5()
Cell2PPselected = PP2.get_PP_for_selected_nuc_corr()
Align.save_mega_alignment_to_file(OutMegFile[:-4]+'Correct1.meg', MEGAseqs_Corrected)
#print (OutMegFile)
Tree_Corrected=OutMegFile[:-4]+'Correct1_rooted.nwk'
Functions3.FastTree(OutMegFile[:-4]+'Correct1.meg','normal')

print('correct FPs and FNs 2')
Cut2=0.7 #PP cut-off
PP2 = PredictCellGenotype('Correct', MEGAseqs_Corrected, Cut2, Tree_Corrected)
MEGAseqs_Corrected_1 = PP2.Correct_error5()
Cell2PPselected = PP2.get_PP_for_selected_nuc_corr()
Align.save_mega_alignment_to_file(OutMegFile[:-4]+'Correct2.meg', MEGAseqs_Corrected_1)
Tree_Corrected=OutMegFile[:-4]+'Correct2_rooted.nwk'
Functions3.FastTree(OutMegFile[:-4]+'Correct2.meg','normal')


Go='n'
if Go=='cla':
   CladeCut=10
   boothreshold=1
   ID2cellLs,Cla2Boo,Anc2DecLs,tree,clade2cellLs=Functions3.GetCladeTip(OutFasNwk,CladeCut,boothreshold)
   Functions3.GetOut(OutFasNwk[:-4]+'_tree.txt',str(tree))    
   print ('add sequence ID')
   STICellLs,Cell2SeqHap=Functions3.ReadFasSeq(STIfas)
   out=['Cell\tTopHapID\n']
   for Cl in clade2cellLs:
    CellLs=clade2cellLs[Cl]
    for Cell in CellLs:
        out.append(Cell+'\t'+Cl+'\n')
   Functions3.GetOut(OutFasNwk[:-4]+'_clone.txt',''.join(out))   
   print ('make consensus clone seq')
   ConsenFas=OutFasNwk[:-4]+'_cloneConsensus.fasta'
   Functions3.makeConsenSeqATGC(clade2cellLs,Cell2SeqHap,ConsenFas,CloneCut1,OutSeq0,0.2)   #refined seqHap makeConsenSeqATGC(clade2cellLs,CellIn2Seq,ConsenFas,MinClone,OutSeq)
   Functions3.FastTree(ConsenFas,'Outgroup') 
elif Go=='ori':
 #print('Compute final PP')
 #PP2 = PredictCellGenotype('Correct', MEGAseqs_Corrected_1, Cut2)
 #MEGAseqs_Corrected_2 = PP2.Correct_error5()
 Cell2PPselected = PP2.get_PP_for_selected_nuc_corr()

 print('clone annotation')
 In0=InFile
 #OutMegFile=In0[:-4]+'_BEAM.meg'

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
 if Cut3>5: Cut3=5
 Tree_Corrected=OutMegFile[:-4]+'Correct2_rooted.nwk'
 Functions3.FastTree(OutMegFile[:-4]+'Correct2.meg','normal')
 Clones = ClusterQuality3(MEGAseqs_Corrected_1,Cut3,Cell2PPls,In,OutMegFile)
 Clones.adjust_cell_genotype1(Tree_Corrected)
 Clones.fill_cloneseq()

#os.remove(OutMegFile[:-4]+'Correct1.meg')
#os.remove(OutMegFile[:-4]+'Correct2.meg')
#
#os.remove('newtree.nwk')
#os.remove('test.nwk')
#os.remove('test1.nwk')
 os.remove(OutMegFile[:-4]+'_withCloneID.meg')  
CSV=glob.glob('All_alignment_PPseq-*.csv')
for i in CSV:
    os.remove(i)
  