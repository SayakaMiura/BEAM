
import os, shutil
from shutil import copy
from alignments.MegaAlignment import MegaAlignment
from ML.MegaML import MegaML
from ML.MegaPP import MegaPP
from ML.MegaAncestor import MegaAncestor
from ML.TreeAnalizer import TreeAnalizer
import random
import scipy.stats

class ClusterQuality3:

    def __init__(self, seqs, num_support_position, Cell2PPls, initial_seq_builder, OutFileName):
       	
        self.cut = num_support_position
        Align = MegaAlignment()
        self.ini_seqs_builder = seqs
        self.CellLs, self.Cell2Seq = Align.name2seq(seqs)
  	
        self.SNVnum = len(self.Cell2Seq[self.CellLs[0]])
      
        self.InMeg = Align.AddNormal(seqs)
        IniCellLs, self.Cell2iniSeq = Align.name2seq(initial_seq_builder)
        self.Cell2PPls = Cell2PPls
        self.out_file_name = OutFileName

    def fill_cloneseq(self):
        Align=MegaAlignment()	

        Cell2BestSeq={}
        
        for Clone in self.Clone2CellLs:
              CellLs=self.Clone2CellLs[Clone]
              DiffNucPosiLs=Align.GetDiffPosi(CellLs, self.Cell2Seq)
         
              Posi2Nuc={}
              for Posi in DiffNucPosiLs:
                   Nuc2PPls={'A':[],'T':[]}
                   for Cell in CellLs:				   
                         Nuc2PPls['A'].append(self.Cell2PPls[Cell][Posi]['A'])
                         Nuc2PPls['T'].append(self.Cell2PPls[Cell][Posi]['T'])
                   TAve=sum(Nuc2PPls['T'])/len(Nuc2PPls['T'])
                   AAve=sum(Nuc2PPls['A'])/len(Nuc2PPls['A'])
                   if TAve>AAve: Posi2Nuc[Posi]='T'
                   else: Posi2Nuc[Posi]='A'
              for Cell in CellLs:
                    CellSeq=self.Cell2Seq['#'+Cell]
              			
                    c=0
                    NewSeq=''
                    while c<self.SNVnum:
                        if CellSeq[c]=='?': NewSeq+='?' 
                        elif (c in Posi2Nuc)!=True: NewSeq+=CellSeq[c]
                        else: NewSeq+=Posi2Nuc[c]
                        c+=1
                    Cell2BestSeq['#'+Cell]=NewSeq	
        self.save_with_cloneID(Cell2BestSeq)
        self.save_without_cloneID(Cell2BestSeq)		
    def adjust_cell_genotype1(self):
        Align = MegaAlignment()
        tree_builder = MegaML()
        tree_analyzer = TreeAnalizer() 		
        status = tree_builder.do_mega_ml(self.InMeg, 'Noresun')
        if status == True:
            tree1 = tree_builder.newick_trees
        else:
            print('failed to run megaML')
        Tree_Rooted1 = tree_analyzer.RootTree_rootBottom(tree1, 'Normal')
	
        InferAncestor = MegaAncestor()
        InferAncestor.alignment_file = self.InMeg
        InferAncestor.input_tree_file = Tree_Rooted1	
     			
        self.ancestor_states, self.offspring2ancestor, cell2code, self.code2cell = InferAncestor.retrieve_ancestor_states()
        ancestor2offspring, self.node2cellclade = InferAncestor.report_anc2dec_lin()	
        for code in self.code2cell:
            if self.code2cell[code].find('Node_') == -1:
                self.node2cellclade[code] = [self.code2cell[code]]	
	
        self.Cellclade_withSupport2SupportCount = self.count_support3()		
        self.Clone2CellLs = self.get_clone3()


    def save_without_cloneID(self, Cell2BestSeq):
        Align = MegaAlignment()
        tree_builder = MegaML()
        tree_analyzer = TreeAnalizer()		
        BestSeq_builder_3 = Align.UpMeg(Cell2BestSeq, self.CellLs)
        Align.save_mega_alignment_to_file(self.out_file_name, BestSeq_builder_3)
		

    def save_with_cloneID(self, CellSeqDic):
        Align = MegaAlignment()
        outSeq_builder = ['#MEGA', '!Title SNVs;', '!Format datatype=dna;', ' ']
        for Clone in self.Clone2CellLs:
            CellLs = self.Clone2CellLs[Clone]
            for Cell in CellLs:
                outSeq_builder += ['#' + Cell + '_{' + Clone + '}', CellSeqDic['#' + Cell]] #######change

        Align.save_mega_alignment_to_file(self.out_file_name[:-4] + '_withCloneID.meg', outSeq_builder)

    def report_cloneID(self):
        return self.Clone2CellLs

  
    def get_clone3(self):

        Clo2CellLs = {}
        self.CloneID2NodeID = {}
        CloID = 1
        AncID=1		
        Cell2CladeID={}		
      	
        for Clade in self.Cellclade_withSupport2SupportCount:
          if self.Cellclade_withSupport2SupportCount[Clade]>=self.cut:		
            CellLs = self.node2cellclade[Clade]	
            for Cell in CellLs:
                if (Cell in Cell2CladeID)!=True: Cell2CladeID[Cell]=[]
                Cell2CladeID[Cell].append(Clade)	
        for Cell in self.Cell2Seq:
            Cell=Cell[1:]		
            if Cell!='Normal' and (Cell in Cell2CladeID)!=True: Cell2CladeID[Cell]=['root']
   
        Cell2CladeIDstr={}
        CladeIDstr2CloType={}
        CladeIDstrLs=[]		
        for Cell in Cell2CladeID:
               CladeID=Cell2CladeID[Cell]
               CladeID.sort()
          
               CladeIDstr='C'.join(map(str, CladeID))
          	   
               Cell2CladeIDstr[Cell]= CladeIDstr
               if CladeIDstrLs.count(CladeIDstr)==0: CladeIDstrLs.append(CladeIDstr)
  
        for Clade in CladeIDstrLs:
		
          Clade0=Clade.split('C')
          if Clade0==['root']: CloTy='root'
          else:		  
            Tip='y'			
            for Clade1 in CladeIDstrLs:
                if Clade1!=Clade:
                    Clade10=Clade1.split('C') 
                    All='y'					
                    for Cl in Clade0:
                        if 	Clade10.count(Cl)==0: All='n'
                    if All=='y': Tip='n'						
            if Tip=='y': CloTy='Tip'
            else: CloTy='Anc'
          CladeIDstr2CloType[Clade]=CloTy
   
        CloCou=1		
        for CladeIDstr in CladeIDstr2CloType:
              CloType=CladeIDstr2CloType[CladeIDstr]
              CellLs=[]
              for Cell in Cell2CladeIDstr:
                   ID=	Cell2CladeIDstr[Cell]
                   if ID==CladeIDstr: CellLs.append(Cell)				   
              Clo2CellLs[CloType+str(CloCou)]=CellLs
              CloCou+=1			  
   
        return Clo2CellLs		
   
    def count_support3(self):
       
        Clade2Support_clean={}
        for code in self.code2cell:
         # print code		
          if self.code2cell[code]!='Normal' and (code in self.offspring2ancestor)==True:		
            if self.code2cell[code].find('Node_') != -1:
                OffSeq=self.ancestor_states[self.code2cell[code]]
            else:
                OffSeq=self.Cell2Seq['#'+self.code2cell[code]]				
            Anc=self.offspring2ancestor[code]
            AncSeq=self.ancestor_states[self.code2cell[Anc]]
            Len=len(AncSeq)
            c=0
            Supp=0			
            while c<Len:
             
                if AncSeq[c].split('\t')[0]=='A' and OffSeq[c].split('\t')[0]=='T': Supp+=1
                c+=1
            Clade2Support_clean[code]=Supp	
  			
        return Clade2Support_clean
  


    def GetOut(self, OutFile, OutIn):
        OutF = open(OutFile, 'w')
        OutF.write(OutIn)
        OutF.close()