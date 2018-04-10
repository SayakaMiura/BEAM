from Bio import Phylo

class MegaAlignment():

    def name2seq(self, align_list):
        self.M = align_list	
        self.clone_order=[]
        self.clone_seq = {}
        Name=''		
        for Seq in self.M:
          Seq=Seq.strip()		
          if Seq!='':		
            if Seq[0]=='#' and Seq.find('#MEGA')==-1 and Seq.find('#mega')==-1:
                if Seq!='#hg19' and Seq!='#Normal':
                    self.clone_order.append(Seq)
                   				
                    self.clone_seq[Seq]=''
                    Name=Seq
                else: Name=''  					
            elif Name!='':
                self.clone_seq[Name] += Seq	
        return self.clone_order, self.clone_seq	
    def name2seq_with_normal(self, align_list):	
        self.M = align_list	
        self.clone_order=[]
        self.clone_seq = {}
        Name=''		
        for Seq in self.M:
          Seq=Seq.strip()		
          if Seq!='':		
            if Seq[0]=='#' and Seq.find('#MEGA')==-1 and Seq.find('#mega')==-1:
                    self.clone_order.append(Seq)
                   				
                    self.clone_seq[Seq]=''
                    Name=Seq					
            elif Name!='':
                self.clone_seq[Name] += Seq	
        return self.clone_order, self.clone_seq		
    def ReadMegSeq(self,Meg): 
      Meg=open(Meg,'r').readlines()
      Read='s'
      out2=''
      NameOrder=[]
      Name2Seq={}
      SeqNum=0
      for i in Meg:
        if i[0]=='#' and i.strip()!='#MEGA' and i.strip()!='#mega' :
          	
            Read='n'
            Name=i.strip()
            NameOrder.append(Name)
            Name2Seq[Name]=''
            SeqNum+=1
        elif Read=='n': Name2Seq[Name]+=i.strip()
        elif Read=='s': out2+=i
      return NameOrder, Name2Seq		
    def ReadFas(self, Meg): 
      Read='s'
      out2=''
      NameOrder=[]
      Cell2Seq={}
      SeqNum=0
      for i in Meg:
        if i[0]=='>' :
          	
            Read='n'
            Name=i.strip()[1:]
            NameOrder.append(Name)
            Cell2Seq[Name]=''
            SeqNum+=1
        elif Read=='n': Cell2Seq[Name]+=i.strip()
    
      return NameOrder, Cell2Seq	
    def find_identical_cellseq(self,Root,Fas):
       IdenCellLs=[]
       NameOrder, Cell2Seq = self.ReadFas(Fas)
       RootSeq=Cell2Seq[Root]
       for Cell in Cell2Seq:
           if Cell!=Root:
                 OtherSeq=Cell2Seq[Cell]		   
                 DiffNum=self.CountDifNum(RootSeq,OtherSeq)				 
                 if DiffNum==0: IdenCellLs.append(Cell)				 
       return IdenCellLs
    def find_identical_cellseq_excmiss(self,Root,Cell2Seq,CutSeqDif):
       IdenCellLs=[]
       RootSeq=Cell2Seq[Root]
       for Cell in Cell2Seq:
           if Cell!=Root:
                 OtherSeq=Cell2Seq[Cell]		   
                 DiffNum=self.CountDifNum_excMiss(RootSeq,OtherSeq)				 
                 if DiffNum<=CutSeqDif: IdenCellLs.append(Cell)				 
       return IdenCellLs	   
    def AddNormal(self, seqs): 
        CellOrder, Cell2Seq = self.name2seq(seqs)	
        Len=len(Cell2Seq[CellOrder[0]])	
        seqs+=['#Normal\n',('A'*Len)+'\n']
        return seqs
    def meg2fas(self, meg_seq):
        Fas=[]
        Cell_order, Cell_seq = self.name2seq_with_normal(meg_seq)
        for Cell in Cell_order:
            Fas += ['>'+Cell[1:], Cell_seq[Cell]]
        return Fas	
    def GetDiffPosi(self, CellLs, Cell2Seq):
         Cell0=CellLs[0]	 
         if Cell0[0]!='#': Cell0='#'+CellLs[0]	
         Len=len(Cell2Seq[Cell0])
         c=0
         DiffPosi=[]		 
         while c<Len:
             NucLs=[]		 
             for Cell in CellLs:
                  if Cell[0]!='#': Cell='#'+Cell 			 
                  if Cell2Seq[Cell][c]!='?': NucLs.append(Cell2Seq[Cell][c])
             NucLs=list(set(NucLs))
             if len(NucLs)>1: DiffPosi.append(c)
             c+=1
         return DiffPosi	
    def GetDiffPosi1(self, CellLs, Cell2Seq):
         Cell0=CellLs[0]	 
         if Cell0[0]!='#': Cell0='#'+CellLs[0]	
         Len=len(Cell2Seq[Cell0])
         c=0
         DiffPosi=[]		 
         while c<Len:
             NucLs=[]		 
             for Cell in CellLs:
                  if Cell[0]!='#': Cell='#'+Cell 			 
                  NucLs.append(Cell2Seq[Cell][c])
             NucLs=list(set(NucLs))
             if len(NucLs)>1: DiffPosi.append(c)
             c+=1
         return DiffPosi
	 
             				  
    def GetSharePosi1(self, sequence, ShareNuc):
        self.name2seq0 = sequence	
        for i in self.name2seq0:
           Name=i
        SharePosi=[]
        Len=len(self.name2seq0[Name])
        c=0
        while c<Len:
            AllMut='y'
            for Ori in self.name2seq0:
              if Ori!='#hg19' and Ori!='#Normal':
               Nuc=self.name2seq0[Ori][c]
               if Nuc!=ShareNuc: AllMut='n'
            if AllMut=='y': 
                 SharePosi.append(c)
            c+=1
        return SharePosi
    def GetSharePosi1_excMis(self, sequence, ShareNuc):
        self.name2seq0 = sequence	
        for i in self.name2seq0:
           Name=i
        SharePosi=[]
        Len=len(self.name2seq0[Name])
        c=0
        while c<Len:
            AllMut='y'
            for Ori in self.name2seq0:
              if Ori!='#hg19' and Ori!='#Normal':
               Nuc=self.name2seq0[Ori][c]
               if Nuc!=ShareNuc and Nuc!='?': AllMut='n'
            if AllMut=='y': 
                 SharePosi.append(c)
            c+=1
        return SharePosi		
		
    def UpMeg(self, Name2Seq0, NameLs):
            if NameLs==[]:
                 for Name in Name2Seq0:
                     NameLs.append(Name)			 
            out=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']
            for Name in NameLs:
                if Name[0]!='#': Name='#'+Name
                out+=[Name,Name2Seq0[Name]]
            return out 
			
   
		
    def CountDifNum(self, Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif	

    def CountDifNum_excMiss(self, Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!='?' and Seq1[c]!='?' and Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif				
			
    def CountAdditionalMut(self, seq1,seq2):
        Len=len(seq1)
        c=0
        DerMut=0		
        while c<Len:
            if seq1[c]=='T' and seq2[c]=='A': DerMut+=1
            c+=1
        return DerMut	
		
    def ModSeq(self, CSeq0,ChangePosi,ChanNuc,Len):					  
                      c=0
                      CutCloSeq=''					  
                      while c<Len:
                          Code1=c in ChangePosi						  
                          if Code1==True: CutCloSeq+=ChanNuc
                          else: CutCloSeq+=CSeq0[c]
                          c+=1
                      return CutCloSeq	
					  
    def GetMutPos(self, Seq):
      TMP=[]
      Len=len(Seq)
      c=0
      while c<Len:
        if Seq[c]=='T': TMP.append(c)
        c+=1
      return TMP 
	
    def get_mega_alignment_string(self, SeqLs):
        result = ''
        for item in SeqLs:#self._mega_seqs:
            result += item + "\n"
        return result
    
    def save_mega_alignment_to_file(self, filename, SeqLs):
        destination = open(filename,'w')
        destination.write(self.get_mega_alignment_string(SeqLs))
        destination.close()        	  