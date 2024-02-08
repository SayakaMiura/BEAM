import os
import sys
from Bio import Phylo
from Bio.Phylo.Consensus import *
from io import StringIO
import numpy as np
import numpy
import glob
from scipy import stats
import pandas as pd
from scipy.stats.distributions import chi2
from scipy.stats import fisher_exact
def ReadFasSeq(pos72):
 StLs=[]
 St2Seq={}
 pos72=open(pos72,'r').readlines()
 for i in pos72:
    if i[0]=='>':
       ID=i.strip()
       StLs.append(ID)
       St2Seq[ID]=''
    else: St2Seq[ID]+=i.strip().replace('a','A').replace('t','T').replace('g','G').replace('c','C')
 return StLs,St2Seq
def ReadMegSeq(Meg): #input is mega alignment file. out is name2seq dictionary and mega head
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
    elif Read=='n': Name2Seq[Name]+=i.strip().replace('a','A').replace('t','T').replace('g','G').replace('c','C')
    elif Read=='s': out2+=i
  return NameOrder, Name2Seq  
def Fas2Meg(FasFile,MegOut):
    SeqLs,Seq2Seq=ReadFasSeq(FasFile)
    out=['#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n']
    for Seq in SeqLs:
        out.append('#'+Seq.replace('>','')+'\n'+Seq2Seq[Seq]+'\n')
    GetOut(MegOut,''.join(out))  

def Nuc2BEAMinWnNoFil(Fas,BEAMin):
    CellLs,Cell2Seq=ReadFasSeq(Fas)
   # IDta=pd.read_csv(Ta,sep='\t')
   # print (IDta)
   # RefSeq=IDta['REF'].to_list()
    Len=len(Cell2Seq[CellLs[0]])
    RefSeq=Cell2Seq['>normal']
  #  if len(Cell2Seq[CellLs[0]])!=len(RefSeq):
  #      print (len(Cell2Seq[CellLs[0]]),len(RefSeq))
  #      open('a','r').readlines()
  #  KeepPosLs=[]
  #  Wild=IDta['Wild'].to_list()  
  #  Mut=IDta['MutTotal'].to_list()  
  #  c=0
  #  KeepPos=[]
  #  while c<Len:
  #     if Wild[c]>=5 and Mut[c]>=5:
  #         KeepPos.append(c)
  #     c+=1
  #  print (KeepPos[:5],len(KeepPos))   
    out=['#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n'] 
  #  outFas=[]
    for Cell in CellLs:
       Seq=Cell2Seq[Cell]
       if len(Seq)!=Len:
           print ('seq length different',len(Seq),Len,Cell,CellLs[0])
           open('a','r').readlines()
       c=0
       NewSeq=''
      # NewATGC=''
       while c < Len:
          if Seq[c]=='?': NewSeq+='?'
          elif Seq[c]!=RefSeq[c]: NewSeq+='T'
          else: NewSeq+='A'
          c+=1
        #  NewATGC+=Seq[c]
        #  c+=1
       if NewSeq.find('T')!=-1: 
           out.append('#'+Cell.replace('>','')+'\n'+NewSeq+'\n')
       #outFas.append(Cell+'\n'+NewATGC+'\n')    
    GetOut(BEAMin,''.join(out))  
   # GetOut(BEAMin[:-4]+'.fasta',''.join(outFas))  
   # GetOut(BEAMin[:-4]+'_Pos.txt','\n'.join(map(str,KeepPos))) 
def Nuc2BEAMinWn(Fas,Ta,BEAMin):
    CellLs,Cell2Seq=ReadFasSeq(Fas)
    IDta=pd.read_csv(Ta,sep='\t')
    print (IDta)
    RefSeq=IDta['REF'].to_list()
    Len=len(RefSeq)
    if len(Cell2Seq[CellLs[0]])!=len(RefSeq):
        print (len(Cell2Seq[CellLs[0]]),len(RefSeq))
        open('a','r').readlines()
    KeepPosLs=[]
    Wild=IDta['Wild'].to_list()  
    Mut=IDta['MutTotal'].to_list()  
    c=0
    KeepPos=[]
    while c<Len:
       if Wild[c]>=5 and Mut[c]>=5:
           KeepPos.append(c)
       c+=1
    print (KeepPos[:5],len(KeepPos))   
    out=['#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n'] 
    outFas=[]
    for Cell in CellLs:
       Seq=Cell2Seq[Cell]
       #c=0
       NewSeq=''
       NewATGC=''
       for c in KeepPos:
          if Seq[c]=='?': NewSeq+='?'
          elif Seq[c]!=RefSeq[c]: NewSeq+='T'
          else: NewSeq+='A'
          NewATGC+=Seq[c]
        #  c+=1
       if NewSeq.find('T')!=-1: 
           out.append('#'+Cell.replace('>','')+'\n'+NewSeq+'\n')
       outFas.append(Cell+'\n'+NewATGC+'\n')    
    GetOut(BEAMin,''.join(out))  
    GetOut(BEAMin[:-4]+'.fasta',''.join(outFas))  
    GetOut(BEAMin[:-4]+'_Pos.txt','\n'.join(map(str,KeepPos))) 
      

      
def Nuc2BEAMin(Fas,Ta,BEAMin):
    CellLs,Cell2Seq=ReadFasSeq(Fas)
    IDta=pd.read_csv(Ta,sep='\t')
    print (IDta)
    RefSeq=IDta['REF'].to_list()
    Len=len(RefSeq)
    if len(Cell2Seq[CellLs[0]])!=len(RefSeq):
        print (len(Cell2Seq[CellLs[0]]),len(RefSeq))
        open('a','r').readlines()
    KeepPosLs=[]
    Wild=IDta['Wild'].to_list()  
    Mut=IDta['Mut'].to_list()  
    c=0
    KeepPos=[]
    while c<Len:
       if Wild[c]>=10 and Mut[c]>=10:
           KeepPos.append(c)
       c+=1
    print (KeepPos[:5],len(KeepPos))   
    out=['#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n'] 
    for Cell in CellLs:
       Seq=Cell2Seq[Cell]
       #c=0
       NewSeq=''
       for c in KeepPos:
          if Seq[c]=='?': NewSeq+='?'
          elif Seq[c]!=RefSeq[c]: NewSeq+='T'
          else: NewSeq+='A'
        #  c+=1
       out.append('#'+Cell.replace('>','')+'\n'+NewSeq+'\n')
    GetOut(BEAMin,''.join(out))  
    GetOut(BEAMin[:-4]+'_Pos.txt','\n'.join(map(str,KeepPos))) 
        
def ReadGV(GV):
    Dec2Anc={}
    NodeMut={}
    Node2In={}
    Edge2In={}	
    GV=open(GV,'r').readlines()#[1:]
    Read='N'
    for i in GV:
        i=i.replace('digraph D {','')	
        if i.find('->')!=-1:
             i0=i.replace(';','').split('->')
             Dec2Anc[i0[1].split('[')[0].strip()]=i0[0].strip()
             Edge2In[i0[0].strip()+'->'+i0[1].split('[')[0].strip()]=i.strip()	
        elif i[0]=='}': Read='D'			
        elif i.strip()!='' and i.find('[')!=-1:
             i0=i.split(' [')
            # print (i0)			 
             Node=i0[0].strip()
             if i0[1].find('\\nMut:')!=-1:			 
                 MutLs=i0[1].split('\\nMut:')[-1].replace('\\n\"]\n','').split(';')
             else: MutLs=[]				 
             MutIn=[]
             for M in MutLs:
                 M=M.strip().split('\\n')			   
                 MutIn+=M	
             NodeMut[Node]=MutIn	
             Node2In[Node]=i		 

	 
        else: pass		
	
    return Dec2Anc,NodeMut,Node2In,Edge2In	
def ReadGV1(GV):
    Dec2Anc={}
    NodeMut={}
    Node2In={}
    Edge2In={}	
    GV=open(GV,'r').readlines()#[1:]
    Read='N'
    for i in GV:
        i=i.replace('digraph D {','')	
        if i.find('->')!=-1:
             i0=i.replace(';','').split('->')
             Dec2Anc[i0[1].split('[')[0].strip()]=i0[0].strip()
             Edge2In[i0[0].strip()+'->'+i0[1].split('[')[0].strip()]=i.strip()	
        elif i[0]=='}': Read='D'			
        elif i.strip()!='' and i.find('[')!=-1:
             i0=i.split(' [')			 
             Node=i0[0].strip()
             MutIn0=i0[1].replace('\"]\n','').split('\\n')[1:]
             MutIn=[]
             for M in MutIn0:
               if M.strip()!='':
                 if M[0]=='B': M='T'+M.replace('Back','').replace(' ','')+'A'
                 elif M[0]=='R': M='Rec A'+M.replace('Rec','').replace(' ','')+'T'
                 else: M='A'+M+'T'
                 MutIn.append(M)
             NodeMut[Node]=MutIn	
             Node2In[Node]=i		 

	 
        else: pass		
	
    return Dec2Anc,NodeMut,Node2In,Edge2In	  
def ReadGVmoa(GV):
    Dec2Anc={}
    NodeMut={}
    Node2In={}
    Edge2In={}	

    GV=open(GV,'r').readlines()#[1:]
    Read='N'
    for i in GV:
        i=i.replace('digraph G {','')	
        if i.find('->')!=-1:
           if i.find('color=\"#ff00005f\"')==-1:
             i0=i.replace(';','').split('->')
             Dec2Anc[i0[1].split('[')[0].strip()]=i0[0].strip()
             Edge2In[i0[0].strip()+'->'+i0[1].split('[')[0].strip()]=i.strip()	
        elif i[0]=='}': Read='D'			
        elif i.strip()!='' and i.find('[')!=-1:
             i0=i.split(' [')			 
             Node=i0[0].strip()
             MutIn=[i0[1].split('\\n')[0].split('label=')[1].replace('\"','').split(',')[0].strip()]
             NodeMut[Node]=MutIn	
             Node2In[Node]=i		 
        else: pass		
	
    return Dec2Anc,NodeMut,Node2In,Edge2In
def ReadNwkWithNodeID(Tree):
    tree = Phylo.read(Tree, "newick")
    C=1
    Name2NodeID={}
    for i in tree.find_clades():
       if i.name==None:
            i.name=C
       Name2NodeID[i.name]=i.comment     
       C+=1
    Tips=tree.get_terminals()
    Dec2Anc={}
    for Tip in Tips:   
       Root2TipLs=tree.get_path(Tip.name)
       D=1
       Len=len(Root2TipLs)
       if Len==1:
           Dec2Anc[Root2TipLs[0].name]='root'       
       while D<Len:
           Dec2Anc[Root2TipLs[D].name]=Root2TipLs[D-1].name
           D+=1
    Dec2Anc[Root2TipLs[0].name]='root'            
    return Dec2Anc,Name2NodeID,Tips	    	
def ReadNodeMap(NodeMap):
    NodeMap=open(NodeMap,'r').readlines()[1:]
    D2A={}
    A2D={}
    N2C={}
    C2N={}
    for i in NodeMap:
        Ls=i.strip().split('\t')
        Line=[]
        for Item in Ls:
            Item=Item.strip().replace(' ','_')
            if Item!='':Line.append(Item)
        N=Line[0]
        C=Line[1]
        N2C[N]=C
        C2N[C]=N
        if Line[3]!='-':
                A2D[C]=[Line[2],Line[3]]
                D2A[Line[2]]=C
                D2A[Line[3]]=C
    return D2A,A2D,N2C,C2N
def ReadCOI(COI):
    COI=open(COI,'r').readlines()
    ParOr=COI[0].strip().split('\t')[1:]
    Len=len(ParOr)
    COI=COI[1:]
    Dec2Anc2COI={}
    for i in COI:
        i=i.split('\t')
        Chi=i[0]
        Dec2Anc2COI[Chi]={}	
        COIs=i[1:]
        c=0
        while c<Len:
           Dec2Anc2COI[Chi][ParOr[c]]=float(COIs[c])
           c+=1
    return Dec2Anc2COI	
def ReadFreTa(Ta):
    Ta=open(Ta,'r').readlines()[1:]
    Pos2Fre={}
    for i in Ta:
        i=i.split('\t')
        Pos2Fre[int(i[0])]=float(i[4])
    return Pos2Fre
def ListCol(File):
  File=open(File,'r').readlines()
  NameOrder,Name2Col=GetHead(File[0])
  File=File[1:]
  Tu2Freq={}
  for Tu in NameOrder:
    Tu2Freq[Tu]=[]
  for i in File:
    i=i.strip().split('\t')
    for Tu in Name2Col:
        Tu2Freq[Tu].append(i[Name2Col[Tu]])
  return Tu2Freq
def GetHead(Head):
    Head=Head.strip().split('\t')
    Len=len(Head)
    c=0
    Name2Col={}
    NameOrder=[]	
    while c<Len:
        Name2Col[Head[c]]=c
        NameOrder.append(Head[c])		
        c+=1
    return NameOrder,Name2Col	
def MOA2GV(MOA):
    GV=MOA[:-4]+'.gv'
    Dec2Anc,NodeInf,Node2In,Edge2In=ReadGV(MOA)
    ID2Lab={}
    for N in Node2In:
       Lab=Node2In[N].split('label=')[-1].split(',')[0].split('\\n')[0]
       ID2Lab[N]=Lab.replace('\"','')	
    out='digraph G {\n'
    for N in Node2In:
       In=Node2In[N]	
       out+=ID2Lab[N]+In[len(N):]
    for E in Edge2In:
      if Edge2In[E].find('label=')!=-1 or E.find('root')!=-1:	
       E=E.split('->')
       out+=ID2Lab[E[0].strip()]+'->'+ID2Lab[E[1].replace(';','').strip()]+';\n'
    out+='}'
    GetOut(GV,out)	
def GetVAFpos(seq_annoFile,VAF):	#Hap+os.sep+'sequence_annotations_Force.txt'
    Pos2Cou={}
    File=open(seq_annoFile,'r').readlines()[1:]
    for i in File:
        i=i.split('\t')
        P=i[0]
        if P not in Pos2Cou: Pos2Cou[P]={'Tot':0,'Var':{}}
        Ref=i[1]
        Alt=i[2]
        Cou=int(i[3])		
        if Ref!=Alt and Ref!='?' and Alt!='?' and Ref!='-' and Alt!='-':
             Pos2Cou[P]['Var'][Alt]=Cou 
        Pos2Cou[P]['Tot']+=Cou
    Pls=[]		
    for P in Pos2Cou:
        Var2Cou=Pos2Cou[P]['Var']
        Tot=Pos2Cou[P]['Tot']		
        for	V in Var2Cou:
            Pro=1.0*Var2Cou[V]/Tot	
            if Pro>VAF: 
                Pls.append(int(P))			
    Pls=list(set(Pls))
    Pls.sort()	
    print (len(Pls))	
    GetOut(seq_annoFile[:-4]+'_Force.txt','\n'.join(map(str,(Pls))))	
		
def GetMostSim(Hap,RefHap):
    DifC2HapLs={}
    Len=len(Hap)
    for Ref in RefHap:
        RefSeq=RefHap[Ref]
        if len(RefSeq)!=Len:
            print (Hap,RefSeq,Len,len(RefSeq))
            open('a','r').readlines()
        else:
           DC=CountDifNum(Hap,RefSeq)
           DifC2HapLs[DC]=DifC2HapLs.get(DC,[])+[Ref.replace('>','')]
    DifCLs=list(DifC2HapLs.keys())
    DifCLs.sort()
    DifC=DifCLs[0]
    Gid=';'.join(DifC2HapLs[DifC])

    return Gid,DifC
def CountDifNum_excMiss(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!='?' and Seq1[c]!='?' and Seq0[c]!='-' and Seq1[c]!='-':			
                    if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif		
def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif
def CountDifPos_from1(Seq0,Seq1):
            Len=len(Seq0)
            Dif=[]
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif.append(c+1)
                c+=1
            return Dif
def makeRaxMLin(OriDic,Qseq,OutFas,Add):
    for Q in Qseq:
        InFas='>Query_'+Q.replace('>','')+'\n'+Qseq[Q]+Add+'\n'
    for O in OriDic:
      InFas+=O.replace('__','_')+'\n'+OriDic[O]+Add+'\n'
    GetOut(OutFas,InFas)

def cleanRaxMLtree(Tree,Out,OutG):
    Len=len(Tree)
    c=0
    Read='s'
    New=''
    while c<Len:
       if Tree[c]=='[':
           Read='r'
       elif Read=='r':
          if Tree[c]==']':
           Read='s'
       elif Read=='s': New+=Tree[c]
       c+=1
    GetOut('Unroot.nwk',New)
    root_tree('Unroot.nwk',OutG)
    GetOut(Out,open('Unroot_rooted.nwk','r').readlines()[0])
    os.remove('Unroot.nwk')
def prune_tree(OriNwk,ExtraLs,OutSeq):
   if os.path.exists('rooted.nwk')==True: os.remove('rooted.nwk')
   trees = list(Phylo.parse(OriNwk, 'newick'))
   for tree in trees:
       tree = tree.root_with_outgroup({'name': OutSeq})
   Phylo.write(trees, 'rooted.nwk', "newick")
   TreeLs=open('rooted.nwk','r').readlines()
   Cou=1
   for Tst in TreeLs:
    tree=Phylo.read(StringIO(Tst.strip()), "newick")
    for tip in ExtraLs:
         tree.prune(tip)
    Phylo.write(tree, OriNwk[:-4]+'_'+str(Cou)+'_prune.nwk','newick')
    Cou+=1

def root_tree(OriNwk,Root):
    Out=OriNwk[:-4]+'_rooted.nwk'
    trees = list(Phylo.parse(OriNwk, 'newick'))
    for tree in trees:
       tree = tree.root_with_outgroup({'name': Root})
    Phylo.write(trees, Out, "newick")
def Getalldec(N,Dec2Anc):
    Anc2Dec=GetAnc2Dec(Dec2Anc)
    DecLs=Anc2Dec.get(N,[])
    DecLsAll=Anc2Dec.get(N,[])	
    Added='y'
    while DecLs!=[]:
        NewDecLs=[]
        for D in DecLs:
            NewDecLs+=Anc2Dec.get(D,[])
 
        DecLsAll+=NewDecLs
        DecLs=NewDecLs
    return DecLsAll	
def GetAnc2Dec(Dec2Anc):
    Anc2Dec={}
    for D in Dec2Anc:
        A=Dec2Anc[D]	 
        Anc2Dec[A]=Anc2Dec.get(A,[])+[D]	
    return Anc2Dec	
def GetTipLs(Dec2Anc):
    Anc2Dec=InvertDic1(Dec2Anc)
    TipLs=[]
    for Dec in Dec2Anc:
        if Dec not in Anc2Dec: TipLs.append(Dec)
    return TipLs        
def ClassifyDif(AncSeq,DecSeq,PosIDOrder):
    ForLs=[]
    BackLs=[]
    Len=len(AncSeq)
    c=0
    while c<Len:
        if AncSeq[c]=='A' and DecSeq[c]=='T': ForLs.append(PosIDOrder[c].strip())
        elif AncSeq[c]=='T' and DecSeq[c]=='A': BackLs.append(PosIDOrder[c].strip())
        c+=1
    return ForLs,BackLs 
def ClassifyDif1(AncSeq,DecSeq,PosIDOrder):
    ForLs=[]
    BackLs=[]
    Len=len(AncSeq)
    print (Len,len(DecSeq))
    c=0
    while c<Len:
        if AncSeq[c] != DecSeq[c]: ForLs.append(AncSeq[c]+PosIDOrder[c].strip()+DecSeq[c])	
        c+=1
    return ForLs,BackLs 	
def MutTree2dot(File):
   Out=File[:-4]+'.gv'

   File=open(File,'r').readlines()[1:]
   Nodes=''
   Paths=''
   for i in File:
       i=i.split('\t')
       
       if len(i)>5: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4]+'\\nCellC:'+i[5].strip()+'\"]\n'	
       else: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4].strip()+'\"]\n'	   
       Paths+=i[1]+'->'+i[0]+'\n'	   

   out='digraph D {\n'+Nodes+Paths+'}\n'
   GetOut(Out,out)
def MutTree2dot1(File):
   Out=File[:-4]+'.gv'

   File=open(File,'r').readlines()[1:]
   Nodes=''
   Paths=''
   for i in File:
       i=i.replace('.','').split('\t')
       
       #if len(i)>5: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4]+'\\nCellC:'+i[5].strip()+'\"]\n'	
       Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'
       MutLs=i[2].split(';')
       In=[]
       c=0
       for Mut in MutLs:
           if c>10: 
               Nodes+=';'.join(In)+'\\n'
               In=[]
               c=0
           
           In.append(Mut)
           c+=1		   
       if In!=[]:    Nodes+=';'.join(In)+'\\n'               		   
       Nodes+='\"]\n'	   
       Paths+=i[1]+'->'+i[0]+'\n'	   

   out='digraph D {\n'+Nodes+Paths+'}\n'
   GetOut(Out,out)
def GetComAnc(NodeLs,Dec2Anc):
    Node2C={}
    for N in NodeLs:
        Node2C[N]=Node2C.get(N,0)+1
        while N in Dec2Anc:
            N=Dec2Anc[N]
            Node2C[N]=Node2C.get(N,0)+1
    Com=[]
    for N in Node2C:			
        if Node2C[N]>1: Com.append(N)
    return Com
def getCOI(TarM,AncMLs,Dec2Anc2COI):
    AncCOI=[]
    for A in AncMLs:
      if A!='':	
       if TarM.replace('Rec','')[1:]!='' and A.replace('Rec','')[1:]!='':	
        if TarM.replace('Rec','')[1:] in Dec2Anc2COI:
           COI0=Dec2Anc2COI[TarM.replace('Rec','')[1:]]
           if A.replace('Rec','')[1:] in COI0:           
               COI=Dec2Anc2COI[TarM.replace('Rec','')[1:]][A.replace('Rec','')[1:]]
               AncCOI.append(A+'='+str(COI))
    return AncCOI	
def PruneTip(GV,Fas,Nwk,Anno):
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    KeepTH=[]
    PruneTH=[]
    PruneTHNwk=[]
    out='digraph D {\n'
    for TH in NodeMut:
      if TH!='root':
        MutLs=NodeMut[TH]
        if MutLs!=[''] and len(MutLs)>0: 
            KeepTH.append(TH)
            out+=Node2In[TH]
        else: 
           if TH[:5]!='Node_': PruneTHNwk.append(TH)
           PruneTH.append(TH)       
      else: out+=Node2In[TH]
    for E in Edge2In:
      if PruneTH.count(E.split('->')[1].strip())==0:
        out+=Edge2In[E]+'\n'
    out+='}\n'    
    GetOut(GV[:-3]+'1.gv',''.join(open(GV,'r').readlines()))    
    KeepTH.append('Outgroup')
    print (PruneTH)
    print ('prune fasta')
    THls,TH2Seq=ReadFasSeq(Fas)
    out=''
    for TH in THls:
        if PruneTH.count(TH.replace('>',''))==0:
            out+=TH+'\n'+TH2Seq[TH]+'\n'
        else: print ('prune',TH)
    GetOut(Fas[:-6]+'COI.fasta',out)
    print ('prune tree')
    prune_tree(Nwk,PruneTHNwk,'Outgroup')
    print ('prune annotation')    
    Anno=open(Anno,'r').readlines()
    out=Anno[0]
    Anno=Anno[1:]
    for i in Anno:
        TH=i.split('\t')[2].strip()
        if PruneTH.count(TH)==0: out+=i        
    GetOut(Nwk.replace(Nwk.split(os.sep)[-1],'')+'TopHapCOI_anno.txt',out)#'TopHapNA__2022_1_FillAllCell_anno1_COI.txt')    
def CleanMutTree(GV,MutIDLs,Nwk):
    OutDir=GV.replace(GV.split(os.sep)[-1],'')
    Out=OutDir+'MutTree.gv'
    print (Out)
    COI=OutDir+'MOAout'+os.sep+'COI_matrix.txt'
    COIout=COI[:-11]+'.txt'
    CloTa=OutDir+'TopHapCOI_anno.txt'#'TopHapNA__2022_1_FillAllCell_anno1_COI.txt'
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    
    Dec2Anc2COI=ReadCOI(COI)
    
    if os.path.exists(MutIDLs)==True:
        MutIDOrder=open(MutIDLs,'r').readlines()
    else: 
        MutLs0=Dec2Anc2COI.keys()
        print (MutLs0)
        #open('a','r').readlines()
        MutLs=[]
        for i in MutLs0:
    
            MutLs.append(int(i[:-1]))
        MutLs.sort()
        MutIDOrder= list(range(0,MutLs[-1]+1))
        print (MutIDOrder)
    COI=open(COI,'r').readlines()
    Head=COI[0].strip().split('\t')
    out=Head[0]
    Head=Head[1:]
    for i in Head:
        Pos=int(i[:-1])
        out+='\t'+ str(MutIDOrder[Pos]).strip()+'_'+i 
    out+='\n'
    COI=COI[1:]
    for i in COI:
        Pos=int(i.split('\t')[0][:-1])
        out+=str(MutIDOrder[Pos]).strip()+'_'+i 
    GetOut(COIout,out)    
    KeepTipLs=[]
    CloTa=open(CloTa,'r').readlines()[1:]
    for i in CloTa:
        i=i.split('\t')[2]
        KeepTipLs.append(i.strip())
    KeepTipLs=list(set(KeepTipLs))
    print ('clone list',KeepTipLs)    
    Anc2Dec=InvertDic1(Dec2Anc)
    TipLs=GetTipLs(Dec2Anc)
    
    PruneLs=[]
    Mut2C={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
            Mut2C[M]=Mut2C.get(M,0)+1
    for Dec in Dec2Anc:
        if Dec[:3]=='All' and KeepTipLs.count(Dec)==0: PruneLs.append(Dec)
        
    for Dec in Dec2Anc:
        if Dec[:5]=='Node_':
            Lins=Anc2Dec[Dec]
            GoodLn=0
            for Lin in Lins:
                AllDecLs=Getalldec(Lin,Dec2Anc)
                AllDecLs.append(Lin)
                Good='n'
                for D in AllDecLs:
                     if D[:3]=='All' and PruneLs.count(D)==0: Good='y'
                if Good=='y': GoodLn+=1   
            if GoodLn<2: PruneLs.append(Dec)
    print ('prune node',PruneLs)   
    NewNode2In={}
    NewEdge2In={}
    PruneLs=[]
    for Node in Node2In:
       In=Node2In[Node]
       if In.find('Mut:')==-1:NewNode2In[Node]=In
       else:    
        MutLs=In.split('Mut:')[-1].replace('\"]','').split(';')
        NewMutLs=[]
        for M0 in MutLs:
          M0=M0.strip().split('\\n')
          for M in M0:
           if M.strip()!='' and M.strip().replace('Rec','')!='':
            Pos=M.strip().replace('Rec','').strip()[1:][:-1]
            print ( Pos,M)
            MID=str(MutIDOrder[int(Pos)]).strip()
            GoodBack='y'
            if Mut2C[M]>1: MID='Rec '+MID #M.find('Rec')!=-1: MID='Rec '+MID
            if M.strip()[-1]=='A':
                 ForFind='n'
                 Anc=Dec2Anc[Node]
                 AncMutLs=NodeMut[Anc]
                 if AncMutLs.count('A'+Pos+'T')!=0: ForFind='y'
                 while Anc in Dec2Anc:
                     Anc=Dec2Anc[Anc]
                     AncMutLs=NodeMut[Anc]
                     if AncMutLs.count('A'+Pos+'T')!=0: ForFind='y'                       
                 MID='Back '+MID
                 if ForFind=='n': GoodBack='n'
            if GoodBack=='y': NewMutLs.append(MID)#.replace(Pos,':'+MID+':').strip())
            else: 
                print ('remove due to no forward mutation',MID)
        if (NewMutLs==[] or NewMutLs==['']) and Node[:3]=='All': PruneLs.append(Node)         
        NewIn=In.split('Mut:')[0]+'\\n'.join(NewMutLs)+'\"]\n'    
        NewNode2In[Node]=NewIn
    for E in Edge2In:
        In=Edge2In[E]
        COIls=In.split('label=\"')[-1].replace('\"];','').split('\\n')
        NewEdge2In[E]=E+';\n'#NewIn    
    out='digraph D {\n'
    for N in NewNode2In:
        out+=NewNode2In[N]
     
    for E in NewEdge2In:
        out+=NewEdge2In[E]+'\n'
    out+='}'
    GetOut(Out,out)	
    os.system('dot -Tpng '+Out+' -o '+Out[:-3]+'.png')   
    print ('pruen list',PruneLs)    
    prune_tree(Nwk,PruneLs,'Outgroup')     
def CleanMutTree1(GV,MutIDLs,CloTa,SNVc): #CleanMutTree1(TopHapMOAgv,MutIDFile,TopHapMOACellAnno) #output '_final.gv' 
    OutDir=GV.replace(GV.split(os.sep)[-1],'')
    Out=GV[:-3]+'_final.gv'#OutDir+'MutTree.gv'
    print (Out)
    COI=OutDir+'MOAout'+os.sep+'COI_matrix.txt'
    COIout=COI[:-11]+'.txt'
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)  
    Dec2Anc2COI=ReadCOI(COI)   
    if os.path.exists(MutIDLs)==True:
        MutIDOrder=open(MutIDLs,'r').readlines()
    else: 
        MutIDOrder= list(range(0,SNVc))
        print (MutIDOrder)
    COI=open(COI,'r').readlines()
    Head=COI[0].strip().split('\t')
    out=Head[0]
    Head=Head[1:]
    for i in Head:
        Pos=int(i[:-1])
        out+='\t'+ str(MutIDOrder[Pos]).strip()+'_'+i 
    out+='\n'
    COI=COI[1:]
    for i in COI:
        Pos=int(i.split('\t')[0][:-1])
        out+=str(MutIDOrder[Pos]).strip()+'_'+i 
    GetOut(COIout,out)    
    KeepTipLs=[]
    CloTa=open(CloTa,'r').readlines()[1:]
    for i in CloTa:
        i=i.split('\t')[1]
        KeepTipLs.append(i.strip())
    KeepTipLs=list(set(KeepTipLs))
    print ('clone list',KeepTipLs)    
    Anc2Dec=InvertDic1(Dec2Anc)
    TipLs=GetTipLs(Dec2Anc)
    
    PruneLs=[]
    Mut2C={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
            Mut2C[M]=Mut2C.get(M,0)+1
    print ('mutID to mut count',Mut2C)        
    for Dec in Dec2Anc:
        if KeepTipLs.count(Dec)==0: PruneLs.append(Dec)
    print ('prune node',PruneLs)   
    NewEdgeIn=[]
    NewNode2Mut={}
    for TN in TipLs:
       if TN not in PruneLs:
           TN1=TN
           NodeID=TN            
           Muts=RenameMut(NodeMut[TN1],Mut2C,MutIDOrder)
           
           while TN1 in Dec2Anc:
              TN1=Dec2Anc[TN1]
              DecC=len(Anc2Dec[TN1])
              if (TN1 in PruneLs or TN1[:3]!='All') and DecC==1:
                  Muts+=RenameMut(NodeMut[TN1],Mut2C,MutIDOrder)
              else: 
                  NewNode2Mut[NodeID]=Muts
                  NewEdgeIn.append(TN1+'->'+NodeID)
                  NodeID=TN1 
                  Muts=RenameMut(NodeMut[TN1],Mut2C,MutIDOrder)
           NewNode2Mut[NodeID]=Muts
    
    out=['digraph D {\n']
    for N in NewNode2Mut:
        if N[:3]!='All': Lab='Node_'+N
        else: Lab=N
        out+=[N+' [label=\"'+Lab+'\\n'+'\\n'.join(NewNode2Mut[N])+'\"]\n']
    NewEdgeIn=list(set(NewEdgeIn)) 
    for E in NewEdgeIn:
        out+=[E+';\n']
    out+=['}']
    GetOut(Out,''.join(out))	 
def RenameMut(MutLs,Mut2C,IDOr):
    NewLs=[]
    for M in MutLs:
        if M[0]=='B' or M[-1]=='A':Back='Back '
        else: Back=''
        if M.find('root')==-1:
           #print (M)
           Mid=int(M.replace('Rec','').replace('Back','').replace('A','').replace('T',''))
           if Mut2C[M]>1: NewLs.append('Rec '+Back+str(IDOr[Mid]).strip())
           else: NewLs.append(Back+str(IDOr[Mid]).strip())
    return NewLs    
def GV2Fas(PosOr0):
    MainFol=PosOr0.replace(PosOr0.split(os.sep)[-1],'')
    COImatrix=MainFol+'TopHapOut'+os.sep+'MOAout'+os.sep+'COI_matrix.txt'
    OutGseq=MainFol+'OutG.fasta' #haplotype
    Nwk=MainFol+'TopHapOut'+os.sep+'TopHap_allMP_AncMP.nwk'
    GV=MainFol+'TopHapOut'+os.sep+'TopHap_MutTree_prune_COI3_ave3_rec_back_COI3_ave3.gv'
    OutF=MainFol+'TopHapOut'+os.sep+'TopHapCOI.fasta'
    CellAnnoTar=MainFol+'TopHapOut'+os.sep+'TopHap*_FillAllCell_anno1.txt'
    PosOr=open(PosOr0,'r').readlines()
    Len=len(PosOr)
    OutGLs,OutG2Seq=ReadFasSeq(OutGseq)
    OutSeq=OutG2Seq[OutGLs[0]]
    print (OutSeq)
    if len(OutSeq)!=Len: open('a','r').readlines()
    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    print (GV, NodeMut)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    for N in Dec2Anc:
        if N.find('Node_')==-1:
            TH='>'+N	
            MLs=NodeMut[N]
            Mut=0
            ForMut=0		
            for M in MLs:
               if M.strip()!='': Mut+=1
               if M.strip()!='': 
                   if M[0]=='T' or M.find('Rec')!=-1: ForMut+=1		   
            if Mut==0: RmTH.append(N)
    		
            else:
               AllMut=MLs
               while N in Dec2Anc:
                   N=Dec2Anc[N]
                   Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
                   if N in NodeMut: 
                       AllMut+=NodeMut[N]
                       print (N,NodeMut[N])
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	
               print (TH,AllMut)               
               for M in AllMut:
                 if M.strip()!='':		   
                   Pos=M.replace('Rec','')[1:][:-1]
                    			   
                   if PosOr.count(Pos)!=0: P=PosOr.index(Pos)
                   else: P=PosOr.index(Pos+'\n')
                   NewSeq[P]=M[-1]
                   if M[0]=='A': NewSeqRB[P]=M[-1]			   
               NewSeq=''.join(NewSeq)	
               NewSeqRB=''.join(NewSeqRB)		   
               print (NewSeq)
               out+=TH+'\n'+NewSeq+'\n'	   
               if ForMut==Mut: 
    		   
                    BackRecTH.append(TH.replace('>',''))	
               else: outRB+=TH+'\n'+NewSeqRB+'\n'	   
    out+=OutGLs[0]+'\n'+OutSeq+'\n'	 
    outRB +=OutGLs[0]+'\n'+OutSeq+'\n'	 
    GetOut(OutF,out)	
    GetOut(OutF[:-6]+'RB.fasta',outRB)
    print ('prune TopHap haplotype',RmTH)
    if os.path.exists(Nwk)==True:
        prune_tree(Nwk,RmTH,OutGLs[0].replace('>',''))
        if os.path.exists(Nwk[:-17]+'COI.nwk')!=True: os.system('megacc -a analyze_user_tree_MP__nucleotide.mao -d '+OutF+' -t '+Nwk[:-4]+'_1_prune.nwk -o '+Nwk[:-17]+'COI.nwk')
        TopHap_DrawMutTree2(OutF,PosOr0,Nwk[:-17]+'COI.nwk','Outgroup')        
        MapCOI3(COImatrix,OutF[:-6]+'_MutTree.gv')        	
        AveCOI3(OutF[:-6]+'_MutTree_COI3.gv')
    CellAnnoLs=glob.glob(CellAnnoTar)           		   
    print (CellAnnoLs)
    for CellAnno in CellAnnoLs:
        Out=CellAnno[:-4]+'_COI.txt'
        CellAnno=open(CellAnno,'r').readlines()[1:]
        out='Group\tSeq\tTopHapID\n'
        for i in CellAnno:
            if RmTH.count(i.split('\t')[2].strip())==0: out+=i	
    
        GetOut(Out,out)	
    RmTH+=BackRecTH
    RmTH=list(set(RmTH))
def GV2Fas1(GV,PosOr0):
    OutF=GV[:-3]+'.fasta'
    OutP=GV[:-3]+'_Pos.txt'
    if os.path.exists(PosOr0)==True:
        PosOr=open(PosOr0,'r').readlines()
    else: 
        PosOr=range(0,99999)
        PosOr=[str(x) for x in PosOr]
    Len=len(PosOr)
    OutSeq='A'*Len    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    print (GV, NodeMut)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    TH2Seq={}
    for N in Dec2Anc:
        
        if N.find('Node_')==-1:
            TH='>'+N	
            MLs=NodeMut[N]
            Mut=0
            ForMut=0		
            for M in MLs:
               if M.strip()!='': Mut+=1
               if M.strip()!='': 
                   if M.find('Back')!=-1 or M.find('Rec')!=-1: ForMut+=1		   
            if Mut==0: RmTH.append(N)
    		
            else:
               AllMut=MLs
               while N in Dec2Anc:
                   N=Dec2Anc[N]
                   Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
                   if N in NodeMut: 
                       AllMut+=NodeMut[N]
                       print (N,NodeMut[N])	
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	              
               for M in AllMut:
                 if M.strip()!='':		   
                   Pos=M.replace('Rec','').replace('Back','').strip().replace('A','').replace('T','').replace('?','')
                    			   
                   if PosOr.count(Pos)!=0: P=PosOr.index(Pos)
                   else: P=PosOr.index(Pos+'\n')
                   NewSeq[P]='T'
                   if M.find('Back')!=-1: NewSeq[P]='A'			   
               NewSeq=''.join(NewSeq)			   
               out+=TH+'\n'+NewSeq+'\n'
               TH2Seq[TH]=NewSeq	   
    print ('extract variable sites')
    MutP ,ExtSeqDic =	ExtractVarPos(TH2Seq) 
    out='>Outgroup\n'+('A'*len(MutP))+'\n'	
    for TH in ExtSeqDic:
        out+=TH+'\n'+ExtSeqDic[TH]+'\n'
    GetOut(OutF,out)
    out=''
    for M in MutP:
        out+= PosOr[M]+'\n'
    GetOut(OutP,out)
  
    return MutP,ExtSeqDic
def GV2Fas2(GV,PosOr0):
    OutF=GV[:-3]+'.fasta'
    OutP=GV[:-3]+'_Pos.txt'
    if os.path.exists(PosOr0)==True:
        PosOr=open(PosOr0,'r').readlines()
    else: 
        PosOr=range(0,99999)
        PosOr=[str(x) for x in PosOr]
    Len=len(PosOr)
    OutSeq='A'*Len    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    print (GV, NodeMut)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    TH2Seq={}
    Seq2THls={}
    for N in Dec2Anc:
        
        if N.find('Node_')==-1:
            TH='>'+N	
            MLs=NodeMut[N]
            Mut=0
            ForMut=0		
            for M in MLs:
               if M.strip()!='': Mut+=1
               if M.strip()!='': 
                   if M.find('Back')!=-1 or M.find('Rec')!=-1: ForMut+=1		   
            Go='y' 
            if Go=='y':            
               AllMut=MLs
               while N in Dec2Anc:
                   N=Dec2Anc[N]
                   Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
                   if N in NodeMut: 
                       AllMut+=NodeMut[N]
                       print (N,NodeMut[N])
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	             
               for M in AllMut:
                 if M.strip()!='':		   
                   Pos=M.replace('Rec','').replace('Back','').strip().replace('A','').replace('T','').replace('?','')
                    			   
                   if PosOr.count(Pos)!=0: P=PosOr.index(Pos)
                   else: P=PosOr.index(Pos+'\n')
                   NewSeq[P]='T'
                   if M.find('Back')!=-1: NewSeq[P]='A'			   
               NewSeq=''.join(NewSeq)			   
               out+=TH+'\n'+NewSeq+'\n'
               TH2Seq[TH]=NewSeq
               Seq2THls[NewSeq]=Seq2THls.get(NewSeq,[])+[TH]               
    print ('extract variable sites')
    MutP ,ExtSeqDic =	ExtractVarPos(TH2Seq) 
    ExpSeq2THls=InvertDic1(ExtSeqDic)
    out='>Outgroup\n'+('A'*len(MutP))+'\n'	
    for Seq in ExpSeq2THls:
        TH=''.join(ExpSeq2THls[Seq])
        out+='>'+TH.replace('>','')+'\n'+Seq+'\n'
    GetOut(OutF,out)
    out=''
    for M in MutP:
        out+= PosOr[M].strip()+'\n'
    GetOut(OutP,out)
  
    return MutP,ExtSeqDic    
def GV2FasMoa(GV,SNVc):
    OutF=GV[:-3]+'.fasta'
    OutP=GV[:-3]+'_Pos.txt'

    PosOr=range(0,SNVc)
    PosOr=[str(x) for x in PosOr]
    Len=len(PosOr)
    OutSeq='A'*Len    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    TH2Seq={}
    for N in Dec2Anc:
               TH='>'+N	
               N1=Dec2Anc[N]
               AllMut=[int(N.replace('\"',''))]
               if N1 !='\"root\"' and N1 !='root':AllMut.append(int(N1.replace('\"','')))
               while N1 in Dec2Anc:
                   N1=Dec2Anc[N1]
                   if N1 !='\"root\"' and N1 !='root':AllMut.append(int(N1.replace('\"','')))
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	             
               for P in AllMut:
                   NewSeq[P]='T'
               NewSeq=''.join(NewSeq)			   
               out+=TH+'\n'+NewSeq+'\n'
               TH2Seq[TH]=NewSeq	   
    print ('extract variable sites')
    MutP ,ExtSeqDic =	ExtractVarPos(TH2Seq) 
    out='>Outgroup\n'+('A'*len(MutP))+'\n'	
    for TH in ExtSeqDic:
        out+=TH+'\n'+ExtSeqDic[TH]+'\n'
    GetOut(OutF,out)

    GetOut(OutP,'\n'.join(map(str,MutP)))
    print (len(MutP))

    return MutP,ExtSeqDic    
def GetMutPos(Seq):
   Len=len(Seq)
   c=0
   Ls=[]
   while c<Len:
      if Seq[c]=='T': Ls.append(c)
      c+=1 
   return Ls       
def ExtractVarPos(CloSeq):
    MutP=[]
    CloneC=0
    for Clo in CloSeq:
        MutP+=GetMutPos(CloSeq[Clo])
        if CloSeq[Clo].find('T')!=-1: CloneC+=1
    MutP=list(set(MutP))
    MutP.sort()
    print (MutP)
    NewSeqDic={}
    for Seq in CloSeq:
        NewSeq=''
        for M in MutP:
            NewSeq+=CloSeq[Seq][M]
        NewSeqDic[Seq]=NewSeq    
    return MutP ,NewSeqDic  
def InvertDic(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=St
 return Hap2ID
def InvertDic1(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=Hap2ID.get(Seq,[])+[St]
 return Hap2ID
def InvertDic2Ls(St2Seq):
 Hap2ID=[]
 for St in St2Seq:
    Seq=St2Seq[St]
    if Hap2ID.count(Seq)==0:	
        Hap2ID.append(Seq)
 	
 return Hap2ID
def InvertDic3(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    SeqLs=St2Seq[St]
    for Seq in SeqLs:	
        Hap2ID[Seq]=Hap2ID.get(Seq,[])+[St]
 return Hap2ID 
def CountDifNum_RmMiss_AT(Seq0,Seq1): #exp, obs
            Len=len(Seq0)
         		
            Dif={'TP':0,'FN':0,'FP':0,'TN':0,'Tot':0}
            c=0
            if len(Seq0)!=len(Seq1):
                print ('skipped',Seq0,Seq1,len(Seq0),len(Seq1))
                return(9999999999999999999999999999999999999)				
            else:				
             while c<Len:
                if Seq0[c]=='A' and Seq1[c]=='T': Dif['FP']+=1
                elif Seq0[c]=='T' and Seq1[c]=='A': Dif['FN']+=1 
                elif Seq0[c]=='T' and Seq1[c]=='T': Dif['TP']+=1  
                elif Seq0[c]=='A' and Seq1[c]=='A': Dif['TN']+=1                
                if Seq0[c]!=Seq1[c] and Seq0[c]!='?' and Seq1[c]!='?': Dif['Tot']+=1
                c+=1
             return Dif 
def CountDifNum_RmMiss(Seq0,Seq1):
            Len=len(Seq0)
         		
            Dif=0
            c=0
            if len(Seq0)!=len(Seq1):
                print ('skipped',Seq0,Seq1,len(Seq0),len(Seq1))
                return(9999999999999999999999999999999999999)				
            else:				
             while c<Len:
                if Seq0[c]!=Seq1[c] and Seq0[c]!='?' and Seq1[c]!='?': Dif+=1
                c+=1
             return Dif
def GetIdenSeq(TarSeq,SeqDic):
    Iden=[]
    for i in SeqDic:
        Seq=SeqDic[i]
        Dif=CountDifNum_RmMiss(TarSeq,Seq)	
      	
        if Dif==0: Iden.append(i)
    return Iden	
def GetIdenSeqBest(TarSeq,SeqDic):
    Iden=[]
    Dif2TarLs={}
    for i in SeqDic:
        Seq=SeqDic[i]
        Dif=CountDifNum_RmMiss(TarSeq,Seq)	
        Dif2TarLs[Dif]=Dif2TarLs.get(Dif,[])+[i]
        #if Dif==0: Iden.append(i)
    DifLs=list(Dif2TarLs.keys())
  #  print (DifLs)
    DifLs.sort()
  #  print (DifLs)
    Iden=Dif2TarLs[DifLs[0]]
  #  open('a','r').readlines()
    return Iden	    
def GetBestMatSeq(THVarSeqDic,CellSeq):  

        Dif2THls={}
        for TH in THVarSeqDic:
            C=CountDifNum_RmMiss(CellSeq,THVarSeqDic[TH])
            Dif2THls[C]=Dif2THls.get(C,[])+[TH]
        DifLs=list(Dif2THls.keys())
        DifLs.sort()
     
        BestLs=Dif2THls[DifLs[0]]
        Dic={}
        for B in BestLs:
            Dic[B]=THVarSeqDic[B]
        return BestLs,Dic      
def GetBest(THVarSeqDic,CellSeq,Var):
        print (THVarSeqDic)
        TP2THls={}        
        for TH in THVarSeqDic:
            Dif=CountDifNum_RmMiss_AT(THVarSeqDic[TH],CellSeq)
            TP2THls[Dif[Var]]=TP2THls.get(Dif[Var],[])+[TH]            
        TPLs=list(TP2THls.keys())
        TPLs.sort()
   
        if Var=='TP' or Var=='TN': BestPos=TPLs[-1]  
        else: BestPos=TPLs[0]  
        BestLs=TP2THls[BestPos]
        return TP2THls,list(set(BestLs)) ,BestPos   
      
def SeqAnnoST(THfas,THclo,STcellFas,MutPFile):
    HitCloLs=[]
    Col2Val=ListCol(THclo)
    CloLs=list(set(Col2Val['Node']+['Outgroup']))
    THls,TH2Seq=ReadFasSeq(THfas)
    out=[]
    for Clo in CloLs:
        out.append('>'+Clo+'\n'+TH2Seq['>'+Clo]+'\n')
    THfas1=THfas[:-6]+'_prune.fasta'
    GetOut(THfas1,''.join(out))
    Col2Val=ListCol(MutPFile)
    MutP=list(map(int,Col2Val['label']))
    print (len(MutP)) 
    SeqNum,LabFile,MOAin1=MakeMOAin2(MutP,'',STcellFas)
    SeqAnno(THfas1,STcellFas[:-6]+'MOAIn.fasta','n')
   
def SeqAnno(THFas,OriFas,PosID):
  
    THCellLs,THVarSeqDic=ReadFasSeq(THFas)    
    OriCellLs,OriVarSeqDic=ReadFasSeq(OriFas)
    
    out=['Cell\tNode\tTP\tFP\tFN\tTN\tMissCount']
    for Cell in OriVarSeqDic:
        CellSeq=OriVarSeqDic[Cell]
        BestTHls,BestTHdic=GetBestMatSeq(THVarSeqDic,CellSeq)
        TP2THls,TPBestLs,TPC=GetBest(BestTHdic,CellSeq,'TP')
        FP2THls,FPBestLs,FPC=GetBest(BestTHdic,CellSeq,'FP')
        FN2THls,FNBestLs,FNC=GetBest(BestTHdic,CellSeq,'FN')
        TN2THls,TNBestLs,TNC=GetBest(BestTHdic,CellSeq,'TN')
        BestT=[]
        TH2FP= InvertDic3(FP2THls) 
        for TB in TPBestLs:
            if TB in FPBestLs: BestT.append(TB)
        if BestT==[]: 
            BestT=TPBestLs            
        TH2FN= InvertDic3(FN2THls)     
        FN2bestTH={}
        for BT in BestT:
           FN=TH2FN[BT][0]
           FN2bestTH[FN]=FN2bestTH.get(FN,[])+[BT]
        FNls=list(FN2bestTH.keys())
        FNls.sort()
        BestT1= FN2bestTH[FNls[0]]
        TH2TN= InvertDic3(TN2THls) 
        TNC=TH2TN[BestT1[0]][0]
        MissC=CellSeq.count('?')    
        if len(BestT1) ==1:     
            if PosID!='y':  out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)]))  
            else: 
                 obsFP=TH2FP[BestT1[0]][0]
                 obsFN=FNls[0]
                 ErrorM=1.0*obsFP/(TPC)
                 ErrorW=1.0*obsFN/(TPC)
             
                 out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)]))             
        else:
            if PosID!='y':  out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)])) 
            else: 
                 obsFP=TH2FP[BestT1[0]][0]
                 obsFN=FNls[0]
                 ErrorM=1.0*obsFP/(TPC)
                 ErrorW=1.0*obsFN/(TPC)
  
                 out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)]))      
  
    GetOut(THFas[:-6]+'_CellAnnoAll.txt','\n'.join(out)) 
def SeqAnnoBestMat(THFas,OriFas):  
    THCellLs,THVarSeqDic=ReadFasSeq(THFas)    
    OriCellLs,OriVarSeqDic=ReadFasSeq(OriFas)
    out=['Cell\tNode\tTP\tFP\tFN\tTN\tMissCount']
    for Cell in OriVarSeqDic:
        Dif2THls={}
        for TH in THVarSeqDic:
            C=CountDifNum_RmMiss(OriVarSeqDic[Cell],THVarSeqDic[TH])
            Dif2THls[C]=Dif2THls.get(C,[])+[TH]
        DifLs=list(Dif2THls.keys())
        DifLs.sort()
  
        BestLs=Dif2THls[DifLs[0]]
     
        if len(BestLs)==1: 
            
            if Cell.replace('>','').replace('#','')!='Normal': 
                  out.append(Cell.replace('>','').replace('#','')+'\t'+BestLs[0].replace('>','').replace('#',''))   
    GetOut(THFas[:-6]+'_CellAnnoAll.txt','\n'.join(out))               
def CellAnnotate(GV,THAnno,MutIDFile,OriFas): 
    if OriFas[-3:]=='meg': OriCellLs,OriCell2Seq=ReadMegSeq(OriFas)
    else: OriCellLs,OriCell2Seq=ReadFasSeq(OriFas)
    
    TH2CellLs=ReadTHAnno(THAnno,1)    
    Cell=TH2CellLs[list(TH2CellLs.keys())[0]][0]
    if OriCellLs.count('>'+Cell)==0 and OriCellLs.count('#'+Cell)==0:
        TH2CellLs=ReadTHAnno(THAnno,0) 
    MutPosLs,THVarSeqDic=GV2Fas1(GV,MutIDFile)

    OriVarSeqDic=MakeVarPosSeq(MutPosLs,OriCell2Seq)
    out=['Seq\tTopHapID']
    outG=['#MEGA\n!Title SNVs;\n!Format datatype=dna;\n']
    for Cell in OriVarSeqDic:
        Dif2THls={}
        for TH in THVarSeqDic:
            C=CountDifNum_RmMiss(OriVarSeqDic[Cell],THVarSeqDic[TH])
            Dif2THls[C]=Dif2THls.get(C,[])+[TH]
        DifLs=list(Dif2THls.keys())
        DifLs.sort()
    
        BestLs=Dif2THls[DifLs[0]]
    
        if len(BestLs)==1: 
            
            if Cell.replace('>','').replace('#','')!='Normal': 
                  outG.append('#'+Cell.replace('>','').replace('#','')+'_{'+BestLs[0].replace('>','').replace('#','')+'}'+'\n'+  OriCell2Seq[Cell])
                  out.append(Cell.replace('>','').replace('#','')+'\t'+BestLs[0].replace('>','').replace('#',''))  
      
    Len=len(OriCell2Seq[Cell])       
    GetOut(GV[:-3]+'_CellAnnoAll.txt','\n'.join(out)) 
    outG.append('#Normal_{Normal}\n'+('A'*Len))    
    GetOut(GV[:-3]+'_GroupMeg.meg','\n'.join(outG))

    os.system('megacc -a distance_estimation_between_grp_avg_nucleotide.mao -d '+GV[:-3]+'_GroupMeg.meg -o '+GV[:-3]+'_GroupMegDist.meg')    
    os.system('megacc -a infer_NJ_nucleotide.mao -d '+GV[:-3]+'_GroupMegDist.meg -o '+GV[:-3]+'_GroupMegDist.nwk')

def AdjustNoude2Cout(Dic,PriLs):
    NewDic={}
    for N in Dic:
       
       if str(N).find(';')!=-1:
          Nls=N.split(';')
          Hit=[]
          for N0 in Nls:
              if '\"'+N0+'\"' in PriLs or N0 in PriLs:
                  Hit.append(N0)
          if len(Hit)==1:                  
         
                  NewDic[int(Hit[0])]=NewDic.get(int(Hit[0]),0)+Dic[N]
               
          elif len(Hit)>1:
                  print ('multiple clones were selected',Hit,Nls,PriLs)          
                  open('a','r').readlines()
       else: NewDic[int(N)]=NewDic.get(int(N),0)+Dic[N]           
  
    return NewDic    

def GetMissPos(SNVc,MutP):
    MissPosLs=[]
    Pos=0
    while Pos<SNVc:
        if Pos not in MutP: MissPosLs.append(Pos)
        Pos+=1 
    return MissPosLs         
def MakeVarPosSeq(PosLs,SeqDic):
    NewDic={}
    for i in SeqDic:
       Seq=SeqDic[i]
       New=''
       for P in PosLs:
          New+=Seq[P]
       NewDic[i]=New
    return (NewDic)

 
def get_observedFrequency(seq_list, target_v1, target_v2):
    observedMM = []
    observedMW = []	
    for Cell in seq_list:
        sequence=seq_list[Cell]

        if sequence[target_v1] =='T' and sequence[target_v2] == 'A':
            observedMW.append(sequence)

        if sequence[target_v1] =='T' and sequence[target_v2] == 'T':
            observedMM.append(sequence)
    return  observedMM, observedMW

def get_expectedFrequency(observedMM, observedMW, beta):
    total = observedMM + observedMW
    expectedMW = total * beta
    expectedMM = total - expectedMW
    return expectedMM, expectedMW

def get_contingencyTable(oMM, oMW, exMM, exMW):
    table = np.array([[oMM, oMW],[exMM, exMW]])
    return table    
                 
def GetSeqWithMut(C2S,PosIncLs):
    NewDic={}
    for C in C2S:
        S=C2S[C]
        Inc='y'
        for P in PosIncLs:
            if S[P]!='T': Inc='n'
        if Inc=='y': NewDic[C]=S 
    return NewDic        
                   
def GetFP(ESeq,OSeq):
    FPs=[]
    Len=len(ESeq)
    c=0
    if len(ESeq)!=len(OSeq): 
        print (len(ESeq),len(OSeq))
        open('a','r').readlines()
    while c<Len:
       if ESeq[c]=='A' and OSeq[c]=='T': FPs.append(c)
       c+=1
    return FPs       
def TopHapAnnotation(TopHap,AllHap):

    Out=TopHap[:-6]+AllHap.split(os.sep)[-1][:-6]+'_anno.txt'
    
    out='Group\tSeq\tTopHapID\n'
    TopLs,TopSeq=ReadFasSeq(TopHap)
    HapLs,HapSeq=ReadFasSeq(AllHap)
    c=0
    for Hap in HapLs:
        HapS=HapSeq[Hap]
        IdenLs=GetIdenSeq(HapS,TopSeq)
  	
        c+=1	
        if len(IdenLs)==1:
     	
            out+='NA\t'+Hap.replace('>','')+'\t'+IdenLs[0].replace('>','')+'\n'
        elif len(IdenLs)>1: pass
    print (Out)
    GetOut(Out,out)  
def TopHapAnnotation1(TopHap,AllHap,Out):

    #Out=TopHap[:-6]+AllHap.split(os.sep)[-1][:-6]+'_anno.txt'
    
    out='Group\tSeq\tTopHapID\n'
    TopLs,TopSeq=ReadFasSeq(TopHap)
    HapLs,HapSeq=ReadFasSeq(AllHap)
    c=0
    for Hap in HapLs:
        HapS=HapSeq[Hap]
      #  print (len(HapS),len(TopSeq[TopLs[0]]))
        #open('a','r').readlines()
        IdenLs=GetIdenSeqBest(HapS,TopSeq)
      #  print (IdenLs)
      #  open('a','r').readlines()
        c+=1	
        if len(IdenLs)==1:
     	
            out+='NA\t'+Hap.replace('>','')+'\t'+IdenLs[0].replace('>','')+'\n'
        elif len(IdenLs)>1: pass
    print (Out)
    GetOut(Out,out) 
def TopHapAnnotation12(TopHap,AllHap,Cell2SeqHap,Out,MinCloSize):

    #Out=TopHap[:-6]+AllHap.split(os.sep)[-1][:-6]+'_anno.txt'
    
    
   TopLs,TopSeq=ReadFasSeq(TopHap)
   ConsenC=len(TopLs)
   if ConsenC==0: return 0
   else:
    HapLs,HapSeq=ReadFasSeq(AllHap)
    Rep='y'
    while Rep=='y':
     out='Group\tSeq\tTopHapID\n'
     THID2SeqID={}
     c=0
     for Hap in HapLs:
        if Hap[0]!='>': Hap='>'+Hap
        HapS=HapSeq[Hap]
      #  print (len(HapS),len(TopSeq[TopLs[0]]))
        #open('a','r').readlines()
        IdenLs=GetIdenSeqBest(HapS,TopSeq)
      #  print (IdenLs)
      #  open('a','r').readlines()
        c+=1	
        if len(IdenLs)==1:
     	
            out+='NA\t'+Hap.replace('>','')+'\t'+IdenLs[0].replace('>','')+'\n'
            THID2SeqID[IdenLs[0].replace('>','')]=THID2SeqID.get(IdenLs[0].replace('>',''),[])+[Hap.replace('>','')]
        elif len(IdenLs)>1: pass
     
     TopLsUp=[]
     TopSeqUp={}
     for THID in THID2SeqID:
         if len(THID2SeqID[THID])<=MinCloSize: print ('a small clone was removed',THID, len(THID2SeqID[THID]))
         else:
            TopLsUp.append('>'+THID)
            TopSeqUp['>'+THID]=TopSeq['>'+THID]            
   #  print ('make consensus again')
   #  ConsenFas=Out[:-4]+'.fasta'
   #  makeConsenSeq(THID2SeqID,Cell2SeqHap,ConsenFas,MinCloSize)
   #  TopLs,TopSeq=ReadFasSeq(ConsenFas) 
     if '>Outgroup' not in TopLsUp: 
         TopLsUp.append('>Outgroup')
         TopSeqUp['>Outgroup']=TopSeq['>Outgroup']         
     if len(TopLsUp)<ConsenC:
        print ('minor clone was removed and repeat')
        ConsenC=len(TopLsUp)
        TopLs=TopLsUp
        TopSeq=TopSeqUp                


     else: Rep='n'   
    FasUp=[]     
    for THID in TopSeq:

        FasUp.append(THID+'\n'+ TopSeq[THID] +'\n')
#    if 'Outgroup' not in THID2SeqID: FasUp.append('>Outgroup\n'+ HapSeq['>Outgroup'] +'\n')           
    print (Out)
    GetOut(Out,out)  
    GetOut(Out[:-4]+'.fasta',''.join(FasUp)) 
    return len(TopSeq)-1    
def TopHapAnnotation2(TopSeq,AllHap):

    #Out=TopHap[:-6]+AllHap.split(os.sep)[-1][:-6]+'_anno.txt'
    
    out='Group\tSeq\tTopHapID\n'
   # TopLs,TopSeq=ReadFasSeq(TopHap)
    HapLs,HapSeq=ReadFasSeq(AllHap)
    c=0
    Clo2CellLs={}
    for Hap in HapLs:
        HapS=HapSeq[Hap]
        IdenLs=GetIdenSeqBest(HapS,TopSeq)
  	
        c+=1	
        if len(IdenLs)==1:
            Clo2CellLs[IdenLs[0].replace('>','')]=Clo2CellLs.get(IdenLs[0].replace('>',''),[])+[Hap.replace('>','')]    	
           # out+='NA\t'+Hap.replace('>','')+'\t'+IdenLs[0].replace('>','')+'\n'
        elif len(IdenLs)>1: pass
   # print (Out)
   # GetOut(Out,out)  
    return Clo2CellLs   
def MakeFullSeq(Fas,PosFile,AnnoFile,SeqLen):
    CloLs,Clo2Seq=ReadFasSeq(Fas)
    PosLs=list(map(int,open(PosFile,'r').readlines()))
    Anno=open(AnnoFile,'r').readlines()[1:]
    Clo2CellLs={}
    for i in Anno:
        i=i.split('\t')
        Clo=i[2].strip()
        Clo2CellLs[Clo]=Clo2CellLs.get(Clo,[])+[i[1]]
    out=[]    
    for Clo in Clo2CellLs:
        if Clo!='Outgroup':
            CellLs=Clo2CellLs[Clo]
            FullSeq=''
            Seq=Clo2Seq['>'+Clo]
            c=0
            p=0
            while c<SeqLen:
                if c not in PosLs: FullSeq+='?'
                else: 
                    FullSeq+=Seq[p]
                    p+=1
                c+=1
            for Cell in CellLs:
                  out.append('>'+Cell+'\n'+FullSeq+'\n')
    GetOut(Fas[:-6]+'_CellFull.fasta',''.join(out))
    return (len(CloLs))    
            

    
    
def makeConsenSeq(clade2cellLs,CellIn2Seq,ConsenFas,MinClone):
  #  CellHapls,CellIn2Seq=ReadMegSeq(BEAMin)
  #  print (CellHapls[:5])
    out=[]
    for Cl in clade2cellLs:
     if Cl !='Outgroup':
      CellLs=clade2cellLs[Cl]

      if MinClone<len(CellLs): 
       print (Cl,len(CellLs))
       Seq=''
       c=0
       Len=len(CellIn2Seq['>'+CellLs[0]])
       while c<Len:
           NucLs=[]
           for Cell in CellLs:
               if  CellIn2Seq['>'+Cell][c]!='?': NucLs.append(CellIn2Seq['>'+Cell][c])
           if NucLs.count('A')>NucLs.count('T') and NucLs.count('A')>1: Seq+='A'
           elif NucLs.count('T')>NucLs.count('A') and NucLs.count('T')>1: Seq+='T'
           else: Seq+='?'
           c+=1
       out.append('>'+Cl+'\n'+Seq+'\n')
      else:
       print ('small clone. removed',Cl,len(CellLs))   
    if out!=[]:   
        out=['>Outgroup\n'+('A'*Len)+'\n']+out       
    GetOut(ConsenFas,''.join(out))  
def makeConsenSeqATGC(clade2cellLs,CellIn2Seq,ConsenFas,MinClone,OutSeq,Cut):
  #  CellHapls,CellIn2Seq=ReadMegSeq(BEAMin)
  #  print (CellHapls[:5])
    out=[]
    for Cl in clade2cellLs:
     if Cl !='Outgroup':
      CellLs=clade2cellLs[Cl]

      if MinClone<=len(CellLs): 
       #print (Cl,len(CellLs))
       Seq=''
       c=0
       Len=len(CellIn2Seq['>'+CellLs[0]])
       while c<Len:
           NucDic={}
           
           for Cell in CellLs:
                Nuc=CellIn2Seq['>'+Cell][c]
                if  CellIn2Seq['>'+Cell][c]!='?': NucDic[Nuc]=NucDic.get(Nuc,0)+1
           Count2Nuc=InvertDic1(NucDic)
           CountLs=list(Count2Nuc.keys())
           CountLs.sort(reverse=True)
          # print (CountLs)           
           if len(Count2Nuc[CountLs[0]])>1:Seq+='?'
           elif len(Count2Nuc[CountLs[0]])==1 and len(CountLs)==1: Seq+=Count2Nuc[CountLs[0]][0]
           else:
              if CountLs[1]<(CountLs[0]*Cut): Seq+=Count2Nuc[CountLs[0]][0]
              else: Seq+='?'

           c+=1
       out.append('>'+Cl+'\n'+Seq+'\n')
      else:
       print ('small clone. removed',Cl,len(CellLs))   
    if out!=[]:   
        out=['>Outgroup\n'+OutSeq+'\n']+out       
    GetOut(ConsenFas,''.join(out))     
def list_internal_nodes_with_two_tips(tree):
    
    internal_nodes_with_two_tips = []
    
    for node in tree.get_nonterminals():
        if len(node.clades) == 2:
            internal_nodes_with_two_tips.append([node.clades[0].name,node.clades[1].name])
    
    return internal_nodes_with_two_tips

def ReadCloAnno(File):
    Ta=open(File,'r').readlines()[1:]
    Clo2Cell={}
    for i in Ta:
        i=i.split('\t')
        Clo=i[2].strip()
        Cell=i[1]
        Clo2Cell[Clo]=Clo2Cell.get(Clo,[])+[Cell]
    return Clo2Cell
  
def MergeClade(CloFas,BEAMinHap,Tree,AnnoOut,OutTree,Cut):	  
    OriCloLs,OriClo2Seq=ReadFasSeq(CloFas)  
    root_tree(Tree,'Outgroup')
    TreeR=Tree[:-4]+'_rooted.nwk'
    print (TreeR)
    tree = Phylo.read(TreeR, "newick")
    Clone2CellLs=ReadCloAnno(AnnoOut)
    Mer='y'
    while Mer=='y':
 # Get internal nodes with two tips
      internal_nodes_with_two_tips = list_internal_nodes_with_two_tips(tree)         
      print (len(internal_nodes_with_two_tips))

      MerLs=[]
      Clo2Seq={}
      Clo2CellLsUp={}
      for Tips in internal_nodes_with_two_tips:
        if Tips[0]!=None and Tips[1]!=None:
            C1=len(Clone2CellLs.get(Tips[0],[]))
            C2=len(Clone2CellLs.get(Tips[1],[]))
            if C1<5 or C2<5:
                if C1<5: 
                   tree.prune(Tips[0])
                if C2<5: 
                   tree.prune(Tips[1])                   
            elif C1<Cut or C2<Cut:
                print ('merge',Tips,C1,C2)
                if C1>C2: 
                   tree.prune(Tips[1])
                 #  Clone2CellLs[Tips[0]]+=Clone2CellLs[Tips[1]]
                 #  Clone2CellLs[Tips[1]]=[]
                   Clo2Seq[Tips[0]]=OriClo2Seq['>'+Tips[0]]
                else: 
                   tree.prune(Tips[0])
                  # Clone2CellLs[Tips[1]]+=Clone2CellLs[Tips[0]]
                  # Clone2CellLs[Tips[0]]=[]
                   Clo2Seq[Tips[1]]=OriClo2Seq['>'+Tips[1]]                   
                MerLs.append(Tips)
        else:
            if Tips[1]!=None:
                   Clo2Seq[Tips[1]]=OriClo2Seq['>'+Tips[1]] 
            if Tips[0]!=None:                   
                   Clo2Seq[Tips[0]]=OriClo2Seq['>'+Tips[0]]                    
      print (len(MerLs)) 
      print ('make consensus clone')
      
      
      print (Clo2Seq.keys(),len(Clo2Seq))      
      if len(MerLs)==0: Mer='n'  
      else:
          print ('minor clones were merged. doing clone annotation...')
          Clone2CellLs=TopHapAnnotation2(Clo2Seq,BEAMinHap)
          print (len(Clone2CellLs))
    Phylo.write(tree, OutTree, "newick")  
    out=['Group\tSeq\tTopHapID\n']
    for Clo in Clone2CellLs:
         Group='_'.join(Clo.split('_')[1:])
         CellLs=Clone2CellLs[Clo]
         if CellLs!=[]:
           for Cell in CellLs:
             out.append(Group+'\t'+Cell+'\t'+Clo.split('_')[0]+'\n')
    GetOut(OutTree[:-4]+'.txt',''.join(out))             
def FillMissBase_consen(Fas,GoodFilLs,SupCut):

    Out=Fas[:-6]+'_FillAllCell.fasta'
    CellLs,Cell2Seq=ReadFasSeq(Fas)
    print ('tot cell c',len(CellLs))
    CellWithMiss=[]
    for Cell in CellLs:
        Seq=Cell2Seq[Cell]
        if Seq.find('?')!=-1 or Seq.find('-')!=-1:
            CellWithMiss.append(Cell)
    print ('cell with missing base',len(CellWithMiss))
    out=''
    Rec=0
    CellC=0
    TotCell=len(CellLs)
    Unique=[]
    print ('filling missing data...')
    for Cell in CellLs:

        CellC+=1	
        Seq=Cell2Seq[Cell]
        if CellWithMiss.count(Cell)==0:
             out+=Cell+'\n'+Cell2Seq[Cell]+'\n'
        else:

           MatchHap2C={}
           for Cell1 in CellLs:
               if Cell!=Cell1:
                    Seq1=Cell2Seq[Cell1]		   
                    DiffNum=CountDifNum_excMiss(Seq,Seq1)
                    if DiffNum==0:
                        MatchHap2C[Seq1]=Cell2Seq.get(Seq1,0)+1

           if MatchHap2C=={}:
               Unique.append(Cell+'\n'+Seq)	   
           Pos2Fill={}
           c=0
           Len=len(Seq)
           Filled=''	   
           UnFillC=0
           FillC=0	   
           while c<Len:
               if Seq[c]=='?' or Seq[c]=='-':
                   Fill={}
                   TotMatch=0			   
                   for Mhap in MatchHap2C:
                       Nuc=Mhap[c]
                       if Nuc!='?' and Nuc!='-':
                             Fill[Nuc]=Fill.get(Nuc,0)+MatchHap2C[Mhap]
                             TotMatch+=MatchHap2C[Mhap]						 

                   Cou2NucLs=InvertDic1(Fill)
                   Cou=list(Cou2NucLs.keys())			   
                   Cou.sort()
                   if len(Cou)==0: 
                        Filled+='?'			   
                        UnFillC+=1  
			
                   else:					
                    BestLs=Cou2NucLs[Cou[-1]]
                    Sup=1.0*Cou[-1]/TotMatch				
			
                    if len(BestLs)==1 and GoodFilLs.count(list(Fill.keys())[0])!=0 and Sup>=SupCut:
		   
                        Filled+=BestLs[0]
                        FillC+=1					
                    else: 
                        Filled+='?'			   
                        UnFillC+=1 	
					
               else: Filled+=Seq[c]						
               c+=1
           if UnFillC==0: 
               out+=Cell+'\n'+Filled+'\n'
               Rec+=1
           elif FillC!=0:
               out+=Cell+'\n'+Filled+'\n'	   
               Rec+=1	
           else: 
               out+=Cell+'\n'+Filled+'\n'	   
               Rec+=1		   
		   
    GetOut(Out,out)   
    GetOut('UniqueSeq.fasta','\n'.join(Unique))    	   
def AnnotateOriID(THanno,cellIDanno):  

    Out=THanno[:-4]+'1.txt'
    cellIDanno=open(cellIDanno,'r').readlines()[1:]
    THID2OriID={}
    for i in cellIDanno:
        i=i.split('\t')
        THID2OriID[i[0]]=i[3]

    THanno=open(THanno,'r').readlines()
    out=THanno[0]
    THanno=THanno[1:]
    for i in THanno:
        i=i.split('\t')
        TH=i[1]
        if TH.split('_')[0] in THID2OriID:	
           Ori=THID2OriID[TH.split('_')[0]]	
           out+=TH+'\t'+Ori+'\t'+i[2].replace('__','_')	

    GetOut(Out,out) 

def GetBestNuc(Dic):
    TMP={}
    for Node in Dic:
        Nuc2Prob=Dic[Node]
        BestP=0
        BestNuc='?'
        for Nuc in Nuc2Prob:
            P=Nuc2Prob[Nuc]
            if P>BestP:
                BestNuc=Nuc
                BestP=P
        TMP[Node]=BestNuc
    return TMP  
def GetCOIID(Mut):	
        Pos=int(Mut[1:][:-1])#+1
        M=Mut[-1]
        MutID=str(Pos)+M
        return MutID
def GetAncLs(MutLs):
    UpID=[]
    for i in MutLs:
      if i!='NA' and i.strip()!='':	
       UpID.append(GetCOIID(i))
    return UpID	
def GetAncMut(Anc,NodeMut,Dec2Anc):
    Ls=[]
    while Anc in Dec2Anc:
       Anc=Dec2Anc[Anc]
       if Anc in NodeMut: Ls+=NodeMut[Anc]
    return Ls	
def GetDecMut(Dec0,NodeMut,Dec2Anc):
   Anc2Dec={}
   for Dec in Dec2Anc:
      Anc=Dec2Anc[Dec]
      Anc2Dec[Anc]=Anc2Dec.get(Anc,[])+[Dec]	  
   Decs=Anc2Dec.get(Dec0,[])
   Add='y'
   while Add=='y':
      AddLs=[]
      Len1=len(Decs)	  
      for i in Decs:
         AddLs+=Anc2Dec.get(i,[])
      AddLs=list(set(AddLs))		 
      Decs+=AddLs	
      Decs=list(set(Decs))
      Len2=len(Decs)	  
      if Len1==Len2: Add='n'
   Ls=[]
   for i in Decs:
       if i in NodeMut: Ls+=NodeMut[i]   
	
   return (Ls)
   

def SumPanda(Ave):
    Ave1=Ave.to_dict()
    Comb={}
    for MeanStd in Ave1:
      
        In=Ave1[MeanStd]
     
        for i in In:
            Comb[i]=Comb.get(i,[])+[In[i]] 
    return (Comb,list(Ave1))

def ReadTHAnno(File,CellPos):	
    File=open(File,'r').readlines()[1:]
    TH2CellLs={}
    for i in File:
        i=i.split('\t')
        TH=i[2].strip()
        TH2CellLs[TH]=TH2CellLs.get(TH,[])+[i[CellPos]]
    return TH2CellLs 
def ReadTHAnno1(File,THPos,CellPos):	
    File=open(File,'r').readlines()[1:]
    TH2CellLs={}
    for i in File:
        i=i.split('\t')
        TH=i[THPos].strip().replace('>','').replace('#','')
        
        TH2CellLs[TH]=TH2CellLs.get(TH,[])+[i[CellPos].replace('>','').replace('#','')]
    return TH2CellLs       
def CompFPFN(CellFas,THAnno):
    THfas=THAnno.replace(THAnno.split(os.sep)[-1],'')+'TopHap.fasta'
    THls,TH2Seq=ReadFasSeq(THfas)
    if CellFas[-4:]=='.meg':CellLs,Cell2Seq=ReadMegSeq(CellFas)
    else: CellLs,Cell2Seq=ReadFasSeq(CellFas)
    TH2CellLs=ReadTHAnno(THAnno,1)
    Cell=TH2CellLs[list(TH2CellLs.keys())[0]][0]
    if CellLs.count('>'+Cell)==0 and CellLs.count('#'+Cell)==0:
        TH2CellLs=ReadTHAnno(THAnno,0)   
    Cout={}
    TH2count={}
    for TH in TH2CellLs:
       TH2count[TH]={}
       THseq=TH2Seq['>'+TH]
       Len=len(THseq)
       Cells= TH2CellLs[TH]
       for Cell in Cells:
          CellSeq=Cell2Seq['>'+Cell]
          if len(CellSeq)!=Len: open('a','r').readlines()
          c=0
          while c<Len: 
             Pair=THseq[c]+CellSeq[c]
             Cout[Pair]=Cout.get(Pair,0)+1 
             TH2count[TH][Pair]=TH2count[TH].get(Pair,0)+1 
             c+=1 
    FP=1.0*Cout.get('AT',0)/(Cout.get('AT',0)+Cout.get('AA',0)+Cout.get('A?',0))
    FN=1.0*Cout.get('TA',0)/(Cout.get('TA',0)+Cout.get('TT',0)+Cout.get('T?',0))       
    return FP,FN ,TH2count   
def GetInter(NodeLs,Dec2Anc):
    ComAnc=GetComAnc(NodeLs,Dec2Anc)

    Inter=[]
    Node2Inter={}
    for N in NodeLs:
        A=Dec2Anc[N]
        Node2Inter[N]=[]
        while ComAnc.count(A)==0: 
           Inter.append(A)
           Node2Inter[N].append(A)#=Node2Inter.get(N,[])+[A]
           if A not in Dec2Anc: A=ComAnc[0]
           else: A=Dec2Anc[A]
    return (Inter,Node2Inter)
def GetNuc(Rec,CellLs,Cell2SeqOri):
    Pos=int(Rec[1:][:-1])
    print (Pos,Rec,CellLs)
    Nucs=[]
    print (len(Cell2SeqOri[list(Cell2SeqOri.keys())[0]]))
    for Cell in CellLs:
      if Cell.replace('>','').replace('#','')!='Normal' and Cell.replace('>','').replace('#','')!='Outgroup':     
        if '>'+Cell in Cell2SeqOri: Nucs.append(Cell2SeqOri['>'+Cell][Pos])
        elif '#'+Cell in Cell2SeqOri: Nucs.append(Cell2SeqOri['#'+Cell][Pos])
        else:
            print (Cell,list(Cell2SeqOri.keys())[:3],len(list(Cell2SeqOri.keys())))
            open('a','r').readlines()        
    return (Nucs)    
def GetDecTH(B,N,Dec2Anc,NodeMut):
    print (B)
    Pos=int(B.strip()[1:][:-1])
    print (B,Pos)
    Anc2Dec=InvertDic1(Dec2Anc) 
    DecLs=Anc2Dec.get(N,[])
    print (N,DecLs)
    GoodDec=[]
    while DecLs!=[]:
        NextDecLs=[]
        for Dec in DecLs:
            MutLs=NodeMut[Dec]
            Good='y'
            for M in MutLs:
                if M=='A'+str(Pos)+'T': Good='n'
            if Good=='y': 
                 print ('h',GoodDec)
                 GoodDec.append(Dec)
                 if Dec in Anc2Dec: NextDecLs+=Anc2Dec[Dec]
            print (Dec,NextDecLs,MutLs,Good)
        DecLs=NextDecLs 
    print ('Res',B,N,GoodDec)        
    return GoodDec    
def GetDecCla(N,Dec2Anc):

    Anc2Dec=InvertDic1(Dec2Anc) 
    DecLs=Anc2Dec.get(N,[])
   
    AllDec=[N]
    if N in Anc2Dec:
        DecLs=[N]
        while DecLs!=[]:
           NextDecLs=[]
           for Dec in DecLs:

                 if Dec in Anc2Dec: NextDecLs+=Anc2Dec[Dec]
        
           DecLs=NextDecLs
           AllDec+=NextDecLs            
       
    return AllDec    

def GetAncTH(N,Dec2Anc):
    Als=[]
    while N in Dec2Anc:
       N=Dec2Anc[N]
       Als.append(N)
    return Als   

 
def CountTarget(DecLs,NodeMut,T):
    C=0
    for D in DecLs:
        MutLs=NodeMut[D]
        for M in MutLs:
            M=M.replace('Rec','').replace('Back','').replace(' ','').replace('A','').replace('T','')
            if int(M)==T:  C+=1
    return C            
def GetFirstConcMut(ConcMutLs,AllAnc):
   
    AllAnc.reverse()    
   
    for i in AllAnc:
        if 'A'+str(i)+'T' in ConcMutLs: 
           First=i
           return i
  #  print ('no concurrent mutation')
    return ''    
           
def CountMut(SeqLs,SeqDic,Pos):
    Dic={}
    for Seq in SeqLs:
        if '>'+Seq in SeqDic: Nuc=SeqDic['>'+Seq][Pos]
        else: Nuc=SeqDic['#'+Seq][Pos]
        Dic[Nuc]=Dic.get(Nuc,0)+1
    return Dic        
                   
    


def GetAncTilFor(Bmut,Bnode,Dec2Anc,NodeMut):
    ForMut='A'+Bmut[1:][:-1]+'T'
  #  print (Bnode,Bmut,ForMut)
    AncLs=[Bnode]	
    Find='n'	
    while Find=='n':
      if Bnode not in Dec2Anc: Find='y'	
      else:	  
       Bnode=Dec2Anc[Bnode]
       AncLs.append(Bnode)	   
       Mls=NodeMut[Bnode]	
        
       if Mls.count(ForMut)!=0: Find='y'
    ForN=Bnode	   
  
    AllDec=Getalldec(ForN,Dec2Anc)
  
    ONls=[]
    for A in AllDec:
       if AncLs.count(A)==0: ONls.append(A)
  	
    return ForN,AncLs,ONls,GetMut(AncLs,NodeMut),GetMut(ONls,NodeMut)	
def GetMut(Nls,NodeMut):
    Mls=[]
    for N in Nls:
        Mls+=NodeMut[N]
    return Mls	
	
def CompMutFre(OriFas): #OriFas_MutFre.txt
    if OriFas[-4:]=='.meg': 
        CellLs,Cell2Seq=ReadMegSeq(OriFas)
        Out=OriFas[:-4]+'_MutFre.txt' 
    else: 
        CellLs,Cell2Seq=ReadFasSeq(OriFas)
        Out=OriFas[:-6]+'_MutFre.txt'  
    Pos2C={}
    Len=len(Cell2Seq[CellLs[0]])
    for Cell in CellLs:
        Seq=Cell2Seq[Cell]
        if Seq.find('T')!=-1:
            c=0
            while c<Len:
                Nuc=Seq[c]
                if c not in Pos2C: Pos2C[c]={}
                Pos2C[c][Nuc]=Pos2C[c].get(Nuc,0)+1
                c+=1
    Pos2Fre={}
    out='Pos from 0\tMut\tWild\tMiss\tMutFre\n'
    c=0
    while c<Len:
        Cou=Pos2C[c]
        M=Cou.get('T',0)
        W=Cou.get('A',0)
        Miss=Cou.get('?',0)
        Tot=M+W
        if Tot==0: Fre=0
        else:
            Fre=1.0*M/Tot 
        Pos2Fre[c]=Fre
        out+=str(c)+'\t'+'\t'.join(map(str,[M,W,Miss,Fre]))+'\n'
        c+=1
    GetOut(Out,out)
   
    return Pos2Fre   
def GetPosMut(MutIDLab,MutIDLs): 
    MutID=MutIDLab.replace('Rec','').replace('Back','').replace('A','').replace('T','').replace('?','').strip()
    if MutIDLs.count(MutID)!=0: Pos=MutIDLs.index(MutID)
    else: Pos= MutIDLs.index(MutID+'\n')
   
    return Pos    

def FillChaMat(File,Min):
    Site2In=ListCol(File)
    Site2Keep={}
    for Site in Site2In:
       if Site!='cellBC':
        In=Site2In[Site]
        Cha2C={}
        for i in In:
           Cha2C[i]=Cha2C.get(i,0)+1
        for Cha in Cha2C:
            if Cha2C[Cha]>=Min: Site2Keep[Site]=Site2Keep.get(Site,[])+[Cha]
            else: print ('remove due to small freq',Site,Cha,Cha2C[Cha])
   # print (Site2Keep)
    SiteOr=list(Site2Keep.keys())
    out=['cellBC\t'+'\t'.join(SiteOr)+'\n']
    Len=len(Site2In['cellBC'])
    c=0
    while c<Len:
       In=[Site2In['cellBC'][c]]
       for Tar in SiteOr:
           Ind=Site2In[Tar][c]
           if Ind in Site2Keep[Tar] or Ind=='-' or Ind=='?': In.append(Ind)
           else: In.append('0')
       c+=1
       out.append('\t'.join(In)+'\n')
    OutF=File[:-4]+'_Fill'+str(Min)+'.txt'
    GetOut(OutF,''.join(out))
    return OutF    
        
def ChaMat2Fas(File):
   
    RmNote=[]
    if os.path.exists(File)!=True:
        print (File,'does not exist')
    else:	
        print (File)
        PosOut=File[:-4]+'_Position.txt'	
        Cha2Val=ListCol(File)	
        ID=File.split('\\')[-1].replace('_character_matrix.alleleThresh.txt','')
    
        File=open(File,'r').readlines()
        ChaOrder0=File[0].strip().split('\t')[1:]
        CellCol=File[0].strip().split('\t')[0]	
        Posout='Character\tState\tposition in Fas\n'
        ChaDic={}
        Pos=0
        ChaOrder=[]	
        for Cha in ChaOrder0:
            ValLs=list(set(Cha2Val[Cha]))
            ValLs.sort()
            MissC=Cha2Val[Cha].count('-')
            MissPeo=1.0*MissC/len(Cha2Val[Cha])
         	
            if MissPeo>=0.3: 
               # print ('remove',Cha,'due to >=30% missing data',MissC,len(ValLs),MissPeo)
                RmNote.append(Cha+'\t'+str(MissC)+'\t'+str(len(Cha2Val[Cha]))+'\t'+str(MissPeo)+'\n')			
            else:
             ChaOrder.append(Cha)		
    			
             ChaDic[Cha]={}	
             StC=0		
             for Val in ValLs:
              if Val!='0' and Val!='-':		
                ChaDic[Cha][Val]=Pos	
                Posout+=Cha+'\t'+Val+'\t'+str(Pos+1)+'\n'
                Pos+=1
                StC+=1			
    		
        GetOut(PosOut,Posout)			
     #   print (ChaOrder)	
    
        File=File[1:]	
        Len=len(Cha2Val[ChaOrder[0]])
        C=0
        Fasout='>Normal\n'+('A'*Pos)+'\n'	
        while C<Len:	
           NewID=Cha2Val[CellCol][C]
           if str(C).find('000')!=-1 and str(C)[-1]=='0': print (NewID,C,Len)	   
    
           Seq=''	
           	   
           for Cha in ChaOrder:
               Cha2Pos=ChaDic[Cha]	   
               State=Cha2Val[Cha][C]
               PosLs=list(Cha2Pos)
               PosLs.sort()
    	   
               if State=='-': Seq+='?'*len(Cha2Pos)	
               elif State=='0': Seq+='A'*len(Cha2Pos)
               else:
                  Find='n'		   
                  for Pos in PosLs:
                      if Pos==State: 
                          Seq+='T'
                          Find='y'					  
                      else: Seq+='A'				  
                  if Find=='n':
    
                      open('A','r').readlines()	
           MissC=Seq.count('?')
           MissPro=1.0*MissC/len(Seq)
           if MissPro>=0.3:
    
               RmNote.append(NewID+'\t'+str(MissC)+'\t'+str(len(Seq))+'\t'+str(MissPro)+'\n')			   
           else:		   
               Fasout+='>'+NewID+'\n'+Seq+'\n'				  
           C+=1	
        	   
        GetOut(PosOut[:-4]+'.fasta',Fasout)
        GetOut(PosOut[:-4]+'_Rm.txt','Site/Cell\tMissing count\tTotal count\tProportion\n'+''.join(RmNote))	
def Fas2Json2(Fa):
 
    CellLs,Cell2Seq=ReadFasSeq(Fa) 
    print ('outgroup sequence assigned: ',CellLs[0])
    OutSeqID=CellLs[0]
    OutSeq=Cell2Seq[OutSeqID]
    Len=len(OutSeq)
    GetOut(Fa[:-3]+'OutSeq.fasta',OutSeq)	
    In=[]
    for Cell in CellLs:
          Seq=Cell2Seq[Cell]
          if Seq==OutSeq:		
            print ('skiped because it is outgroup sequence: ',Cell,len(Seq))
          else:		
            VarLs=[]
            C=0
            while C<Len:
                if Seq[C]!=OutSeq[C]: VarLs+=[str(C),'\"'+Seq[C]+'\"']
                C+=1
            In.append('{'+'\"V\":['+','.join(VarLs)+'],\"I\":[\"20220101\",\"'+Cell.replace('>','')+'\",\"NA\",\"'+Cell.replace('>','')+'\",\"'+Cell.replace('>','')+'\"]}')
    out='{\"sequences\":['+','.join(In)+']}'
    GetOut(Fa[:-3]+'.json',out)	#count from 0
def Fas2Json3(Fa,OutSeq):
 
    CellLs,Cell2Seq=ReadFasSeq(Fa) 
   # print ('outgroup sequence assigned: ',CellLs[0])
    #OutSeqID=CellLs[0]
   # OutSeq=Cell2Seq[OutSeqID]
    Len=len(OutSeq)
    GetOut(Fa[:-6]+'_OutSeq.fasta',OutSeq)	
    In=[]
    for Cell in CellLs:
            Seq=Cell2Seq[Cell]
      #    if Seq==OutSeq:		
      #      print ('skiped because it is outgroup sequence: ',Cell,len(Seq))
      #    else:		
            VarLs=[]
            C=0
            while C<Len:
                if Seq[C]!=OutSeq[C]: VarLs+=[str(C),'\"'+Seq[C]+'\"']
                C+=1
            In.append('{'+'\"V\":['+','.join(VarLs)+'],\"I\":[\"20220101\",\"'+Cell.replace('>','')+'\",\"NA\",\"'+Cell.replace('>','')+'\",\"'+Cell.replace('>','')+'\"]}')
    out='{\"sequences\":['+','.join(In)+']}'
    GetOut(Fa[:-6]+'.json',out)	#count from 0    	   
    
def MakeAncSeq(Ta,Ref):

    NodeMap=Ta.replace('.csv','_nodeMap.txt')

    Out=Ta.replace('.csv','WithAnc.meg')
    NodeMap=open(NodeMap,'r').readlines()
    NodeMap=NodeMap[1:]
    Node2Label={}
    RefNodeID=''
    Dec2Anc={}
    Label2Count={}
    Label2Seq={}
    Label2Cell={}
    
    for i in NodeMap:
        i=i.strip().split('\t')
        ii=[]
        for Item in i:
            Item=Item.strip()
            if Item!='': ii.append(Item)
        Label=ii[0].strip().replace(' ','_')
        Node=ii[1].strip()
        Decs=[ii[2].strip(),ii[3].strip()]
        if Label[:5]=='Node_': 
            for Dec in Decs:
                Dec2Anc[Dec]=Node
        Label2Cell['Node_'+Node]=Label			
        Label2Count[Label]=0
        Label2Seq['Node_'+Node]=''
        if Label==Ref or Label=='Normal': RefNodeID=Node
        Node2Label[Node]=Label

    RefAnc=Dec2Anc[RefNodeID]
    Code=RefAnc in Dec2Anc
    RmNode=Node2Label[RefAnc]
    if Code==True: print ('Error: Root is incorrect. Fix the nwk file and redo ancestor analysis', Ta)
    else:
        AddNode=[]
        Ta=open(Ta,'r').readlines()
        Head=Ta[0].strip().split(',')
        Node2Col={}
        Node2Prob={}
        c=0
        Len=len(Head)
        while c<Len:
            i=Head[c].strip()
            Code=i in Label2Seq
            if Code==True: Node2Col[i]=c
            c+=1
        Ta=Ta[1:]
    
        PosiP=1
        for i in Ta:
            i=i.strip().split(',')
            PosiC=int(i[0])
            Nuc=i[1]
            if PosiC!=PosiP:
                PosiP=PosiC
                Node2BestNuc=GetBestNuc(Node2Prob)
                for Node in Node2BestNuc:
                    Label2Seq[Node]+=Node2BestNuc[Node]
                Node2Prob={}
            for Node in Node2Col:
                Prob=i[Node2Col[Node]]
                Code=Node in Node2Prob
                if Code!=True: Node2Prob[Node]={}
                Node2Prob[Node][Nuc]=float(Prob)
    
    Node2BestNuc=GetBestNuc(Node2Prob)
    for Node in Node2BestNuc:
                    Label2Seq[Node]+=Node2BestNuc[Node]

    out=''
  
    for Label in Label2Seq:

     if Label2Seq[Label]!='':
          out+='#'+Label+'CellID'+Label2Cell[Label]+'\n'+Label2Seq[Label]+'\n'

    GetOut(Out,out)     

def CleanFas(Fas,OutFas,Cut):
    CellLs,Cell2Seq=ReadFasSeq(Fas)
    out=[]
    for Cell in CellLs:
        Seq=Cell2Seq[Cell]
        MissC=1.0*Seq.count('?')/len(Seq)
        TC=Seq.count('T')
        if MissC>=Cut:pass
          #  print ('too many missing',Cell,MissC)
        elif TC==0:
             print ('No mutation',Cell)        
        else: out.append(Cell+'\n'+Seq+'\n')
    print ('original',len(CellLs),'pass',len(out))    
    out.append('>Outgroup\n'+('A'*len(Seq))+'\n')    
    GetOut(OutFas,''.join(out))        
# Define a function to traverse the tree iteratively and extract tips and branch lengths for each clade
def traverse_tree_iteratively(tree,BraCut):
    clade_info_dict = {}  # Dictionary to store clades and their respective tips and branch lengths

    # Initialize a stack for iterative traversal
    stack = [(tree,None)]  # Each item in the stack is a tuple (node, parent_branch_length)
  #  print (stack.pop())
  #  open('a','r').readlines()
    AllInter=[]
    while stack:
        node, parent_branch_length = stack.pop()

        if node.is_leaf(): pass
        else:
            # If it's an internal node, traverse its children and update clade info
          #  print (node.name,node.children)
            ChiLs=node.children
            c=0
            while ChiLs!=[]:
               # print ('child list',ChiLs)
                Up=[]
                for child in ChiLs:
                    if child.is_leaf(): pass
                    else: Up+=child.children
                if Up==[]: ChiLs=[]
                else: ChiLs=Up  
                AllInter+=Up         
    print (len(AllInter))  
    NodeID=1   
    NodeID2TipLs={}    
    for node in AllInter:
        if node.is_leaf(): pass
        else:         
    #    node, parent_branch_length = stack.pop()

        # If it's a leaf node (tip), add its name and branch length to the respective clade
    #    if node.is_leaf():
    #        clade_info_dict[node.name] = {'tips': [node.name], 'branch_length': parent_branch_length}
    #    else:
            # If it's an internal node, traverse its children and update clade info
          # print (node.name,node.dist,node.children)
           if node.dist>=BraCut:
           #open('a','r').readlines()
            ChiLs=node.children
            TipLs=[]
            c=0
            while ChiLs!=[]:
               # print ('child list',ChiLs)
                Up=[]
                for child in ChiLs:
                   # print (child,child.is_leaf())
                    if child.is_leaf():TipLs.append(child.name)
                    else: Up+=child.children
                if Up==[]: ChiLs=[]
                else: ChiLs=Up  
               # print (ChiLs)
               # if c==2:open('a','r').readlines()  
                c+=1                
                # Push the child node and its branch length to the stack for further traversal
               # stack.append((child, child.dist + (parent_branch_length or 0)))
            #print (len(TipLs),TipLs[:5])
           # if len(TipLs)>=CladeCut:
            NodeID2TipLs['Clade'+str(NodeID)]=TipLs
            NodeID+=1
    print (len(NodeID2TipLs))        
   # open('a','r').readlines()    

    return NodeID2TipLs
def get_descendant_ids(node):
    if node.is_terminal():
        return [node.name]
    else:
        descendant_ids = []
        for child in node.clades:
            descendant_ids.extend(get_descendant_ids(child))
        print (descendant_ids)
        open('a','r').readlines()        
        return descendant_ids    
def find_clades_with_n_tips(tree, n, threshold):
    clades = {}
    Cla2Boo={}
    Tipls=[]
    Anc2DecLs={}
   # clade_id_mapping = {}
    clade_id_counter = 1
    for clade in tree.find_clades(order='preorder'):
       # clade_id_mapping[clade] = 'Clade'+str(clade_id_counter)
       if clade.is_terminal(): pass
       else:       
        clade.name='Clade'+str(clade_id_counter)
        clade_id_counter += 1 
  #  print (tree)
  #  open('ar','r').readlines()    
    for clade in tree.find_clades():
        if clade.is_terminal(): Tipls.append(clade.name)
        else:
            descendant_ids = []
            for child in clade.clades:
                #descendant_ids.extend(get_descendant_ids(child)) 
                descendant_ids.append(child.name)
          #  print (clade.clades) 
          #  print (clade.name,descendant_ids) 
          #  open('a','r').readlines()            
            Anc2DecLs[clade.name] =  descendant_ids             
        if len(clade.get_terminals()) >= n and clade.confidence and clade.confidence >= threshold:
            TipInfo=clade.get_terminals()#.append(clade)
            Tips=[]
            for i in TipInfo:
                Tips.append(i.name)
            clades[clade.name]=Tips
            Cla2Boo[clade.name]=clade.confidence
    
  #  print (tree)      
    return clades,Cla2Boo,Anc2DecLs,tree#,Tipls
def GetAllDecs(Anc2DecLs):
    Anc2AllDecLs={}
    for Anc in Anc2DecLs:
        All=Anc2DecLs[Anc]
        DecLs=Anc2DecLs[Anc]
        while DecLs!=[]:
            NewDecs=[]
            for D in DecLs:
                  NewDecs+=Anc2DecLs.get(D,[])
            All+=NewDecs 
            DecLs=NewDecs
           # print (Anc,All,DecLs)
        Anc2AllDecLs[Anc]=All 
    return (Anc2AllDecLs)        
def GetCladeTip(nwk_file_path,n_tips,boo_threshold):
  #  nwk_file_path = 'E:\\Desktop\\STI\\Sim1\\out\\G7snv_1000_1_1000SNV_0.0FP_0.2FN_0.95Miss_10000cell_10clone_1_ATGC.fasta\\preds_all_samples1_rooted.nwk'  # Replace with the path to your Newick file
  #  n_tips = 10
  #  boo_threshold = 0.8
    try:
        tree = Phylo.read(nwk_file_path, 'newick')
    except FileNotFoundError:
        print(f"Error: File not found at {nwk_file_path}")
        return
    except Exception as e:
        print(f"Error reading the Newick file: {e}")
        return

    ID2cellLs,Cla2Boo,Anc2DecLs,tree = find_clades_with_n_tips(tree, n_tips, boo_threshold)
    Anc2AllDecLs=GetAllDecs(Anc2DecLs)
    ID2cellFilLs={}
    for ID in ID2cellLs:
        CellLs=ID2cellLs[ID]
        DecLs=Anc2AllDecLs[ID]
       # print (ID,DecLs)
       # open('a','r').readlines()
        RmLs=[]
        for D in DecLs:
            if D in ID2cellLs:
               #  print (ID,'found dec',D)
                 RmLs+=ID2cellLs[D]
        FilLs=[]
       # print (len(RmLs))
        for Cell in CellLs:
             if Cell not in RmLs:
                  FilLs.append(Cell)
       # print (len(CellLs),len(FilLs))
        if len(FilLs)>=10:
            ID2cellFilLs[ID]= FilLs   
       # print (FilLs)        
        
   # print (len(ID2cellFilLs),ID2cellFilLs.keys())
   # open('a','r').readlines()
  #  if not clades_with_n_tips:
  #      print(f"No clades found with {n_tips} tips.")
  #  else:
  #      print(f"Clades with {n_tips} tips:")
  #      ID2cellLs={}
  #      for idx, clade in enumerate(clades_with_n_tips, start=1):
  #          print(f"{idx}. {clade}")
  #          TipLs=clades_with_n_tips[clade]
  #          for Tip in TipLs:
  #              print (Tip.name)
  #              ID2cellLs['Cla'+str(idx)]=ID2cellLs.get('Cla'+str(idx),[])+[Tip.name]
  #      print (len(clades_with_n_tips))  
  #  print (ID2cellLs,Cla2Boo,Tipls)  
    return ID2cellLs,Cla2Boo,Anc2DecLs,tree,ID2cellFilLs#,Tipls

def GetClade(Nwk,OutGroup,BraCut,CladeCut):
    # Load the phylogenetic tree (replace 'tree_file.xml' with your file)
    print (Nwk)
    import ete3
    tree = ete3.Tree(Nwk)
    outgroup_node = tree&OutGroup
    tree.set_outgroup(outgroup_node)
    tree.write(outfile=Nwk[:-4]+'_rooted.nwk', format=1) 
    clade_info_dict = traverse_tree_iteratively(tree,BraCut)
    all_tips = tree.get_leaves()
    Clade2TipC={}
    for cl in clade_info_dict:
        TipC=len(clade_info_dict[cl])
        Clade2TipC[cl]=TipC
    TipC2ClaLs=InvertDic1(Clade2TipC)
  #  print (TipC2ClaLs)    
    TipCls=list(TipC2ClaLs.keys())
    TipCls.sort()
  #  print (TipCls)
    Done=[]
    clade_info_dict_up={}
    for TipC in TipCls:
        Clals=TipC2ClaLs[TipC]
        if TipC<CladeCut: 
            for Cla in Clals:        
                Done+=clade_info_dict[Cla]
        else:
           # print (TipC,Clals)
            for Cla in Clals:
                 TipLs=clade_info_dict[Cla]
                 UpLs=[]
                 for Tip in TipLs:
                     if Tip not in Done: UpLs.append(Tip)
                 Done+=TipLs
                 if len(UpLs)>=CladeCut: 
                      clade_info_dict_up[Cla]=UpLs
                    #  print (Cla,len(UpLs))                      
    print (len(clade_info_dict_up),clade_info_dict_up.keys())     
    return clade_info_dict_up    
def FastTree(Fas,OutID):                                   
    if Fas[-4:]=='.meg':
        CellLs,Cell2Seq=ReadMegSeq(Fas)
        out=[]
        for Cell in CellLs:
          if Cell.replace('#','')!= OutID:     
            Seq=Cell2Seq[Cell]
            out.append('>'+Cell.replace('#','')+'\n'+Seq+'\n')
        out.append('>'+OutID+'\n'+('A'*len(Seq))+'\n')
        Fas=Fas[:-4]+'.fasta'
        GetOut(Fas,''.join(out))
    if os.path.exists(Fas[:-6]+'.nwk')!=True:        
        CMP='FastTree.exe -nt \"'+Fas+'\" >\"'+Fas[:-6]+'.nwk\"'	
        os.system(CMP)	
        print (CMP)
    root_tree(Fas[:-6]+'.nwk',OutID)
    Rtree=Fas[:-6]+'_rooted.nwk'
      
    
def GetOut(OutFile,outIn):
 OutF=open(OutFile,'w')
 OutF.write(outIn)
 OutF.close()
def GetFasSeq(Dic, OutF):
   out=[]
   for i in Dic:
       S=Dic[i]
       if i[0]!='>': i='>'+i.replace('#','')
       out.append(i+'\n'+S+'\n')
   GetOut(OutF,''.join(out))    
