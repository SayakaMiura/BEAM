from alignments.MegaAlignment import MegaAlignment
from Bio import Phylo
from ML.MegaML import MegaML
from ML.MegaAncestor import MegaAncestor
class TreeAnalizer:

    def find_branch_with_mutation(self, ancestor_states, offspring2ancestor, code2cell):
        SNV2branch={}
        Len=len(ancestor_states[list(ancestor_states.keys())[0]])
        Posi=0
     #   print(list(ancestor_states.keys()),code2cell)		
        while Posi<Len:
            MutBra=''
            for Dec in offspring2ancestor:
              if (Dec in code2cell)==True:			
               if code2cell[Dec]!='Normal':			
                   Anc=offspring2ancestor[Dec]
                   if (code2cell[Dec] in ancestor_states)!=True: Seq=self.CellExtant2Seq['#'+code2cell[Dec]]
                   else: Seq=ancestor_states[code2cell[Dec]]	
                #   print Seq[Posi].split('\t')[0], ancestor_states[code2cell[Anc]][Posi].split('\t')[0]				   
                   if Seq[Posi].split('\t')[0]=='T' and ancestor_states[code2cell[Anc]][Posi].split('\t')[0]=='A':
                         MutBra=Dec	
            SNV2branch[Posi]=MutBra
            Posi+=1
        return SNV2branch 			
		
    def RootTree(self, OriNwk, Root):
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()
        trees = list(Phylo.parse('test.nwk', 'newick'))	
        for tree in trees:
             tree = tree.root_with_outgroup({'name': Root})		 
        Phylo.write(trees, 'newtree.nwk', "newick")
        Tree=open('newtree.nwk','r').readlines()[0].strip()
        Len=len(Tree)
        Posi=Tree.find(','+Root+':0.00000')
        PosRev=-1*(Len-Posi)
        
        LastBraLen=''
        Rm=''
        while Tree[PosRev]!=':':
            LastBraLen+=Tree[PosRev]
            PosRev=PosRev-1
        BraLen= LastBraLen[::-1]
        NewTree='('+Root+':'+BraLen+Tree[2:].replace('):'+BraLen+Root+':0.00000','')+'\n'	
        return NewTree		



    def RootTree_rootBottom(self, OriNwk, RootTaxa):	
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse('test.nwk', 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': RootTaxa})
        Phylo.write(trees, 'test1.nwk', "newick")	

        return open('test1.nwk','r').readlines()[0]	
  
