#from alignments.AncestralState import AncestralState
#from alignments.AncestralStatesList import AncestralStatesList
import os.path
import traceback

class nodeMapParser(object):

    """
        Parses ancestral states text files that are generated by the MEGA/MEGA-CC
        software as a part of Maximum Parsimony tree construction
        
        On success, call get_ancestral_states to a list of AncestralState objects
        generated from the ancestral states text file
    """
    def __init__(self):
        self._input_file_name = ''
     #   self._ancestral_states = AncestralStatesList()
        self._messages = []
        
    @property
    def input_file_name(self):
        return self._input_file_name
    
    @input_file_name.setter
    def input_file_name(self, value):    
        self._input_file_name = value
        
    @property
    def messages(self):
        return self._messages
    
    def parse(self):
        try:
            print 'parsing nodeMap file...'
            result = False
            if os.path.isfile(self._input_file_name) == False:
                IOError('ancestral states file not found')
            input_file = open(self._input_file_name, 'r')
            lines = input_file.readlines()
            input_file.close()
          #  self._set_filename(self._input_file_name)#lines[1])
          #  self._set_newick_string(lines[4])
           # index = 7
            #while index < len(lines):
             #   self._parse_line(lines[index])
              #  index += 1
           # result = (self._ancestral_states.num_nodes > 3) 
            self._parse_line(lines)		   
            if len(lines)>2:#result == True:
                pass#print 'successfully parsed ancestral states file: ' + str(self._ancestral_states.num_nodes) + ' taxa'
            else:
                print 'failed to parse ancestral states file'
            return result
        except Exception as e:
            traceback.print_exc()           
            self._messages.append(str(e))
            return False            
        
            
    def _parse_line(self, NodeMap):
        self._Dec2Anc={}
        self.A2D={}
        self.A2lin={}		
        self.Cell2Code={}
        self.Code2Cell={}
        for i in NodeMap:
            Ls=i.strip().split('\t')
            Line=[]
            for Item in Ls:
                Item=Item.strip()
                if Item!='':Line.append(Item)
          #  print Line			
           # for ii in Line:
            N=Line[0]
            C=Line[1]
            self.Cell2Code[N]=C
            self.Code2Cell[C]=N
            if Line[3]!='-':
                    self.A2D[C]=[Line[2],Line[3]]
                    self._Dec2Anc[Line[2]]=C
                    self._Dec2Anc[Line[3]]=C
        for A in self.A2D:
            Dlin=[]
            NewD=self.A2D[A]
            while NewD!=[]:
                NewDup=[]
                for D in NewD:
                   # print D				
                    if self.Code2Cell.has_key(D)==True: 
                          if self.Code2Cell[D].find('Node_')==-1:Dlin.append(self.Code2Cell[D])
                          else: NewDup+=self.A2D[D]						  
                    elif self.A2D.has_key(D)==True:NewDup+=self.A2D[D]
                    
                NewD=NewDup
           	self.A2lin[A]=Dlin				
    
    def get_anc2dec_lin(self):
        	
        return self.A2D, self.A2lin	

    def get_nodeMap_states(self):
        return self._Dec2Anc, self.Code2Cell, self.Cell2Code
    