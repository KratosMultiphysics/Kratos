from numpy import *

class LsDynaSet:
    def __init__(self):
        self.node_ids = []
        self.sid = -1
        self.sid = -1


class LsDynaMeshConverter:
    '''very simple function to parse a Gid mesh file and extract nodal coordinates and tetrahedra connectivity'''
    def __init__(self,inputfilename, outputfilename):
        try:
            self.input_file = open( inputfilename)
        except:
            raise Exception("attempting to open an inexisting file")
        self.out_file = open( outputfilename, 'w')
        
        self.coordinates = []
        self.NodeIds = []
        
        self.TetraIds = []
        self.TetraConnectivity = []   
        self.TetraProp = []
        self.TriangleIds = []
        self.TriangleConnectivity = []
        self.TriProp = []
        self.sets_in_file = []
                
    def Read(self):
        
        
        if self.input_file:
            
            self.input_file.seek(0)
            for line in self.input_file:
                if line.find("*MESH_VOLUME_NODE") != -1:
                    self.ReadNodes()
                    
            self.input_file.seek(0)
            for line in self.input_file:        
                if line.find("*MESH_VOLUME_ELEMENT") != -1:
                    self.ReadElements()
                
            #read all sets
            ##count sets in file
            self.input_file.seek(0)
            number_of_sets = 0
            line_counter = 0
            begin_of_set = []
            for line in self.input_file:
                if line.find("*ICFD_SET_NODE_LIST") != -1:
                    number_of_sets+=1
                    begin_of_set.append(line_counter)
                line_counter += 1
            print("sets found = ",number_of_sets)
            print("begin of set ",begin_of_set)

            
            
            for i in range(0,number_of_sets):
                #print("current beginning at line ",begin_of_set[i])
                self.input_file.seek(0)
                line_counter = 0
                for line in self.input_file:
                    line_counter+=1
                    if(line_counter >= begin_of_set[i]):
                        if line.find("*ICFD_SET_NODE_LIST") != -1:
                            myset = self.ReadSet()
                            self.sets_in_file.append(myset)
                            break
                    
    def ReadWords(self,line):
        i  = line.find("//")
        if i != -1:
            return line[:i].split() 
        return line.split() 
         
               
    def ReadNodes(self):
        counter = 0
        for line in self.input_file:
            if line.find("*") == -1:
                if line.find("$") != -1: #line of comments
                    pass
                else:
                    counter += 1
                    words = self.ReadWords(line)
                    self.NodeIds.append( int(words[0]) )
                    coordinates = [float(words[1]), float(words[2]), float(words[3])]
                    self.coordinates.append( coordinates )
            else:
                break
        print("read ",counter, " Nodes in total")
 
        
    def ReadElements(self):
        tetra_counter = 0
        tri_counter = 0
        for line in self.input_file:
            if line.find("*") == -1:
                if line.find("$") != -1: #line of comments
                    pass
                else:                
                    words = self.ReadWords(line)
                    #print (words)
                    el_id, prop_id, connectivity = int(words[0]), int(words[1]), [int(x) for  x in words[2:]]
                    
                    connectivity_lenght = 1
                    for i in range(1,len(connectivity)):
                        if connectivity[i] != connectivity[i-1]:
                            connectivity_lenght+=1
                        else:
                            break
                    
                    active_connectivity = [int(x) for  x in connectivity[0:(connectivity_lenght-2)]]
                    
                    #print(active_connectivity)
                    if(connectivity_lenght == 4): #tetrahedra
                        tetra_counter += 1
                        self.TetraIds.append( el_id )
                        self.TetraConnectivity.append(active_connectivity)
                        self.TetraProp.append(prop_id)
                        
                    elif(connectivity_lenght == 3):
                        tri_counter += 1
                        self.TriangleIds.append( el_id )
                        self.TriangleConnectivity.append(active_connectivity)
                        self.TriProp.append(prop_id)
                    else:
                        print("last line was ",words)
                        print(active_connectivity)
                        print(connectivity)
                        raise Exception("sorry the connectivity lenght was unexpected : ",connectivity_lenght)
            else:
                break
        print("read ",tetra_counter, " Tetras in total")
        print("read ",tri_counter, " Triangles in total")
        
    def ReadSet(self):
        newset = LsDynaSet()
        
        for line in self.input_file:
            if line.find("*") == -1:
                if line.find("$") != -1: #line of comments
                    pass
                else:
                    words = self.ReadWords(line)
                    if(len(words) == 2):
                        #print(words)
                        sid = int(words[0])
                        pid = int(words[1])
                        newset.sid = sid
                        newset.pid = pid
                    else:
                        #print(words)
                        for word in words:
                            myid = int(word)
                            if(myid != 0):
                                newset.node_ids.append(myid)
            else:
                break
                 
        return newset
  
    def WritePropertyBlock(self):
        #count the properties which are unique
        unique_properties = set(self.TetraProp)
        print("properties found =",unique_properties)
        
        for property_id in unique_properties:
            self.out_file.write("Begin Properties " + str(property_id) + "\n")
            self.out_file.write("End Properties\n")


    def WriteNodes(self):
        self.out_file.write("Begin Nodes\n")
        
        for i in range(0,len(self.NodeIds) ):
            coords = self.coordinates[i]
            self.out_file.write( str(self.NodeIds[i]) + " " + str(float(coords[0])) + " " + str(float(coords[1])) + " " + str(float(coords[2])) + "\n" )
        
        self.out_file.write("End Nodes\n")
        
    def WriteTetras(self, elementname):
        self.out_file.write("Begin Elements " + elementname +"\n")
        
        for i in range(0,len(self.TetraIds)):
            connectivity = self.TetraConnectivity[i]
            self.out_file.write( str(self.TetraIds[i]) + " " + str(self.TetraProp[i]) + " ")
            for word in connectivity:
                self.out_file.write( str(word) + " " )
            self.out_file.write(  "\n" )
        
        self.out_file.write("End Elements\n")
        

    def WriteSet(self,myset):
        self.out_file.write("Begin Mesh "+str(myset.sid)+"\n")
        self.out_file.write("Begin MeshNodes\n")
        
        for i in range(0,len(myset.node_ids) ):
            coords = self.coordinates[i]
            self.out_file.write( str(myset.node_ids[i]) + "\n" )
        
        self.out_file.write("End MeshNodes\n")
        self.out_file.write("End Mesh\n")
 
    def WriteMeshes(self):
        for myset in self.sets_in_file:
            self.WriteSet(myset)

                   
                   