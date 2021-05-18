from numpy import *




class GidMeshConverter:
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
        self.TriangleIds = []
        self.TriangleConnectivity = []
        
    def ReadAndWriteResult(self,resfile,varname):
        tmp = open( resfile)
        
        for line in tmp:
                if line.find(varname) != -1:
                    ids,values = self.ReadResult(tmp,varname)
                                       
                    self.out_file.write("Begin NodalData "+varname+"\n")
                    
                    for i in range(0,len(ids) ):
                        self.out_file.write( str(ids[i]) + " 0 " + str(values[i]) + "\n" )
                    
                    self.out_file.write("End NodalData \n")
        
                    break
        
        
                    
        
        tmp.close()
        
    def Read(self):
        if self.input_file:
            for line in self.input_file:
                if line.find("Coordinates") != -1:
                    self.ReadNodes()
                if line.find("Elements") != -1:
                    self.ReadElements()
                                                            
    def ReadWords(self,line):
        i  = line.find("//")
        if i != -1:
            return line[:i].split() 
        return line.split() 
         
               
    def ReadNodes(self):
        counter = 0
        for line in self.input_file:
            if line.find("End") == -1:
                counter += 1
                words = self.ReadWords(line)
                self.NodeIds.append( int(words[0]) )
                coordinates = [float(words[1]), float(words[2]), float(words[3])]
                self.coordinates.append( coordinates )
            else:
                break
        print("read ",counter, " Nodes in total")
 
    def ReadResult(self,inputfile,varname):
        counter = 0
        ids = []
        values = []
        for line in inputfile:
            if line.find("End") == -1:
                if line.find("Values") == -1:
                    counter += 1
                    words = self.ReadWords(line)
                    ids.append( int(words[0]) )
                    values.append( float(words[1]))
                else:
                    pass
            else:
                break
            
        print("read ",counter, " values in total")
        
        return ids,values
        
    def ReadElements(self):
        tetra_counter = 0
        tri_counter = 0
        for line in self.input_file:
            if line.find("End") == -1:
                
                words = self.ReadWords(line)
                #print (words)
                el_id, connectivity = int(words[0]), [int(x) for  x in words[1:]]
                
                if(len(connectivity) == 4):
                    tetra_counter += 1
                    self.TetraIds.append( el_id )
                    self.TetraConnectivity.append(connectivity)
                elif(len(connectivity) == 3):
                    tri_counter += 1
                    self.TriangleIds.append( el_id )
                    self.TriangleConnectivity.append(connectivity)
        print("read ",tetra_counter, " Tetras in total")
        print("read ",tri_counter, " Triangles in total")
  
    def WritePropertyBlock(self,property_id):
        self.out_file.write("Begin Properties " + str(property_id) + "\n")
               
        self.out_file.write("End Properties\n")


    def WriteNodes(self):
        self.out_file.write("Begin Nodes\n")
        
        for i in range(0,len(self.NodeIds) ):
            coords = self.coordinates[i]
            self.out_file.write( str(self.NodeIds[i]) + " " + str(float(coords[0])) + " " + str(float(coords[1])) + " " + str(float(coords[2])) + "\n" )
        
        self.out_file.write("End Nodes\n")
        
    def WriteTetras(self,property_id, elementname):
        self.out_file.write("Begin Elements " + elementname +"\n")
        
        for i in range(0,len(self.TetraIds)):
            connectivity = self.TetraConnectivity[i]
            self.out_file.write( str(self.TetraIds[i]) + " " + str(property_id) + " ")
            for word in connectivity:
                self.out_file.write( str(word) + " " )
            self.out_file.write(  "\n" )
        
        self.out_file.write("End Elements\n")


              

                   
                   