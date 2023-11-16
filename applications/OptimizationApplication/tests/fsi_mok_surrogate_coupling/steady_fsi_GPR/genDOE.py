import os, sys, shutil, copy
import numpy as np
from SU2.util import ordered_bunch 
from SU2.util.ordered_dict import OrderedDict

inf = 1.0e20

# ----------------------------------------------------------------------
#  DOE Class
# ----------------------------------------------------------------------

class genDOE(ordered_bunch):
    """ config = SU2.io.Config(filename="")
        Parameters can be accessed by item or attribute
        ie: config['MESH_FILENAME'] or config.MESH_FILENAME

        Methods:
            read()       - read from a config file
            write()      - write to a config file (requires existing file)
    """
    _mesh = {}

    def __init__(self,*args,**kwarg):

        # initialize ordered bunch
        super(genDOE,self).__init__(*args,**kwarg)

    def readMesh(self):
        """ reads from a mesh file """
        try:
            msh = 'fluid'
            self._mesh['fluid'] = read_mesh(self.mesh['path']['fluid'],self.mesh['meshFiles']['fluid'],msh)
            msh = 'solid'
            self._mesh['solid'] = read_mesh(self.mesh['path']['solid'],self.mesh['meshFiles']['solid'],msh)
        except IOError:
            print('Could not find mesh file for: %s' % msh)
        except:
            print('Unexpected error: ', sys.exc_info()[0])
            raise 

    def write(self):
        pathTrain = os.path.join(self.outDir, 'train')
        if not os.path.isdir(pathTrain):
            os.mkdir(pathTrain)
        write_list(pathTrain, self.nDOE['train'], self.doeListString)
        write_doe(pathTrain, self.nDOE['train'], self.doeListString, self.velProfile, self.inletMarker, self.mesh['path']['fluid'], self.velMul['train'])

        pathValid = os.path.join(self.outDir, 'valid')
        if not os.path.isdir(pathValid):
            os.mkdir(pathValid)
        write_list(pathValid, self.nDOE['valid'], self.doeListString)
        write_doe(pathValid, self.nDOE['valid'], self.doeListString, self.velProfile, self.inletMarker, self.mesh['path']['fluid'], self.velMul['valid'])

    def __getattr__(self,k):
        try:
            return super(genDOE,self).__getattr__(k)
        except AttributeError:
            raise AttributeError('Config parameter not found')

    def __getitem__(self,k):
        try:
            return super(genDOE,self).__getitem__(k)
        except KeyError:
            raise KeyError('Config parameter not found: %s' % k)


    def local_files(self):
        """ removes path prefix from all *_FILENAME params
        """
        for key, value in self.items():
            if key.split('_')[-1] == 'FILENAME':
                self[key] = os.path.basename(value)

    def __repr__(self):
        #return '<Config> %s' % self._filename
        return self.__str__()

    def __str__(self):
        output = 'Config: %s' % self._filename
        for k,v in self.items():
            output +=  '\n    %s= %s' % (k,v)
        return output
#: class Config


# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def read_mesh(path,filename,msh):
    """ reads a mesh file """

    input_file = open(path + '/' +filename)


    # process each line
    while 1:
        # read the line
        line = input_file.readline()
        if not line:break
        # remove line returns
        line = line.strip('\r\n')
        if (len(line) == 0):continue

        if 'NDIME' in line:
            line = line.split('=')
            dims = int(line[1])
            continue
        
        if 'NPOIN' in line:
            # split across equals sign
            line = line.split('=')
            nNodes = int(line[1])
            coords = np.empty((nNodes,dims))
            i = 0
            while i<nNodes:
                # read the line
                line = input_file.readline()
                try: 
                    if not line:
                        raise ImportError
                except:print('Mesh file incomplete...no coordinate info'); raise
                # remove line returns
                line = line.strip('\r\n')
                if (len(line) == 0):continue

                line = line.split()
                try:float(line[0])
                except:print('Coordinates info incomplete'); raise
                coords[i] = [line[0], line[1]]

                i += 1

            coordsFile = open(path + '/coordinates_{}.txt'.format(msh), 'w')
            np.savetxt(coordsFile, coords)
            coordsFile.close()

            #continue
                
        if ('MARKER_TAG' in line):
            line = line.strip('\r\n').split('= ')
            markerName = line[1]
            line = input_file.readline()
            # remove line returns
            line = line.strip('\r\n')
            if 'MARKER_ELEMS' in line:
                line = line.split('=')
                nEle = int(line[1])
                markerNodes = np.array([], int)
                i = 0
                while i<nEle:
                     # read the line
                    line = input_file.readline()
                    # remove line returns
                    line = line.strip('\r\n').split()
                    markerNodes = np.append(markerNodes,[int(line[1]), int(line[2])],0)
                    i +=1
                markerNodes = np.unique(markerNodes)
                markerNodesFile = open(path + ('/{}.txt'.format(markerName)), 'w')
                #markerNodesFile.write(msh + '\n')
                np.savetxt(markerNodesFile, markerNodes, '%d')
                markerNodesFile.close()

    return

#: def read_mesh()

# -------------------------------------------------------------------
#  Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def write_list(outDir, nDOE, string):
    
    doeList = open(outDir + '/listDOE.txt','w')
    for i in range(nDOE):
        str = outDir + '/' + string + '_{}'.format(i) + '_0' + '.dat'
        doeList.write(str + "\n")
    doeList.close()

#: def write_list()

def write_doe(outDir,nDOE, string, velProf, marker, meshPath, velMul):
    if nDOE == 1:
        velmul = velMul
    else:
        velmul = np.linspace(velMul[0], velMul[1], num=nDOE, endpoint=True)
    for i in range(nDOE):
        coords = np.loadtxt(meshPath + '/coordinates_fluid.txt', dtype=float)
        markerNodes = np.loadtxt(meshPath + '/{}.txt'.format(marker), dtype=int, skiprows=0)
        velFunc = velProf(min(coords.T[1]), max(coords.T[1]))
        str = string + '_{}'.format(i) + '_0'
        doe = open(outDir + '/' + str + '.dat','w')
        doe.write('NMARK= 1\n')
        doe.write('MARKER_TAG= {}\n'.format(marker))
        doe.write('NROW= {}\n'.format(len(markerNodes)))
        doe.write('NCOL= {}\n'.format((2*(np.shape(coords)[1]) + 2)))
        coordinates = coords[markerNodes]
        velMag = np.empty([len(markerNodes)])
        temperature = np.empty([len(markerNodes)])
        xVelVec = np.empty([len(markerNodes)])
        yVelVec = np.empty([len(markerNodes)])
        for j in range(len(velMag)):
            velMag[j] = velFunc(coordinates[j][1]) * velmul[i]
            temperature[j] = 1.0
            xVelVec[j] = 1.0
            yVelVec[j] = 0.0
        doeData = np.vstack((coordinates.T,temperature.T,velMag.T,xVelVec.T,yVelVec.T))
        np.savetxt(doe,doeData.T,fmt='%f')
        doe.close()

