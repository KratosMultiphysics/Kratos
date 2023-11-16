import numpy as np
import math, os
import genSampling, genDOE

def velProf(y1, y2):
    rhs = np.array([0.0, 1.0, 0.5])
    m = (y1 + y2)/2
    A = np.array([[y1**2, y1, 1.0], [y2**2, y2, 1.0], [m**2, m, 1.0]])
    b = np.linalg.solve(A,rhs)

    return lambda y: b[0]*(y**2) + b[1]*y + b[2]


doeOptions = {'mesh'        : {'path'      : {'fluid' : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/mesh/fluid',
                                              'solid' : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/mesh/solid'},
                               'meshFiles' : {'fluid' : 'mesh_channel.su2',
                                              'solid' : 'mesh_cantilever.su2'}},
              'nDOE'        : {'train' : 20,
                               'valid' : 1},
              'inletMarker' : 'inlet',
              'velProfile'  : velProf,
              'doeListString' : 'inlet_velocity',
              'velMul'      : {'train' : [10.0,20.0],
                               'valid' : [16.75]},
              'outDir'      : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/DOE'
}

samplingOptions = {'pathDOE'    : {'train' : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/DOE/train',
                                   'valid' : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/DOE/valid'},
                   'cfgOpts'    : {'cfgPath'  : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/configFiles',
                                   'cfgFiles' : {'fsi'   : 'config_fsi_steady.cfg',
                                                 'fluid' : 'config_channel.cfg',
                                                 'solid' : 'config_cantilever.cfg'},
                                   'parVar'   : {'name' : 'INLET_FILENAME',
                                                 'cfgFileType' : 'fluid'}},
                   'outDir'     : {'train' : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/samples/train',
                                   'valid' : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/samples/valid'} 
}

surrogateOptions = {'nodeInfo'              : {'fluid'              : {'input'  : (os.path.join(doeOptions['mesh']['path']['fluid'],'inlet.txt'),
                                                                                   os.path.join(doeOptions['mesh']['path']['solid'],'feabound.txt')),
                                                                       'output' : os.path.join(doeOptions['mesh']['path']['fluid'],'flowbound.txt')},

                                               'solid'              : {'input'  : os.path.join(doeOptions['mesh']['path']['fluid'],'flowbound.txt'),
                                                                       'output' : os.path.join(doeOptions['mesh']['path']['solid'],'feabound.txt')}},

                    'samples'               : {'path'               : {'train'  : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/samples/train',
                                                                       'valid'  : '/home/virtro/Desktop/Master_thesis/steady_fsi_GPR/samples/valid'},
                                               'fileSuffix'         : 'restart_'},

                    'ioVars'                : {'fluid'              : {'input'  : ('Velocity_y', 'Displacement_x', 'Displacement_y'),
                                                                       'output' : 'Pressure'},
                                               'solid'              : {'input'  : 'Pressure',
                                                                       'output' : ('Displacement_x', 'Displacement_y')}},
                    'interpolationMethod'   : 'interpolation.Gaussian',  # interpolation.Kriging, interpolation.TPS
                    'interpolationOptions'  : {'augmentation'       : 1,
                                               'regularization'     : True,
                                               'lowFidelityModel'   : None,
                                               'tuneModus'          : 'ML',  # None, 'CV', 'RandomSearch'
                                               'tuneOptions'        : None,
                                               'optimizationOptions': None,},
                    'outputDir'             : '/home/virtro/Desktop/Tutorials-master/multiphysics/steady_fsi',
}

genDOE = genDOE.genDOE(doeOptions)
genDOE.readMesh()
genDOE.write()

genSamples = genSampling.genSampling(samplingOptions)
genSamples.read()
genSamples.write()
genSamples.run()

#genSurrogate = genSurrogate.genSurrogate(surrogateOptions)
#genSurrogate.snaps()







