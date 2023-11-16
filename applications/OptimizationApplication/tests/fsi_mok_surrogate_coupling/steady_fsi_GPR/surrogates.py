import os
import genSurrogate


surrogateOptions = {'ioData'               :{'train'  : '/home/virtro/Kratos/bin/mpi-gcc/applications/OptimizationApplication/tests/fsi_mok/surrogate_io/train.xlsx',
                                             'valid'  : '/home/virtro/Kratos/bin/mpi-gcc/applications/OptimizationApplication/tests/fsi_mok/surrogate_io/valid.xlsx'
                                             },

                    'fluid'                 : {'input'  : ('YOUNG_MODULUS','DISPLACEMENT_X','DISPLACEMENT_Y'),
                                               'output' : ('PRESSURE',)},
                    'solid'                 : {'input'  : ('PRESSURE',),
                                               'output' : ('DISPLACEMENT_X','DISPLACEMENT_Y')},
                    'interpolationMethod'   : 'interpolation.TPS',  # interpolation.Kriging, interpolation.TPS
                    'interpolationOptionsFluid'  : {'augmentation'       : -1,
                                               'regularization'     : True,
                                               'lowFidelityModel'   : None,
                                               'scale'              : False,
                                               'tuneModus'          : None,  # None, 'CV', 'RandomSearch'
                                               'tuneOptions'        : None,
                                               'optimizationOptions': {'methodOptions' : {'popsize' : 5}}},
                    'interpolationOptionsSolid'  : {'augmentation'       : 1,
                                               'regularization'     : True,
                                               'lowFidelityModel'   : None,
                                               'scale'              : False,
                                               'tuneModus'          : None,  # None, 'CV', 'RandomSearch'
                                               'tuneOptions'        : None,
                                               'optimizationOptions': {'methodOptions' : {'popsize' : 5}}},#, 'disp' : True}}},
                    'outputDir'             : '/home/virtro/Desktop/Tutorials-master/multiphysics/steady_fsi'
}

genSurrogate = genSurrogate.genSurrogate(surrogateOptions)
genSurrogate.snaps()
genSurrogate.surrogates()
genSurrogate.predict()
