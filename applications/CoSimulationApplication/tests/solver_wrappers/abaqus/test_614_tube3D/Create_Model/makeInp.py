from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import os
mdb = Mdb(pathName='CSM_Time0.cae')
tubeModel = mdb.ModelFromInputFile(name='Model-1',inputFileName='Base.inp')
tubeMaterial = tubeModel.Material(name = 'Material');
tubeMaterial.Elastic(table= ((300000.0, 0.3), ))
tubeMaterial.Density(table=((1200.0, ), ))
tubeAssembly = tubeModel.rootAssembly
tubeInstance = tubeAssembly.instances['PART-1-1']
tubePart = tubeModel.parts['PART-1']
tubeModel.HomogeneousShellSection(integrationRule=SIMPSON, material='Material', name='ShellSection', nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, thickness=0.001, thicknessModulus=None, useDensity=OFF)
tubePart.SectionAssignment(offset=0.5, region=Region(elements=tubePart.elements), sectionName='ShellSection')
tubeAssembly.regenerate()
step1 = tubeModel.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=0.0001, nlgeom=ON, maxNumInc=1, haftol=1, initialInc=0.0001, minInc=0.0001, maxInc=0.0001, amplitude=RAMP, noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED, application=QUASI_STATIC)
step1.Restart(frequency = 99999, overlay  = ON)
tubeModel.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', region=tubeAssembly.sets['LEFTFIXED'], u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
tubeModel.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', region=tubeAssembly.sets['RIGHTFIXED'], u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
movingSurface0 = tubeAssembly.Surface(name='MOVINGSURFACE0', side2Elements=tubeInstance.elements)
tubeAssembly.regenerate()
tubeModel.Pressure(createStepName='Step-1', distributionType=USER_DEFINED, field='', magnitude=1.0, name='Load-1', region=movingSurface0)
tubeModel.SurfaceTraction(createStepName='Step-1', directionVector=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), distributionType=USER_DEFINED, field='', localCsys=None, magnitude=1.0, name='Load-2', region=movingSurface0, traction=GENERAL)
tubeAssembly.Set(name='WALLOUTSIDE', nodes=tubeInstance.nodes)
tubeModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1', region=tubeAssembly.sets['WALLOUTSIDE'], variables=('COORD', 'U'))
tubeModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
tubeModel.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1', variables=PRESELECT)
jobName = 'CSM_Time0'
tubeJob = mdb.Job(name = jobName, model = 'Model-1', description = 'Tube')
tubeJob.writeInput()
