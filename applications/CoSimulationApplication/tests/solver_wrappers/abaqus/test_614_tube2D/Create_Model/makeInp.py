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
from makeSurface import *
import os
mdb = Mdb(pathName='CSM_Time0.cae')
tubeModel = mdb.ModelFromInputFile(name='Model-1',inputFileName='Base.inp')
tubeMaterial = tubeModel.Material(name = 'Material');
tubeMaterial.Elastic(table= ((300000.0, 0.3), ))
tubeMaterial.Density(table=((1200.0, ), ))
tubeAssembly = tubeModel.rootAssembly
tubeInstance = tubeAssembly.instances['PART-1-1']
tubePart = tubeModel.parts['PART-1']
tubePart.setValues(space = AXISYMMETRIC, type = DEFORMABLE_BODY)
tubeModel.HomogeneousSolidSection(material='Material', name='TubeSection', thickness=1.0)
tubePart.SectionAssignment(offset=0.0, region=Region(elements=tubePart.elements), sectionName='TubeSection')
step1 = tubeModel.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=0.0001, nlgeom=ON, maxNumInc=1, haftol=1, initialInc=0.0001, minInc=0.0001, maxInc=0.0001, amplitude=RAMP, noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED, application=QUASI_STATIC)
step1.Restart(frequency = 99999, overlay  = ON)
movingSurface0 = SurfaceFromNodeSet(tubeAssembly, tubeInstance, 'BEAMINSIDEMOVING', 'MOVINGSURFACE0')
tubeModel.Pressure(name = 'DistributedPressure', createStepName = 'Step-1', distributionType = USER_DEFINED, field = '', magnitude = 1, region=movingSurface0)
tubeModel.SurfaceTraction(name = 'DistributedShear', createStepName = 'Step-1', region = movingSurface0, magnitude = 1, traction = GENERAL, directionVector = ((0,0,0), (1,0,0)), distributionType = USER_DEFINED)
tubeModel.DisplacementBC(name = 'FixedEnds', createStepName = 'Step-1', region = tubeAssembly.sets['BEAMINSIDEFIXED'], u1 = 0, u2 = 0, ur3 = UNSET)
tubeModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1', region=tubeAssembly.sets['BEAMINSIDEMOVING'], variables=('COORD', 'U'))
tubeModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
tubeModel.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1', variables=PRESELECT)
jobName = 'CSM_Time0'
tubeJob = mdb.Job(name = jobName, model = 'Model-1', description = 'Tube')
tubeJob.writeInput()
