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

os.system('Convert Mesh.neu CPE8R')
mdb = Mdb(pathName='Project.cae')
beamModel = mdb.ModelFromInputFile(name='Model-1',inputFileName='Mesh.inp')
beamMaterial = beamModel.Material(name = 'Material');
beamMaterial.Elastic(table= ((10000000.0, 0.49), ))
beamMaterial.Density(table=((1100.0, ), ))
beamAssembly = beamModel.rootAssembly
beamInstance = beamAssembly.instances['PART-1-1']
beamPart = beamModel.parts['PART-1']
beamModel.HomogeneousSolidSection(material='Material', name='BeamSection', thickness=1.0)
beamPart.SectionAssignment(offset=0.0, region=Region(elements=beamPart.elements), sectionName='BeamSection')
step1 = beamModel.ImplicitDynamicsStep(name = 'Step-1', previous = 'Initial', timePeriod = 0.01, nlgeom = ON, maxNumInc = 10000, haftol = 0.1, initialInc = (0.01/1000.0), amplitude = RAMP, minInc = (0.01/1000000000.0), maxInc=0.009, application=MODERATE_DISSIPATION)
step1.Restart(frequency = 99999, overlay  = ON)
movingSurface0 = SurfaceFromNodeSet(beamAssembly, beamInstance, 'BEAMINSIDEMOVING0', 'MOVINGSURFACE0')
movingSurface1 = SurfaceFromNodeSet(beamAssembly, beamInstance, 'BEAMINSIDEMOVING1', 'MOVINGSURFACE1')
movingSurface2 = SurfaceFromNodeSet(beamAssembly, beamInstance, 'BEAMINSIDEMOVING2', 'MOVINGSURFACE2')
beamModel.Pressure(name = 'DistributedPressure0', createStepName = 'Step-1', distributionType = USER_DEFINED, field = '', magnitude = 1.0, region = movingSurface0)
beamModel.Pressure(name = 'DistributedPressure1', createStepName = 'Step-1', distributionType = USER_DEFINED, field = '', magnitude = 1.0, region = movingSurface1)
beamModel.Pressure(name = 'DistributedPressure2', createStepName = 'Step-1', distributionType = USER_DEFINED, field = '', magnitude = 1.0, region = movingSurface2)
beamModel.SurfaceTraction(name = 'DistributedShear0', createStepName = 'Step-1', region = movingSurface0, magnitude = 1, traction = GENERAL, directionVector = ((0,0,0), (1,0,0)), distributionType = USER_DEFINED)
beamModel.SurfaceTraction(name = 'DistributedShear1', createStepName = 'Step-1', region = movingSurface1, magnitude = 1, traction = GENERAL, directionVector = ((0,0,0), (1,0,0)), distributionType = USER_DEFINED)
beamModel.SurfaceTraction(name = 'DistributedShear2', createStepName = 'Step-1', region = movingSurface2, magnitude = 1, traction = GENERAL, directionVector = ((0,0,0), (1,0,0)), distributionType = USER_DEFINED)
if 9.81 > 0.0:
	beamModel.Gravity(name = 'Gravity', createStepName = 'Step-1', comp2 = -9.81)
beamModel.DisplacementBC(name = 'FixedTopEnd', createStepName = 'Step-1', region = beamAssembly.sets['BEAMFIXED'], u1 = 0, u2 = 0, ur3 = 0)
beamModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1', variables=('COORD', 'U'))
beamModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
beamModel.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1', variables=PRESELECT)
jobName = 'Base'
beamJob = mdb.Job(name = jobName, model = 'Model-1', description = 'Beam')
beamJob.writeInput()
