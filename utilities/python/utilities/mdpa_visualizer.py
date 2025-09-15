#!/usr/bin/env python
import sys

if len(sys.argv) != 2:
    raise RuntimeError("Please provide a mdpa name to create vtk visualization.")

mdpa_name = sys.argv[1]
if mdpa_name.endswith(".mdpa"):
    mdpa_name = mdpa_name[:-5]

import KratosMultiphysics as Kratos
try:
    import KratosMultiphysics.StructuralMechanicsApplication
except:
    pass
try:
    import KratosMultiphysics.FluidDynamicsApplication
except:
    pass
try:
    from KratosMultiphysics.HDF5Application.create_xdmf_file import WriteMultifileTemporalAnalysisToXdmf
except:
    pass

from KratosMultiphysics.vtk_output_process import VtkOutputProcess
from hdf5_utilities import GetHDF5File, OutputModelPartToHDF5, OutputNodalResultsToHDF5

model = Kratos.Model()
model_part = model.CreateModelPart(mdpa_name)

Kratos.ModelPartIO(mdpa_name).ReadModelPart(model_part)

default_settings = Kratos.Parameters("""{"model_part_name":"", "save_output_files_in_folder": false}""")
default_settings["model_part_name"].SetString(mdpa_name)
vtk_output_process = VtkOutputProcess(model, default_settings)
vtk_output_process.PrintOutput()

h5_file = GetHDF5File(mdpa_name + ".h5", "truncate")
OutputModelPartToHDF5(model_part, h5_file)
OutputNodalResultsToHDF5(model_part, h5_file, [])
del h5_file
WriteMultifileTemporalAnalysisToXdmf(mdpa_name + ".h5", "/ModelData", "/ResultsData")


