# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import KratosMultiphysics.kratos_utilities as kratos_utilities
if kratos_utilities.CheckIfApplicationsAvailable("MeshingApplication"):
    import KratosMultiphysics.MeshingApplication as MeshingApplication
    missing_meshing_dependencies = True
else:
    missing_meshing_dependencies = False

def CreateRemeshingProcess(main_model_part, adaptative_remesh_parameters):
    if main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
        remeshing_process = MeshingApplication.MmgProcess2D(main_model_part, adaptative_remesh_parameters["remeshing_parameters"])
    else:
        is_surface = False
        for elem in main_model_part.Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            break
        if is_surface:
            remeshing_process = MeshingApplication.MmgProcess3DSurfaces(main_model_part, adaptative_remesh_parameters["remeshing_parameters"])
        else:
            remeshing_process = MeshingApplication.MmgProcess3D(main_model_part, adaptative_remesh_parameters["remeshing_parameters"])

    return remeshing_process

def CreateMetricProcess(main_model_part, adaptative_remesh_parameters):
    if main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
        metric_process = MeshingApplication.MetricErrorProcess2D(main_model_part, adaptative_remesh_parameters["metric_error_parameters"])
    else:
        metric_process = MeshingApplication.MetricErrorProcess3D(main_model_part, adaptative_remesh_parameters["metric_error_parameters"])

    return metric_process
