from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class AuxiliarMethodsAdaptiveRemeshing(object):
    """
    This class is an auxiliar script when using adaptative remeshing put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, analysis):
        """This function is the constructor of the class

            Keyword arguments:
            self It signifies an instance of a class.
            analysis The AnalysisStage to be computed
        """
        # Saves the analysis
        self.analysis = analysis

    def AdaptativeRemeshingDetectBoundary(self):
        """This function detects the boundary to preserve pure node BC

            Keyword arguments:
            self It signifies an instance of a class.
        """
        solver = self.analysis._GetSolver()
        computing_model_part = solver.GetComputingModelPart()
        if not self.analysis.process_remesh:
            convergence_criteria = solver.get_convergence_criterion()
            convergence_criteria.Initialize(computing_model_part)

        # Ensuring to have conditions on the BC before remesh
        is_surface = False
        for elem in computing_model_part.Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            break

        if not is_surface:
            list_model_parts = []
            # We need to detect the conditions in the boundary conditions
            if self.analysis.project_parameters.Has("constraints_process_list"):
                constraints_process_list = self.analysis.project_parameters["constraints_process_list"]
                for i in range(0,constraints_process_list.size()):
                    item = constraints_process_list[i]
                    list_model_parts.append(item["Parameters"]["model_part_name"].GetString())
            skin_detection_parameters = KratosMultiphysics.Parameters("""
            {
                "list_model_parts_to_assign_conditions" : []
            }
            """)
            root_model_part_name = computing_model_part.GetRootModelPart().Name
            for name in list_model_parts:
                name = name.replace(root_model_part_name + ".", "")
                skin_detection_parameters["list_model_parts_to_assign_conditions"].Append(name)

            if computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                detect_skin = KratosMultiphysics.SkinDetectionProcess2D(computing_model_part, skin_detection_parameters)
            else:
                detect_skin = KratosMultiphysics.SkinDetectionProcess3D(computing_model_part, skin_detection_parameters)
            detect_skin.Execute()
        solver.SetEchoLevel(self.analysis.echo_level)
