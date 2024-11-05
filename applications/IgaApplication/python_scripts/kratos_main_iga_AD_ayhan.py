# Importing Kratos
from weakref import ref
import KratosMultiphysics
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

#Quadrature_points
import sys
import os
from get_gauss_integration_points_regular_element import evaluate_face_integration_points
from get_gauss_integration_points_regular_element import evaluate_shape_functions_and_derivatives
from catmull_clark_subdivision_constant_boundaries import catmull_clark_subdivision
from catmull_clark_subdivision_constant_boundaries import order_control_points_in_faces
class StructuralMechanicsAnalysis(AnalysisStage):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        solver_settings = project_parameters["solver_settings"]

        if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            raise Exception("StructuralMechanicsAnalysis: " + '"domain_size" defined both in "problem_data" and "solver_settings"!')

        if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            raise Exception("StructuralMechanicsAnalysis: " + '"model_part_name" defined both in problem_data" and "solver_settings"!')

        if solver_settings.Has("time_stepping") and project_parameters["problem_data"].Has("time_Step"):
            raise Exception("StructuralMechanicsAnalysis: " + '"time_stepping" defined both in "problem_data" and "solver_settings"!')

        if not solver_settings.Has("time_stepping"):
            raise Exception("StructuralMechanicsAnalysis: Using the old way to pass the time_step, this was removed!")

        if not solver_settings.Has("domain_size"):
            raise Exception("StructuralMechanicsAnalysis: Using the old way to pass the domain_size, this was removed!")

        # Detect is a contact problem
        # NOTE: We have a special treatment for contact problems due to the way the convergence info is printed (in a table). Not doing this will provoque that the table is discontinous (and not fancy and eye-candy)
        solver_settings = project_parameters["solver_settings"]
        self.contact_problem = solver_settings.Has("contact_settings") or solver_settings.Has("mpc_contact_settings")

        super().__init__(model, project_parameters)

    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """

        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel() #1
        self._ModelersPrepareGeometryModel() #2

        ######################
        ## Refinement step# ##

        ########################## Extract data from input ###################################
        if os.path.isfile("physics.iga.json"):
            physics_file = open("physics.iga.json",'r')
            PhysicsParameters = KratosMultiphysics.Parameters(physics_file.read())

        surface_id = 0
        sub_model_part = ""

        # Assign material to model_parts (if Materials.json exists)
        if os.path.isfile("materials.json"):
            materials_file = open("materials.json",'r')
            MaterialParameters = KratosMultiphysics.Parameters(materials_file.read())
        
        # create property for shell elements
        shell_properties = self.model.GetModelPart("IgaModelPart").GetProperties()[1]
        shell_properties.SetValue(KratosMultiphysics.THICKNESS, MaterialParameters["properties"][0]["Material"]["Variables"]["THICKNESS"].GetDouble())
        shell_properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, MaterialParameters["properties"][0]["Material"]["Variables"]["YOUNG_MODULUS"].GetDouble())
        shell_properties.SetValue(KratosMultiphysics.POISSON_RATIO, MaterialParameters["properties"][0]["Material"]["Variables"]["POISSON_RATIO"].GetDouble())
        shell_properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw())

        # create property for shell elements
        load_properties = self.model.GetModelPart("IgaModelPart").GetProperties()[2]

        sub_model_part_bc = []
        sub_model_part_bc_param = []
        sub_model_part_bc_type = []
        
        for item in PhysicsParameters["element_condition_list"]:
            if item["geometry_type"].GetString() == "GeometrySurface" and item["parameters"]["type"].GetString() == "element":
                surface_id = item["brep_ids"][0].GetInt()
                sub_model_part = item["iga_model_part"].GetString()
            if item["geometry_type"].GetString() == "GeometrySurface" and item["parameters"]["type"].GetString() == "condition":
                surface_id_cond = item["brep_ids"][0].GetInt()
                sub_model_part_cond = item["iga_model_part"].GetString()
            if item["geometry_type"].GetString() == "GeometrySurfaceNodes" or item["geometry_type"].GetString() == "GeometrySurfaceVariationNodes":
                sub_model_part_bc.append(item["iga_model_part"].GetString())
                local_parameters = []
                local_parameters.append(item["parameters"]["local_parameters"][0].GetDouble())
                local_parameters.append(item["parameters"]["local_parameters"][1].GetDouble())
                sub_model_part_bc_type.append(item["geometry_type"].GetString())
                sub_model_part_bc_param.append(local_parameters)

        surface = self.model.GetModelPart("IgaModelPart").GetGeometry(surface_id) 
        nurbs_surface = surface.GetGeometryPart(KratosMultiphysics.Geometry.BACKGROUND_GEOMETRY_INDEX)
        print('Nurbs_surface: \n',nurbs_surface)

        ########################### Subdivision stuff ################################
    	#TO DO
        rows = 3
        cols = 3
        grid_size = (rows, cols)
        nGp = 2
        weight = 0.25
        faces = [
            [0, 1, 4, 3],       # Face 1
            [1, 2, 5, 4],       # Face 2
            [3, 4, 7, 6],       # Face 3
            [4, 5, 8, 7],       # Face 4
        ]

        control_points = []
        for i in range(0, len(nurbs_surface)):
            control_point = [nurbs_surface[i].X, nurbs_surface[i].Y, nurbs_surface[i].Z]
            control_points.append(control_point)
        print('Control points of the geometry:', control_points)
        print('\n')

        # print(len(nurbs_surface))
        # print(nurbs_surface)

        refined_control_points, refined_faces = catmull_clark_subdivision(control_points, faces)

        refined_rows = 5
        refined_cols = 5
        refined_grid_size = (refined_rows, refined_cols)

        #order control points in faces
        refined_faces = order_control_points_in_faces(refined_control_points, refined_faces)

        print('Refined Control Points:\n',refined_control_points)
        print('\n')
        print('Refined Faces:\n',refined_faces)

        # Recreate nodes in model part to ensure correct assignment of dofs
        node_id = (self.model.GetModelPart("IgaModelPart").Nodes)[len(self.model.GetModelPart("IgaModelPart").Nodes)].Id

        nodes = KratosMultiphysics.NodesVector()
        for i, vertex in enumerate(refined_control_points):
            # print(node_id+i+1)
            node = self.model.GetModelPart("IgaModelPart").CreateNewNode(node_id+i+1, vertex[0], vertex[1], vertex[2])
            nodes.append(node)

        KratosMultiphysics.CreateQuadraturePointsUtility.SetInternals(surface, nodes)

        
        ######################
        self.model.GetModelPart("IgaModelPart").CreateSubModelPart(sub_model_part)
        self.model.GetModelPart("IgaModelPart").CreateSubModelPart(sub_model_part_cond)



        # self._ModelersSetupModelPart() #3 skip iga_modeler
        # ########################################################################################
        # ########################## Remove elements #############################################
        # self.model.GetModelPart("IgaModelPart").Elements.clear()
        # self.model.GetModelPart("IgaModelPart").GetSubModelPart(sub_model_part).Elements.clear()

        # ########################## Remove conditions #############################################
        # self.model.GetModelPart("IgaModelPart").Conditions.clear()
        # self.model.GetModelPart("IgaModelPart").GetSubModelPart(sub_model_part_cond).Conditions.clear()
        # print(self.model.GetModelPart("IgaModelPart"))
        ########################################################################################

        # Here GaussIntegrationPoints are in the parametric domain of the element:
        #GaussIntegrationPoints_parametric, GaussIntegrationPoints_physical, element_control_points = evaluate_face_integration_points(faces, control_points, rows, cols, nGp, weight)
        GaussIntegrationPoints_parametric, GaussIntegrationPoints_physical, element_control_points = evaluate_face_integration_points(refined_faces, refined_control_points, refined_rows, refined_cols, nGp, weight)
        print('Gauss Integration Points in Parametric Domain:\n',GaussIntegrationPoints_parametric)
        print('\n')
        print('Gauss Integration Points in Physical Domain:\n',GaussIntegrationPoints_physical)
        print('\n')
        print('Element Control Points:\n',element_control_points)
        print('\n')

       
        
        element_id = 1
        element_id_2 = 1
        print(len(element_control_points))
        
        for i in range(0, len(element_control_points)):
        
            print("Element number: ", element_id)
            print('\n')
            _element_control_point = element_control_points[element_id-1]
            element_control_point = KratosMultiphysics.Vector(len(_element_control_point))

            for j in range(len(_element_control_point)):
                element_control_point[j] = _element_control_point[j]
                self.model.GetModelPart("IgaModelPart").GetSubModelPart(sub_model_part).AddNode(self.model.GetModelPart("IgaModelPart").Nodes[_element_control_point[j]+node_id+1], 0 )
                self.model.GetModelPart("IgaModelPart").GetSubModelPart(sub_model_part_cond).AddNode(self.model.GetModelPart("IgaModelPart").Nodes[_element_control_point[j]+node_id+1], 0 )
            

            print("Element control point for the selected element: ", element_control_point)
            print('\n')
            ## check only one element!
            integration_points = []
            integration_points_paramentric = []
        
            for j in range(len(GaussIntegrationPoints_parametric[element_id-1])):
                integration_points_paramentric.append([GaussIntegrationPoints_parametric[element_id-1][j][0],GaussIntegrationPoints_parametric[element_id-1][j][1]
                                    ,GaussIntegrationPoints_parametric[element_id-1][j][2],GaussIntegrationPoints_parametric[element_id-1][j][3]])
                
            print("Integration points parametric: ", integration_points_paramentric)
            
            for j in range(len(GaussIntegrationPoints_physical[element_id-1])):
                integration_points.append([GaussIntegrationPoints_physical[element_id-1][j][0],GaussIntegrationPoints_physical[element_id-1][j][1]
                                    ,GaussIntegrationPoints_physical[element_id-1][j][2],GaussIntegrationPoints_physical[element_id-1][j][3]])
            print("Integration points: ", integration_points)

            index_print = 1
            # compute shape functions Subdivision
            N = []
            dN_dxi = []
            dN_deta = []
            d2N_dxi2 = []
            d2N_deta2 = []
            d2N_dxi_deta = []
            d2N_deta_dxi = []
            
            for integration_point in integration_points_paramentric:
                print(len(element_control_points[i]))
                _N, _dN_dxi, _dN_deta, _d2N_dxi2, _d2N_deta2, _d2N_dxi_deta, _d2N_deta_dxi = evaluate_shape_functions_and_derivatives(integration_point[0], integration_point[1], len(element_control_points[i]))
                print(f"Shape functions N for Point {index_print}: ", _N) #2D
                print('\n')
                print(f"Shape function derivatives w.r.t to xi for Point {index_print}: ", _dN_dxi) # Now 2D
                print('\n')
                print(f"Shape function derivatives w.r.t to eta for Point {index_print}: ", _dN_deta) # Now 2D
                print('\n')
                index_print+=1
                N.append(_N)
                dN_dxi.append(_dN_dxi)
                dN_deta.append(_dN_deta)
                d2N_dxi2.append(_d2N_dxi2)
                d2N_deta2.append(_d2N_deta2)
                d2N_dxi_deta.append(_d2N_dxi_deta) #_d2N_dxi_deta = _d2N_deta_dxi
                d2N_deta_dxi.append(_d2N_deta_dxi)
                #CreateQuadraturePoint utility
            
            print('N\n: ', N)
            print('\n')
            print('dN_dxi\n: ', dN_dxi)
            print('\n')
            print('dN_deta\n: ', dN_deta)
            print('\n')
            print('d2N_dxi_deta\n: ', d2N_dxi_deta)
            print('\n')
            print('d2N_deta_dxi\n: ', d2N_deta_dxi)
            print('\n')
            
            '''
            # Convert shape functions and their derivatives to Kratos Matrix format
            N1 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape functions N for an element (N)
            N2 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape function derivatives w.r.t to xi for an element (dN_dxi)
            N3 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape function derivatives w.r.t to eta for an element (dN_deta)

            for j in range(len(element_control_points[element_id-1])):
                N1[0,j] = N[0][j]
                N2[0,j] = dN_dxi[0][j]
                N3[0,j] = dN_deta[0][j]
            

            print('N1\n: ', N1)
            print('\n')
            print('N2\n: ', N2)
            print('\n')
            print('N3\n: ', N3)
            print('\n')
            '''

            # create quadrature_point_geometries
            for j in range(4): # Iterate over each gauss point within the element
                index = int(j)
                print(index)
                # Convert shape functions and their derivatives to Kratos Matrix format
                N1 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape functions N for an element (N)
                N2 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape function derivatives w.r.t to xi for an element (dN_dxi)
                N3 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape function derivatives w.r.t to eta for an element (dN_deta)

                N4 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape function derivatives w.r.t to xi2 for an element (d2N_dxi2)
                N5 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape function derivatives w.r.t to xieta for an element (d2N_dxieta)
                N6 = KratosMultiphysics.Matrix(1,len(element_control_points[element_id-1])) # Shape function derivatives w.r.t to eta2 for an element (d2N_deta2)

                print(i)
                print(element_control_points[element_id-1])#
                for k in range(len(element_control_points[element_id-1])):
                    N1[0,k] = N[index][k]
                    N2[0,k] = dN_dxi[index][k]
                    N3[0,k] = dN_deta[index][k]
                    N4[0,k] = d2N_dxi2[index][k]
                    N5[0,k] = d2N_dxi_deta[index][k]
                    N6[0,k] = d2N_deta2[index][k]

                print('N1\n: ', N1)
                print('\n')
                print('N2\n: ', N2)
                print('\n')
                print('N3\n: ', N3)
                print(integration_points_paramentric[index])
                print(element_control_point)                                                

                quadrature_point_geometry = KratosMultiphysics.GeometriesVector()
                quadrature_point_geometry = KratosMultiphysics.CreateQuadraturePointsUtility.CreateQuadraturePoint(surface,
                                                                                                                    3,
                                                                                                                    2,
                                                                                                                    integration_points_paramentric[index],
                                                                                                                    #index,
                                                                                                                    element_control_point, 
                                                                                                                    N1, 
                                                                                                                    N2, 
                                                                                                                    N3, N4, N5, N6) #control_points, N, dN_dxi, dN_deta
   
                self.model.GetModelPart("IgaModelPart").GetSubModelPart(sub_model_part).CreateNewElement('Shell3pElement', element_id_2, quadrature_point_geometry, shell_properties)

                self.model.GetModelPart("IgaModelPart").GetSubModelPart(sub_model_part_cond).CreateNewCondition('LoadCondition', element_id_2, quadrature_point_geometry, load_properties)

                element_id_2 += 1
            element_id += 1
        

        # Strong boundary conditions
        for i in range(0, len(sub_model_part_bc)): 
            self.model.GetModelPart("IgaModelPart").CreateSubModelPart(sub_model_part_bc[i])

            ##
            u_start = 0
            u_end = refined_rows
            v_start = 0
            v_end = refined_cols

            if sub_model_part_bc_type[i] == "GeometrySurfaceNodes":
                if sub_model_part_bc_param[i][0] >= 0:
                    u_start = int(sub_model_part_bc_param[i][0]) * (refined_rows - 1)
                    u_end = int(sub_model_part_bc_param[i][0]) * (refined_rows - 1) + 1
                if sub_model_part_bc_param[i][1] >= 0:
                    v_start = int(sub_model_part_bc_param[i][1]) * (refined_cols - 1)
                    v_end = int(sub_model_part_bc_param[i][1]) * (refined_cols - 1) + 1
            elif sub_model_part_bc_type[i] == "GeometrySurfaceVariationNodes":
                if sub_model_part_bc_param[i][0] == 0:
                    u_start = 1
                    u_end = 2
                if sub_model_part_bc_param[i][0] == 1:
                    u_start = refined_rows - 2
                    u_end = refined_rows - 1
                if sub_model_part_bc_param[i][1] == 0:
                    v_start = 1
                    v_end = 2
                if sub_model_part_bc_param[i][1] == 1:
                    v_start = refined_cols - 2
                    v_end = refined_cols - 1

            for m in range(u_start, u_end):
                for n in range(v_start, v_end):
                    # print(self.model.GetModelPart("IgaModelPart").Nodes[(m + n * refined_rows) + node_id + 1])
                    self.model.GetModelPart("IgaModelPart").GetSubModelPart(sub_model_part_bc[i]).AddNode(self.model.GetModelPart("IgaModelPart").Nodes[(m + n * refined_rows) + node_id + 1])
                    

        # print(self.model.GetModelPart("IgaModelPart"))
        # exit()

        ########################################################################################
        ########################################################################################
        

        self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        self._AnalysisStage__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()
        self.Check()

        self.ModifyAfterSolverInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time

        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

        # In case of contact problem
        if self.contact_problem:
            self._GetSolver().SetEchoLevel(self.echo_level)
            # To avoid many prints
            if self.echo_level == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """

        # In case of contact problem
        if self.contact_problem:
            # First we check if one of the output processes will print output in this step this is done to save computation in case none of them will print
            is_output_step = False
            for output_process in self._GetListOfOutputProcesses():
                if output_process.IsOutputStep():
                    is_output_step = True
                    break

            if is_output_step:
                # Informing the output will be created
                KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
                KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "TIME: ", self.time)

        # Creating output
        super().OutputSolutionStep()

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return structural_solvers.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "json_output_process",
                "json_check_process", "check_analytic_results_process", "contact_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        info_msg  = "Using the old way to create the processes, this was removed!\n"
                        info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                        info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                        info_msg += "for a description of the new format"
                        raise Exception("StructuralMechanicsAnalysis: " + info_msg)

            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this was removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                raise Exception("StructuralMechanicsAnalysis: " + info_msg)
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _GetSimulationName(self):
        return "::[KSM Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
