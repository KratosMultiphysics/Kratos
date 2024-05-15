# Import Python modules
import math

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.MPMApplication as KratosMechanics
from KratosMultiphysics.MPMApplication import python_solvers_wrapper_mpm
from KratosMultiphysics.MPMApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

have_external_solvers = KratosUtilities.CheckIfApplicationsAvailable("LinearSolversApplication")

@KratosUnittest.skipUnless(have_external_solvers, "Missing required application: LinearSolversApplication")
class ManufacturedSolutionTestMPM(KratosUnittest.TestCase):
    def testManufacturedSolutionMPM(self):
        self.runTest()

    def setUp(self):
        self.print_output   = False
        self.print_type     = "Vtk" # Available fields: "GiD" "Vtk" "Both"
        self.print_convergence_plot = False
        self.work_folder    = "manufactured_solution_test_mpm"
        self.settings       = "ManufacturedSolutionTestParametersMPM.json"
        self.meshes_list    = [ "manufactured_solution_ref0",
                                "manufactured_solution_ref1",
                                "manufactured_solution_ref2",
                                "manufactured_solution_ref3"];
        self.stabilization_list = ["asgs"]; # Available list: "ppp","osgs","asgs"


    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            for filename in self.meshes_list:
                KratosUtilities.DeleteFileIfExisting(filename + '.time')

    def runTest(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.settings, 'r') as parameter_file:
                self.OriginalProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())


            with open('Results_convergence_test.txt', 'w') as f:
                f.write('Results convergence test \n')

                ## Solve the problem for each one of the stabilizations
                for stab_name in  self.stabilization_list:
                    h = []
                    err_p = []
                    err_v = []

                    mesh_0_characteristic_size = 0.2
                    den = 1

                    ## Solve the manufactured solution problem for each one of the refinements
                    for mesh_name in self.meshes_list:
                        ## Solve the problem imposing the previously obtained values
                        CaseProjectParameters = self.OriginalProjectParameters.Clone()
                        MPMProblem = ManufacturedSolutionProblem(CaseProjectParameters, mesh_name, stab_name, self.print_output, self.print_type)
                        MPMProblem.Run()

                        ## Compute the obtained solution error
                        h.append(mesh_0_characteristic_size/den)
                        err_p.append(MPMProblem.ComputePressureErrorNorm())
                        err_v.append(MPMProblem.ComputeDisplacementErrorNorm())
                        den *= 2

                        print('----------------------------------------------------')
                        print('Results for the stabilization', stab_name, '\n')
                        print('----------------------------------------------------')
                        print('Element size (h): ', h)
                        print('L2 displacement error: ', err_v)
                        print('L2 pressure error: ', err_p)

                    f.write('\n')
                    f.write('--------------------------------------------------- \n')
                    f.write('Results for the stabilization ' + str(stab_name) + '\n')
                    f.write('--------------------------------------------------- \n')
                    f.write('Element size (h): ' + str(h) + '\n')
                    f.write('L2 displacement error: ' + str(err_v) + '\n')
                    f.write('L2 pressure error: ' + str(err_p) + '\n')

                    if stab_name=="ppp":
                        err_v_ppp=err_v
                        err_p_ppp=err_p
                    elif stab_name=="asgs":
                        err_v_asgs=err_v
                        err_p_asgs=err_p
                    elif stab_name=="osgs":
                        err_v_osgs=err_v
                        err_p_osgs=err_p

            f.close()

            ## Convergence plot print
            if (self.print_convergence_plot == True):
                ## Plot the convergence graphs
                import matplotlib.pyplot as plt
                h_1_v = []
                h_2_v = []
                h_3_v = []
                h_1_p = []
                h_2_p = []
                h_3_p = []

                den =1
                for i in range(0,len(err_p)):
                    h_1_v.append(err_v[0]/den)
                    h_2_v.append(err_v[0]/den**2)
                    h_3_v.append(err_v[0]/den**3)
                    h_1_p.append(err_p[0]/den)
                    h_2_p.append(err_p[0]/den**2)
                    h_3_p.append(err_p[0]/den**3)
                    den *= 2

                plt.figure(1)
                plt.rc('text', usetex=True)
                plt.rc('font', family='serif')
                if "ppp" in self.stabilization_list:
                    plt.loglog(h, err_v_ppp, '-x', color = 'r', label = 'ppp')
                if "asgs" in self.stabilization_list:
                    plt.loglog(h, err_v_asgs, '-x', color = 'b', label = 'asgs')
                if "osgs" in self.stabilization_list:
                    plt.loglog(h, err_v_osgs, '-o', color = 'g', label = 'osgs')

                plt.loglog(h, h_1_v, '--', color = 'k', label = 'Slope 1')
                plt.loglog(h, h_2_v, ':' , color = 'k', label = 'Slope 2')

                plt.ylabel(r'$L^2$ displacement error')
                plt.xlabel('Element size (h)')
                plt.legend(loc=4, ncol=1)
                plt.tight_layout()
                plt.show
                plt.savefig('l2_norm_convergence_disp.png')

                plt.figure(2)

                # Plot pressure
                plt.rc('text', usetex=True)
                plt.rc('font', family='serif')
                if "ppp" in self.stabilization_list:
                    plt.loglog(h, err_p_ppp, '-x', color = 'r', label = 'ppp')
                if "asgs" in self.stabilization_list:
                    plt.loglog(h, err_p_asgs, '-x', color = 'b', label = 'asgs')
                if "osgs" in self.stabilization_list:
                    plt.loglog(h, err_p_osgs, '-o', color = 'g', label = 'osgs')

                plt.loglog(h, h_1_p, '--', color = 'k', label = 'Slope 1')
                plt.loglog(h, h_2_p, ':' , color = 'k', label = 'Slope 2')

                plt.ylabel(r'$L^2$ pressure error')
                plt.xlabel('Element size (h)')
                plt.legend(loc=4, ncol=1)
                plt.tight_layout()
                plt.show
                plt.savefig('l2_norm_convergence_pressure.png')


                # Check obtained solution
            expected_displacement_errors =  [0.13811251114366857, 0.10989112694618389, 0.04668841129297075, 0.014368072149072656]
            expected_pressure_errors = [0.8966431547459842, 0.5662128821749041, 0.2215006730602814, 0.06640894005793938]

            for i in range(len(self.meshes_list)):
                self.assertAlmostEqual(err_v[i], expected_displacement_errors[i])
                self.assertAlmostEqual(err_p[i], expected_pressure_errors[i])


class ManufacturedSolutionProblem(ParticleMechanicsAnalysis):

    def __init__(self, project_parameters, input_file_name, stab_name, print_output, print_type):
        self.print_output = print_output
        self.print_type = print_type
        self.input_file_name = input_file_name
        self.stab_name = stab_name
        self.project_parameters = project_parameters
        self.model = KratosMultiphysics.Model()

        ## Set the current mesh case problem info
        self.project_parameters["problem_data"]["problem_name"].SetString(self.input_file_name)
        self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(self.input_file_name+"_Body")
        self.project_parameters["solver_settings"]["grid_model_import_settings"]["input_filename"].SetString(self.input_file_name+"_Grid")
        ## Set the stabilization employed
        self.project_parameters["solver_settings"]["stabilization"].SetString(self.stab_name)

        ## If required, set up the GiD I/O
        if self.print_output:
            if (self.print_type == "GiD"):
                self._AddGiDOutput()
            elif (self.print_type == "Vtk"):
                self._AddVtkOutput()

        ## Note that the base particle mechanics analysis constructor is called after the creation of the model and parameters customization
        super().__init__(self.model, self.project_parameters)

    def ModifyAfterSolverInitialize(self):
        super().ModifyAfterSolverInitialize()

        ### Calculate NODAL_AREA
        nodal_area_process = KratosMultiphysics.CalculateNonHistoricalNodalAreaProcess(self._GetSolver().GetComputingModelPart())
        nodal_area_process.Execute()

        #### Fix the pressure in one node (bottom left corner)
        for node in self._GetGridModelPart().Nodes:
            if ((node.X<0.0001) and (node.Y<0.0001)):
                analytical_press = self._ComputeNodalPressureManufacturedSolution(node)
                node.Fix(KratosMultiphysics.PRESSURE)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, analytical_press)

        ### Initialize the buffer with the analytical solution
        for i_buff in range(self._GetSolver().GetComputingModelPart().GetBufferSize()):
            self._SetManufacturedSolutionValues(buffer_position = i_buff, fix = False)

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        ### Apply manufactured solution BCs
        # Set the manufactured solution source terms
        self._SetManufacturedSolutionValues(fix=True, set_only_boundaries=True)
        self._SetManufacturedSolutionSourceValues()

    ### We enhance the class with these two methods to calculate the error norms
    def ComputeDisplacementErrorNorm(self):
        err_v = 0
        deno = 0
        model_part = self._GetGridModelPart()
        for node in model_part.Nodes:
            weight = node.GetValue(KratosMultiphysics.NODAL_AREA)
            disp_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
            disp_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)

            analytical_disp = self._ComputeNodalDisplacementManufacturedSolution(node)

            err_x = abs(analytical_disp[0] - disp_x)
            err_y = abs(analytical_disp[1] - disp_y)
            err_node = err_x**2 + err_y**2
            err_v += weight*err_node
            deno  += weight*(abs(analytical_disp[0])**2 + abs(analytical_disp[1])**2)

        denom      = math.sqrt(deno)
        error_disp = math.sqrt(err_v)/denom

        return error_disp # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    def ComputePressureErrorNorm(self):
        err_p = 0
        deno = 0
        model_part = self._GetGridModelPart()
        for node in model_part.Nodes:
            weight = node.GetValue(KratosMultiphysics.NODAL_AREA)
            pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
            analytical_pres = self._ComputeNodalPressureManufacturedSolution(node)
            err_p += weight*(abs(analytical_pres - pres)**2)
            deno  += weight*(abs(analytical_pres)**2)

        denom =  math.sqrt(deno)
        error_press =  math.sqrt(err_p)/denom

        return error_press # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    ### Internal methods required for the manufactured solution calculation

    def _SetManufacturedSolutionValues(self, buffer_position = 0, fix=True, set_only_boundaries=False):
        ### Set the analytical solution for the manufactured solution computation
        time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        root_model_part = self._GetGridModelPart()

        for node in root_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto1").Nodes:
            disp = self._ComputeNodalDisplacementManufacturedSolution(node)
            press = self._ComputeNodalPressureManufacturedSolution(node)

            if fix:

                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.PRESSURE)

            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, buffer_position, disp[0])
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, buffer_position, disp[1])
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, buffer_position, press)

    def _GetGridModelPart(self):
        if not self.model.HasModelPart("Background_Grid"):
            raise Exception("The GridModelPart was not created yet!")
        return self.model.GetModelPart("Background_Grid")

    def _GetBodyModelPart(self):
        if not self.model.HasModelPart("Initial_MPM_Material.Parts_Material_domain_Material"):
            raise Exception("The MPM_Material was not created yet!")
        return self.model.GetModelPart("Initial_MPM_Material.Parts_Material_domain_Material")

    def _SetManufacturedSolutionSourceValues(self):
        ### Set the body force as source term
        model_part = self._GetGridModelPart()

        for node in model_part.Nodes:
            f_exac = self._ComputeNodalSourceTermManufacturedSolution(node)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, f_exac[0])  # Set the x-component body force field
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, f_exac[1])  # Set the y-component body force field

            disp = self._ComputeNodalDisplacementManufacturedSolution(node)
            press = self._ComputeNodalPressureManufacturedSolution(node)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, disp[0])
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, disp[1])
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, press)

    def _ComputeNodalSourceTermManufacturedSolution(self, node):
        x=node.X
        y=node.Y
        pi=math.pi
        k=1e-2

        E = self._GetBodyModelPart().Properties[1].GetValue(KratosMultiphysics.YOUNG_MODULUS)
        nu = self._GetBodyModelPart().Properties[1].GetValue(KratosMultiphysics.POISSON_RATIO)
        mu=E/(2*(1+nu))

        Deviatoric_x = -(2*k*mu*math.exp(x + y)*(4*k*x**2*math.exp(x + y) + 4*k*y**2*math.exp(x + y) + 8*k*x*math.exp(x + y) + 8*k*y*math.exp(x + y) + 8*k*x*y*math.exp(x + y) - 3)*(x**2 + 2*x*y + 4*x + y**2 + 4*y + 2))/3
        Deviatoric_y = -(2*k*mu*math.exp(x + y)*(4*k*x**2*math.exp(x + y) + 4*k*y**2*math.exp(x + y) + 8*k*x*math.exp(x + y) + 8*k*y*math.exp(x + y) + 8*k*x*y*math.exp(x + y) + 3)*(x**2 + 2*x*y + 4*x + y**2 + 4*y + 2))/3

        Volumetric_x = 2*k*mu*pi*math.exp(x + y)*math.cos(2*pi*y)*math.sin(2*pi*x)*(x**2 + 2*x*y + 2*x + y**2 + 2*y) - 2*mu*pi*math.cos(2*pi*x)*math.sin(2*pi*y)*(k*math.exp(x + y)*(x + y)**2 + k*math.exp(x + y)*(2*x + 2*y) - 1)
        Volumetric_y = 2*mu*pi*math.cos(2*pi*y)*math.sin(2*pi*x)*(k*math.exp(x + y)*(x + y)**2 + k*math.exp(x + y)*(2*x + 2*y) + 1) - 2*k*mu*pi*math.exp(x + y)*math.cos(2*pi*x)*math.sin(2*pi*y)*(x**2 + 2*x*y + 2*x + y**2 + 2*y)

        fx =  - Deviatoric_x - Volumetric_x
        fy =  - Deviatoric_y - Volumetric_y

        return [fx, fy]

    def _ComputeNodalDisplacementManufacturedSolution(self, node):
        k=1e-2
        ux = k*((node.X+node.Y)**2)*math.exp(node.X+node.Y)
        uy = -k*((node.X+node.Y)**2)*math.exp(node.X+node.Y)

        return [ux, uy]

    def _ComputeNodalPressureManufacturedSolution(self, node):
        E = self._GetBodyModelPart().Properties[1].GetValue(KratosMultiphysics.YOUNG_MODULUS)
        nu = self._GetBodyModelPart().Properties[1].GetValue(KratosMultiphysics.POISSON_RATIO)
        mu=E/(2*(1+nu))
        press = mu*math.sin(2*math.pi*node.X)*math.sin(2*math.pi*node.Y)

        return press

    def _AddGiDOutput(self):
        gid_output_settings_Body = KratosMultiphysics.Parameters("""{
            "python_module" : "particle_gid_output_process",
            "kratos_module" : "KratosMultiphysics.MPMApplication",
            "process_name"  : "ParticleMPMGiDOutputProcess",
            "Parameters"    : {
                "model_part_name"        : "MPM_Material",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 1,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "gauss_point_results"         : ["MP_VELOCITY","MP_DISPLACEMENT"],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                },
                "output_name": []
            }
        }""")
        gid_output_settings_Grid = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "Parameters"    : {
                "model_part_name"        : "Background_Grid",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 1,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["DISPLACEMENT","PRESSURE"],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                },
                "output_name": []
            }
        }""")
        gid_output_settings_Body["Parameters"]["output_name"].SetString(self.input_file_name+"_Body")
        gid_output_settings_Grid["Parameters"]["output_name"].SetString(self.input_file_name+"_Grid")
        self.project_parameters["output_processes"]["gid_output_processes"].Append(gid_output_settings_Body)
        self.project_parameters["output_processes"]["gid_output_processes"].Append(gid_output_settings_Grid)



    def _AddVtkOutput(self):
        vtk_output_settings_Body=KratosMultiphysics.Parameters("""{
            "python_module" : "particle_vtk_output_process",
            "kratos_module" : "KratosMultiphysics.MPMApplication",
            "process_name"  : "ParticleMPMVTKOutputProcess",
            "Parameters"    : {
                "model_part_name"             : "MPM_Material",
                "output_control_type"         : "step",
                "output_interval"             : 1,
                "file_format"                 : "ascii",
                "output_precision"            : 7,
                "folder_name"                 : [],
                "output_sub_model_parts"      : false,
                "save_output_files_in_folder" : true,
                "gauss_point_results"         : ["MP_VELOCITY","MP_DISPLACEMENT","MP_PRESSURE"]
            }
        }""")
        vtk_output_settings_Grid=KratosMultiphysics.Parameters("""{
            "python_module" : "vtk_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "VtkOutputProcess",
            "Parameters"    : {
                "model_part_name"                    : "Background_Grid",
                "output_control_type"                : "step",
                "output_interval"                    : 1,
                "file_format"                        : "ascii",
                "output_precision"                   : 7,
                "output_sub_model_parts"             : false,
                "save_output_files_in_folder"        : true,
                "output_path"                        : [],
                "nodal_solution_step_data_variables" : ["DISPLACEMENT","PRESSURE","BODY_FORCE"]
            }
        }""")
        vtk_output_settings_Body["Parameters"]["folder_name"].SetString("VtkPost_"+self.input_file_name)
        vtk_output_settings_Grid["Parameters"]["output_path"].SetString("VtkPost_"+self.input_file_name)
        self.project_parameters["output_processes"]["vtk_output_processes"].Append(vtk_output_settings_Grid)
        self.project_parameters["output_processes"]["vtk_output_processes"].Append(vtk_output_settings_Body)


if __name__ == '__main__':
    KratosUnittest.main()

