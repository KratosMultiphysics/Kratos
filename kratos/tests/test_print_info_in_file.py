import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM
import KratosMultiphysics.print_info_in_file_process as PrintProcess
import KratosMultiphysics.kratos_utilities as KratosUtils
import os

if KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication as SMA

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestPrintNodalInfoInFile(KratosUnittest.TestCase):
    """ In this text we check that nodal values are printed
    correctly, we check DISPLACEMENT and REACTION
    """
    def test_nodal_print(self):
        current_path = os.path.abspath(os.getcwd())

        current_model = KM.Model()
        model_part = current_model.CreateModelPart("test")
        properties = KM.Properties(1)
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        model_part.ProcessInfo = KM.ProcessInfo()
        model_part.ProcessInfo[KM.STEP] = 1

        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        nodes = [node1, node2]

        geom = KM.Tetrahedra3D4(node1,node2,node3,node4)

        submodel_1 = model_part.CreateSubModelPart("to_plot")
        submodel_1.AddNodes([1,2])

        submodel_2 = model_part.CreateSubModelPart("to_plot_2")
        submodel_2.AddNodes([3,4])

        node1.SetSolutionStepValue(KM.DISPLACEMENT, [1.1,0.0,0.0])
        node2.SetSolutionStepValue(KM.DISPLACEMENT, [2.6,0.0,0.0])
        node3.SetSolutionStepValue(KM.DISPLACEMENT, [3.0,0.1,0.0])
        node4.SetSolutionStepValue(KM.DISPLACEMENT, [4.3,0.0,0.3])

        settings = KM.Parameters("""{
            "model_part_name"                    : "test.to_plot",
            "variable_name"                      : "DISPLACEMENT",
            "results_type"                       : "nodal_historical",
            "file_name"                          : "info_file.txt",
            "output_control_type"                : "step",
            "erase_previous_info"                : true,
            "output_interval"                    : 1,
            "sum_results_from_multiple_entities" : true,
            "write_buffer_size"                  : 1,
            "output_path"                        : ""}""")

        process = PrintProcess.PrintInfoInFileProcess(current_model, settings)
        process.PrintOutput()

        ref_file_name = os.path.abspath(settings["file_name"].GetString())

        if not os.path.isfile(ref_file_name):
            err_msg  = 'The specified reference file name "'
            err_msg += ref_file_name
            err_msg += '" is not valid!'
            raise Exception(err_msg)

        with open(ref_file_name, "r") as plot_file:
            contents = plot_file.readlines()
            for line in contents:
                if line[0] != "#":
                    numbers = line.split("\t")
                    if len(numbers) > 2:
                        if float(numbers[0]) != 0.0 or float(numbers[1]) != 3.7 or float(numbers[2]) != 0.0 or float(numbers[3]) != 0.0:
                            raise Exception("The print does not give the expected result...")

        # now we test the second submodel
        settings["model_part_name"].SetString("test.to_plot_2")
        process = PrintProcess.PrintInfoInFileProcess(current_model, settings)
        process.PrintOutput()

        with open(ref_file_name, "r") as plot_file:
            contents = plot_file.readlines()
            for line in contents:
                if line[0] != "#":
                    numbers = line.split("\t")
                    if len(numbers) > 2:
                        if float(numbers[0]) != 0.0 or float(numbers[1]) != 7.3 or float(numbers[2]) != 0.1 or float(numbers[3]) != 0.3:
                            raise Exception("The print does not give the expected result...")
        process.ExecuteFinalize()
        os.remove(ref_file_name)


class TestPrintElementalInfoInFile(KratosUnittest.TestCase):
    """ In this text we check that elemental values are printed
    correctly, we check CAUCHY_STRESS_VECTOR and VON_MISES_STRESS
    """
    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def test_elemental_values(self):
        dim = 2
        current_model = KM.Model()
        self.mp = current_model.CreateModelPart("MainModelPart")
        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = dim

        self.mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.mp.AddNodalSolutionStepVariable(KM.REACTION)
        self.mp.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)

        #define properties
        self.mp.GetProperties()[1].SetValue(KM.YOUNG_MODULUS, 210e9)
        self.mp.GetProperties()[1].SetValue(KM.POISSON_RATIO, 0.3)
        self.mp.GetProperties()[1].SetValue(KM.THICKNESS, 1.0)
        self.mp.GetProperties()[1].SetValue(KM.DENSITY, 1.0)

        # create nodes
        self.mp.CreateNewNode(1, 0.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(2, 1.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(3, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(4, 0.00000, 1.00000, 0.00000)

        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X,self.mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y,self.mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z,self.mp)

        #create a submodelpart for boundary conditions
        bcs = self.mp.CreateSubModelPart("FixedEdgeNodes")
        bcs.AddNodes([1, 4])

        bcmn = self.mp.CreateSubModelPart("MovingNodes")
        bcmn.AddNodes([2,3])

        g = [0,0,0]
        self.mp.GetProperties()[1].SetValue(KM.VOLUME_ACCELERATION,g)

        cl = SMA.LinearElasticPlaneStrain2DLaw()
        self.mp.GetProperties()[1].SetValue(KM.CONSTITUTIVE_LAW, cl)

        #create Element
        elem = self.mp.CreateNewElement("SmallDisplacementElement2D4N", 1, [1, 2, 3, 4], self.mp.GetProperties()[1])
        elem_mp = self.mp.CreateSubModelPart("Element")
        elem_mp.AddElements([1])
        elem.Initialize(self.mp.ProcessInfo)

        # We fix the displ in some nodes
        bcs = self.mp.GetSubModelPart("FixedEdgeNodes")
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_X, 0.0, bcs.Nodes)
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_Y, 0.0, bcs.Nodes)
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_Z, 0.0, bcs.Nodes)

        # KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, bcs.Nodes)
        # KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, bcs.Nodes)

        bcmn = self.mp.GetSubModelPart("MovingNodes")
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_X, 0.01, bcmn.Nodes)
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_Y, 0.0, bcmn.Nodes)
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_Z, 0.0, bcmn.Nodes)

        settings = KM.Parameters("""{
            "model_part_name"                    : "MainModelPart.Element",
            "variable_name"                      : "VON_MISES_STRESS",
            "results_type"                       : "elemental",
            "file_name"                          : "elemental_von.txt",
            "output_control_type"                : "step",
            "erase_previous_info"                : true,
            "output_interval"                    : 1,
            "sum_results_from_multiple_entities" : false,
            "write_buffer_size"                  : 1,
            "integration_point_number" : 0,
            "output_path"                        : ""
            }
        """)

        process = PrintProcess.PrintInfoInFileProcess(current_model, settings)
        process.PrintOutput()

        ref_file_name = os.path.abspath(settings["file_name"].GetString())

        if not os.path.isfile(ref_file_name):
            err_msg  = 'The specified reference file name "'
            err_msg += ref_file_name
            err_msg += '" is not valid!'
            raise Exception(err_msg)

        with open(ref_file_name, "r") as plot_file:
            contents = plot_file.readlines()
            for line in contents:
                if line[0] != "#":
                    numbers = line.split("\t")
                    if len(numbers) > 2:
                        if float(numbers[0]) != 0.0 or float(numbers[1]) != 2.4565e+09:
                            raise Exception("The print does not give the expected result...")
        process.ExecuteFinalize()
        os.remove(ref_file_name)


        settings = KM.Parameters("""{
            "model_part_name"                    : "MainModelPart.Element",
            "variable_name"                      : "GREEN_LAGRANGE_STRAIN_VECTOR",
            "results_type"                       : "elemental",
            "file_name"                          : "elemental_strain.txt",
            "output_control_type"                : "step",
            "erase_previous_info"                : true,
            "output_interval"                    : 1,
            "sum_results_from_multiple_entities" : false,
            "write_buffer_size"                  : 1,
            "integration_point_number" : 0,
            "output_path"                        : ""
            }
        """)

        process = PrintProcess.PrintInfoInFileProcess(current_model, settings)
        process.PrintOutput()

        ref_file_name = os.path.abspath(settings["file_name"].GetString())

        if not os.path.isfile(ref_file_name):
            err_msg  = 'The specified reference file name "'
            err_msg += ref_file_name
            err_msg += '" is not valid!'
            raise Exception(err_msg)

        with open(ref_file_name, "r") as plot_file:
            contents = plot_file.readlines()
            for line in contents:
                if line[0] != "#":
                    numbers = line.split("\t")
                    print(numbers)
                    if len(numbers) > 2:
                        if float(numbers[0]) != 0.0 or float(numbers[1]) != 1.0E-2 or float(numbers[2]) != 0.0 or float(numbers[3]) != 0.0:
                            raise Exception("The print does not give the expected result...")
        process.ExecuteFinalize()
        os.remove(ref_file_name)


if __name__ == '__main__':
    KratosUnittest.main()
