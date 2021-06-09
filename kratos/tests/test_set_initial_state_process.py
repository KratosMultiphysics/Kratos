import KratosMultiphysics as km
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.set_initial_state_process as set_initial_state_process


class TestSetInitialStateProcess(KratosUnittest.TestCase):
    def test_set_initial_state_process(self):
        #ExecuteBasicVTKoutputProcessCheck()
        km.Logger.GetDefaultOutput().SetSeverity(km.Logger.Severity.WARNING)
        current_model = km.Model()
        mp = current_model.CreateModelPart("Main")
        SetupModelPart3D(mp)

        parameters = km.Parameters(
            """{
            "Parameters" : {
                "model_part_name": "Main",
                "imposed_strain_multiplier": "0.0 * t",
                "imposed_strain": [1,0,0,0,0,0],
                "imposed_stress_multiplier": "0.0 * t",
                "imposed_stress": [1,0,0,0,0,0],
                "imposed_deformation_gradient": [[1,0,0], [0,1,0], [0,0,1]]
            }
        }"""
        )

        #process =  vtk_output_process.Factory(parameters, current_model)
        process =  set_initial_state_process.Factory(parameters, current_model)

        process.ExecuteInitialize()
        time = 0.0
        step = 0
        while (time <= 1.0):
            time += 0.2
            step += 1
            mp.ProcessInfo[km.STEP] += 1
            SetupSolution(mp)
            #process.ExecuteInitializeSolutionStep()
            mp.CloneTimeStep(time)
            process.ExecuteFinalizeSolutionStep()



def SetupModelPart3D(mp):
    mp.AddNodalSolutionStepVariable(km.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(km.VELOCITY)
    mp.AddNodalSolutionStepVariable(km.PRESSURE)

    # Create nodes
    mp.CreateNewNode(1, 0.0 , 1.0 , 1.0)
    mp.CreateNewNode(2, 0.0 , 1.0 , 0.0)
    mp.CreateNewNode(3, 0.0 , 0.0 , 1.0)
    mp.CreateNewNode(4, 1.0 , 1.0 , 1.0)
    mp.CreateNewNode(5, 0.0 , 0.0 , 0.0)
    mp.CreateNewNode(6, 1.0 , 1.0 , 0.0)
    mp.CreateNewNode(7, 1.0 , 0.0 , 1.0)
    mp.CreateNewNode(8, 1.0 , 0.0 , 0.0)
    mp.CreateNewNode(9, 2.0 , 1.0 , 1.0)
    mp.CreateNewNode(10, 2.0 , 1.0 , 0.0)
    mp.CreateNewNode(11, 2.0 , 0.0 , 1.0)
    mp.CreateNewNode(12, 2.0 , 0.0 , 0.0)
    mp.CreateNewNode(13, 0.0 , 0.0 , 2.0)
    mp.CreateNewNode(14, 1.0 , 0.0 , 2.0)
    mp.CreateNewNode(15, 1.0 , 1.0 , 2.0)

    # Create elements
    mp.CreateNewElement("Element3D4N", 1, [12, 10, 8, 9], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 2, [4, 6, 9, 7], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 3, [11, 7, 9, 8], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 4, [5, 3, 8, 6], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 5, [4, 6, 7, 3], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 6, [2, 3, 5, 6], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 7, [10, 9, 6, 8], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 8, [7, 8, 3, 6], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 9, [7, 8, 6, 9], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 10, [4, 1, 6, 3], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 11, [9, 12, 11, 8], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D4N", 12, [3, 2, 1, 6], mp.GetProperties()[1])
    mp.CreateNewElement("Element3D6N", 13, [3, 7, 4, 13, 14, 15], mp.GetProperties()[1])

    # Create a submodelpart for boundary conditions
    bcs = mp.CreateSubModelPart("FixedEdgeNodes")
    bcs.AddNodes([1, 2, 5])

    bcmn = mp.CreateSubModelPart("MovingNodes")
    bcmn.AddNodes([13, 14, 15])


def SetupSolution(mp):
    time = mp.ProcessInfo[km.TIME] + 0.158
    step = mp.ProcessInfo[km.STEP]

    for node in mp.Nodes:
        node.SetSolutionStepValue(km.DISPLACEMENT,0,[node.X*time,node.Y,node.Z*step])
        node.SetSolutionStepValue(km.VELOCITY,0,[2*node.X,2*node.Y,2*node.Z])
        node.SetSolutionStepValue(km.PRESSURE,0,node.X*time*step)

    for i_elem, elem in enumerate(mp.Elements):
        elem.SetValue(km.DETERMINANT, [i_elem*0.189,time,time*step])

if __name__ == '__main__':
    KratosUnittest.main()
