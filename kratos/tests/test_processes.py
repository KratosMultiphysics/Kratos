from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *
import math

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName


class TestProcesses(KratosUnittest.TestCase):

    def test_assign_processes(self):
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        #reset all data
        for node in model_part.Nodes:
            node.Free(DISPLACEMENT_X)
            node.Free(DISPLACEMENT_Y)
            node.Free(DISPLACEMENT_Z)
            node.Free(VELOCITY_X)
            node.Free(VELOCITY_Y)
            node.Free(VELOCITY_Z)
            node.SetSolutionStepValue(DENSITY,0,0.0)
            node.SetSolutionStepValue(VISCOSITY,0,0.0)
            node.SetSolutionStepValue(DISPLACEMENT_X,0,0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y,0,0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Z,0,0.0)
            node.SetSolutionStepValue(VELOCITY_X,0,0.0)
            node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
            node.SetSolutionStepValue(VELOCITY_Z,0,0.0)

        settings = Parameters(
            """
            {
                "process_list" : [
                    {
                        "python_module"   : "assign_scalar_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableProcess",
                        "Parameters"            : {
                            "model_part_name" : "Main",
                            "variable_name"   : "VISCOSITY",
                            "interval"        : [0.0, 10.0],
                            "constrained"		  : true,
                            "value"      : "x+100.0*y*t**2"
                        }
                    },
                    {
                        "python_module"   : "assign_scalar_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableProcess",
                        "Parameters"            : {
                            "model_part_name" : "Main",
                            "variable_name"   : "DENSITY",
                            "value"      : "x*x+y*y+z*z+t"
                        }
                    },
                    {
                        "python_module"   : "assign_scalar_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableProcess",
                        "Parameters"            : {
                            "model_part_name" : "Main",
                            "variable_name"   : "DISPLACEMENT_X",
                            "interval"        : [0.0, 5.0],
                            "constrained"		  : true,
                            "value"      : "sqrt(x**2+y**2)*t",
                            "local_axes"               :{
                                "origin" : [0.0, 0.0, 0.0],
                                "axes"  : [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
                            }
                        }
                    },
                    {
                        "python_module"   : "assign_vector_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignVectorVariableProcess",
                        "Parameters"            : {
                                "model_part_name"      : "Main",
                                "variable_name"        : "DISPLACEMENT",
                                "interval"             : [11.0, 15.0],
                                "value"                : [10.0, null, "t"],
                                "local_axes"           : {}
                            }
                    },
                    {
                        "python_module"   : "assign_vector_by_direction_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignVectorByDirectionProcess",
                        "Parameters"            : {
                                "model_part_name"      : "Main",
                                "variable_name"        : "VELOCITY",
                                "interval"             : [11.0, 15.0],
                                "modulus"              : 10.0,
                                "constrained"          : false,
                                "direction"            : [1.0, 0.0, 0.0],
                                "local_axes"           : {}
                            }
                    },
                    {
                        "python_module"   : "assign_vector_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignVectorVariableProcess",
                        "Parameters"            : {
                                "model_part_name"      : "Main",
                                "variable_name"        : "DISPLACEMENT",
                                "interval"             : [20.0, 24.0],
                                "constrained"          : false,
                                "value"                : [10.0, null, "t"],
                                "local_axes"           : {}
                            }
                    },
                    {
                        "python_module"   : "assign_vector_by_direction_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignVectorByDirectionProcess",
                        "Parameters"            : {
                                "model_part_name"      : "Main",
                                "variable_name"        : "VELOCITY",
                                "interval"             : [20.0, 24.0],
                                "modulus"              : "sin(x*pi*t)",
                                "constrained"          : false,
                                "direction"            : [0.0, 1.0, 0.0],
                                "local_axes"           : {}
                            }
                    },
                    {
                        "python_module"   : "assign_vector_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignVectorProcess",
                        "Parameters"            : {
                                "model_part_name"      : "Main",
                                "variable_name"        : "DISPLACEMENT",
                                "interval"             : [25.0, "End"],
                                "constrained"          : [true,true,false],
                                "value"                : [null, "x+y*t", "t"],
                                "local_axes"           : {}
                            }
                    },
                    {
                        "python_module"   : "assign_vector_by_direction_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignVectorByDirectionProcess",
                        "Parameters"            : {
                                "model_part_name"      : "Main",
                                "variable_name"        : "VELOCITY",
                                "interval"             : [25.0, "End"],
                                "modulus"              : "sqrt(abs(x*y))",
                                "constrained"          : true,
                                "direction"            : [0.0, 1.0, 1.0],
                                "local_axes"           : {}
                            }
                    }
                ]
                }
            """
            )

        Model = {"Main":model_part}

        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( settings["process_list"] )

        for node in model_part.Nodes:
            self.assertFalse(node.IsFixed(DISPLACEMENT_X))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Y))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Z))

        ############################################################
        ##time = 3 - both within the active interval
        model_part.CloneTimeStep(3.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), math.sqrt(node.X**2+node.Y**2)*t)
            self.assertEqual(node.GetSolutionStepValue(DENSITY), node.X**2+node.Y**2+node.Z**2+t)
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), node.X+100.0*node.Y*t**2)
            self.assertTrue(node.IsFixed(DENSITY))
            self.assertTrue(node.IsFixed(VISCOSITY))
            self.assertTrue(node.IsFixed(DISPLACEMENT_X))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Y))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Z))

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        for node in model_part.Nodes:
            self.assertFalse(node.IsFixed(DENSITY))
            self.assertFalse(node.IsFixed(VISCOSITY))
            self.assertFalse(node.IsFixed(DISPLACEMENT_X))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Y))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Z))

        ############################################################
        ##time = 3 - DISPLACEMENT_X is not in the active interval
        model_part.CloneTimeStep(6.0)

        for node in model_part.Nodes:
            self.assertFalse(node.IsFixed(DISPLACEMENT_X))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Y))
            self.assertFalse(node.IsFixed(DISPLACEMENT_Z))

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), math.sqrt(node.X**2+node.Y**2)*3.0) ##still the old value
            self.assertEqual(node.GetSolutionStepValue(DENSITY), node.X**2+node.Y**2+node.Z**2+t)
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), node.X+100.0*node.Y*t**2)
            self.assertTrue(node.IsFixed(DENSITY))
            self.assertTrue(node.IsFixed(VISCOSITY))
            self.assertFalse(node.IsFixed(DISPLACEMENT_X)) #it is left unfixed at the end of the previous interval

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        for node in model_part.Nodes:
            self.assertFalse(node.IsFixed(DENSITY))
            self.assertFalse(node.IsFixed(VISCOSITY))
            self.assertFalse(node.IsFixed(DISPLACEMENT_X))

        ############################################################
        ##time = 12 - DISPLACEMENT applied as a vector. x,z components fixed, y component not imposed
        ##time = 12 - VELOCITY applied as a vector by componentes. All components free. x component is not zero.
        model_part.CloneTimeStep(12.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 10.0)
            self.assertTrue(node.IsFixed(DISPLACEMENT_X))

            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 0.0) #not applied!!
            self.assertFalse(node.IsFixed(DISPLACEMENT_Y))

            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), 12.0)
            self.assertTrue(node.IsFixed(DISPLACEMENT_Z))

            self.assertEqual(node.GetSolutionStepValue(VELOCITY_X), 10.0)
            self.assertFalse(node.IsFixed(VELOCITY_X))

            self.assertEqual(node.GetSolutionStepValue(VELOCITY_Y), 0.0)
            self.assertFalse(node.IsFixed(VELOCITY_Y))

            self.assertEqual(node.GetSolutionStepValue(VELOCITY_Z), 0.0)
            self.assertFalse(node.IsFixed(VELOCITY_Z))

        #print("**********************************************")
        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        ############################################################
        ##time >= 20 - DISPLACEMENT applied as a vector. x,z components fixed, y component not imposed
        ##time >= 20 - VELOCITY applied as a vector by componentes. All components free. y component is not zero.
        model_part.CloneTimeStep(20.1)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        #print("Checking time = ", t)
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 10.0)
            self.assertFalse(node.IsFixed(DISPLACEMENT_X))

            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 0.0) #not applied!!
            self.assertFalse(node.IsFixed(DISPLACEMENT_Y))

            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), t)
            self.assertFalse(node.IsFixed(DISPLACEMENT_Z))

            self.assertEqual(node.GetSolutionStepValue(VELOCITY_X), 0.0)
            self.assertFalse(node.IsFixed(VELOCITY_X))

            self.assertEqual(node.GetSolutionStepValue(VELOCITY_Y), math.sin(node.X*math.pi*t))
            self.assertFalse(node.IsFixed(VELOCITY_Y))

            self.assertEqual(node.GetSolutionStepValue(VELOCITY_Z), 0.0)
            self.assertFalse(node.IsFixed(VELOCITY_Z))

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        ############################################################
        ##time >= 25 - DISPLACEMENT applied as a vector. x,z components fixed, y component not imposed
        ##time >= 25 - VELOCITY applied as a vector by componentes. All components fixed. y and z components are not zero.
        model_part.CloneTimeStep(26.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        #print("Checking time = ", t)
        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 10.0) #previous value
            self.assertFalse(node.IsFixed(DISPLACEMENT_X)) #not fixed since set as null

            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), node.X+node.Y*t) #not applied!!
            self.assertTrue(node.IsFixed(DISPLACEMENT_Y)) #set to true

            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), t)
            self.assertFalse(node.IsFixed(DISPLACEMENT_Z))

            self.assertEqual(node.GetSolutionStepValue(VELOCITY_X), 0.0)
            self.assertTrue(node.IsFixed(VELOCITY_X))

            self.assertAlmostEqual(node.GetSolutionStepValue(VELOCITY_Y), (math.sqrt(abs(node.X*node.Y)))/math.sqrt(2))
            self.assertTrue(node.IsFixed(VELOCITY_Y))

            self.assertAlmostEqual(node.GetSolutionStepValue(VELOCITY_Z), (math.sqrt(abs(node.X*node.Y)))/math.sqrt(2))
            self.assertTrue(node.IsFixed(VELOCITY_Z))

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

    def test_rotated_system(self):
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        #note that y and z are inverted in the rotated system
        settings = Parameters(
            """
            {
                "process_list" : [
                    {
                        "python_module"   : "assign_scalar_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableProcess",
                        "Parameters"            : {
                            "model_part_name" : "Main",
                            "variable_name"   : "VISCOSITY",
                            "interval"        : [0.0, 10.0],
                            "constrained"     : false,
                            "value"      : "x+100.0*y*t**2",
                            "local_axes"               :{
                                "origin" : [10.0, 0.0, 0.0],
                                "axes"  : [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0] ]
                            }
                        }
                    }
                    ]
                }
            """
            )

        Model = {"Main":model_part}

        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( settings["process_list"] )

        model_part.CloneTimeStep(3.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        ##verify the result
        t = model_part.ProcessInfo[TIME]
        for node in model_part.Nodes:
            x = node.X - 10.0
            y = node.Z
            z = node.Y
            self.assertEqual(node.GetSolutionStepValue(VISCOSITY), x+100.0*y*t**2)
            self.assertFalse(node.IsFixed(VISCOSITY))

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()


    def test_assign_scalar_value_to_conditions(self):
        model_part = ModelPart("Main")
        model_part_io = ModelPartIO(GetFilePath("test_processes"))
        model_part_io.ReadModelPart(model_part)

        settings = Parameters(
            """
            {
                "process_list" : [
                    {
                        "python_module"   : "assign_scalar_variable_to_conditions_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableToConditionsProcess",
                        "Parameters"            : {
                            "model_part_name":"Main",
                            "variable_name": "PRESSURE",
                            "value" : 15.0
                        }
                    },
                    {
                        "python_module"   : "assign_scalar_variable_to_conditions_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableToConditionsProcess",
                        "Parameters"            : {
                            "model_part_name":"Main",
                            "variable_name": "VISCOSITY",
                            "value" : 2
                        }
                    }
                    ]
                }
            """
            )

        Model = {"Main":model_part}

        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( settings["process_list"] )

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        for cond in model_part.Conditions:
            self.assertEqual(cond.GetValue(PRESSURE), 15.0)
            self.assertEqual(cond.GetValue(VISCOSITY), 2)


    def test_assign_scalar_field_to_conditions(self):
        model_part = ModelPart("Main")
        model_part_io = ModelPartIO(GetFilePath("test_processes"))
        model_part_io.ReadModelPart(model_part)

        settings = Parameters(
            """
            {
                "process_list" : [
                    {
                        "python_module"   : "assign_scalar_variable_to_conditions_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableToConditionsProcess",
                        "Parameters"            : {
                            "model_part_name":"Main",
                            "variable_name": "INITIAL_STRAIN",
                            "value" : "x+y*t+z"
                        }
                    }
                    ]
                }
            """
            )

        Model = {"Main":model_part}

        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( settings["process_list"] )

        model_part.CloneTimeStep(5.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        t = model_part.ProcessInfo[TIME]
        for cond in model_part.Conditions:
            v = cond.GetValue(INITIAL_STRAIN)

            i = 0
            for node in cond.GetNodes():
                self.assertEqual(v[i],node.X+node.Y*t+node.Z)
                i=i+1
                
    def test_assign_scalar_field_component_to_conditions(self):
        model_part = ModelPart("Main")
        model_part_io = ModelPartIO(GetFilePath("test_processes"))
        model_part_io.ReadModelPart(model_part)

        settings = Parameters(
            """
            {
                "process_list" : [
                    {
                        "python_module"   : "assign_scalar_variable_to_conditions_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignScalarVariableToConditionsProcess",
                        "Parameters"            : {
                            "model_part_name":"Main",
                            "variable_name": "DISPLACEMENT_X",
                            "value" : "t"
                        }
                    }
                    ]
                }
            """
            )

        Model = {"Main":model_part}

        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( settings["process_list"] )

        model_part.CloneTimeStep(5.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        t = model_part.ProcessInfo[TIME]
        for cond in model_part.Conditions:
            v = cond.GetValue(DISPLACEMENT)
            self.assertEqual(v[0],t)

    def test_find_nodal_h_process(self):
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(NODAL_H)
        model_part_io = ModelPartIO(GetFilePath("test_processes"))
        model_part_io.ReadModelPart(model_part)

        FindNodalHProcess(model_part).Execute();

        for i in range(1,len(model_part.Nodes)):
            self.assertEqual(model_part.GetNode(i).GetSolutionStepValue(NODAL_H), 0.25)
        self.assertEqual(model_part.GetNode(len(model_part.Nodes)).GetSolutionStepValue(NODAL_H), 0.5)

    def test_assign_acceleration_to_nodes(self):
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(ACCELERATION)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)

        model_part_io = ModelPartIO(GetFilePath("test_processes"))
        model_part_io.ReadModelPart(model_part)

        settings = Parameters(
            """
            {
                "process_list" : [
                    {
                        "python_module"   : "assign_time_derivative_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"          : "AssignTimeDerivativeProcess",
                        "Parameters"            : {
                            "model_part_name":"Main",
                            "variable_name" : "ACCELERATION",
                            "variable_to_be_solved_for" : "DISPLACEMENT",
                            "value" : ["t",null,"z"],
                            "interval" : [3.0,4.0]
                        }
                    }
                    ]
                }
            """
            )

        Model = {"Main":model_part}

        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( settings["process_list"] )

        ################### here we are within the interval
        model_part.CloneTimeStep(3.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()


        for node in model_part.Nodes:
            self.assertEqual(node.IsFixed(ACCELERATION_X), True)
            self.assertEqual(node.IsFixed(ACCELERATION_Y), False)
            self.assertEqual(node.IsFixed(ACCELERATION_Z), True)
            self.assertEqual(node.IsFixed(DISPLACEMENT_X), True)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Y), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Z), True)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_X), 3.0) #t = 3.0
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Z), node.Z)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 0.0) #displacements remain unmodified, they will be assigned by the scheme
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), 0.0)

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        for node in model_part.Nodes:
            self.assertEqual(node.IsFixed(ACCELERATION_X), False)
            self.assertEqual(node.IsFixed(ACCELERATION_Y), False)
            self.assertEqual(node.IsFixed(ACCELERATION_Z), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_X), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Y), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Z), False)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_X), 3.0) #t = 3.0
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Z), node.Z)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 0.0) #displacements remain unmodified, they will be assigned by the scheme
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), 0.0)

        ################### here we are outside of the interval - values do not change but everything is free
        model_part.CloneTimeStep(8.0)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertEqual(node.IsFixed(ACCELERATION_X), False)
            self.assertEqual(node.IsFixed(ACCELERATION_Y), False)
            self.assertEqual(node.IsFixed(ACCELERATION_Z), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_X), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Y), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Z), False)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_X), 3.0) #t = 3.0
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Z), node.Z)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 0.0) #displacements remain unmodified, they will be assigned by the scheme
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), 0.0)

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        for node in model_part.Nodes:
            self.assertEqual(node.IsFixed(ACCELERATION_X), False)
            self.assertEqual(node.IsFixed(ACCELERATION_Y), False)
            self.assertEqual(node.IsFixed(ACCELERATION_Z), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_X), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Y), False)
            self.assertEqual(node.IsFixed(DISPLACEMENT_Z), False)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_X), 3.0) #t = 3.0
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(ACCELERATION_Z), node.Z)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_X), 0.0) #displacements remain unmodified, they will be assigned by the scheme
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Y), 0.0)
            self.assertEqual(node.GetSolutionStepValue(DISPLACEMENT_Z), 0.0)
if __name__ == '__main__':
    KratosUnittest.main()
