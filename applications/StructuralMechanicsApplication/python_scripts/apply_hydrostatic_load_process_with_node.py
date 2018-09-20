from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import sys
from math import *


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return ApplyHydrostaticLoadProcess(Model, settings["Parameters"])


class ApplyHydrostaticLoadProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
                    {
                    "mesh_id"                       : 0,
                    "properties_id"                 : 0,
                    "main_model_part_name"          : "Structure",
                    "model_part_name"               : "SurfacePressure3D_hemisphere",
                    "specific_weight"               : "200*t",
                    "interval"                      : [0.0, 1e30],
                    "local_axes"                    : {},
                    "fluid_volume"                  : 1.5,
                    "nodal_element_id"              : 0,
                    "plane_normal"                  : [0.0,0.0,1.0],
                    "initial_free_surface_radius"   : 0.0              
                
            }
            """
                                                         )

        # assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # here i do a trick, since i want to allow "specific_weight" to be a string or a double specific_weight
        if(settings.Has("specific_weight")):
            if(settings["specific_weight"].IsString()):
                default_settings["specific_weight"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(
            "SPECIFIC_WEIGHT")

        if(type(self.variable) != KratosMultiphysics.DoubleVariable):
            msg = "Error in ApplyHydrostaticLoadProcess. Variable type of variable specific_weight is incorrect . Must be a scalar "
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.main_model_part = Model[settings["main_model_part_name"].GetString(
        )]
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        self.specific_weight_is_numeric = False

        if settings["specific_weight"].IsNumber():
            self.specific_weight_is_numeric = True
            self.specific_weight = settings["specific_weight"].GetDouble()
        else:
            self.function_string = settings["specific_weight"].GetString()
            self.aux_function = KratosMultiphysics.PythonGenericFunctionUtility(
                self.function_string, settings["local_axes"])

            if(self.aux_function.DependsOnSpace()):
                raise RuntimeError(
                    "Specific density varying in space is not yet implemented")

        self.step_is_active = False

        self.fluid_volume = settings["fluid_volume"].GetDouble()
        nodal_element_id = settings["nodal_element_id"].GetInt()
        print("####NODAL ELEMENT ID :: ", nodal_element_id)

        if (nodal_element_id == 0):
            raise RuntimeError(
                "The elment Id of NodalFluidConcentratedElement3D3N has to be given as the argument")

        n_x = settings["plane_normal"].GetVector()[0]
        n_y = settings["plane_normal"].GetVector()[1]
        n_z = settings["plane_normal"].GetVector()[2]
        self.plane_normal = [n_x, n_y, n_z]
        self.initial_free_surface_radius = settings["initial_free_surface_radius"].GetDouble(
        )
        properties_id = settings["properties_id"].GetInt()

        self.properties = self.main_model_part.GetProperties()[properties_id]

        for cond in self.model_part.Conditions:
            cond.Properties = self.properties

        self.nodal_elem = self.main_model_part.GetElement(nodal_element_id)
        self.nodal_elem.Properties = self.properties
        self.free_surface_centre = [self.nodal_elem.GetNode(
            0).X, self.nodal_elem.GetNode(0).Y, self.nodal_elem.GetNode(0).Z]

        print("####### AfterConstructor######################")
        print("NODAL ELEMENT")
        print(self.nodal_elem)
        print("FREE SURFACE CENRE")
        print(self.free_surface_centre)

    def ExecuteInitialize(self):

        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_RADIUS, self.initial_free_surface_radius)
        self.properties.SetValue(
            StructuralMechanicsApplication.FLUID_VOLUME, self.fluid_volume)
        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_NORMAL, self.plane_normal)

        self.VolumeCalcUtilty = StructuralMechanicsApplication.VolumeCalculationUnderPlaneUtility(
            self.free_surface_centre, self.initial_free_surface_radius, self.plane_normal)
        self.VolumeCalcUtilty.UpdatePositionOfPlaneBasedOnTargetVolume(
            self.model_part, self.fluid_volume, 1E-6, 20)
        self.free_surface_centre = self.VolumeCalcUtilty.GetCentre()

        self.nodal_elem.GetNode(0).X = self.free_surface_centre[0]
        self.nodal_elem.GetNode(0).Y = self.free_surface_centre[1]
        self.nodal_elem.GetNode(0).Z = self.free_surface_centre[2]

        self.nodal_elem.GetNode(0).X0 = self.free_surface_centre[0]
        self.nodal_elem.GetNode(0).Y0 = self.free_surface_centre[1]
        self.nodal_elem.GetNode(0).Z0 = self.free_surface_centre[2]

        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_AREA, self.VolumeCalcUtilty.GetIntersectedArea())
        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_CENTRE, self.free_surface_centre)

        avg_nodes = 10
        avg_elems = 10

        neighbhor_finder = KratosMultiphysics.FindNodalNeighboursProcess(
            self.main_model_part, avg_elems, avg_nodes)

        neighbhor_finder.Execute()

        self.nodal_elem.SetValue(
            StructuralMechanicsApplication.STIFFNESS_SCALING, 1.0)
        self.nodal_elem.SetValue(
            StructuralMechanicsApplication.IS_DYNAMIC, False)

        master_creator = StructuralMechanicsApplication.BuildMasterConditionsForHydrostaticLoadingProcess(
            self.model_part, self.nodal_elem)

        master_creator.Execute()

        print("####### AfterConstructor######################")
        print("NODAL ELEMENT")
        print(self.nodal_elem)
        print("FREE SURFACE CENRE")
        print(self.free_surface_centre)

    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.interval.IsInInterval(current_time)):

            self.step_is_active = True

            if self.specific_weight_is_numeric:
                self.properties.SetValue(
                    StructuralMechanicsApplication.SPECIFIC_WEIGHT, self.specific_weight)
            else:
                if self.aux_function.DependsOnSpace() == False:  # depends on time only
                    self.specific_weight = self.aux_function.CallFunction(
                        0.0, 0.0, 0.0, current_time)
                    self.properties.SetValue(
                        StructuralMechanicsApplication.SPECIFIC_WEIGHT, self.specific_weight)

                # most general case - space varying function (possibly also time varying)
                else:
                    raise RuntimeError(
                        "Specific density varying in space is not yet implemented")

    def ExecuteFinalizeSolutionStep(self):

        if (self.step_is_active):

            centre = self.properties.GetValue(
                StructuralMechanicsApplication.FREE_SURFACE_CENTRE)
            specific_wt = self.properties.GetValue(
                StructuralMechanicsApplication.SPECIFIC_WEIGHT)

            self.VolumeCalcUtilty.SetPlaneParameters(
                centre, self.free_surface_radius, self.plane_normal)
            vol = self.VolumeCalcUtilty.CalculateVolume(self.main_model_part)
            area = self.VolumeCalcUtilty.GetIntersectedArea()

            print("#######################################")
            print("Specific_wt ::", specific_wt)
            print("Free Surface Centre :: ", centre)
            print("Fluid volume :: ", vol)
            print("Free surface area :: ", area)
            print("#######################################")

            print("#######################################")

            print("Nodal Surface Centre X :: ", self.nodal_elem.GetNode(0).X)
            print("Nodal Surface Centre Y :: ", self.nodal_elem.GetNode(0).Y)
            print("Nodal Surface Centre Z :: ", self.nodal_elem.GetNode(0).Z)

            print("#######################################")

            self.step_is_active = False
