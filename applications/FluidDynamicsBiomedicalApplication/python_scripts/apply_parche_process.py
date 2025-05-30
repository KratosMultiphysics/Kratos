import math
import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

import numpy as np

import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

import KratosMultiphysics.FluidDynamicsBiomedicalApplication as KratosBio

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyParcheProcess(Model, settings["Parameters"])


class ApplyParcheProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"               : 0,
            "model_part_name"       : "",
            "value"                 : {},
            "thickness"             : 0.0,
            "resistance_multiplier" : 0.0,
            "max_porosity"          : 0.0,
            "initial_time"          : 0.0,
            "ramp_exponent"         : 0.0,
            "interval"              : [0.0,"End"]
        }
        """)

        # Trick: allows "value" to be a double, a string or a table value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")
            elif settings["value"].IsNumber():
                default_settings["value"].SetDouble(0.0)
        else:
            err_msg = "Provided settings have no 'value'. This needs to be provided."
            raise Exception(err_msg)

        settings.ValidateAndAssignDefaults(default_settings)

        # Set a Kratos parameters suitable for the core processes to set the PRESSURE
        model_settings = settings.Clone()

        # Create a copy of the PRESSURE settings to set the EXTERNAL_PRESSURE

        # Check the core processes input data
        if (model_settings["model_part_name"].GetString() == ""):
            raise Exception("Empty parche model part name. Set a valid model part name.")
        elif (model_settings["value"].IsString()):
            if (model_settings["value"].GetString == ""):
                raise Exception("Porosity function sting is empty.")

#        self.hydrostatic_outlet = settings["hydrostatic_outlet"].GetBool()
#        self.h_top = settings["h_top"].GetDouble()
        self.thickness      = settings["thickness"].GetDouble()
        self.Cd             = settings["resistance_multiplier"].GetDouble()
        self.max_porosity   = settings["max_porosity"].GetDouble()
        self.initial_time   = settings["initial_time"].GetDouble()
        self.ramp_exponent  = settings["ramp_exponent"].GetDouble()

        if (self.ramp_exponent < 0):
            raise Exception("The exponent must be positive for a ramp!")

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.parche_model_part = Model[model_settings["model_part_name"].GetString()]

#        self.structure_model_part = self._GetSolver().GetStructureComputingModelPart()
#        self.parche_smp = self.structure_model_part.GetSubModelPart("DISPLACEMENT_Parche")

        self.one_forth = 1.0/4.0


    def ExecuteInitializeSolutionStep(self):

        t = self._GetSolver().GetFluidComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]

        phi0 = self.max_porosity
        t0 = self.initial_time
        esp = self.ramp_exponent


        if t > t0:

            porosity = 0

            self.FindApproximatingVolume(self)

            # print(normal)

            # print(cdiff)

            i = 0
            phi = phi0*(1 - math.exp(-esp*(t-t0)))

            for element in self._GetSolver().GetFluidComputingModelPart().Elements:

                center = element.GetGeometry().Center()
                dist2center = (center[0] - self.pos_average[0])*(center[0] - self.pos_average[0]) + (center[1] - self.pos_average[1])*(center[1] - self.pos_average[1]) + (center[2] - self.pos_average[2])*(center[2] - pos_average[2])
                porosity = 0
                if (center[2] - self.cdiff < self.a*center[0] + self.b*center[1] + self.c < center[2] + self.cdiff and dist2center < self.dist2max*1.5):
                    i = i+1
                    nodes_el = element.GetNodes()
                    vg = np.array([0.0, 0.0, 0.0])
                    for node in nodes_el:
                        vg0 = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                        vg += vg0
                    vg *= self.one_forth

                    norm_vg = math.sqrt(vg[0]*vg[0] + vg[1]*vg[1] + vg[2]*vg[2])

                    porosity = 0.5*self.Cd*phi/(1.0 - phi)*norm_vg
                element.SetValue(KratosFluid.RESISTANCE, porosity)
            print(i,phi)




    # Private methods section
    def FindApproximatingVolume(self):

        lhs = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        rhs = np.array([0.0, 0.0, 0.0])

        self.pos_average = np.array([0.0, 0.0, 0.0])

        for node in self.parche_smp.Nodes:
            self.pos_average += np.array([node.X, node.Y, node.Z])
        self.pos_average /= self.parche_smp.NumberOfNodes()

        # print(pos_average)

        self.dist2max = 0

        x0 = self.pos_average[0]
        y0 = self.pos_average[1]
        z0 = self.pos_average[2]

        for node in self.parche_smp.Nodes:

            xi = node.X
            yi = node.Y
            zi = node.Z

            dist2 = (xi - x0)*(xi - x0) + (yi - y0)*(yi - y0) + (zi - z0)*(zi - z0)
            if dist2 > self.dist2max:
                self.dist2max = dist2
                #print(dist2max)

            lhs += np.array([[xi*xi, xi*yi, xi],
                                [xi*yi, yi*yi, yi],
                                [ xi,    yi,   1 ]])

            rhs += np.array([xi*zi, yi*zi, zi])

        plane_coefficients = np.linalg.solve(lhs, rhs)

        self.a = plane_coefficients[0]
        self.b = plane_coefficients[1]
        self.c = plane_coefficients[2]

        self.normal = np.array([self.a ,self.b, -1.0])/np.sqrt(self.a*self.a + self.b*self.b + 1)
        self.cdiff = np.abs(self.thickness/self.normal[2])


