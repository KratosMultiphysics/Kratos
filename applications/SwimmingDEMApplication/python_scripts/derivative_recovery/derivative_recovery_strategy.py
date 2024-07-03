# importing the Kratos Library
import sys

from . import recoverer
from . import standard_recoverer
from . import zhang_guo_recoverer
from . import L2_projection_recoverer
from . import pouliot_2012_recoverer
from . import pouliot_2012_edge_recoverer
from . import lagrangian_recoverer

import weakref

class DerivativeRecoveryStrategy:
    def __init__(self, parameters, fluid_model_part, custom_functions_tool = None):
        self.parameters = parameters
        self.fluid_model_part = fluid_model_part
        if custom_functions_tool is not None:
            self.custom_functions_tool = weakref.proxy(custom_functions_tool)
        else:
            self.custom_functions_tool = None

        self.store_full_gradient = self.parameters["store_full_gradient_option"].GetBool()
        self.mat_deriv_type = self.parameters["material_acceleration_calculation_type"].GetInt()
        self.laplacian_type = self.parameters["laplacian_calculation_type"].GetInt()
        self.vorticity_type = self.parameters["vorticity_calculation_type"].GetInt()
        self.pressure_grad_type = self.parameters["pressure_grad_recovery_type"].GetInt()
        self.fluid_fraction_grad_type = self.parameters["coupling"]["backward_coupling"]["fluid_fraction_grad_type"].GetInt()

        self.do_pre_recovery = False
        self.must_reconstruct_gradient = self.laplacian_type in {0, 3, 4, 5, 6} and self.mat_deriv_type in {2,3,4}

        # if self.parameters["custom_fluid"]["fluid_already_calculated"].GetBool(): # the fluid has been calculated before, and the derivatives fed to the fluid_model_part
        #     self.pre_computed_derivatives = self.parameters["load_derivatives"].GetBool()
        # else:
        self.pre_computed_derivatives = False

        self.mat_deriv_tool = self.GetMatDerivTool()
        self.laplacian_tool = self.GetLaplacianTool()
        self.vorticity_tool = self.GetVorticityTool()
        self.velocity_grad_tool = self.GetVelocityGradTool()
        self.pressure_grad_tool = self.GetPressureGradTool()
        self.fluid_fraction_grad_tool = self.GetFluidFractionGradTool()

    def GetMatDerivTool(self):
        if self.pre_computed_derivatives:
            return recoverer.MaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        elif self.mat_deriv_type == 0:
            return recoverer.EmptyGradientRecoverer(self.parameters, self.fluid_model_part)
        elif self.mat_deriv_type == 1:
            return standard_recoverer.StandardMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        elif self.mat_deriv_type == 2:
            if self.laplacian_type == 2:
                return zhang_guo_recoverer.ZhangGuoMaterialAccelerationAndLaplacianRecoverer(self.parameters, self.fluid_model_part)
            else:
                return zhang_guo_recoverer.ZhangGuoMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        elif self.mat_deriv_type == 3:
            return L2_projection_recoverer.L2ProjectionDirectMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        elif self.mat_deriv_type == 4:
            return L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        elif self.mat_deriv_type == 5:
            return L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        elif self.mat_deriv_type == 6:
            return pouliot_2012_edge_recoverer.Pouliot2012EdgeMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
            #return pouliot_2012_recoverer.Pouliot2012MaterialAccelerationRecoverer(self.parameters, self.fluid_model_part, self.do_pre_recovery)
        elif self.mat_deriv_type == 7:
            if self.store_full_gradient:
                return zhang_guo_recoverer.ZhangGuoMaterialAccelerationAndLaplacianRecoverer(self.parameters, self.fluid_model_part)
            else:
                raise Exception('The value of material_acceleration_calculation_type is ' + str(self.mat_deriv_type) + ' , which is only compatible with store_full_gradient = 1. Instead it is store_full_gradient = ' + str(self.store_full_gradient) + '.')
        elif self.mat_deriv_type == 8:
            return lagrangian_recoverer.LagrangianMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        else:
            raise Exception('The value of material_acceleration_calculation_type is ' + str(self.mat_deriv_type) + ' , which does not correspond to any valid option.')

    def GetLaplacianTool(self):
        if self.laplacian_type == 0 or self.pre_computed_derivatives:
            return recoverer.EmptyLaplacianRecoverer(self.parameters, self.fluid_model_part)
        elif self.laplacian_type == 1:
            return standard_recoverer.StandardLaplacianRecoverer(self.parameters, self.fluid_model_part)
        elif self.laplacian_type == 2:
            # if self.mat_deriv_type == 2:
            #     return zhang_guo_recoverer.ZhangGuoMaterialAccelerationAndLaplacianRecoverer(self.parameters, self.fluid_model_part)
            # else:
            return zhang_guo_recoverer.ZhangGuoDirectLaplacianRecoverer(self.parameters, self.fluid_model_part)
        elif self.laplacian_type in {3, 4, 5, 6}:
            if self.store_full_gradient:
                return L2_projection_recoverer.L2ProjectionLaplacianRecoverer(self.parameters, self.fluid_model_part)
            else:
                raise Exception('The value of laplacian_calculation_type is ' + str(self.laplacian_type) + ' , which is only compatible with store_full_gradient = 1. Instead it is store_full_gradient = ' + str(self.store_full_gradient) + '.')
        elif self.laplacian_type == 7:
            if self.store_full_gradient:
                return zhang_guo_recoverer.ZhangGuoMaterialAccelerationAndLaplacianRecoverer(self.parameters, self.fluid_model_part)
            else:
                raise Exception('The value of laplacian_calculation_type is ' + str(self.mat_deriv_type) + ' , which is only compatible with store_full_gradient = 1. Instead it is store_full_gradient = ' + str(self.store_full_gradient) + '.')
        else:
            raise Exception('The value of laplacian_calculation_type is ' + str(self.laplacian_type) + ' , which does not correspond to any valid option.')

    def GetVorticityTool(self):
        #if self.pre_computed_derivatives:
        return recoverer.VorticityRecoverer(self.parameters, self.fluid_model_part)
        if self.vorticity_type == 0:
            return recoverer.EmptyVorticityRecoverer(self.parameters, self.fluid_model_part)
        elif self.vorticity_type == 5: # NOT WORKING YET WITH VECTORS
            return recoverer.EmptyVorticityRecoverer(self.parameters, self.fluid_model_part)
        else:
            raise Exception('The value of vorticity_calculation_type is ' + str(self.vorticity_type) + ' , which does not correspond to any valid option.')

    def GetVelocityGradTool(self):
        if self.pre_computed_derivatives:
            return recoverer.EmptyGradientRecoverer(self.parameters, self.fluid_model_part)
        if self.must_reconstruct_gradient:
            if self.mat_deriv_type == 2:
                return zhang_guo_recoverer.ZhangGuoGradientRecoverer(self.parameters, self.fluid_model_part)
            elif self.mat_deriv_type == 6:
                return pouliot_2012_edge_recoverer.Pouliot2012EdgeGradientRecoverer(self.parameters, self.fluid_model_part)
                # return pouliot_2012_recoverer.Pouliot2012GradientRecoverer(self.parameters, self.fluid_model_part)
            else:
                return L2_projection_recoverer.L2ProjectionGradientRecoverer(self.parameters, self.fluid_model_part)
        else:
            return recoverer.EmptyGradientRecoverer(self.parameters, self.fluid_model_part)

    def GetPressureGradTool(self):
        if self.pressure_grad_type == 0:
            return recoverer.EmptyGradientRecoverer(self.parameters, self.fluid_model_part)
        elif self.pressure_grad_type == 1:
            return standard_recoverer.StandardGradientRecoverer(self.parameters, self.fluid_model_part)
        elif self.pressure_grad_type == 2:
            return zhang_guo_recoverer.ZhangGuoGradientRecoverer(self.parameters, self.fluid_model_part)
        else:
            raise Exception('The value of pressure_grad_recovery_type is ' + str(self.pressure_grad_type) + ' , which does not correspond to any valid option.')

    def GetFluidFractionGradTool(self):
        if self.fluid_fraction_grad_type == 0:
            return recoverer.EmptyGradientRecoverer(self.parameters, self.fluid_model_part)
        elif self.fluid_fraction_grad_type == 1:
            return standard_recoverer.StandardGradientRecoverer(self.parameters, self.fluid_model_part)
        elif self.fluid_fraction_grad_type == 2:
            return zhang_guo_recoverer.ZhangGuoGradientRecoverer(self.parameters, self.fluid_model_part)
        elif self.fluid_fraction_grad_type == 3:
            return L2_projection_recoverer.L2ProjectionFluidFractionGradientRecoverer(self.parameters, self.fluid_model_part)
        else:
            raise Exception('The value of fluid_fraction_grad_type is ' + str(self.fluid_fraction_grad_type) + ' , which does not correspond to any valid option.')

    def Recover(self):
        # Some of the following may be empty, and some may do the work of others for efficiency.
        self.mat_deriv_tool.RecoverMaterialAcceleration()
        # standard = standard_recoverer.StandardMaterialAccelerationRecoverer(self.parameters, self.fluid_model_part)
        # standard.cplusplus_recovery_tool.SmoothVectorField(self.fluid_model_part, MATERIAL_ACCELERATION, VELOCITY_COMPONENT_GRADIENT)
        self.velocity_grad_tool.RecoverGradientOfVelocity()
        self.pressure_grad_tool.RecoverPressureGradient()
        self.vorticity_tool.RecoverVorticityFromGradient()
        self.laplacian_tool.RecoverVelocityLaplacian()

    def RecoverFluidFractionGradient(self):
        self.fluid_fraction_grad_tool.RecoverFluidFractionGradient()
