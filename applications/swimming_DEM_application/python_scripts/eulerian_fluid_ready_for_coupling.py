from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

class Solution(object):

    def __init__(self, model):
        self.model = model
        import ProjectParameters as pp
        self.pp = pp
        self.fluid_model_part = self.model.CreateModelPart("FluidPart")

    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

    def Initialize(self):
        #self.fluid_model_part = self.all_model_parts.Get("FluidPart')
        self.SetFluidSolverModule()
        self.AddFluidVariables()
        self.ReadFluidModelPart()
        self.SetFluidBufferSizeAndAddDofs()
        self.SetFluidSolver()
        self.fluid_solver.Initialize()
        self.ActivateTurbulenceModel()

    def SetFluidSolverModule(self):
        self.SolverSettings = self.pp.FluidSolverConfiguration
        self.solver_module = import_solver(self.SolverSettings)

    def SetFluidSolver(self):
        self.fluid_solver = self.solver_module.CreateSolver(self.fluid_model_part, self.SolverSettings)

    def AddFluidVariables(self):
        # caution with breaking up this block (memory allocation)! {
        self.AddFluidVariablesByFluidSolver()
        self.AddFluidVariablesBySwimmingDEMAlgorithm()

    def AddFluidVariablesByFluidSolver(self):

        if "REACTION" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

        self.solver_module.AddVariables(self.fluid_model_part, self.SolverSettings)

    def AddFluidVariablesBySwimmingDEMAlgorithm(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.pp.fluid_vars)

    def SetFluidBufferSizeAndAddDofs(self):
        self.fluid_model_part.SetBufferSize(3)
        self.solver_module.AddDofs(self.fluid_model_part, self.SolverSettings)

    def ReadFluidModelPart(self):
        os.chdir(self.main_path)
        model_part_io_fluid = ModelPartIO(self.pp.problem_name)
        model_part_io_fluid.ReadModelPart(self.fluid_model_part)

    def ActivateTurbulenceModel(self):

        for element in self.fluid_model_part.Elements:
            element.SetValue(C_SMAGORINSKY, 0.0)

        if self.pp.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras":
            # apply the initial turbulent viscosity on all of the nodes
            turb_visc = self.pp.FluidSolverConfiguration.TurbulentViscosity
            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
                visc = node.GetSolutionStepValue(VISCOSITY)
                node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
                if node.IsFixed(VELOCITY_X) and node.GetSolutionStepValue(VELOCITY_X, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Y) and node.GetSolutionStepValue(VELOCITY_Y, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Z) and node.GetSolutionStepValue(VELOCITY_Z, 0) != 0.0:
                    node.Fix(TURBULENT_VISCOSITY)

            # select nodes on the wall
            self.fluid_solver.wall_nodes = []
            for i in self.SolverSettings.SA_wall_group_ids:
                # get the nodes of the wall for SA.
                nodes = self.fluid_model_part.GetNodes(i)
                for node in nodes:
                    self.fluid_solver.wall_nodes.append(node)
                    node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
                    node.Fix(TURBULENT_VISCOSITY)


    def Finalize(self):
        pass



if __name__ == "__main__":
    model = Model()
    Solution(model).Run()
