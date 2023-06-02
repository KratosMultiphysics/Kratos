import KratosMultiphysics

def GetResidualOutputAnalysisClass(cls):
    class ResidualOutputAnalysis(cls):

        def __init__(self, model,project_parameters):
            super().__init__(model,project_parameters)

            self.residuals_list = []

        def ModifyInitialGeometry(self):
            super().ModifyInitialGeometry()

            self.main_model_part = self.model.GetModelPart("Structure")


        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

            self.strategy = self._GetSolver()._GetSolutionStrategy()
            self.buildsol = self._GetSolver()._GetBuilderAndSolver()
            self.scheme = self._GetSolver()._GetScheme()

            A = self.strategy.GetSystemMatrix()
            b = self.strategy.GetSystemVector()

            space = KratosMultiphysics.UblasSparseSpace()

            space.SetToZeroMatrix(A)
            space.SetToZeroVector(b)

            # Create dummy vector
            xD = space.CreateEmptyVectorPointer()
            space.ResizeVector( xD, space.Size1(A) )
            space.SetToZeroVector(xD)

            self.buildsol.Build(self.scheme, self.main_model_part, A, b)

            # Apply constraints
            self.buildsol.ApplyConstraints(self.scheme, self.main_model_part, A, b)
            # Apply boundary conditions
            self.buildsol.ApplyDirichletConditions(self.scheme, self.main_model_part, A, xD, b)
            
            self.residuals_list.append(b)

        def GetResiduals(self):
            return self.residuals_list

    return ResidualOutputAnalysis