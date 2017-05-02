from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.AdjointFluidApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestCase(KratosUnittest.TestCase):

    def addVariables(self,model_part):
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(ACCELERATION)
        model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
        model_part.AddNodalSolutionStepVariable(PRESSURE)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(BODY_FORCE)
        model_part.AddNodalSolutionStepVariable(REACTION)
        model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
        model_part.AddNodalSolutionStepVariable(ADJOINT_VELOCITY)
        model_part.AddNodalSolutionStepVariable(ADJOINT_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(ADJOINT_PRESSURE)
        model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)
        model_part.AddNodalSolutionStepVariable(NORMAL_SENSITIVITY)

    def assignConstantSolutionStepValue(self,model_part,var,step,val):
        for node in model_part.Nodes:
            node.SetSolutionStepValue(var,step,val)

    def getBufferSize(self):
        return 2

    def getDomainSize(self):
        return 2

    def getAlphaBossak(self):
        return -0.3

    def getDeltaTime(self):
        return 1.0

    def createModelPart(self):
        # create test model part
        model_part = ModelPart("test")
        self.addVariables(model_part)
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        prop = model_part.GetProperties()[0]
        model_part.CreateNewElement("VMS2D3N", 1, [1, 2, 3], prop)
        buffer_size = self.getBufferSize()
        model_part.SetBufferSize(buffer_size)
        # fluid properties
        for step in range(buffer_size):
            self.assignConstantSolutionStepValue(model_part,DENSITY,step,1.0)
            self.assignConstantSolutionStepValue(model_part,VISCOSITY,step,1.0e-5)
        # dofs
        for node in model_part.Nodes:
            node.AddDof(VELOCITY_X,REACTION_X)
            node.AddDof(VELOCITY_Y,REACTION_Y)
            node.AddDof(VELOCITY_Z,REACTION_Z)
            node.AddDof(PRESSURE,REACTION_WATER_PRESSURE)
        # bcs
        model_part.GetNode(1).Fix(VELOCITY_X)
        model_part.GetNode(1).Fix(VELOCITY_Y)
        model_part.GetNode(2).Fix(VELOCITY_X)
        model_part.GetNode(2).Fix(VELOCITY_Y)
        model_part.GetNode(3).Fix(VELOCITY_X)
        model_part.GetNode(3).Fix(VELOCITY_Y)
        model_part.GetNode(2).Fix(PRESSURE)
        model_part.GetNode(3).SetSolutionStepValue(VELOCITY_X,0,1.0)
        model_part.GetNode(3).SetSolutionStepValue(VELOCITY_X,1,1.0)
        model_part.GetNode(3).SetSolutionStepValue(VELOCITY_Y,0,1.0)
        model_part.GetNode(3).SetSolutionStepValue(VELOCITY_Y,1,1.0)
        # process_info
        model_part.ProcessInfo[TIME] = self.getDeltaTime()
        model_part.ProcessInfo[DELTA_TIME] = self.getDeltaTime()
        model_part.ProcessInfo[DOMAIN_SIZE] = self.getDomainSize()
        model_part.ProcessInfo[OSS_SWITCH] = 0
        model_part.ProcessInfo[DYNAMIC_TAU] = 1.0
        # sub model parts
        model_part.CreateSubModelPart("structure")
        model_part.GetSubModelPart("structure").AddNode(model_part.GetNode(1),0)
        model_part.GetSubModelPart("structure").AddNode(model_part.GetNode(3),0)
        # flags
        for node in model_part.GetSubModelPart("structure").Nodes:
            node.Set(STRUCTURE,True)
        return model_part

    def createAdjointModelPart(self):
        # create test model part
        model_part = ModelPart("test")
        self.addVariables(model_part)
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        prop = model_part.GetProperties()[0]
        model_part.CreateNewElement("VMSAdjointElement2D", 1, [1, 2, 3], prop)
        buffer_size = self.getBufferSize()
        model_part.SetBufferSize(buffer_size)
        # fluid properties
        for step in range(buffer_size):
            self.assignConstantSolutionStepValue(model_part,DENSITY,step,1.0)
            self.assignConstantSolutionStepValue(model_part,VISCOSITY,step,1.0e-5)
        # dofs
        for node in model_part.Nodes:
            node.AddDof(ADJOINT_VELOCITY_X)
            node.AddDof(ADJOINT_VELOCITY_Y)
            node.AddDof(ADJOINT_VELOCITY_Z)
            node.AddDof(ADJOINT_PRESSURE)
        # bcs
        model_part.GetNode(1).Fix(ADJOINT_VELOCITY_X)
        model_part.GetNode(1).Fix(ADJOINT_VELOCITY_Y)
        model_part.GetNode(2).Fix(ADJOINT_VELOCITY_X)
        model_part.GetNode(2).Fix(ADJOINT_VELOCITY_Y)
        model_part.GetNode(3).Fix(ADJOINT_VELOCITY_X)
        model_part.GetNode(3).Fix(ADJOINT_VELOCITY_Y)
        model_part.GetNode(2).Fix(ADJOINT_PRESSURE)
        # process_info
        model_part.ProcessInfo[TIME] = self.getDeltaTime()
        model_part.ProcessInfo[DELTA_TIME] =-self.getDeltaTime()
        model_part.ProcessInfo[DOMAIN_SIZE] = self.getDomainSize()
        model_part.ProcessInfo[OSS_SWITCH] = 0
        model_part.ProcessInfo[DYNAMIC_TAU] = 1.0
        # sub model parts
        model_part.CreateSubModelPart("structure")
        model_part.GetSubModelPart("structure").AddNode(model_part.GetNode(1),0)
        model_part.GetSubModelPart("structure").AddNode(model_part.GetNode(3),0)
        model_part.CreateSubModelPart("boundary")
        model_part.GetSubModelPart("boundary").AddNode(model_part.GetNode(1),0)
        # flags
        for node in model_part.GetSubModelPart("structure").Nodes:
            node.Set(STRUCTURE,True)
        for node in model_part.GetSubModelPart("boundary").Nodes:
            node.Set(BOUNDARY,True)
        return model_part

    def calculateDrag(self,model_part):
        # force in x-direction on STRUCTURE nodes
        drag = 0.0
        for node in model_part.GetSubModelPart("structure").Nodes:
            drag = drag - node.GetSolutionStepValue(REACTION_X)
        print('drag = ', drag)
        return drag

    def getDragFlagVector(self):
        DragFlagVector = Vector(9)
        for i in range(9):
            DragFlagVector[i] = 0.0
        i = 0
        for node in self.adjoint_model_part.Nodes:
            if node.Is(STRUCTURE):
                DragFlagVector[i] = 1.0 # x-direction
            i = i + 3
        return DragFlagVector

    def invert2x2(self,m):
        invm = Matrix(2,2)
        idet = 1.0 / (m[0,0] * m[1,1] - m[1,0] * m[0,1])
        invm[0,0] = idet * m[1,1]
        invm[0,1] =-idet * m[0,1]
        invm[1,0] =-idet * m[1,0]
        invm[1,1] = idet * m[0,0]
        return invm

    def solve2x2(self,LHS,RHS):
        invLHS = self.invert2x2(LHS)
        x = invLHS * RHS
        return x

    def transposeMatrix(self,m):
        mt = Matrix(m.Size2(),m.Size1())
        for i in range(m.Size1()):
            for j in range(m.Size2()):
                mt[j,i] = m[i,j]
        return mt

    def copyNodalSolutionStepValues(self, origin_model_part, destination_model_part, var):
        if var == ACCELERATION:
            # for consistency with adjoint element, the weighted acceleration is stored in the current
            # time step. for this test case, the acceleration is always zero so this has no effect.
            alpha_bossak = self.getAlphaBossak()
            for origin_node in origin_model_part.Nodes:
                destination_node = destination_model_part.GetNode(origin_node.Id)
                val = (1.0 - alpha_bossak) * origin_node.GetSolutionStepValue(var,0) + alpha_bossak * origin_node.GetSolutionStepValue(var,1)
                destination_node.SetSolutionStepValue(var,0,val)
        else:
            for origin_node in origin_model_part.Nodes:
                destination_node = destination_model_part.GetNode(origin_node.Id)
                val = origin_node.GetSolutionStepValue(var,0)
                destination_node.SetSolutionStepValue(var,0,val)

    def calculateAdjointFromElement(self,adjoint_element):
        ft = [2,8]
        alpha_bossak = self.getAlphaBossak()
        gamma = 0.5 - alpha_bossak
        delta_time = self.model_part.ProcessInfo[DELTA_TIME]
        ElemMassMatrix0 = adjoint_element.Calculate(MASS_MATRIX_0,self.adjoint_model_part.ProcessInfo)
        ElemAdjointMatrix1 = adjoint_element.Calculate(ADJOINT_MATRIX_1,self.adjoint_model_part.ProcessInfo)
        ElemAdjointSystemMatrix = Matrix(9,9)
        for i in range(9):
            for j in range(9):
                ElemAdjointSystemMatrix[i,j] = ElemAdjointMatrix1[i,j] - (1.0 - alpha_bossak) * ElemMassMatrix0[j,i] / (gamma * delta_time)
        AdjointSystemMatrix = Matrix(2,2) # system for the 2 free pressure dofs
        for i in range(2):
            for j in range(2):
                AdjointSystemMatrix[i,j] = ElemAdjointSystemMatrix[ft[i],ft[j]]
        InvAdjointSystemMatrix = self.invert2x2(AdjointSystemMatrix)
        DragFlagVector = self.getDragFlagVector()
        ElemGradientOfDrag = ElemAdjointSystemMatrix * DragFlagVector
        MinusGradientOfDrag = Vector(2)
        for i in range(2):
            MinusGradientOfDrag[i] = -ElemGradientOfDrag[ft[i]]
        Adjoint = InvAdjointSystemMatrix * MinusGradientOfDrag
        return Adjoint

    def setUp(self):
        self.model_part = self.createModelPart()
        self.adjoint_model_part = self.createAdjointModelPart()
        domain_size = self.getDomainSize()
        alpha_bossak = self.getAlphaBossak()
        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            alpha_bossak,0,domain_size)
        self.linear_solver = SkylineLUFactorizationSolver()
        self.builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        # linear strategy is used when all velocity dofs are fixed
        self.solver = ResidualBasedLinearStrategy(self.model_part,
                                                  self.time_scheme,
                                                  self.linear_solver,
                                                  self.builder_and_solver,
                                                  True, False, False, False)

        #self.convergence_criteria = VelPrCriteria(1e-3,1e-5,1e-3,1e-5)
        #solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,time_scheme,linear_solver,convergence_criteria,10,True,False,False)
        self.solver.SetEchoLevel(0)
        self.solver.Check()

        objective_params = Parameters("""{
            "structure_model_part_name": "structure",
            "drag_direction": [1.0, 0.0, 0.0]
        }""")
        scheme_params = Parameters("""{
            "boundary_model_part_name": "boundary",
            "alpha_bossak": -0.3,
            "adjoint_start_time": 0.0,
            "adjoint_end_time": 1.0
        }""")
        self.objective_function = DragObjectiveFunction2D(objective_params)
        self.adjoint_time_scheme = AdjointBossakScheme(scheme_params, self.objective_function)
        self.adjoint_builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.adjoint_solver = ResidualBasedLinearStrategy(self.adjoint_model_part,
                                                          self.adjoint_time_scheme,
                                                          self.linear_solver,
                                                          self.adjoint_builder_and_solver,
                                                          False, False, False, False)

        self.adjoint_solver.SetEchoLevel(0)
        self.adjoint_solver.Check()
        self.adjoint_solver.Initialize()

    def test_PrimalGradient(self):
        # this is an intermediate test, useful for debugging.
        # the goal is to verify that the primal gradient obtained from the adjoint element
        # is the same as the primal gradient estimated from the fluid solver.

        alpha_bossak = self.getAlphaBossak()
        gamma = 0.5 - alpha_bossak
        delta_time = self.model_part.ProcessInfo[DELTA_TIME]
        ft = [2,8] # local indices for element matrices corresponding to free pressure dofs
        h = 0.00000001 # finite difference step size

        # solve
        self.solver.Solve()
        # copy solution to adjoint model part
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,VELOCITY)
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,PRESSURE)
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,ACCELERATION)

        # do finite differencing of pressure
        node1 = self.model_part.GetNode(1)
        # unperturbed
        P1 = self.model_part.GetNode(1).GetSolutionStepValue(PRESSURE,0)
        P3 = self.model_part.GetNode(3).GetSolutionStepValue(PRESSURE,0)
        # x-perturbation
        x = node1.X
        node1.X = x + h
        self.solver.Solve()
        FDdP1dX1 = (self.model_part.GetNode(1).GetSolutionStepValue(PRESSURE,0) - P1) / h
        FDdP3dX1 = (self.model_part.GetNode(3).GetSolutionStepValue(PRESSURE,0) - P3) / h
        node1.X = x
        # y-perturbation
        y = node1.Y
        node1.Y = y + h
        self.solver.Solve()
        FDdP1dY1 = (self.model_part.GetNode(1).GetSolutionStepValue(PRESSURE,0) - P1) / h
        FDdP3dY1 = (self.model_part.GetNode(3).GetSolutionStepValue(PRESSURE,0) - P3) / h
        node1.Y = y

        # use adjoint element
        adjoint_element = self.adjoint_model_part.GetElement(1)
        ElemMassMatrix0 = adjoint_element.Calculate(MASS_MATRIX_0,self.adjoint_model_part.ProcessInfo)
        ElemAdjointMatrix1 = adjoint_element.Calculate(ADJOINT_MATRIX_1,self.adjoint_model_part.ProcessInfo)
        ElemAdjointMatrix1 = self.transposeMatrix(ElemAdjointMatrix1)
        LHS = Matrix(2,2) # \partial f / \partial(p1, p3)
        for i in range(2):
            for j in range(2):
                LHS[i,j] = ElemAdjointMatrix1[ft[i],ft[j]] - (1.0 - alpha_bossak) * ElemMassMatrix0[ft[i],ft[j]] / (gamma * delta_time)
        RHS = Vector(2)
        ElemShapeDerivativeMatrix1 = adjoint_element.Calculate(SHAPE_DERIVATIVE_MATRIX_1,self.adjoint_model_part.ProcessInfo)
        ElemShapeDerivativeMatrix1 = self.transposeMatrix(ElemShapeDerivativeMatrix1)
        RHS[0] =-ElemShapeDerivativeMatrix1[ft[0],0] # -\partial f / \partial x1
        RHS[1] =-ElemShapeDerivativeMatrix1[ft[1],0]
        x = self.solve2x2(LHS,RHS)
        dP1dX1 = x[0]
        dP3dX1 = x[1]
        RHS[0] =-ElemShapeDerivativeMatrix1[ft[0],1] # -\partial f / \partial y1
        RHS[1] =-ElemShapeDerivativeMatrix1[ft[1],1]
        x = self.solve2x2(LHS,RHS)
        dP1dY1 = x[0]
        dP3dY1 = x[1]

        # compare analytical pressure gradient with finite difference
        self.assertNotEqual(dP1dX1, 0.0) # prevent silent failure
        self.assertAlmostEqual(dP1dX1, FDdP1dX1,6)
        self.assertNotEqual(dP1dY1, 0.0) 
        self.assertAlmostEqual(dP1dY1, FDdP1dY1,6)
        self.assertNotEqual(dP3dX1, 0.0) 
        self.assertAlmostEqual(dP3dX1, FDdP3dX1,6)
        self.assertNotEqual(dP3dY1, 0.0) 
        self.assertAlmostEqual(dP3dY1, FDdP3dY1,6)

    def test_ElementSensitivity(self):
        # this is an intermediate test, useful for debugging.
        # the goal is to verify that the sensitivity obtained from the adjoint element
        # is the same as the sensitivity estimated from the fluid solver. it directly
        # tests the adjoint element without dependencies on scheme, strategy, etc.
        # test_PrimalGradient should work before this test is debugged.

        ft = [2,8] # local indices for element matrices corresponding to free pressure dofs
        h = 0.00000001 # finite difference step size

        # solve
        self.solver.Solve()
        # copy solution to adjoint model part
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,VELOCITY)
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,PRESSURE)
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,ACCELERATION)

         # do finite differencing of drag
        FDSensitivity = Vector(2)
        node1 = self.model_part.GetNode(1)
        # unperturbed
        drag0 = self.calculateDrag(self.model_part)
        # x-perturbation
        x = node1.X
        node1.X = x + h
        self.solver.Solve()
        node1.X = x
        FDSensitivity[0] = (self.calculateDrag(self.model_part) - drag0) / h
        # y-perturbation
        y = node1.Y
        node1.Y = y + h
        self.solver.Solve()
        node1.Y = y
        FDSensitivity[1] = (self.calculateDrag(self.model_part) - drag0) / h
        
        # use adjoint element
        adjoint_element = self.adjoint_model_part.GetElement(1)
        Adjoint = self.calculateAdjointFromElement(adjoint_element)
        DragFlagVector = self.getDragFlagVector()
        ElemShapeDerivativeMatrix1 = adjoint_element.Calculate(SHAPE_DERIVATIVE_MATRIX_1,self.adjoint_model_part.ProcessInfo)
        ShapeDerivativeMatrix1 = Matrix(2,2)
        for i in range(2):
            for j in range(2):
                ShapeDerivativeMatrix1[i,j] = ElemShapeDerivativeMatrix1[i,ft[j]]
        ElemPartialGradientOfObjective = ElemShapeDerivativeMatrix1 * DragFlagVector
        PartialGradientOfObjective = Vector(2)
        for i in range(2):
            PartialGradientOfObjective[i] = ElemPartialGradientOfObjective[i]
        Sensitivity = PartialGradientOfObjective + (ShapeDerivativeMatrix1 * Adjoint)
        
        # compare analytical sensitivity with finite difference
        for i in range(2):
            self.assertNotEqual(Sensitivity[i], 0.0) # prevent silent failure
            self.assertAlmostEqual(Sensitivity[i], FDSensitivity[i])

    def test_AdjointBossakDragScheme(self):
        # the goal of this test is to verify that the scheme and strategy correctly solve
        # the adjoint equation of the first time step.
        # test_ElementSensitivity should work before this test is debugged.

        # solve
        self.solver.Solve()
        # copy solution to adjoint model part
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,VELOCITY)
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,PRESSURE)
        self.copyNodalSolutionStepValues(self.model_part,self.adjoint_model_part,ACCELERATION)

        # use AdjointBossakDragScheme
        self.adjoint_solver.Initialize()
        self.adjoint_solver.Solve()
        Adjoint = Vector(2)
        Adjoint[0] = self.adjoint_model_part.GetNode(1).GetSolutionStepValue(ADJOINT_PRESSURE,0)
        Adjoint[1] = self.adjoint_model_part.GetNode(3).GetSolutionStepValue(ADJOINT_PRESSURE,0)
        
        # compare with element calculation (independent of scheme)
        adjoint_element = self.adjoint_model_part.GetElement(1)
        ElemAdjoint = self.calculateAdjointFromElement(adjoint_element)

        for i in range(2):
            self.assertNotEqual(Adjoint[i], 0.0)
            self.assertAlmostEqual(Adjoint[i], ElemAdjoint[i])

        Sensitivity = Vector(2)
        Sensitivity[0] = self.adjoint_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_X,0)
        Sensitivity[1] = self.adjoint_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_Y,0)

         # do finite differencing of drag
        h = 0.00000001 # finite difference step size
        FDSensitivity = Vector(2)
        node1 = self.model_part.GetNode(1)
        # unperturbed
        drag0 = self.calculateDrag(self.model_part)
        # x-perturbation
        x = node1.X
        node1.X = x + h
        self.solver.Solve()
        node1.X = x
        FDSensitivity[0] = (self.calculateDrag(self.model_part) - drag0) / h
        # y-perturbation
        y = node1.Y
        node1.Y = y + h
        self.solver.Solve()
        node1.Y = y
        FDSensitivity[1] = (self.calculateDrag(self.model_part) - drag0) / h
        
        # compare analytical sensitivity with finite difference
        for i in range(2):
            self.assertNotEqual(Sensitivity[i], 0.0)
            self.assertAlmostEqual(Sensitivity[i], FDSensitivity[i])

if __name__ == '__main__':
    KratosUnittest.main()
