#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(POROSITY)

#    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(DIAMETER)

    print "variables for the edgebased incompressible fluid solver added correctly"

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);


class EdgeBasedLevelSetSolver:
    
    def __init__(self,model_part,domain_size,body_force,viscosity,density):

        print "entered in EdgeBasedLevelSetSolver python constructor"
        #data of the problem
        self.model_part = model_part
        self.domain_size = domain_size
        self.body_force = body_force
        self.density = density
        self.viscosity = viscosity

        self.use_mass_correction = True

        self.redistance_frequency = 5;
        self.step = 0;

        self.extrapolation_layers = 5

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()

        #definition of the solvers
#        pDiagPrecond = DiagonalPreconditioner()
#        self.pressure_linear_solver =  CGSolver(1e-3, 5000,pDiagPrecond)
#        self.pressure_linear_solver =  CGSolver(1e-3, 5000)

        pDiagPrecond = DiagonalPreconditioner()
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
#        self.pressure_linear_solver =  BICGSTABSolver(1e-6, 5000)

        ##initializing the press proj to -body_force
        press_proj_init = Vector(3)
        press_proj_init[0] = body_force[0]*density
        press_proj_init[1] = body_force[1]*density
        press_proj_init[2] = body_force[2]*density
        for node in self.model_part.Nodes:
            eps = node.GetSolutionStepValue(POROSITY)
            node.SetSolutionStepValue(PRESS_PROJ,0,press_proj_init*eps)
        print "entered in EdgeBasedLevelSetSolver initialize"



    def Initialize(self):
        print "entered in EdgeBasedLevelSetSolver python constructor"
        #build the edge data structure
        if(self.domain_size == 2):
            self.matrix_container = MatrixContainer2D()
        else:
            self.matrix_container = MatrixContainer3D()

        self.matrix_container.ConstructCSRVector(self.model_part)
        self.matrix_container.BuildCSRData(self.model_part)

        ##for 3D problems we need to evaluate the condition's neighbours
        if(self.domain_size == 3):
            self.condition_neighbours_finder = FindConditionsNeighboursProcess(self.model_part,self.domain_size,10)
            self.condition_neighbours_finder.Execute()

        ##constructing the solver
        if(self.domain_size == 2):
            self.distance_utils = SignedDistanceCalculationUtils2D()
            self.fluid_solver = EdgeBasedLevelSet2D(self.matrix_container,self.model_part,self.viscosity,self.density,self.body_force,self.use_mass_correction)
        else:
            self.distance_utils = SignedDistanceCalculationUtils3D()
#            self.distance_utils = SignedDistanceCalculationBinBased3D()
            self.fluid_solver = EdgeBasedLevelSet3D(self.matrix_container,self.model_part,self.viscosity,self.density,self.body_force,self.use_mass_correction)

#
        self.max_edge_size = self.distance_utils.FindMaximumEdgeSize(self.model_part)
        self.distance_size = self.max_edge_size * 3.0;
        print "###################### max distance = ",self.distance_size;
        
#        self.reorder = True
#        self.distance_tools = BodyDistanceCalculationUtils()

#TO BE DELETED SOON
##        print "ASSIGNING DIAMETER FROM SOLVER: TO BE REMOVED SOOOOOON!!!!!!!!!!!!!!!!!!!!"
        for node in self.model_part.Nodes:
##            if(node.GetSolutionStepValue(POROSITY) != 1.0):
            node.SetSolutionStepValue(DIAMETER,0, 0.01)

        self.fluid_solver.Initialize()

#        for node in self.model_part.Nodes:
#            dist = node.GetSolutionStepValue(DISTANCE)
#            node.SetSolutionStepValue(DISTANCE,1,dist)
#        self.Redistance()


        print "**********************************************"
        print "finished EdgeBasedLevelSetSolver initialize"

    ################################################################
    ################################################################
#    def CalculateDistances(self):
#        if(self.domain_size == 2):
#            self.distance_tools.CalculateDistances2D(self.model_part.Elements,DISTANCE, self.reorder);
#        else:
#            self.distance_tools.CalculateDistances3D(self.model_part.Elements,DISTANCE, self.reorder);

    ################################################################
    ################################################################
    def Redistance(self):
        self.distance_utils.CalculateDistances(self.model_part,DISTANCE,self.distance_size)
#        self.distance_utils.CalculateDistances(self.model_part,DISTANCE,self.distance_size)


#        self.fluid_solver.MarkInternalNodes()
#        self.CalculateDistances()
#        self.fluid_solver.MarkExternalAndMixedNodes()
#        self.fluid_solver.ChangeSignToDistance()
#        self.CalculateDistances()
#        self.fluid_solver.ChangeSignToDistance()
        

    ################################################################
    ################################################################
    def FluidOnlySolve(self):
        print "entered in EdgeBasedLevelSetSolver fluid only solve"
        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)

        (self.fluid_solver).SolveStep1();
        (self.fluid_solver).SolveStep2(self.pressure_linear_solver);
        (self.fluid_solver).SolveStep3();

        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
        print "finished EdgeBasedLevelSetSolver fluid only solve"
   
    ################################################################
    ################################################################
    def Solve(self):

        (self.fluid_solver).UpdateFixedVelocityValues()

        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)

        ##convect levelset function
       # self.convection_solver.Solve();
        (self.fluid_solver).ConvectDistance()

        ##solve fluid
        (self.fluid_solver).SolveStep1();
        (self.fluid_solver).SolveStep2(self.pressure_linear_solver);
        (self.fluid_solver).SolveStep3();

        if(self.step == self.redistance_frequency):
            self.Redistance()
            self.step = 0
            print "redistance was executed"
        self.step += 1

    ################################################################
#    ################################################################
#    def Solve(self):
#        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
#
#        ##convect levelset function
#        (self.fluid_solver).ConvectDistance()
#
#        convection_success = (self.fluid_solver).CheckDistanceConvection()
#        if(convection_success == False):
#            #time step reduction is needed
#            print "############### distance convection failed!! ###############"
#            return False
#            errrrrrrrrr
#
#        ##solve fluid
#        (self.fluid_solver).SolveStep1();
#        (self.fluid_solver).SolveStep2(self.pressure_linear_solver);
#        (self.fluid_solver).SolveStep3();
#
#        if(self.step == self.redistance_frequency):
#            self.Redistance()
#            self.step = 0
#            print "redistance was executed"
#        self.step += 1
#
#        return True

    ################################################################
    ################################################################
    def EstimateTimeStep(self,safety_factor,max_Dt):
        dt = (self.fluid_solver).ComputeTimeStep(safety_factor,max_Dt);

        if(dt > max_Dt):
            dt = max_Dt

        print dt

        return dt

    ################################################################
    ################################################################
    def CalculateInitialPressureDistribution(self):
        ##prova!
        dt_aux = 1e-6
        self.model_part.ProcessInfo.SetValue(DELTA_TIME,dt_aux)
        (self.fluid_solver).SolveStep1();
        aaa = Vector(3)
        aaa[0] = self.body_force[0]*dt_aux
        aaa[1] = self.body_force[1]*dt_aux
        aaa[2] = self.body_force[2]*dt_aux
        for node in self.model_part.Nodes:
            if(node.IsFixed(VELOCITY_X) == False):
                node.SetSolutionStepValue(VELOCITY,0,aaa)

#        for node in self.model_part.Nodes:
#            print node.GetSolutionStepValue(VELOCITY)

        (self.fluid_solver).SolveStep2(self.pressure_linear_solver);

        zero = Vector(3); zero[0]=0.0; zero[1]=0.0; zero[2]=0.0;
        for node in self.model_part.Nodes:
            if(node.IsFixed(VELOCITY_X) == False):
                node.SetSolutionStepValue(VELOCITY,0,zero)
        self.model_part.ProcessInfo.SetValue(DELTA_TIME,0.0)



        
        

