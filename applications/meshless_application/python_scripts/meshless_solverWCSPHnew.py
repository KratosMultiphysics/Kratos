#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.MeshlessApplication import *


def AddVariables(model_part):  #this way we only need one command to add all the variables to our problem
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(RADIUS);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(SEARCH_RADIUS);
    model_part.AddNodalSolutionStepVariable(DENSITY_NORM_PARAM);
    model_part.AddNodalSolutionStepVariable(DIV_OF_VEL);
    model_part.AddNodalSolutionStepVariable(XSPH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(DENS_VARIATION);
    model_part.AddNodalSolutionStepVariable(EFFECTIVE_RADIUS);
    model_part.AddNodalSolutionStepVariable(BODYFORCE_ACC);
    model_part.AddNodalSolutionStepVariable(VISCOUS_ACC);
    model_part.AddNodalSolutionStepVariable(PRESSURE_ACC);
    model_part.AddNodalSolutionStepVariable(BOUNDARY_ACC);

    model_part.AddNodalSolutionStepVariable(DUMMY_NORMALIZE_RHS);
    model_part.AddNodalSolutionStepVariable(DUMMY_APPLY_XSPH);
    model_part.AddNodalSolutionStepVariable(DUMMY_BOUNDARY_PRESSURES);

    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE );


    model_part.AddNodalSolutionStepVariable(ERASE_FLAG );
    model_part.AddNodalSolutionStepVariable(IS_WET );
    model_part.AddNodalSolutionStepVariable(INI_PRESSURE );
    model_part.AddNodalSolutionStepVariable(OLD_VEL );
    model_part.AddNodalSolutionStepVariable(RHS );
    model_part.AddNodalSolutionStepVariable(DELTA_TIME );
    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY);





def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);

    print "variables for the Meshless solver added correctly"

class MeshlessSolverWCSPHnew:
    #######################################################################
    def __init__(self,model_part,domain_size,box_corner1,box_corner2):
        self.model_part = model_part
        self.domain_size = domain_size
        #~ self.soundspeed = soundspeed
        #~ self.ReferenceDensity = ReferenceDensity
        #~ self.ini_spacing = ini_spacing
        #~ self.boundary_force_choice = boundary_force_choice
        #~ self.state_equation_choice = state_equation_choice
        #~ self.gama = gama
        #~ self.kernel_choice_for_density = kernel_choice_for_density;
        #~ self.kernel_choice_for_velocity = kernel_choice_for_velocity;
        #~ self.kernel_choice_for_viscosity = kernel_choice_for_viscosity;
        #~ self.kernel_choice_for_pressure = kernel_choice_for_pressure;

        ##saving the limits of the box (all the nodes external to this will be erased)
        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2


      ##  self.EraseNodes = NodeAndElementEraseProcess(model_part);

        self.SPHTimeIntegrator = SPHTimeIntegrator()


         #######################################################################
    def Initialize(self):
        #creating the solution strategy
        print "Initializing Meshless solver solver"
        (self.SPHTimeIntegrator).Initialize(self.model_part)


           #######################################################################
    def MarkOutsideNodes(self):
        #creating the solution strategy
        print "Initializing Meshless solver"
        (self.SPHTimeIntegrator).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );

      #######################################################################
    def EstimateDeltaTime(self,time_step):
		print "Checking the time step"

		dt = (self.SPHTimeIntegrator).EstimateDeltaTime(time_step,self.model_part)

		return dt

     #######################################################################
    def Solve(self):
		print "just before solve"
		(self.SPHTimeIntegrator).Solve(self.model_part,self.box_corner1,self.box_corner2,(self.model_part).Nodes )

		print "SOLVED!!"


 #######################################################################
    def OutputToGID(self,gid_io,time):
		print "Writing outut to GID"

		#writing mesh
		gid_io.InitializeMesh( time );
		gid_io.WriteNodeMesh((self.model_part).GetMesh());
		gid_io.WriteMesh((self.model_part).GetMesh());
		gid_io.FinalizeMesh();
		gid_io.InitializeResults(time, (self.model_part).GetMesh());
		gid_io.WriteNodalResults(VELOCITY, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(PRESSURE_ACC, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(VISCOUS_ACC, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(BOUNDARY_ACC, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(BODYFORCE_ACC, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(PRESSURE, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(IS_STRUCTURE, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(DENSITY_NORM_PARAM, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(VISCOSITY, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(DENSITY, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(NODAL_MASS, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(IS_FREE_SURFACE, (self.model_part).Nodes, time, 0);
		gid_io.WriteNodalResults(DENS_VARIATION, (self.model_part).Nodes, time, 0);

		gid_io.Flush()
		gid_io.FinalizeResults()
		print "Results Written to GID"

