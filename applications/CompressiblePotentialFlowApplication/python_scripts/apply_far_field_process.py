import KratosMultiphysics
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "inlet_phi": 1.0,
                "velocity_infinity": [1.0,0.0,0],
                "density_infinity"  : 1.225,
                "gamma"                 : 1.4,
                "speed_sound_infinity"  : 340,
                "AOAdeg" : 0
            }  """ );
        
            
        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.model_part             = Model[settings["model_part_name"].GetString()]
        self.inlet_phi              = settings["inlet_phi"].GetDouble()
        self.velocity_infinity      = KratosMultiphysics.Vector(3)#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0]   = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1]   = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2]   = settings["velocity_infinity"][2].GetDouble()            
        self.density_infinity       = settings["density_infinity"].GetDouble()
        #self.density_infinity      = settings["density_infinity"].GetDouble() #TODO: must read this from the properties
        self.gamma                  = settings["gamma"].GetDouble()
        self.speed_sound_infinity   = settings["speed_sound_infinity"].GetDouble()
        self.AOAdeg                 = settings["AOAdeg"].GetDouble()
        

        #convert angle to radians
        self.AOArad = self.AOAdeg*pi/180
        
        #rotate velocity vector according to the angle of attack
        self.velocity_infinity[0]   = settings["velocity_infinity"][0].GetDouble()*cos(self.AOArad)
        self.velocity_infinity[2]   = settings["velocity_infinity"][0].GetDouble()*sin(self.AOArad)
        self.velocity_infinity[1]   = settings["velocity_infinity"][2].GetDouble()
                
        from numpy import linalg as LA
        velocity_norm = LA.norm(self.velocity_infinity)
        
        self.mach                   = velocity_norm/self.speed_sound_infinity
        
        print('\nAngle of attack   =  \t'   , self.AOAdeg)
        print('\nvelocity_infinity =  \t'   , self.velocity_infinity)
        print('\nvelocity_norm     =  \t'   , velocity_norm)
        print('\nmach =               \t'   , self.mach)
        print('\ndensity_infinity =   \t'   , self.density_infinity)
        print('\nspeed_sound_infinity =\t'  , self.speed_sound_infinity)
        print('\ngamma =              \t'   , self.gamma)
        
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.VELOCITY,       self.velocity_infinity)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DENSITY,        self.density_infinity)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.GAMMA,          self.gamma)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.SOUND_VELOCITY, self.speed_sound_infinity)
        
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,0)
        
        print(self.model_part.ProcessInfo.GetValue(KratosMultiphysics.IS_RESTARTED))
        
        if(self.model_part.ProcessInfo.GetValue(KratosMultiphysics.IS_RESTARTED)):
            print('\nComputing Compressible Flow\n')
        else:
            print('\nComputing Incompressible Flow\n')
        
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part,self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])   

    def Execute(self):
        for cond in self.model_part.Conditions:
            cond.SetValue(KratosMultiphysics.VELOCITY, self.velocity_infinity)

        #select the first node
        for node in self.model_part.Nodes:
            node1 = node
            break
        
        #find the node with the minimal x
        x0 = node1.X
        y0 = node1.Y
        z0 = node1.Z

        pos = 1e30
        for node in self.model_part.Nodes:
            #node.Set(KratosMultiphysics.STRUCTURE)#setting structure fla
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0
            
            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
            
            if(tmp < pos):
                pos = tmp

        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0
            
            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
            
            if(tmp < pos+1e-9):
                node.Fix(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
                node.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,self.inlet_phi)
        
    def ExecuteInitializeSolutionStep(self):
        self.Execute()
        
        