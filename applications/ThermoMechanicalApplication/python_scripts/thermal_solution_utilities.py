from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *

class ThermalSolutionUtilities:
    def __init__(self, model_part,ProjectParameters):
        self.model_part = model_part
        self.ProjectParameters = ProjectParameters
        
        import detect_solidus_liquidus_temperature
        detect_solidus_liquidus_temperature.DetectSolidusLiquidusTemperature(model_part)
        
        #settings for pure convection of temperature
        self.temperature_convection_settings = ConvectionDiffusionSettings()
        self.temperature_convection_settings.SetUnknownVariable(TEMPERATURE)
        self.temperature_convection_settings.SetConvectionVariable(VELOCITY)
        self.temperature_convection_settings.SetMeshVelocityVariable(MESH_VELOCITY)
        #self.temperature_convection_settings.SetDensityVariable(TEMPERATURE)
        
    def InitializeThermalSimulation(self):
        pass
    
    def BeforeThermalSolution(self,max_distance):
        pass
      
        #mark as inactive the temperature elements and condition with distance greater than zero
        #for elem in self.model_part.Elements:
            #npos = 0
            #for node in elem.GetNodes():
                #if(node.GetSolutionStepValue(DISTANCE) > 0.9*max_distance):
                    #npos = npos + 1    
            #if(npos > 0):
                #elem.Set(ACTIVE,False)
            #else:
                #elem.Set(ACTIVE,True)
            
        #for cond in self.model_part.Elements:
            #npos = 0
            #for node in cond.GetNodes():
                #if(node.GetSolutionStepValue(DISTANCE) > 0.9*max_distance):
                    #npos = npos + 1                  
            #if(npos > 0):
                #cond.Set(ACTIVE,False)                
            #else:
                #cond.Set(ACTIVE,True) 
                
                
        #now detect if there is a node has just been switched on
        #latent_heat = self.model_part.ProcessInfo[LATENT_HEAT]        
        #specific_heat_table = self.model_part.GetTable(2)
        #liquidus = self.model_part.ProcessInfo[FLUID_TEMPERATURE]
        #solidus = self.model_part.ProcessInfo[SOLID_TEMPERATURE]
        #ambient = self.ProjectParameters.AMBIENT_TEMPERATURE
        ##air_temperature = 0.5*(solidus+liquidus) #
        #air_temperature = self.ProjectParameters.FLUID_TEMPERATURE #liquidus #0.5*(solidus + liquidus)
        
        ##speficif_heat = specific_heat_table.GetValue(air_temperature) #self.ProjectParameters.FLUID_TEMPERATURE)
        
        #for node in self.model_part.Nodes:
            #dist = node.GetSolutionStepValue(DISTANCE)
            #if(dist > 0):   
                #T = air_temperature #self.ProjectParameters.FLUID_TEMPERATURE
                ##H = speficif_heat*T + latent_heat
                #for i in range(0,3):
                    #node.SetSolutionStepValue(TEMPERATURE,i,T) ##assign to the past the current temperature
                    ##node.SetSolutionStepValue(ENTHALPY,i,H)
                
    def AfterThermalSolution(self):
        pass
