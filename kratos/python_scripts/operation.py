#this class is designed as the base class to add "custom operations"
#it is designed so that only the constructor can be modified and all of the other 
#functions should retain the same interface
class Operation:
  
    #model_part -> model_part to which the operation will be applied
    #groups -> python array containing the ids of the groups to which the operation should be applied
    #group_container -> container of the groups. Can be queried for the nodes that correspond to a given group Id
    #echo_level -> level of expected echo for the operation: echo_level=0 implies no echo
    def __init__(self,model_part,groups,group_container,echo_level=0):
	self.model_part = model_part
	self.groups = groups
	self.group_container = group_container
	self.echo_level = echo_level
	
    def PrintInfo(self):
       raise Exception ("Method PrintInfo Should be implemented on all derived classes")

    #this function is designed for being called at the beginning of the computations
    #right after reading the model and the groups 
    def ExecuteInitialize(self):
      if(self.echo_level > 0):
	print "Finished ExecuteInitialize for Class",self.PrintInfo()
            
    #this function is designed for being execute once before the solution loop but after all of the 
    #solvers where built
    def ExecuteBeforeSolutionLoop(self):
      if(self.echo_level > 0):
	print "Finished ExecuteBeforeSolutionLoop for Class",self.PrintInfo()
      
    #this function will be executed at every time step BEFORE performing the solve phase
    def ExecuteInitializeSolutionStep(self):
      if(self.echo_level > 0):
	print "Finished ExecuteInitializeSolutionStep for Class",self.PrintInfo()
      
    #this function will be executed at every time step AFTER performing the solve phase
    def ExecuteFinalizeSolutionStep(self):
      if(self.echo_level > 0):
	print "Finished ExecuteFinalizeSolutionStep for Class",self.PrintInfo()

    #this function will be executed at every time step BEFORE  writing the output
    def ExecuteBeforeOutputStep(self):
      if(self.echo_level > 0):
	print "Finished ExecuteBeforeOutputStep for Class",self.PrintInfo()
      
    #this function will be executed at every time step AFTER writing the output
    def ExecuteAfterOutputStep(self):
      if(self.echo_level > 0):
	print "Finished ExecuteAfterOutputStep for Class",self.PrintInfo()
      
    #this function is designed for being called at the beginning of the computations
    #right after reading the model and the groups 
    def ExecuteFinalize(self):      
      if(self.echo_level > 0):
	print "Finished ExecuteFinalize for Class",self.PrintInfo()
    

	
  