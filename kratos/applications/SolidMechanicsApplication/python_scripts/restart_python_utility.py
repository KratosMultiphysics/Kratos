#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()

class RestartUtility:
    #######################################################################
    def __init__(self,model_part,problem_path,problem_name,restart_step):

        self.model_part = model_part
        
        #set restart flags
        self.load_restart_flag = False
        self.save_restart_flag = False

        #set problem name
        self.problem_name = problem_name

        #set problem path
        self.problem_path = problem_path
        
        #set restart step
        self.restart_step = restart_step

        #set restart path
        self.restart_path = self.problem_path + "/" + self.problem_name + "_" + str(restart_step)
       
        #set serializer 
        self.serializer = Serializer(self.restart_path,SERIALIZER_NO_TRACE)
        #self.serializer = Serializer(self.restart_path,SERIALIZER_TRACE_ERROR)

        #save restart variables
        self.step_to_restart  = 1
        self.restart_interval = 10

        #time variables
        self.start_steps  = 0
        self.steps_number = 0
        self.time_step    = 0
        self.start_time   = 0

        self.step_i       = 0
        self.step_n       = 0
        self.current_step = 0

    #######################################################################
    def Initialize(self,load_restart,save_restart,restart_interval,main_solver,list_files,rotation_dofs):

         if(load_restart == "True"):
             self.load_restart_flag = True
             self.serializer.Load("ModelPart",self.model_part,rotation_dofs)
             self.CleanPosteriorFiles(list_files)
         else:
             self.load_restart_flag = False
             main_solver.AddVariables(self.model_part,rotation_dofs)
             #self.CleanPreviousFiles(list_files) #remove these files in the problem type script
            
         if(save_restart == "True"):
             self.save_restart_flag = True
             self.restart_interval  = restart_interval
             if(self.load_restart_flag == True):
                 #re-write the restart if it was loaded
                 self.step_to_restart = 0
             else:
                 self.step_to_restart = 1

    #######################################################################   
    def CleanPreviousFiles(self,list_files):
        print "Start: -remove previous problem files-"
        #remove previous results:
        filelist1 = [ f for f in os.listdir(self.problem_path) if f.endswith(".bin") ]
        for f in filelist1:
            os.remove(f)

        filelist1_2 = [ f for f in os.listdir(self.problem_path) if f.endswith(".res") ]
        for f in filelist1_2:
            os.remove(f)

        filelist1_3 = [ f for f in os.listdir(self.problem_path) if f.endswith(".msh") ]
        for f in filelist1_3:
            os.remove(f)

            
        #remove previous restart files:
        filelist3 = [ f for f in os.listdir(self.problem_path) if f.endswith(".rest") ]
        for f in filelist3:
            os.remove(f)
            
        #remove previous graph files:
        filelist4 = [ f for f in os.listdir(self.problem_path) if f.endswith(".png") ]
        for f in filelist4:
            os.remove(f)


        #remove previous list files:
        list_files.RemoveListFiles();

    #######################################################################   
    def CleanPosteriorFiles(self,list_files):
        print "Restart: -remove restart step posterior problem files-"

        #remove previous results after restart:
        total_files = 0;
        filelist1 = []
        filelist2 = []
        for f in os.listdir(self.problem_path):
            if(f.endswith(".post.bin")):
                total_files = total_files + 1
                
        total_files = total_files - self.restart_step + 250  #arbitrary number to ensure to remove all rest files

        for rfile in range(0,total_files):
            for f in os.listdir(self.problem_path):
                if(f.endswith("_"+str(self.restart_step+rfile+1)+".post.bin")):
                    filelist1.append(f);
                if(f.endswith("_"+str(self.restart_step+rfile+1)+".png")):
                    filelist2.append(f);
                                
                                
        for f in filelist1:
            os.remove(f)

        for f in filelist2:
            os.remove(f)


        #remove previous restart files:
        counter = 1;
        total_files = 0;
        filelist2 = []
        for f in os.listdir(self.problem_path):
            if(f.endswith(".post.bin")):
                total_files = total_files + 1
                                                
        total_files = total_files + 100  #arbitrary number to ensure to remove all rest files
                                                
        for rfile in range(0,total_files):
            for f in os.listdir(self.problem_path):
                if(f.endswith("_"+str(self.restart_step+rfile+1)+".rest")):
                    filelist2.append(f);
                                                            
        for f in filelist2:
            os.remove(f)

        #remove previous list files and rebuild it:
        list_files.RemoveListFiles();
        list_files.BuildListFiles();
    

    #######################################################################   
    def PrintRestartFile(self,write_id):

        if( self.save_restart_flag == True and self.step_to_restart == self.restart_interval ): 
            self.restart_path= self.problem_path + "/" + self.problem_name + "_" + str(write_id-1)
            #self.serializer = Serializer(self.restart_path,SERIALIZER_TRACE_ALL);
            #self.serializer = Serializer(self.restart_path,SERIALIZER_TRACE_ERROR);
            self.serializer = Serializer(self.restart_path,SERIALIZER_NO_TRACE);        
            self.serializer.Save("ModelPart",self.model_part);
            self.serializer = 0;
            print " RESTART PRINTED :", write_id-1;
            self.step_to_restart = 1;
        else:
            self.step_to_restart = self.step_to_restart + 1;


    #######################################################################   
    def StartModelRead(self,buffer_size,domain_size,problem_type,rotation_dofs,main_solver,modeler):
        
        if(self.load_restart_flag == False):
            #reading the model 
            self.model_part_io = ModelPartIO(self.problem_name)
            self.model_part_io.ReadModelPart(self.model_part)

            #Note: 
            #the buffer size should be set once the mesh is read for the first time

            #set the buffer size
            self.model_part.SetBufferSize(buffer_size)

            #set the degrees of freedom
            main_solver.AddDofs(self.model_part,problem_type,rotation_dofs)

            #set the constitutive law
            import constitutive_law_python_utility as constitutive_law_utils
            constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.model_part,domain_size);
            constitutive_law.Initialize();

            #set mesh searches and modeler
            modeler.InitializeDomains();

            #find nodal h
            modeler.SearchNodalH();

        else:

            #set mesh searches and modeler
            modeler.InitializeDomains();
            
            # Note:
            # master elements and neighbours not loaded in the restart
            # they are initialized and the search is performed again

        #to ensure the nodal variables initialization in buffer
        self.start_steps = buffer_size - 1 


    #######################################################################   
    def StartModelRead(self,buffer_size,domain_size,problem_type,rotation_dofs,main_solver):
        
        if(self.load_restart_flag == False):
            #reading the model 
            self.model_part_io = ModelPartIO(self.problem_name)
            self.model_part_io.ReadModelPart(self.model_part)

            #Note: 
            #the buffer size should be set once the mesh is read for the first time

            #set the buffer size
            self.model_part.SetBufferSize(buffer_size)

            #set the degrees of freedom
            main_solver.AddDofs(self.model_part,problem_type,rotation_dofs)

            #set the constitutive law
            import constitutive_law_python_utility as constitutive_law_utils
            constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.model_part,domain_size);
            constitutive_law.Initialize();

 
        #to ensure the nodal variables initialization in buffer
        self.start_steps = buffer_size - 1 


    #######################################################################   
    def SetTimeVariables(self,time_step,steps_number):

        self.step_i       = 0
        self.step_n       = steps_number + self.start_steps + 1

        self.time_step    = time_step
        self.steps_number = steps_number
        self.time_start   = self.start_steps * self.time_step
 
    #######################################################################   
    def InitializeTimeIntegration(self,gid_print,modeler,graph_plot,conditions):

        #initialize time variables
        write_id          = 0
        self.current_step = 0
              

        if(self.load_restart_flag == True):
            self.step_i       = self.model_part.ProcessInfo[TIME_STEPS]+1;
            self.time_step    = self.model_part.ProcessInfo[DELTA_TIME];
            write_id          = self.model_part.ProcessInfo[WRITE_ID];
            self.current_step = self.step_i - self.step_n
        else:
            self.model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = self.time_step;
            conditions.Initialize();
            self.step_i       = 0

            
        #initialize print results variables for time integration
        gid_print.Initialize(self.time_start,self.step_n,self.current_step,write_id)

        #initialize mesh modeling variables for time integration
        modeler.Initialize(self.current_step,self.current_step)

        #initialize graph plot variables for time integration
        graph_plot.Initialize(self.current_step)
    
    #######################################################################   
    def InitializeTimeIntegration(self,gid_print,conditions):

        #initialize time variables
        write_id          = 0
        self.current_step = 0
              

        if(self.load_restart_flag == True):
            self.step_i       = self.model_part.ProcessInfo[TIME_STEPS]+1;
            self.time_step    = self.model_part.ProcessInfo[DELTA_TIME];
            write_id          = self.model_part.ProcessInfo[WRITE_ID];
            self.current_step = self.step_i - self.step_n
        else:
            self.model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = self.time_step;
            conditions.Initialize();
            self.step_i       = 0

            
        #initialize print results variables for time integration
        gid_print.Initialize(self.time_start,self.step_n,self.current_step,write_id)

    

    #######################################################################   
    def TimeStep(self):
        return self.time_step
       
    #######################################################################
