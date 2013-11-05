#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()
import os

class RestartUtility:
    #######################################################################
    def __init__(self,model_part,problem_path,problem_name):

        self.model_part = model_part
        
        #set restart flags
        self.load_restart_flag = False
        self.save_restart_flag = False

        #set problem name
        self.problem_name = problem_name

        #set problem path
        self.problem_path = problem_path
             
        #set serializer flag
        self.serializer_flag = "SERIALIZER_NO_TRACE"
        #self.serializer_flag = "SERIALIZER_TRACE_ERROR"
         

    #######################################################################
    def Load(self,restart_time):

        self.load_restart_flag = True
        restart_path = os.path.join(self.problem_path,self.problem_name + "_" + str(restart_time))
        serializer = Serializer(restart_path,self.serializer_flag)
        serializer.Load("ModelPart",self.model_part)
 
    #######################################################################   
    def CleanPreviousFiles(self,list_files):
        print "Start: -remove previous problem files-"
        #remove previous results:
        filelist1 = [ f for f in os.listdir(self.problem_path) if f.endswith("post.bin") ]
        for f in filelist1:
            os.remove(f)

        filelist1_2 = [ f for f in os.listdir(self.problem_path) if f.endswith("post.res") ]
        for f in filelist1_2:
            os.remove(f)

        filelist1_3 = [ f for f in os.listdir(self.problem_path) if f.endswith("post.msh") ]
        for f in filelist1_3:
            os.remove(f)
            
        #remove previous restart files:
        filelist3 = [ f for f in os.listdir(self.problem_path) if f.endswith(".rest") ]
        for f in filelist3:
	    try:
            	os.remove(f)
            except WindowsError:
                pass
                      
        #remove previous graph files:
        filelist4 = [ f for f in os.listdir(self.problem_path) if f.endswith(".png") ]
        for f in filelist4:
            os.remove(f)

        #remove previous list files:
        list_files.RemoveListFiles();

    #######################################################################   
    def CleanPosteriorFiles(self,time_step,restart_time,list_files):
        print "Restart: -remove restart step posterior problem files-"

        #remove previous results after restart:
        total_files = 0;
        filelist1 = []
        for f in os.listdir(self.problem_path):
            if(f.endswith(".post.bin")):
                total_files = total_files + 1
                
        total_files = total_files + 100  #arbitrary number to ensure to remove all files

        for rfile in range(0,total_files):
            for f in os.listdir(self.problem_path):
                step_time = restart_time + rfile * time_step
                if(f.endswith("_"+str(step_time)+".post.bin")):
                    filelist1.append(f);                           
                                
        for f in filelist1:
            os.remove(f)

        #remove previous graphs after restart:
        total_files = 0;
        filelist2 = []
        for f in os.listdir(self.problem_path):
            if(f.endswith(".png")):
                total_files = total_files + 1
                
        total_files = total_files + 100  #arbitrary number to ensure to remove all files

        for rfile in range(0,total_files):
            for f in os.listdir(self.problem_path):
                step_time = restart_time + rfile * time_step
                if(f.endswith("_"+str(step_time)+".png")):
                    filelist2.append(f);
                                
        for f in filelist2:
            os.remove(f)


        #remove previous restart files:
        counter = 1;
        total_files = 0;
        filelist3 = []
        for f in os.listdir(self.problem_path):
            if(f.endswith(".rest")):
                total_files = total_files + 1
                                                
        total_files = total_files + 100  #arbitrary number to ensure to remove all rest files


        for rfile in range(0,total_files):
            for f in os.listdir(self.problem_path):
                step_time = restart_time + rfile * time_step
                if(f.endswith("_"+str(step_time)+".rest")):
                    filelist3.append(f);
                                                            
        for f in filelist3:
            os.remove(f)

        #remove previous list files and rebuild it:
        list_files.RemoveListFiles();
        list_files.BuildListFiles();
    

    #######################################################################   
    def Save(self,current_time,current_step):

        restart_path = os.path.join(self.problem_path,self.problem_name + "_" + str(current_time))
        serializer = Serializer(restart_path,self.serializer_flag)
        serializer.Save("ModelPart",self.model_part);
        print " RESTART PRINTED :", current_step;
 
    #######################################################################
