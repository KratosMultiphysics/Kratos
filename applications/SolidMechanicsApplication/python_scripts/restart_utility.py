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
    def CleanPreviousFileType(self,file_ending_type):

        filelist = [ f for f in os.listdir(self.problem_path) if f.endswith(file_ending_type) ]

        for f in filelist:
            try:
                os.remove(f)
            except WindowsError:
                pass

    #######################################################################   
    def CleanPosteriorFileType(self,file_ending_type):

        total_files = 0;
        filelist = []
        for f in os.listdir(self.problem_path):
            if(f.endswith(file_ending_type)):
                total_files = total_files + 1
                
        total_files = total_files + 100  #arbitrary number to ensure to remove all files

        for rfile in range(0,total_files):
            for f in os.listdir(self.problem_path):
                step_time = restart_time + rfile
                if(f.endswith("_"+str(step_time)+file_ending_type)):
                    filelist.append(f);                           
                                
        for f in filelist:
            try:
                os.remove(f)
            except WindowsError:
                pass
 
    #######################################################################   
    def CleanPreviousFiles(self):

        print "Start: -remove previous problem files-"

        #remove previous results:
        file_ending_type= ".post.bin"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type= ".post.res"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type= ".post.msh"
        self.CleanPreviousFileType(file_ending_type)
    
        #remove previous graph files:
        file_ending_type= ".png"
        self.CleanPreviousFileType(file_ending_type)

        #remove previous restart files:
        file_ending_type= ".rest"
        self.CleanPreviousFileType(file_ending_type)
                              

    #######################################################################   
    def CleanPosteriorFiles(self,restart_time):
        print "Restart: -remove restart step posterior problem files-"

        #remove previous results after restart:
        file_ending_type= ".post.bin"
        self.CleanPosteriorFileType(file_ending_type)

        file_ending_type= ".post.res"
        self.CleanPosteriorFileType(file_ending_type)

        file_ending_type= ".post.msh"
        self.CleanPosteriorFileType(file_ending_type)

        #remove previous graphs after restart:

        file_ending_type= ".png"
        self.CleanPosteriorFileType(file_ending_type)

        #remove previous restart files:
        file_ending_type= ".rest"
        self.CleanPosteriorFileType(file_ending_type)
 

    #######################################################################   
    def Save(self,current_time,current_step):

        restart_path = os.path.join(self.problem_path,self.problem_name + "_" + str(current_time))
        serializer = Serializer(restart_path,self.serializer_flag)
        serializer.Save("ModelPart",self.model_part);
        print " RESTART PRINTED :", current_step;
 
    #######################################################################
