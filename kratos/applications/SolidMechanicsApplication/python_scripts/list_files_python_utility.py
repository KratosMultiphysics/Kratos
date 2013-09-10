#Class to manage list files and result files
from KratosMultiphysics import *
CheckForPreviousImport()

class ListFilesUtility:
    #######################################################################
    def __init__(self,problem_path,problem_name,print_lists):

        #set problem name
        self.problem_name = problem_name

        #set problem path
        self.problem_path = problem_path

        #set print list
        if(print_lists == "True"):
            self.print_lists = True
        else:
            self.print_lists = False

        #initialize lists
        self.header_in_list = []
        self.listprint      = [] 
        
    #######################################################################
    def Initialize(self,file_list):
        
        self.file_list = file_list

        for lfile in file_list:
            self.listprint.append(1);
            self.header_in_list.append(True);

    #######################################################################   
    def BuildListFiles(self):
        
        if(self.print_lists == True):
            total_files = 0;
            for f in os.listdir(self.problem_path):
                if(f.endswith(".post.bin")):
                    total_files = total_files + 1
     
            for rfile in range(0,total_files):
                for f in os.listdir(problem_path):
                    if(f.endswith("_"+str(rfile)+".post.bin")):
                        problempath= self.problem_path + "/" + self.problem_name + "_1.post.lst"
                        if(os.path.exists(problempath) == False):
                            listfile = open(problempath,"a")
                            problemname = "Multiple\n"
                            listfile.write(problemname)
                            problemname = self.problem_name + "_" + str(rfile) + ".post.bin\n" 
                            listfile.write(problemname)
                            listfile.close()
                            file_set = True
                        else:
                            listfile = open(problempath,"a")
                            problemname = self.problem_name + "_" + str(rfile) + ".post.bin\n" 
                            listfile.write(problemname)
                            listfile.close()
                            file_set = True
            
                        num_list_files = len(self.file_list) 
                        for lfile in range(0,num_list_files):
                            if( self.file_list[lfile] == self.listprint[lfile] ):
                                problempath= self.problem_path + "/" + self.problem_name + "_"+ str(self.file_list[lfile]) + ".post.lst"
                                if(os.path.exists(problempath) == False):
                                    listfile = open(problempath,"a")
                                    problemname = "Multiple\n" 
                                    listfile.write(problemname)
                                    self.header_in_list[lfile] = False
                                    problemname = self.problem_name + "_" + str(rfile) + ".post.bin\n" 
                                    listfile.write(problemname)
                                    listfile.close()
                                    self.listprint[lfile] = 1
                                else:
                                    listfile = open(problempath,"a")
                                    problemname = self.problem_name + "_" + str(rfile) + ".post.bin\n" 
                                    listfile.write(problemname)
                                    listfile.close()
                                    self.listprint[lfile] = 1
                            else:
                                self.listprint[lfile] = self.listprint[lfile]+1
            
        
        
    #######################################################################   
    def RemoveListFiles(self):

        #remove previous list files:
        filelist2 = [ f for f in os.listdir(self.problem_path) if f.endswith(".lst") ]
        for f in filelist2:
            os.remove(f)

    #######################################################################   
    def PrintListFiles(self,current_step):

         #print list files:
         if(self.print_lists == True):
             problempath= self.problem_path + "/" + self.problem_name + "_1.post.lst"
             if(os.path.exists(problempath) == False):
                 listfile = open(problempath,"a")
                 #if(current_step == 0 and general_variables.LoadRestart == "False"):
                 problemname = "Multiple\n" 
                 listfile.write(problemname)
                 #endif
                 problemname = self.problem_name + "_" + str(current_step) + ".post.bin\n" 
                 listfile.write(problemname)
                 listfile.close()
             else:
                 listfile = open(problempath,"a")
                 problemname = self.problem_name + "_" + str(current_step) + ".post.bin\n" 
                 listfile.write(problemname)
                 listfile.close()
      
             num_list_files = len(self.file_list) 
             for lfile in range(0,num_list_files):
                 if( general_variables.file_list[lfile] == listprint[lfile] ):
                     problempath= self.problem_path + "/" + self.problem_name + "_"+ str(self.file_list[lfile]) + ".post.lst"
                     if(os.path.exists(problempath) == False):
                         listfile = open(problempath,"a")
                         #if(header_in_list[lfile] == True and general_variables.LoadRestart == "False"):
                         problemname = "Multiple\n" 
                         listfile.write(problemname)
                         self.header_in_list[lfile] = False
                         #endif
                         problemname = self.problem_name + "_" + str(current_step) + ".post.bin\n" 
                         listfile.write(problemname)
                         listfile.close()
                         self.listprint[lfile] = 1
                     else:
                         listfile = open(problempath,"a")
                         problemname = self.problem_name + "_" + str(current_step) + ".post.bin\n" 
                         listfile.write(problemname)
                         listfile.close()
                         self.listprint[lfile] = 1
                 else:
                     self.listprint[lfile] = self.listprint[lfile]+1

     #######################################################################   
