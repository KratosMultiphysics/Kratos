from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FemToDemApplication  import *
import shutil 
CheckForPreviousImport()

def Wait():
    input("Press Something")

#============================================================================================================================
class AdaptiveMeshRefinementUtility:

    def __init__(self, ProjectParameters, starting_time, solver_constructor, constitutive_law_utility,
     gid_output_utility, conditions_util, ProblemPath):
        
        ## Parameters and utilities
        self.ProjectParameters = ProjectParameters
        self.solver_constructor = solver_constructor
        self.constitutive_law_utility = constitutive_law_utility
        self.gid_output_utility = gid_output_utility
        
        self.conditions_util = conditions_util
        self.problem_path    = os.getcwd()
        self.AMR_files_path  = os.path.join(self.problem_path, "AMR_Files")
        
        self.n_refinements = 0
        self.last_refinement_id = 1 #Same as the initial current_id
        
        ## Time operations initialization
        self.ending_time = ProjectParameters["problem_data"]["end_time" ].GetDouble()
        self.delta_time  = ProjectParameters["problem_data"]["time_step"].GetDouble()
        
        # set AMR frequency
        self.amr_frequency = ProjectParameters["AMR_data"]["refinement_frequency"].GetDouble()
        if(self.amr_frequency < self.delta_time):
            self.amr_frequency = self.delta_time

        # set time counter
        self.time_counter = starting_time + self.amr_frequency

        # set time operation tolerance
        self.tolerance = self.delta_time * 1e-10;

#============================================================================================================================
    def Initialize(self):
        
        self.gid_path = self.ProjectParameters["AMR_data"]["gid_path"].GetString()
        AMR_info_path = os.path.join(self.problem_path,"AMR_info.txt")

        activate_AMR = True
        
        if(self.amr_frequency > self.ending_time):
            activate_AMR = False
        
        # creates the AMR_Files folder or roemove if exists
        if not os.path.isdir(self.AMR_files_path):
            os.makedirs(str(self.AMR_files_path))
        else:  
            shutil.rmtree(str(self.AMR_files_path), ignore_errors = True)
            os.makedirs(str(self.AMR_files_path))

        if os.path.isfile(str(AMR_info_path)):
            shutil.rmtree(str(AMR_info_path), ignore_errors = True)
        
        return activate_AMR
        
#============================================================================================================================
    def CheckAMR(self, current_time):
        
        refine = False
        last_mesh = False
        
        if( (current_time + self.tolerance >= self.time_counter) and (current_time + self.tolerance < self.ending_time) ):
            self.time_counter = self.time_counter + self.amr_frequency
            refine = True
        elif(current_time + self.tolerance >= self.ending_time):
            last_mesh = True

        return refine, last_mesh
        
#============================================================================================================================       
    def Execute(self, model_part, main_step_solver, gid_output_util, current_time, current_id):
        
        ## Previous definitions ---------------------------------------------------------------------------------------------
        
        problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()
        #GidOutputConfiguration = self.ProjectParameters.GidOutputConfiguration
        output_mode = self.ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString()
        output_multiple_files = self.ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["MultiFileFlag"].GetString()
        plane_state = self.ProjectParameters["AMR_data"]["plane_state"].GetString()
        mesh_optimality_criteria = self.ProjectParameters["AMR_data"]["mesh_optimality_criteria"].GetString()
        permissible_error = self.ProjectParameters["AMR_data"]["permissible_error"].GetDouble()
        Mapping_Procedure = self.ProjectParameters["AMR_data"]["Mapping_Procedure"].GetString()
        #self.main_model_part = model_part
        ## Finalize previous post results -----------------------------------------------------------------------------------
        #print("dentro de execute1", os.getcwd())
        #i = 69
        #print(str("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.msh "))
        #Wait()
        #gid_output_util.finalize_results()
        """
        if(output_mode=="GiD_PostBinary"):
            if(output_multiple_files=="MultipleFiles"):
                for i in range(self.last_refinement_id, current_id + 1):
                    # os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.bin")
                    src = os.path.join(self.problem_path, str(problem_name) + str(i) + ".post.bin")
                    dst = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) +"_step_" + str(i) + ".post.bin")
                    shutil.move(src, dst)
            else:
                # os.system("move "+str(self.problem_path)+"/"+str(problem_name)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.bin")
                src = os.path.join(self.problem_path, str(problem_name) + ".post.bin ")
                dst = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + ".post.bin")
                shutil.move(src, dst)

        else:
            if(output_multiple_files=="MultipleFiles"):
                for i in range(self.last_refinement_id, current_id + 1):
                    #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.msh")
                    #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.res")
                    src_mesh = os.path.join(self.problem_path, str(problem_name) + "_" + str(i) + ".post.msh")
                    src_res  = os.path.join(self.problem_path, str(problem_name) + "_" + str(i) + ".post.res")
                    dst_mesh = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_"+  str(self.n_refinements) + "_step_" + str(i) + ".post.msh")
                    dst_res  = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + "_step_" + str(i) + ".post.res")

                    shutil.move(src_mesh, dst_mesh)
                    shutil.move(src_res , dst_res)

            else:
                #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.msh")
                #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.res")

                src_mesh = os.path.join(self.problem_path, str(problem_name) + "_" + str(self.last_refinement_id) +".post.msh")
                src_res  = os.path.join(self.problem_path, str(problem_name) + "_" + str(self.last_refinement_id) +".post.res" )
                dst_mesh = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + ".post.msh")
                dst_res  = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + ".post.res")

                shutil.move(src_mesh, dst_mesh)
                shutil.move(src_res , dst_res)
            """
        #print("antes de copy  ")   #cornejo
        #Wait()        
        #os.system("copy "+str(self.problem_path)+"/"+str(problem_name)+".mdpa "+str(self.AMR_files_path)+"/"+str(problem_name)+"_mesh_"+str(self.n_refinements)+".mdpa")
        src_mdpa = os.path.join(self.problem_path, str(problem_name) +".mdpa")
        dst_mdpa = os.path.join(self.AMR_files_path, str(problem_name) + "_mesh_" + str(self.n_refinements) + ".mdpa")
        shutil.copy(src_mdpa, dst_mdpa)


        ## MESH LOOP --------------------------------------------------------------------------------------------------------
        
        #print("dentro de execute2")   #cornejo
        #Wait()
        
        print("----------------------------------")
        print("      START MESH REFINEMENT       ")
        print("----------------------------------")
        
        mesh_convergence = False
        max_num_iter = 5
        iteration_number = 0
        
        while(mesh_convergence == False and iteration_number < max_num_iter ):
            
            iteration_number = iteration_number + 1

            if(iteration_number == 1):
                AdaptiveMeshRefinementProcess(model_part,plane_state,
                                                         problem_name,
                                                         self.problem_path,
                                                         mesh_optimality_criteria,
                                                         permissible_error,
                                                         self.n_refinements).Execute()
                #Wait()

                # Move the posts of the amr to the AMR_folder
                src_mesh = os.path.join(self.problem_path, str(problem_name)   + "_AMR_parameters.post.msh")
                src_res  = os.path.join(self.problem_path, str(problem_name)   + "_AMR_parameters.post.res")
                dst_mesh = os.path.join(self.AMR_files_path, str(problem_name) + "_AMR_parameters_mesh_" + str(self.n_refinements) + ".post.msh")
                dst_res  = os.path.join(self.AMR_files_path, str(problem_name) + "_AMR_parameters_mesh_" + str(self.n_refinements) + ".post.res")

                shutil.move(src_mesh, dst_mesh)
                shutil.move(src_res, dst_res)

            else:
                AdaptiveMeshRefinementProcess(model_part, plane_state,
                                                          problem_name,
                                                          self.problem_path,
                                                          mesh_optimality_criteria,
                                                          permissible_error,
                                                          self.n_refinements).ExecuteAfterOutputStep()
            
            #print("despues execute after output step")
           # Wait()
            

            # GID GENERATE THE NEW MESH BASED ON THE BACKGROUND MESH  ----> TODO
            #os.system("cd " + str(self.gid_path) + " && ./gid -b "+ str(self.problem_path)+"/"+str(problem_name)+".bch -n")
            #os.system("cd && mv "+str(self.problem_path)+"/"+str(problem_name)+".dat "+str(self.problem_path)+"/"+str(problem_name)+".mdpa && rm "+str(self.problem_path)+"/"+str(problem_name)+"-1.dat")
            
            # Execute .bch file with GiD 
            os.system("cd " + str(self.gid_path[:-8]) + " && gid -b " + os.path.join(self.problem_path, str(problem_name) + ".bch -n"))
            #print("despues de gid commands gid path: ", self.gid_path[:-8])
            #Wait()            
            #shutil.move(os.path.join(str(self.problem_path), str(problem_name) + ".dat" ) , os.path.join(str(self.problem_path), str(problem_name) + ".mdpa"))
            #shutil.rmtree(os.path.join(str(self.problem_path), str(problem_name)+"-1.dat"))
            
            
            ## Finalize previous mesh ---------------------------------------------------------------------------------------
            
            #Finalize previous solver
            #main_step_solver.Finalize()
            #main_step_solver.FinalizeSolutionStep() # cornejo
            
            #Save previous model_part
            model_part_old = model_part
            
            #------------------------------------------------------------------->> Aqui estamos
            ## Generate new Model Part --------------------------------------------------------------------------------------

            # Definition of model part
            #model_part = ModelPart("SolidDomain")
            delta_time = self.ProjectParameters["problem_data"]["time_step" ].GetDouble()
            step = model_part_old.ProcessInfo[STEP]

            #current_time += delta_time

            model_part = ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
            model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())
            model_part.ProcessInfo.SetValue(DELTA_TIME,  delta_time)
            model_part.ProcessInfo.SetValue(TIME, current_time)  # curent or current+delta_t ??
            model_part.ProcessInfo.SetValue(STEP, step)


            #print("time", model_part.ProcessInfo[TIME])
            #Wait()
            ###TODO replace this "model" for real one once available in kratos core
            self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : model_part}


            # we delete 3 added parameters to avoid error --> TODO
            self.ProjectParameters["solver_settings"].RemoveValue("damp_factor_m")
            self.ProjectParameters["solver_settings"].RemoveValue("dynamic_factor")
            self.ProjectParameters["solver_settings"].RemoveValue("stabilization_factor")


            #construct the solver (main setting methods are located in the solver_module)
            solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
            main_step_solver   = solver_module.CreateSolver(model_part, self.ProjectParameters["solver_settings"])


            #print("delta time:" , model_part.ProcessInfo[DELTA_TIME])
            #Wait()

            # Add variables (always before importing the model part)
            main_step_solver.AddVariables()

            # Read model_part (note: the buffer_size is set here) (restart is read here)
            main_step_solver.ImportModelPart()

            # Add dofs (always after importing the model part)
            if((model_part.ProcessInfo).Has(IS_RESTARTED)):
                if(model_part.ProcessInfo[IS_RESTARTED] == False):
                    main_step_solver.AddDofs()
            else:
                main_step_solver.AddDofs()

            # Add materials (assign material to model_parts if Materials.json exists)
            #self.AddMaterials(model_part)

            #test
            model_part.ProcessInfo.SetValue(STEP, step)
            dyn_mode  = model_part_old.ProcessInfo[IS_DYNAMIC]
            model_part.ProcessInfo.SetValue(IS_DYNAMIC, dyn_mode)
            #print("stepppp1 :", model_part.ProcessInfo[STEP])
            #Wait()
            # Add processes -> en main
            #self.model_processes = self.AddProcesses(model_part)
            #self.model_processes.ExecuteInitialize()

            #### START SOLUTION ####

            #self.computing_model_part = main_step_solver.GetComputingModelPart()

            ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
            #main_step_solver.Initialize()

            #neighbour_elemental_finder =  FindElementalNeighboursProcess(model_part, 2, 5)
            #neighbour_elemental_finder.Execute()

            #model_part.ProcessInfo[STEP] += 1
            #model_part.CloneTimeStep(current_time) 


            # Definition of output utility
            #gid_output_util = self.gid_output_utility.GidOutputUtility(self.ProjectParameters,
                                                                                 #problem_name,
                                                                                 #current_time,
                                                                                 #self.ending_time,
                                                                                # self.delta_time)
            
            
            #model_part.ProcessInfo[MESH_REFINED] = 1

            ## Mapping of variables -----------------------------------------------------------------------------------------
            
            MappingVariablesProcess(model_part_old, model_part, "Constant", Mapping_Procedure).Execute()

            ## Erase old Model Part -----------------------------------------------------------------------------------------
            model_part_old = None
            
            mesh_convergence = True # only for testing TODO
            
        if(mesh_convergence == True):
            print("NEW MESH CONVERGED AFTER ", iteration_number," ITERATIONS")
        else:
            print("### WARNING: NO MESH CONVERGED AFTER ", iteration_number, " ITERATIONS ###")
        
        ## Saving files of new mesh -----------------------------------------------------------------------------------------
        src = os.path.join(str(self.problem_path),str(problem_name)+".bgm ")
        dst = os.path.join(str(self.AMR_files_path), str(problem_name)+"_AMR_" + str(self.n_refinements+1) + ".bgm")
        shutil.copy(src, dst)
        
        self.last_refinement_id = current_id + 1
        self.n_refinements = self.n_refinements + 1
        
        return model_part, main_step_solver



#============================================================================================================================
    def Finalize(self, model_part, current_id):
        
        ## Previous definitions ---------------------------------------------------------------------------------------------
        
        problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()
        #GidOutputConfiguration = self.ProjectParameters.GidOutputConfiguration
        output_mode  = self.ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString()
        output_multiple_files  = self.ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["MultiFileFlag"].GetString()
        plane_state = self.ProjectParameters["AMR_data"]["plane_state"].GetString()
        mesh_optimality_criteria = self.ProjectParameters["AMR_data"]["mesh_optimality_criteria"].GetString()
        permissible_error = self.ProjectParameters["AMR_data"]["permissible_error"].GetDouble()

        ## Finalize previous post results -----------------------------------------------------------------------------------
        """
        if(output_mode == "GiD_PostBinary"):
            if(output_multiple_files=="MultipleFiles"):
                for i in range(self.last_refinement_id, current_id + 1):
                    #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.bin")
                    src = os.path.join(str(self.problem_path), str(problem_name) + "_" + str(i) + ".post.bin")
                    dst = os.path.join(str(self.AMR_files_path), str(problem_name) + "_results_mesh_" + str(self.n_refinements) + "_step_" + str(i) + ".post.bin")
                    shutil.move(src, dst)
            else:
                #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.bin")
                src = os.path.join(str(self.problem_path), str(problem_name) +  ".post.bin")
                dst = os.path.join(str(self.AMR_files_path), str(problem_name) + "_results_mesh_" + str(self.n_refinements) + ".post.bin")
                shutil.move(src, dst)
        else:
            if(output_multiple_files == "MultipleFiles"):
                for i in range(self.last_refinement_id, current_id + 1):
                    #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.msh")
                    #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.res")
                    mshsrc = os.path.join(str(self.problem_path), str(problem_name)+"_"+str(i)+".post.msh")
                    ressrc = os.path.join(str(self.problem_path), str(problem_name)+"_"+str(i)+".post.res")
                    mshdst = os.path.join(str(self.AMR_files_path), str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.msh")
                    resdst = os.path.join(str(self.AMR_files_path), str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.res")
                    shutil.move(mshsrc, mshdst)
                    shutil.move(ressrc, resdst)

            else:
                #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.msh")
                #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.res")
                mshsrc = os.path.join(str(self.problem_path), str(problem_name) + "_" + str(self.last_refinement_id)+".post.msh")
                ressrc = os.path.join(str(self.problem_path), str(problem_name) + "_" + str(self.last_refinement_id)+".post.res")
                mshdst = os.path.join(str(self.AMR_files_path), str(problem_name) + "_results_mesh_" + str(self.n_refinements)+".post.msh")
                resdst = os.path.join(str(self.AMR_files_path), str(problem_name) + "_results_mesh_" + str(self.n_refinements)+".post.res")
                shutil.move(mshsrc, mshdst)
                shutil.move(ressrc, resdst)

        
        #os.system("cp "+str(self.problem_path)+"/"+str(problem_name)+".mdpa "+str(self.AMR_files_path)+"/"+str(problem_name)+"_mesh_"+str(self.n_refinements)+".mdpa")
        src = os.path.join(str(self.problem_path), str(problem_name) + ".mdpa")
        dst = os.path.join(str(self.AMR_files_path), str(problem_name) + "_mesh_" + str(self.n_refinements) + ".mdpa")
        shutil.copy(src, dst)
        """
        ## Compute and save info of last mesh -------------------------------------------------------------------------------
        
        AdaptiveMeshRefinementProcess(model_part, 
                                      plane_state,
                                      problem_name,
                                      self.problem_path,
                                      mesh_optimality_criteria,
                                      permissible_error,
                                      self.n_refinements).ExecuteFinalize()
        
        #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_AMR_parameters.post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_AMR_parameters_mesh_"+str(self.n_refinements)+".post.msh")
        #os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_AMR_parameters.post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_AMR_parameters_mesh_"+str(self.n_refinements)+".post.res")

        srcmsh = os.path.join(str(self.problem_path),  str(problem_name) + "_AMR_parameters.post.msh")
        srcres = os.path.join(str(self.problem_path),  str(problem_name) + "_AMR_parameters.post.res")
        dstmsh = os.path.join(str(self.AMR_files_path), str(problem_name) + "_AMR_parameters_mesh_" + str(self.n_refinements) + ".post.msh")
        dstres = os.path.join(str(self.AMR_files_path), str(problem_name) + "_AMR_parameters_mesh_" + str(self.n_refinements) + ".post.res")
        shutil.move(srcmsh, dstmsh)
        shutil.move(srcres, dstres)

        #os.system("mv "+str(self.problem_path)+"/AMR_info.txt "+str(self.AMR_files_path))
        src = os.path.join(str(self.problem_path),   "AMR_info.txt")
        dst = os.path.join(str(self.AMR_files_path), "AMR_info.txt")
        shutil.move(src, dst)


# copy of the main ============================================================================================================================        

    def AddMaterials(self, main_model_part):

        # Assign material to model_parts (if Materials.json exists)
        import process_factory

        if os.path.isfile("Materials.json"):
            materials_file = open("Materials.json",'r')
            MaterialParameters = Parameters(materials_file.read())

            if(MaterialParameters.Has("material_models_list")):

                ## Get the list of the model_part's in the object Model
                for i in range(self.ProjectParameters["solver_settings"]["problem_domain_sub_model_part_list"].size()):
                    part_name = self.ProjectParameters["solver_settings"]["problem_domain_sub_model_part_list"][i].GetString()
                    if(main_model_part.HasSubModelPart(part_name) ):
                        self.Model.update({part_name: main_model_part.GetSubModelPart(part_name)})


                assign_materials_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( MaterialParameters["material_models_list"] )

            for process in assign_materials_processes:
                process.Execute()
        else:
            print(" No Materials.json found ")

#============================================================================================================================

    def AddProcesses(self, main_model_part):

        # Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
        ## Get the list of the submodel part in the object Model
        for i in range(self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            if(main_model_part.HasSubModelPart(part_name) ):
                self.Model.update({part_name: main_model_part.GetSubModelPart(part_name)})
        
        # Obtain the list of the processes to be applied
        import process_handler

        process_parameters = Parameters("{}") 
        process_parameters.AddValue("echo_level", self.ProjectParameters["problem_data"]["echo_level"])
        process_parameters.AddValue("constraints_process_list", self.ProjectParameters["constraints_process_list"])
        process_parameters.AddValue("loads_process_list", self.ProjectParameters["loads_process_list"])
        if( self.ProjectParameters.Has("problem_process_list") ):
            process_parameters.AddValue("problem_process_list", self.ProjectParameters["problem_process_list"])
        if( self.ProjectParameters.Has("output_process_list") ):
            process_parameters.AddValue("output_process_list", self.ProjectParameters["output_process_list"])

        return (process_handler.ProcessHandler(self.Model, process_parameters))