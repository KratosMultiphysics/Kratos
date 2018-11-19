from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

#Base class to develop other solvers
class ModelManager(object):
    """The base class for solid mechanic model build process.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs.

    """
    def __init__(self, Model, custom_settings):

        default_settings = KratosMultiphysics.Parameters("""
        {
           "model_name": "solid_domain",
           "dimension": 3,
           "bodies_list": [],
           "domain_parts_list": [],
           "processes_parts_list": [],
           "output_model_part": "output_domain",
           "solving_model_part": "computing_domain",
           "composite_solving_parts": [],
           "input_file_settings": {
                "type" : "mdpa",
                "name" : "unknown_name",
                "label": 0
           },
           "variables":[]
        }
        """)

        # attention dofs mover to solid_solver
        if(custom_settings.Has("dofs")):
            custom_settings.RemoveValue("dofs")
            print(" WARNING: [MODEL_MANAGER] dofs moved to SolidSolver")

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings["input_file_settings"].ValidateAndAssignDefaults(default_settings["input_file_settings"])

        # Set void model
        self.model = Model
        self.main_model_part = self._create_main_model_part()

        # Process Info
        self.process_info = self.main_model_part.ProcessInfo

        # Variables settings
        self.nodal_variables = []

        # Composite solving parts
        self.transfer_solving_parts = []

    ########

    #
    def ExecuteInitialize(self):
        self.ImportModel()

    #
    def ExecuteBeforeSolutionLoop(self):
        pass
    #
    def ExecuteInitializeSolutionStep(self):
        if( self._domain_parts_updated() ):
            self._update_composite_solving_parts()
            # print(" UPDATE_SOLVING_PARTS Initialize")
    #
    def ExecuteFinalizeSolutionStep(self):
        pass
    #
    def ExecuteBeforeOutputStep(self):
        if( self._domain_parts_updated() ):
            self._update_composite_solving_parts()
            # print(" UPDATE_SOLVING_PARTS Before output")
    #
    def ExecuteAfterOutputStep(self):
        if( self._domain_parts_updated() ):
            self._update_composite_solving_parts()
            # print(" UPDATE_SOLVING_PARTS After output")

    ########

    #
    def ImportModel(self):

        self._add_variables()

        #print(self._class_prefix()+" Importing model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["input_file_settings"]["name"].GetString()

        if(self.settings["input_file_settings"]["type"].GetString() == "mdpa"):
            # Import model part from mdpa file.
            print(self._class_prefix()+" Reading file: "+ input_filename + ".mdpa")
            #print("   " + os.path.join(problem_path, input_filename) + ".mdpa ")
            sys.stdout.flush()

            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SPACE_DIMENSION, self.settings["dimension"].GetInt())
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.settings["dimension"].GetInt()) # Legacy

            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()

            # Somewhere must ask if you want to clean previous files
            self._clean_previous_result_files()

        elif(self.settings["input_file_settings"]["type"].GetString() == "rest"):
            # Import model part from restart file.
            restart_path = os.path.join(problem_path, self.settings["input_file_settings"]["name"].GetString() + "__" + str(self.settings["input_file_settings"]["label"].GetInt() ) )
            if(os.path.exists(restart_path+".rest") == False):
                raise Exception("Restart file not found: " + restart_path + ".rest")
            print("   Loading Restart file: ", restart_path + ".rest ")
            # set serializer flag
            serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

            serializer = KratosMultiphysics.FileSerializer(restart_path, serializer_flag)
            serializer.Load(self.main_model_part.Name, self.main_model_part)

            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
            self.main_model_part.ProcessInfo[KratosSolid.RESTART_STEP_TIME] = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

            #I use it to rebuild the contact conditions.
            load_step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] +1
            self.main_model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step
            # print("   Finished loading model part from restart file ")

            print( self.main_model_part )

            self._build_composite_solving_parts()

        else:
            raise Exception("Other input options are not yet implemented.")


        dofs = self.main_model_part.NumberOfNodes() * self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]
        #print (self._class_prefix()+" Finished importing model part")
        print (self._class_prefix()+" Model Ready (DOFs:"+str(dofs)+")")

    #
    def ExportModel(self):
        name_out_file = self.settings["input_file_settings"]["name"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        # Model part writing
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    #
    def CleanModel(self):
        self._clean_body_parts()

    ########

    def GetProcessInfo(self):
        return self.process_info

    def GetModel(self):
        return self.model

    def GetMainModelPart(self):
        return self.main_model_part

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["solving_model_part"].GetString())

    def GetOutputModelPart(self):
        #return self.main_model_part.GetSubModelPart(self.settings["output_model_part"].GetString())
        return self.main_model_part.GetSubModelPart(self.settings["solving_model_part"].GetString())

    def SaveRestart(self):
        pass #one should write the restart file here

    def SetVariables(self, variables):
        self.nodal_variables = self.nodal_variables + variables


    #### Model manager internal methods ####

    def _create_main_model_part(self):
        # Defining the model_part
        main_model_part = self.model.CreateModelPart(self.settings["model_name"].GetString())
        return main_model_part

    def _create_sub_model_part(self, part_name):
        self.main_model_part.CreateSubModelPart(part_name)

    def _add_variables(self):

        self._set_input_variables()

        self.nodal_variables = list(set(self.nodal_variables))

        self.nodal_variables = [self.nodal_variables[i] for i in range(0,len(self.nodal_variables)) if self.nodal_variables[i] != 'NOT_DEFINED']
        self.nodal_variables.sort()

        print(" Variables :",self.nodal_variables)

        for variable in self.nodal_variables:
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))
            #print(" Added variable ", KratosMultiphysics.KratosGlobals.GetVariable(variable),"(",variable,")")

        print(self._class_prefix()+" General Variables ADDED")


    def _set_input_variables(self):
        variables_list = self.settings["variables"]
        for i in range(0, variables_list.size() ):
            self.nodal_variables.append(variables_list[i].GetString())

    #
    def _execute_after_reading(self):

        # Build bodies
        if self._has_bodies():
            self._build_bodies()

        # Build solving model parts
        self._build_solving_model_part()

        # Build composite solving model parts
        self._build_composite_solving_parts()
        self._update_composite_solving_parts()

        # Build output model part
        #self._create_sub_model_part(self.settings["output_model_part"].GetString())

    #
    def _build_bodies(self):

        #construct body model parts:
        solid_body_model_parts = []
        fluid_body_model_parts = []
        rigid_body_model_parts = []

        void_flags = []

        bodies_list = self.settings["bodies_list"]
        for i in range(bodies_list.size()):
            #create body model part
            body_model_part_name = bodies_list[i]["body_name"].GetString()
            self.main_model_part.CreateSubModelPart(body_model_part_name)
            body_model_part = self.main_model_part.GetSubModelPart(body_model_part_name)

            print(self._class_prefix()+" Body Created: "+body_model_part_name)
            body_model_part.ProcessInfo = self.main_model_part.ProcessInfo
            body_model_part.Properties  = self.main_model_part.Properties

            #build body from their parts
            body_parts_name_list = bodies_list[i]["parts_list"]
            body_parts_list = []
            for j in range(body_parts_name_list.size()):
                body_parts_list.append(self.main_model_part.GetSubModelPart(body_parts_name_list[j].GetString()))

            body_model_part_type = bodies_list[i]["body_type"].GetString()

            for part in body_parts_list:
                entity_type = "Nodes"
                if (body_model_part_type=="Fluid"):
                    part.Set(KratosMultiphysics.FLUID)
                    assign_flags = [KratosMultiphysics.FLUID]
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                    transfer_process.Execute()
                elif (body_model_part_type=="Solid"):
                    part.Set(KratosMultiphysics.SOLID)
                    assign_flags = [KratosMultiphysics.SOLID]
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                    transfer_process.Execute()
                elif (body_model_part_type=="Rigid"):
                    part.Set(KratosMultiphysics.RIGID)
                    assign_flags = [KratosMultiphysics.RIGID]
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                    transfer_process.Execute()

                entity_type = "Elements"
                transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type)
                transfer_process.Execute()
                entity_type = "Conditions"
                transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type)
                transfer_process.Execute()

            if( body_model_part_type == "Solid" ):
                body_model_part.Set(KratosMultiphysics.SOLID)
                solid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
            if( body_model_part_type == "Fluid" ):
                body_model_part.Set(KratosMultiphysics.FLUID)
                fluid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
            if( body_model_part_type == "Rigid" ):
                body_model_part.Set(KratosMultiphysics.RIGID)
                rigid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))

    #
    def _build_solving_model_part(self):

        # The solving_model_part is labeled 'KratosMultiphysics.ACTIVE' flag (in order to recover it)
        self._create_sub_model_part(self.settings["solving_model_part"].GetString())

        solving_model_part_name    = self.settings["solving_model_part"].GetString()
        domain_model_part_names    = self.settings["domain_parts_list"]
        processes_model_part_names = self.settings["processes_parts_list"]

        fluid_parts = False
        solid_parts = False
        domain_parts = []
        for i in range(domain_model_part_names.size()):
            domain_part = self.main_model_part.GetSubModelPart(domain_model_part_names[i].GetString())
            if( domain_part.Is(KratosMultiphysics.FLUID) ):
                fluid_parts = True
            elif( domain_part.Is(KratosMultiphysics.SOLID) ):
                solid_parts = True

            domain_parts.append(domain_part)


        processes_parts = []
        for i in range(processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(processes_model_part_names[i].GetString()))

        solving_model_part = self.main_model_part.GetSubModelPart(solving_model_part_name)
        solving_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        solving_model_part.Properties  = self.main_model_part.Properties

        #set flag to identify the fluid/solid body parts in the computing domain
        if( solid_parts ):
            solving_model_part.Set(KratosMultiphysics.SOLID)
        if( fluid_parts ):
            solving_model_part.Set(KratosMultiphysics.FLUID)

        #set flag to identify the computing model part
        solving_model_part.Set(KratosMultiphysics.ACTIVE)

        entity_type = "Nodes"
        transfer_process = KratosSolid.TransferEntitiesProcess(solving_model_part,self.main_model_part,entity_type)
        transfer_process.Execute()

        for part in domain_parts:
            entity_type = "Elements"
            transfer_process = KratosSolid.TransferEntitiesProcess(solving_model_part,part,entity_type)
            transfer_process.Execute()

        for part in processes_parts:
            part.Set(KratosMultiphysics.BOUNDARY)
            entity_type = "Conditions"
            #condition flags as BOUNDARY or CONTACT are reserved to composite or contact conditions (do not set it here)
            transfer_process = KratosSolid.TransferEntitiesProcess(solving_model_part,part,entity_type)
            transfer_process.Execute()

    #
    def _build_composite_solving_parts(self):

        print(self._class_prefix()+" Composite Solving Parts")
        solving_parts = self.settings["composite_solving_parts"]
        for i in range(0,solving_parts.size()):
            print(self._class_prefix()+" Build Part: "+solving_parts[i]["model_part_name"].GetString())
            solving_part_transfer = KratosSolid.TransferSolvingModelPartProcess(self.main_model_part,solving_parts[i])
            self.transfer_solving_parts.append(solving_part_transfer)

    #
    def _update_composite_solving_parts(self):
        self.current_update_time = self.process_info[KratosMultiphysics.TIME]
        print(self._class_prefix()+" Update Solving Parts")
        for transfer in self.transfer_solving_parts:
            transfer.Execute()

    #
    def _clean_body_parts(self):

        #delete body parts: (materials have to be already assigned)
        if( self._has_bodies() ):
            bodies_list = self.settings["bodies_list"]
            for i in range(bodies_list.size()):
                #get body parts
                body_parts_name_list = bodies_list[i]["parts_list"]
                for j in range(body_parts_name_list.size()):
                    self.main_model_part.RemoveSubModelPart(body_parts_name_list[j].GetString())
                    #print(self._class_prefix()+" Body Part Removed: "+ body_parts_name_list[j].GetString())

    #
    def _domain_parts_updated(self):
        update_time = False
        if not self._is_not_restarted():
            if self.process_info.Has(KratosSolid.RESTART_STEP_TIME):
                update_time = self._check_current_time_step(self.process_info[KratosSolid.RESTART_STEP_TIME])
                #print(" RESTART_STEP_TIME ",self.process_info[KratosSolid.RESTART_STEP_TIME], update_time)

        if not update_time and self.process_info.Has(KratosSolid.MESHING_STEP_TIME):
            update_time = self._check_previous_time_step(self.process_info[KratosSolid.MESHING_STEP_TIME])
            #print(" MESHING_STEP_TIME ",self.process_info[KratosSolid.MESHING_STEP_TIME], update_time)


        if not update_time and self.process_info.Has(KratosSolid.CONTACT_STEP_TIME):
            update_time = self._check_previous_time_step(self.process_info[KratosSolid.CONTACT_STEP_TIME])
            #print(" CONTACT_STEP_TIME ",self.process_info[KratosSolid.CONTACT_STEP_TIME], update_time)

        if update_time:
            update_time  = not self._check_current_time_step(self.current_update_time)

        return update_time
    #
    def _check_current_time_step(self, step_time):
        current_time  = self.process_info[KratosMultiphysics.TIME]
        delta_time    = self.process_info[KratosMultiphysics.DELTA_TIME]
        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001

        if( step_time > current_time-tolerance and step_time < current_time+tolerance ):
            return True
        else:
            return False
    #
    def _check_previous_time_step(self, step_time):
        current_time  = self.process_info[KratosMultiphysics.TIME]
        delta_time    = self.process_info[KratosMultiphysics.DELTA_TIME]
        previous_time = current_time - delta_time

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001

        if( step_time > previous_time-tolerance and step_time < previous_time+tolerance ):
            return True
        else:
            return False
    #
    def _is_not_restarted(self):
        if self.process_info.Has(KratosMultiphysics.IS_RESTARTED):
            if self.process_info[KratosMultiphysics.IS_RESTARTED]:
                return False
            else:
                return True
        else:
            return True
    #
    def _has_bodies(self):
        if( self.settings.Has("bodies_list") ):
            if( self.settings["bodies_list"].size() > 0 ):
                return True
        return False

    #
    def _clean_previous_result_files(self):

        file_endings = [".post.bin",".post.msh",".post.res",".post.lst",".post.csv",".rest"]
        problem_path = os.getcwd()
        for file_end in file_endings:
            filelist = [f for f in os.listdir(problem_path) if f.endswith(file_end)]

            for f in filelist:
                try:
                    os.remove(f)
                except OSError:
                    pass
    #
    @classmethod
    def _class_prefix(self):
        header = "::[---Model_Manager---]::"
        return header
