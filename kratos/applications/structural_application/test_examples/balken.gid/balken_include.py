##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
##### note: KRATOS is released under LGPL                    #####
##################################################################
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
kernel = Kernel()   #defining kernel
##################################################################
##################################################################
class Model:
        def __init__( self, problem_name, path ):
                #setting the domain size for the problem to be solved
                self.domain_size = 3
                ##################################################################
                ## DEFINE MODELPART ##############################################
                ##################################################################
                self.model_part = ModelPart("ekate_simulation")
                self.path = path
                self.problem_name = problem_name
                ##################################################################
                ## DEFINE SOLVER #################################################
                ##################################################################
                # reading simulation parameters
                number_of_time_steps = 1
                self.analysis_parameters = []
                # content of analysis_parameters:
                # perform_contact_analysis_flag
                # penalty value for normal contact
                # maximum number of uzawa iterations
                # friction coefficient
                # penalty value for frictional contact
                # contact_double_check_flag
                # contact_ramp_penalties_flag
                # maximum penalty value for normal contact
                # ramp criterion for normal contact
                # ramp factor for normal contact
                # maximum penalty value for frictional contact
                # ramp criterion for frictional contact
                # ramp factor for frictional contact
                perform_contact_analysis_flag = True
                # performing contact analysis: reading contact parameters
                penalty =        1e+10
                maxuzawa = 25
                friction =            0
                frictionpenalty =        1e+05
                contact_double_check_flag = False
                contact_ramp_penalties_flag = False
                maxpenalty = penalty
                rampcriterion = 0.0
                rampfactor = 0.0
                fricmaxpenalty = penalty
                fricrampcriterion = 0.0
                fricrampfactor = 0.0
                self.analysis_parameters.append(perform_contact_analysis_flag)
                self.analysis_parameters.append(penalty)
                self.analysis_parameters.append(maxuzawa)
                self.analysis_parameters.append(friction)
                self.analysis_parameters.append(frictionpenalty)
                self.analysis_parameters.append(contact_double_check_flag)
                self.analysis_parameters.append(contact_ramp_penalties_flag)
                self.analysis_parameters.append(maxpenalty)
                self.analysis_parameters.append(rampcriterion)
                self.analysis_parameters.append(rampfactor)
                self.analysis_parameters.append(fricmaxpenalty)
                self.analysis_parameters.append(fricrampcriterion)
                self.analysis_parameters.append(fricrampfactor)
		self.analysis_parameters.append(False)
                
                abs_tol =        1e-06
                rel_tol =       0.0001
                
                ## generating solver
                import structural_solver_advanced
                self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
                structural_solver_advanced.AddVariables( self.model_part )

                ##################################################################
                ## READ MODELPART ################################################
                ##################################################################
                #reading a model
                write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
                write_elements = WriteConditionsFlag.WriteElementsOnly
                post_mode = GiDPostMode.GiD_PostBinary
                multi_file_flag = MultiFileFlag.SingleFile
                self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
                self.gid_io.ReadModelPart(self.model_part)
                self.meshWritten = False
                ## READ DEACTIVATION FILE ########################################
                self.deac_file = open(self.path+self.problem_name+".deac",'r' )
                self.activation_flags = [0]
                for line in self.deac_file:
                        val_set = line.split(' ')
                        elem_num = int(val_set[0])
                        act_level = int(val_set[1])
                        self.activation_flags.append(act_level)
                print "input data read OK"
                #print "+++++++++++++++++++++++++++++++++++++++"
                #for node in self.model_part.Nodes:
                #        print node
                #print "+++++++++++++++++++++++++++++++++++++++"
                
                #the buffer size should be set up here after the mesh is read for the first time
                self.model_part.SetBufferSize(2)

                ##################################################################
                ## ADD DOFS ######################################################
                ##################################################################                
                structural_solver_advanced.AddDofs( self.model_part )

                ##################################################################
                ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
                ##################################################################
                #defining linear solver
                #plinear_solver = SkylineLUFactorizationSolver()
                plinear_solver = MKLPardisoSolver()
                self.solver.structure_linear_solver = plinear_solver
                self.solver.Initialize()
                (self.solver.solver).SetEchoLevel(2);

                ##################################################################
                ## INITIALISE RESTART UTILITY ####################################
                ##################################################################
                #restart_utility= RestartUtility( self.problem_name )
                
        def SetUpActivationLevels( self, model_part, activation_list ):
                for element in self.model_part.Elements:
                        element.SetValue(ACTIVATION_LEVEL, activation_list[element.Id])

#        def write_restart_file( self, time ):
#                print("------------> restart file written for time step: "+str(time))
#                self.restart_utility.ChangeFileName(problem_name+str(time))
#                self.restart_utility.StoreNodalVariables(model_part)
#                self.restart_utility.StoreInSituStress(model_part)
#                self.restart_utility.StoreConstitutiveLawVariables(model_part)
#
#        def restart_time_step( self, time, Dt ):
#                print("############ time step solution has to be restarted ############")
#                time = time-Dt
#                model_part.CloneTimeStep(time)
#                for step in range(1,11):
#                        time = time+ Dt/10.0
#                        model_part.CloneTimeStep(time)
#                        #####################################################################################################
#                        model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, True )
#                        model_part.ProcessInfo.SetValue( FIRST_TIME_STEP, False )
#                        #####################################################################################################
#                        solver.Solve()
#                        print("~~~~~~~~~~~~~~ RESTARTED STEP ( DT= "+str(Dt/10.0)+" / Step= "+str(step)+" ) ~~~~~~~~~~~~~~")
#                print("############ restart finished ############")
#
#        def write_to_file( self, time ):
#                for i in range(0, len(self.layer_nodes_sets['top'])):
#                        settlements.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetSolutionStepValue(DISPLACEMENT_Z))+"\n")
#        for i in range(0, len(layer_nodes_sets['side'])):
#                pressure_air.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(AIR_PRESSURE))+"\n")
#                pressure_water.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(WATER_PRESSURE))+"\n")
#

        def ApplyInsituWaterPressure( self, free_node_list_water, free_node_list_air, z_zero):
                gravity_z =          -50;
                water_density= 1000.0;
                for i in range(1, len(self.model_part.Nodes)+1):
                        if self.model_part.Nodes[i].HasDofFor(WATER_PRESSURE):
                                if (self.model_part.Nodes[i].IsFixed(WATER_PRESSURE)==0):                
                                        water_pressure= 0.0
                                        self.model_part.Nodes[i].SetSolutionStepValue(WATER_PRESSURE, water_pressure)
                                        self.model_part.Nodes[i].SetSolutionStepValue(WATER_PRESSURE_EINS, water_pressure)
                                        self.model_part.Nodes[i].SetSolutionStepValue(WATER_PRESSURE_NULL, water_pressure)
                                        self.model_part.Nodes[i].Fix(WATER_PRESSURE)
                                        free_node_list_water.append(i)
                        if self.model_part.Nodes[i].HasDofFor(AIR_PRESSURE):
                                if (self.model_part.Nodes[i].IsFixed(AIR_PRESSURE)==0):
                                        self.model_part.Nodes[i].SetSolutionStepValue(AIR_PRESSURE, 0.0)
                                        self.model_part.Nodes[i].SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
                                        self.model_part.Nodes[i].SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)
                                        self.model_part.Nodes[i].Fix(AIR_PRESSURE)
                                        free_node_list_air.append(i)

        def FreePressureNodes(self,free_node_list_water, free_node_list_air):
                for item in free_node_list_water:
                        self.model_part.Nodes[item].Free(WATER_PRESSURE)
                for item in free_node_list_air:
                        self.model_part.Nodes[item].Free(AIR_PRESSURE)
                        
        def WriteMaterialParameters( self, time, indices ):
                self.gid_io.OpenResultFile( self.path+self.problem_name, GiDPostMode.GiD_PostBinary)
                for index in indices:
                        self.gid_io.SuperPrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, index)
                self.gid_io.CloseResultFile()
                
        def WriteOutput( self, time ):
                if( self.meshWritten == False ):
                        self.gid_io.InitializeMesh( 0.0 )
                        mesh = self.model_part.GetMesh()
                        self.gid_io.WriteMesh( mesh )
                        self.meshWritten = True
                        self.gid_io.FinalizeMesh()
                self.gid_io.InitializeResults( 0.0, self.model_part.GetMesh() )
                print("write nodal displacements")
                self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.Nodes, time, 0)
                self.gid_io.FinalizeResults()
                                
        def InitializeModel( self ):
                ##################################################################
                ## INITIALISE CONSTITUTIVE LAWS ##################################
                ##################################################################
                #set material parameters
                params1 = Vector(7)
                params1[0] =      2.1e+11 #Young's Modulus
                params1[1] =          0.3 #Poisson Ratio
                params1[2] = 0.0                          #Internal Friction Angle
                params1[3] = 0.0                          #Cohesion
                params1[4] = 0.0                          #Compressive Strength 
                params1[5] = 0.0                          #Tensile Strength
                params1[6] = 0.0                          #RMR Variance
                self.model_part.Properties[1].SetValue(MATERIAL_PARAMETERS, params1 )
                params2 = Vector(7)
                params2[0] =        3e+10 #Young's Modulus
                params2[1] =          0.3 #Poisson Ratio
                params2[2] = 0.0                          #Internal Friction Angle
                params2[3] = 0.0                          #Cohesion
                params2[4] = 0.0                          #Compressive Strength 
                params2[5] = 0.0                          #Tensile Strength
                params2[6] = 0.0                          #RMR Variance
                self.model_part.Properties[2].SetValue(MATERIAL_PARAMETERS, params2 )
                self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
                print "Linear elastic model selected"
                self.model_part.Properties[2].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
                print "Linear elastic model selected"
                ##################################################################
                ## STORE LAYER SETS ##############################################
                ##################################################################
                ## ELEMENTS on layers ############################################
                self.layer_sets = {}
                layer_elements_list = [
                1 ,
                2 ,
                3 ,
                4 ,
                5 ,
                6 ,
                7 ,
                8 ,
                9 ,
                10 ,
                11 ,
                12 ,
                13 ,
                14 ,
                15 ,
                16 ,
                17 ,
                18 ,
                19 ,
                20 ,
                21 ,
                22 ,
                23 ,
                24 ,
                25 ,
                ]
                self.layer_sets['Layer0'] = layer_elements_list
                ## ELEMENTS on inner boundaries ##################################
                self.inner_boundary_elements = [
                ]
                ## NODES on layers ###############################################
                self.layer_nodes_sets = {}
                layer_nodes_list = [
            1 ,
            2 ,
            3 ,
            4 ,
            5 ,
            6 ,
            7 ,
            8 ,
            9 ,
            10 ,
            11 ,
            12 ,
            13 ,
            14 ,
            15 ,
            16 ,
            17 ,
            18 ,
            19 ,
            20 ,
            21 ,
            22 ,
            23 ,
            24 ,
            25 ,
            26 ,
            27 ,
            28 ,
            29 ,
            30 ,
            31 ,
            32 ,
            33 ,
            34 ,
            35 ,
            36 ,
            37 ,
            38 ,
            39 ,
            40 ,
            41 ,
            42 ,
            43 ,
            44 ,
            45 ,
            46 ,
            47 ,
            48 ,
            49 ,
            50 ,
            51 ,
            52 ,
            53 ,
            54 ,
            55 ,
            56 ,
            57 ,
            58 ,
            59 ,
            60 ,
            61 ,
            62 ,
            63 ,
            64 ,
            65 ,
            66 ,
            67 ,
            68 ,
            69 ,
            70 ,
            71 ,
            72 ,
            73 ,
            74 ,
            75 ,
            76 ,
            77 ,
            78 ,
            79 ,
            80 ,
            81 ,
            82 ,
            83 ,
            84 ,
            85 ,
            86 ,
            87 ,
            88 ,
            89 ,
            90 ,
            91 ,
            92 ,
            93 ,
            94 ,
            95 ,
            96 ,
            97 ,
            98 ,
            99 ,
            100 ,
            101 ,
            102 ,
            103 ,
            104 ,
            105 ,
            106 ,
            107 ,
            108 ,
            109 ,
            110 ,
            111 ,
            112 ,
            113 ,
            114 ,
            115 ,
            116 ,
            117 ,
            118 ,
            119 ,
            120 ,
            121 ,
            122 ,
            123 ,
            124 ,
            125 ,
            126 ,
            127 ,
            128 ,
            129 ,
            130 ,
            131 ,
            132 ,
            133 ,
            134 ,
            135 ,
            136 ,
            137 ,
            138 ,
            139 ,
            140 ,
            141 ,
            142 ,
            143 ,
            144 ,
            145 ,
            146 ,
            147 ,
            148 ,
            149 ,
            150 ,
            151 ,
            152 ,
            153 ,
            154 ,
            155 ,
            156 ,
            157 ,
            158 ,
            159 ,
            160 ,
            161 ,
            162 ,
            163 ,
            164 ,
            165 ,
            166 ,
            167 ,
            168 ,
            169 ,
            170 ,
            171 ,
            172 ,
            173 ,
            174 ,
            175 ,
            176 ,
            177 ,
            178 ,
            179 ,
            180 ,
            181 ,
            182 ,
            183 ,
            184 ,
            185 ,
            186 ,
            187 ,
            188 ,
            189 ,
            190 ,
            191 ,
            192 ,
            193 ,
            194 ,
            195 ,
            196 ,
            197 ,
            198 ,
            199 ,
            200 ,
            201 ,
            202 ,
            203 ,
            204 ,
            205 ,
            206 ,
            207 ,
            208 ,
            209 ,
            210 ,
            211 ,
            212 ,
            213 ,
            214 ,
            215 ,
            216 ,
            217 ,
            218 ,
            219 ,
            220 ,
            221 ,
            222 ,
            223 ,
            224 ,
            225 ,
            226 ,
            227 ,
            228 ,
            229 ,
            230 ,
            231 ,
            232 ,
            233 ,
            234 ,
            235 ,
            236 ,
            237 ,
            238 ,
            239 ,
            240 ,
            241 ,
            242 ,
            243 ,
            244 ,
            245 ,
            246 ,
            247 ,
            248 ,
            249 ,
            250 ,
            251 ,
            252 ,
            253 ,
            254 ,
            255 ,
            256 ,
            257 ,
            258 ,
            259 ,
            260 ,
            261 ,
            262 ,
            263 ,
            264 ,
            265 ,
            266 ,
            267 ,
            268 ,
            269 ,
            270 ,
            271 ,
            272 ,
            273 ,
            274 ,
            275 ,
            276 ,
            277 ,
            278 ,
            279 ,
            280 ,
            281 ,
            282 ,
            283 ,
            284 ,
            285 ,
            286 ,
            287 ,
            288 ,
            289 ,
            290 ,
            291 ,
            292 ,
            293 ,
            294 ,
            295 ,
            296 ,
            297 ,
            298 ,
            299 ,
            300 ,
            301 ,
            302 ,
            303 ,
            304 ,
            305 ,
            306 ,
            307 ,
            308 ,
            309 ,
            310 ,
            311 ,
            312 ,
            313 ,
            314 ,
            315 ,
            316 ,
            317 ,
            318 ,
            319 ,
            320 ,
            321 ,
            322 ,
            323 ,
            324 ,
            325 ,
            326 ,
            327 ,
            328 ,
            329 ,
            330 ,
            331 ,
            332 ,
            333 ,
            334 ,
            335 ,
            336 ,
            337 ,
            338 ,
            339 ,
            340 ,
            341 ,
            342 ,
            343 ,
            344 ,
            345 ,
            346 ,
            347 ,
            348 ,
            349 ,
            350 ,
            351 ,
            352 ,
            353 ,
            354 ,
            355 ,
            356 ,
            357 ,
            358 ,
            359 ,
            360 ,
            361 ,
            362 ,
            363 ,
            364 ,
            365 ,
            366 ,
            367 ,
            368 ,
            369 ,
            370 ,
            371 ,
            372 ,
            373 ,
            374 ,
            375 ,
            376 ,
            377 ,
            378 ,
            379 ,
            380 ,
            381 ,
            382 ,
            383 ,
            384 ,
            385 ,
            386 ,
            387 ,
            388 ,
            389 ,
            390 ,
            391 ,
            392 ,
            393 ,
            394 ,
            395 ,
            396 ,
            397 ,
            398 ,
            399 ,
            400 ,
            401 ,
            402 ,
            403 ,
            404 ,
            405 ,
            406 ,
            407 ,
            408 ,
            409 ,
            410 ,
            411 ,
            412 ,
            413 ,
            414 ,
            415 ,
            416 ,
            417 ,
            418 ,
            419 ,
            420 ,
            421 ,
            422 ,
            423 ,
            424 ,
            425 ,
            426 ,
            427 ,
            428 ,
            429 ,
            430 ,
            431 ,
            432 ,
            433 ,
            434 ,
            435 ,
            436 ,
            437 ,
            438 ,
            439 ,
            440 ,
            441 ,
            442 ,
            443 ,
            444 ,
            445 ,
            446 ,
            447 ,
            448 ,
            449 ,
            450 ,
            451 ,
            452 ,
            453 ,
            454 ,
            455 ,
            456 ,
            457 ,
            458 ,
            459 ,
            460 ,
            461 ,
            462 ,
            463 ,
            464 ,
            465 ,
            466 ,
            467 ,
            468 ,
                ]
                self.layer_nodes_sets['Layer0'] = layer_nodes_list
                ## CONTACT MASTER NODES ##########################################
                self.contact_master_nodes = [
                6 ,
                8 ,
                15 ,
                22 ,
                25 ,
                33 ,
                37 ,
                41 ,
                49 ,
                55 ,
                56 ,
                67 ,
                73 ,
                74 ,
                79 ,
                91 ,
                92 ,
                97 ,
                109 ,
                110 ,
                115 ,
                127 ,
                128 ,
                133 ,
                145 ,
                146 ,
                151 ,
                163 ,
                164 ,
                169 ,
                181 ,
                182 ,
                187 ,
                199 ,
                200 ,
                205 ,
                217 ,
                218 ,
                223 ,
                231 ,
                233 ,
                241 ,
                244 ,
                245 ,
                250 ,
                262 ,
                263 ,
                268 ,
                280 ,
                281 ,
                286 ,
                298 ,
                299 ,
                304 ,
                316 ,
                317 ,
                322 ,
                334 ,
                335 ,
                340 ,
                352 ,
                353 ,
                358 ,
                370 ,
                371 ,
                376 ,
                388 ,
                389 ,
                394 ,
                406 ,
                407 ,
                412 ,
                424 ,
                425 ,
                430 ,
                442 ,
                443 ,
                448 ,
                456 ,
                458 ,
                466 ,
                ]
                ## CONTACT SLAVE NODES ###########################################
                self.contact_slave_nodes = [
                3 ,
                7 ,
                12 ,
                20 ,
                21 ,
                32 ,
                40 ,
                45 ,
                52 ,
                62 ,
                66 ,
                72 ,
                83 ,
                87 ,
                90 ,
                102 ,
                105 ,
                108 ,
                120 ,
                123 ,
                126 ,
                138 ,
                141 ,
                144 ,
                156 ,
                159 ,
                162 ,
                174 ,
                177 ,
                180 ,
                192 ,
                195 ,
                198 ,
                210 ,
                213 ,
                216 ,
                228 ,
                232 ,
                239 ,
                253 ,
                257 ,
                261 ,
                273 ,
                276 ,
                279 ,
                291 ,
                294 ,
                297 ,
                309 ,
                312 ,
                315 ,
                327 ,
                330 ,
                333 ,
                345 ,
                348 ,
                351 ,
                363 ,
                366 ,
                369 ,
                381 ,
                384 ,
                387 ,
                399 ,
                402 ,
                405 ,
                417 ,
                420 ,
                423 ,
                435 ,
                438 ,
                441 ,
                453 ,
                457 ,
                464 ,
                ]
                ## INNER BOUNDARY NODES ##########################################
                self.inner_boundary_nodes = [
                ]
                ##################################################################
                print "layer sets stored"                
                ##################################################################
                ## STORE NODES ON GROUND SURFACE #################################
                ##################################################################
                self.top_surface_nodes = []
                print "nodes on ground surface stored"
                ##################################################################
                ## ACTIVATION ####################################################
                ##################################################################
                self.deac = DeactivationUtility()
                self.SetUpActivationLevels( self.model_part, self.activation_flags )
                self.deac.Initialize( self.model_part )
                print "activation utility initialized"
                self.SetCalculateInSituStress( False )
                
        def FinalizeModel( self ):
                self.gid_io.CloseResultFile()
                
        def SetCalculateInSituStress( self, calculation_flag ):
                self.insitu_stress_flag = calculation_flag

        def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
                self.deac.Reactivate( self.model_part, from_reac, to_reac )
                self.deac.Deactivate( self.model_part, from_deac, to_deac )
                self.model_part.CloneTimeStep(time)
                self.model_part.ProcessInfo.SetValue( CALCULATE_INSITU_STRESS, self.insitu_stress_flag )
                self.solver.Solve()
##################################################################
