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
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.EkateAuxiliaryApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ExternalConstitutiveLawsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *
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
        number_of_time_steps = 10
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
        perform_contact_analysis_flag = False
        penalty = 0.0
        maxuzawa = 0.0
        friction = 0.0
        frictionpenalty = 0.0
        contact_double_check_flag = False
        contact_ramp_penalties_flag = False
        maxpenalty = 0.0
        rampcriterion = 0.0
        rampfactor = 0.0
        fricmaxpenalty = 0.0
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
        #PrintSparsityInfoFlag
        self.analysis_parameters.append(False)
        
        abs_tol =        1e-06
        #rel_tol =       0.0001
        rel_tol = 1e-10
        
        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        #import ekate_solver_parallel
        #self.solver = ekate_solver_parallel.EkateSolver( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
        #ekate_solver_parallel.AddVariables( self.model_part )
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False

        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
        print "input data read OK"
        #print "+++++++++++++++++++++++++++++++++++++++"
        #for node in self.model_part.Nodes:
        #    print node
        #print "+++++++++++++++++++++++++++++++++++++++"
        
        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################        
        for node in self.model_part.Nodes:
            node.AddDof( WATER_PRESSURE )
        structural_solver_advanced.AddDofs( self.model_part )
        #ekate_solver_parallel.AddDofs( self.model_part )

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = MKLPardisoSolver()
        #plinear_solver = ParallelMKLPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2);

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )
        
    def SetUpActivationLevels( self, model_part, activation_list, cond_activation_list ):
        for element in self.model_part.Elements:
            element.SetValue(ACTIVATION_LEVEL, activation_list[element.Id])
        for condition in self.model_part.Conditions:
            if( not (condition.GetValue(IS_TYING_MASTER) or condition.GetValue(IS_CONTACT_MASTER) ) ):
                condition.SetValue(ACTIVATION_LEVEL, activation_list[cond_activation_list[condition.Id-1]])

#    def write_restart_file( self, time ):
#        print("------------> restart file written for time step: "+str(time))
#        self.restart_utility.ChangeFileName(problem_name+str(time))
#        self.restart_utility.StoreNodalVariables(model_part)
#        self.restart_utility.StoreInSituStress(model_part)
#        self.restart_utility.StoreConstitutiveLawVariables(model_part)
#
#    def restart_time_step( self, time, Dt ):
#        print("############ time step solution has to be restarted ############")
#        time = time-Dt
#        model_part.CloneTimeStep(time)
#        for step in range(1,11):
#            time = time+ Dt/10.0
#            model_part.CloneTimeStep(time)
#            #####################################################################################################
#            model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, True )
#            model_part.ProcessInfo.SetValue( FIRST_TIME_STEP, False )
#            #####################################################################################################
#            solver.Solve()
#            print("~~~~~~~~~~~~~~ RESTARTED STEP ( DT= "+str(Dt/10.0)+" / Step= "+str(step)+" ) ~~~~~~~~~~~~~~")
#        print("############ restart finished ############")
#
#    def write_to_file( self, time ):
#        for i in range(0, len(self.layer_nodes_sets['top'])):
#            settlements.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetSolutionStepValue(DISPLACEMENT_Z))+"\n")
#    for i in range(0, len(layer_nodes_sets['side'])):
#        pressure_air.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(AIR_PRESSURE))+"\n")
#        pressure_water.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(WATER_PRESSURE))+"\n")
#

#    def FixPressureNodes( self, free_node_list_water, free_node_list_air):
#        for i in range(1, len(self.model_part.Nodes)+1):
#            if i in self.model_part.Nodes:
#                i_node = self.model_part.Nodes[i]
#                if i_node.HasDofFor(WATER_PRESSURE):
#                    if (i_node.IsFixed(WATER_PRESSURE)==0):        
#                        i_node.Fix(WATER_PRESSURE)
#                        free_node_list_water.append(i)
#                if i_node.HasDofFor(AIR_PRESSURE):
#                    if (i_node.IsFixed(AIR_PRESSURE)==0):
#                        i_node.Fix(AIR_PRESSURE)
#                        free_node_list_air.append(i)                

    def FixPressureNodes( self, free_node_list_water, free_node_list_air):
        for node in self.model_part.Nodes:
            if (node.IsFixed(WATER_PRESSURE)==0):                
                node.Fix(WATER_PRESSURE)
                free_node_list_water.append(node)
            if (node.IsFixed(AIR_PRESSURE)==0):
                node.Fix(AIR_PRESSURE)
                free_node_list_air.append(node)                                

    def ApplyInsituWaterPressure( self, free_node_list_water, free_node_list_air, z_zero, gravity_z):
        water_density=1000.0;
        for node in self.model_part.Nodes:                              
            water_pressure= water_density*gravity_z*(z_zero-(node.Z-node.GetSolutionStepValue(DISPLACEMENT_Z,0)))
            node.SetSolutionStepValue(WATER_PRESSURE, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_EINS, water_pressure)
            node.SetSolutionStepValue(WATER_PRESSURE_NULL, water_pressure)
        for node in self.model_part.Nodes:              
            node.SetSolutionStepValue(AIR_PRESSURE, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_EINS, 0.0)
            node.SetSolutionStepValue(AIR_PRESSURE_NULL, 0.0)

#        def FreePressureNodes(self,free_node_list_water, free_node_list_air):
#                for item in free_node_list_water:
#                        self.model_part.Nodes[item].Free(WATER_PRESSURE)
#                        #item.Free(WATER_PRESSURE)
#                for item in free_node_list_air:
#                        self.model_part.Nodes[item].Free(AIR_PRESSURE)
#                        #item.Free(AIR_PRESSURE)

    def FreePressureNodes(self,free_node_list_water, free_node_list_air):
        for item in free_node_list_water:
            #self.model_part.Nodes[item].Free(WATER_PRESSURE)
            item.Free(WATER_PRESSURE)
        for item in free_node_list_air:
            #self.model_part.Nodes[item].Free(AIR_PRESSURE)
            item.Free(AIR_PRESSURE)
            
    def WriteMaterialParameters( self, time, indices ):
        self.gid_io.OpenResultFile( self.path+self.problem_name, GiDPostMode.GiD_PostBinary)
        #self.gid_io.ChangeOutputName( self.path+self.problem_name +str(time), GiDPostMode.GiD_PostBinary )
        for index in indices:
            self.gid_io.SuperPrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, index)
        self.gid_io.CloseResultFile()

    def WriteMonitoringSectionResults( self, time ):
        outfile = open("step_"+str(time)+".dat",'w')
        outfile.write("ekate result file for step "+str(time)+"\n")
        outfile.close()
        
    def WriteOutput( self, time ):
                self.gid_io.InitializeMesh( time )
                mesh = self.model_part.GetMesh()
                #self.gid_io.WriteNodeMesh( mesh )
                self.gid_io.WriteMesh( mesh )
                print("mesh written...")
                self.gid_io.FinalizeMesh()
                self.gid_io.InitializeResults( time, self.model_part.GetMesh() )
                #self.gid_io.PrintOnGaussPoints(PLASTICITY_INDICATOR, self.model_part, time, 0)
                #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 0)
                #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 1)
                #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 2)
                #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 3)
                #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 4)
                #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 5)
                #self.gid_io.PrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part, time, 6)
                print("write nodal displacements")
                self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.Nodes, time, 0)
                self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
                self.gid_io.PrintOnGaussPoints(WATER_PRESSURE, self.model_part, time)
                #self.gid_io.WriteNodalResults(WATER_PRESSURE, self.model_part.Nodes, time, 0)
                self.gid_io.PrintOnGaussPoints(EXCESS_PORE_WATER_PRESSURE, self.model_part, time)
                self.gid_io.PrintOnGaussPoints(SATURATION, self.model_part, time)
                #self.gid_io.PrintOnGaussPoints(CONTACT_PENETRATION, self.model_part, time)
                #self.gid_io.PrintOnGaussPoints(NORMAL, self.model_part, time, 0)
                self.gid_io.FinalizeResults()
                
    def InitializeModel( self, linear_elastic ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        self.model_part.Properties[1].SetValue(DENSITY,         2000 )        
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,      1.5e+09 )        
        self.model_part.Properties[1].SetValue(POISSON_RATIO,          0.3 )
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
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
        ]
        self.layer_sets['Layer0'] = layer_elements_list
        layer_elements_list = [
        ]
        self.layer_sets['surface'] = layer_elements_list
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
        310 ,
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
        469 ,
        470 ,
        471 ,
        472 ,
        473 ,
        474 ,
        475 ,
        476 ,
        477 ,
        478 ,
        479 ,
        480 ,
        481 ,
        482 ,
        483 ,
        484 ,
        485 ,
        486 ,
        487 ,
        488 ,
        489 ,
        490 ,
        491 ,
        492 ,
        493 ,
        494 ,
        495 ,
        498 ,
        499 ,
        500 ,
        501 ,
        502 ,
        503 ,
        504 ,
        505 ,
        506 ,
        507 ,
        508 ,
        509 ,
        510 ,
        511 ,
        512 ,
        513 ,
        514 ,
        515 ,
        516 ,
        517 ,
        518 ,
        519 ,
        520 ,
        521 ,
        522 ,
        523 ,
        524 ,
        525 ,
        526 ,
        527 ,
        528 ,
        529 ,
        530 ,
        531 ,
        532 ,
        533 ,
        534 ,
        535 ,
        536 ,
        538 ,
        540 ,
        541 ,
        542 ,
        543 ,
        544 ,
        545 ,
        546 ,
        547 ,
        548 ,
        549 ,
        550 ,
        551 ,
        552 ,
        553 ,
        556 ,
        557 ,
        558 ,
        559 ,
        560 ,
        561 ,
        562 ,
        563 ,
        564 ,
        565 ,
        566 ,
        567 ,
        568 ,
        569 ,
        570 ,
        571 ,
        572 ,
        573 ,
        574 ,
        575 ,
        576 ,
        577 ,
        578 ,
        579 ,
        580 ,
        581 ,
        582 ,
        583 ,
        584 ,
        585 ,
        586 ,
        587 ,
        588 ,
        589 ,
        590 ,
        591 ,
        594 ,
        595 ,
        596 ,
        597 ,
        598 ,
        599 ,
        600 ,
        601 ,
        602 ,
        603 ,
        604 ,
        605 ,
        606 ,
        607 ,
        608 ,
        609 ,
        610 ,
        611 ,
        612 ,
        613 ,
        614 ,
        615 ,
        616 ,
        617 ,
        618 ,
        619 ,
        620 ,
        621 ,
        622 ,
        623 ,
        624 ,
        625 ,
        626 ,
        627 ,
        628 ,
        629 ,
        630 ,
        631 ,
        632 ,
        633 ,
        634 ,
        635 ,
        636 ,
        637 ,
        638 ,
        639 ,
        640 ,
        641 ,
        642 ,
        643 ,
        644 ,
        645 ,
        647 ,
        649 ,
        652 ,
        653 ,
        654 ,
        655 ,
        656 ,
        657 ,
        658 ,
        659 ,
        660 ,
        661 ,
        662 ,
        663 ,
        664 ,
        665 ,
        666 ,
        667 ,
        668 ,
        669 ,
        670 ,
        671 ,
        672 ,
        673 ,
        674 ,
        675 ,
        676 ,
        677 ,
        678 ,
        679 ,
        680 ,
        681 ,
        682 ,
        683 ,
        684 ,
        685 ,
        686 ,
        687 ,
        690 ,
        691 ,
        692 ,
        693 ,
        694 ,
        695 ,
        696 ,
        697 ,
        698 ,
        699 ,
        700 ,
        701 ,
        702 ,
        703 ,
        704 ,
        705 ,
        706 ,
        707 ,
        708 ,
        709 ,
        710 ,
        711 ,
        712 ,
        713 ,
        714 ,
        716 ,
        717 ,
        718 ,
        719 ,
        720 ,
        721 ,
        722 ,
        723 ,
        724 ,
        725 ,
        726 ,
        727 ,
        728 ,
        729 ,
        730 ,
        731 ,
        732 ,
        733 ,
        734 ,
        735 ,
        736 ,
        737 ,
        738 ,
        739 ,
        740 ,
        741 ,
        742 ,
        743 ,
        744 ,
        745 ,
        746 ,
        747 ,
        748 ,
        749 ,
        750 ,
        751 ,
        752 ,
        753 ,
        754 ,
        757 ,
        758 ,
        759 ,
        760 ,
        761 ,
        762 ,
        763 ,
        764 ,
        767 ,
        768 ,
        769 ,
        770 ,
        771 ,
        772 ,
        773 ,
        774 ,
        775 ,
        776 ,
        777 ,
        778 ,
        779 ,
        780 ,
        781 ,
        782 ,
        783 ,
        784 ,
        785 ,
        786 ,
        787 ,
        788 ,
        789 ,
        790 ,
        791 ,
        794 ,
        795 ,
        796 ,
        797 ,
        798 ,
        799 ,
        800 ,
        801 ,
        802 ,
        805 ,
        806 ,
        807 ,
        808 ,
        809 ,
        810 ,
        811 ,
        812 ,
        813 ,
        814 ,
        815 ,
        816 ,
        817 ,
        818 ,
        819 ,
        820 ,
        823 ,
        824 ,
        825 ,
        826 ,
        827 ,
        828 ,
        829 ,
        830 ,
        831 ,
        832 ,
        833 ,
        834 ,
        835 ,
        836 ,
        837 ,
        838 ,
        839 ,
        840 ,
        841 ,
        842 ,
        843 ,
        844 ,
        845 ,
        846 ,
        847 ,
        848 ,
        851 ,
        852 ,
        853 ,
        854 ,
        855 ,
        856 ,
        857 ,
        858 ,
        859 ,
        860 ,
        861 ,
        862 ,
        863 ,
        864 ,
        865 ,
        868 ,
        869 ,
        870 ,
        871 ,
        872 ,
        873 ,
        874 ,
        875 ,
        876 ,
        877 ,
        878 ,
        879 ,
        ]
        self.layer_nodes_sets['Layer0'] = layer_nodes_list
        layer_nodes_list = [
        291 ,
        309 ,
        311 ,
        360 ,
        361 ,
        384 ,
        385 ,
        424 ,
        446 ,
        447 ,
        496 ,
        497 ,
        537 ,
        539 ,
        554 ,
        555 ,
        592 ,
        593 ,
        646 ,
        648 ,
        650 ,
        651 ,
        688 ,
        689 ,
        715 ,
        755 ,
        756 ,
        765 ,
        766 ,
        792 ,
        793 ,
        803 ,
        804 ,
        821 ,
        822 ,
        849 ,
        850 ,
        866 ,
        867 ,
        880 ,
        ]
        self.layer_nodes_sets['surface'] = layer_nodes_list
        ## CONTACT MASTER NODES ##########################################
        self.contact_master_nodes = [
        ]
        ## CONTACT SLAVE NODES ###########################################
        self.contact_slave_nodes = [
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
        ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
        ##################################################################
        self.node_groups = {}
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        #self.SetUpActivationLevels( self.model_part, self.activation_flags, self.cond_activation_flags )
        self.deac.Initialize( self.model_part )
        print "activation utility initialized"
        ##################################################################
        ## MESH TYING ####################################################
        ##################################################################
        #self.mesh_tying_utility= MeshTyingUtility()
        ##self.mesh_tying_utility.InitializeMeshTyingUtilityLagrange(self.model_part)
        #self.mesh_tying_utility.InitializeMeshTyingUtility(self.model_part)
        #print "mesh-tying utility successfully initialized"
        #self.model_part.Check( self.model_part.ProcessInfo )
        print "model successfully initialized"

        
    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()
        
    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
##################################################################
