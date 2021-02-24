"""
 --- apply_mass_local_conservation_process.py ---- Fri, May 22, 2020 11:27:27 AM ----
  Altair Manufacturing Solver

  Author: ddiez --- Maintained by: ddiez
  Copyright: Altair Engineering, Inc. 2015 - 2020
 ************************************************************
"""

from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory

if KratosMultiphysics.DataCommunicator.GetDefault().IsDistributed():
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as trilinos_linear_solver_factory



elementId=[568,505,700]
def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalMassConservationCheckProcess( Model, settings["Parameters"] )

class ApplyLocalMassConservationCheckProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                    : "fluid_model_part",
            "mass_correction_settings"           : {
                "echo_level" : 1
            },
            "write_to_log_file"                      : true,
            "log_file_name"                      : "mass_conservation.log",
            "convector_settings"                 : {},
            "check_volume_loss_at_the_end"       : false,
            "min_dt_factor"                      : 1e-5,
            "max_dt_factor"                      : 1.0,
            "correct_backwards"                  : true,
            "maximum_iterations"                 : 5
        }""" )

        settings.ValidateAndAssignDefaults(default_parameters)
        convector_settings = KratosMultiphysics.Parameters( """{
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "max_cfl" : 1.0,
            "cross_wind_stabilization_factor" : 0.7,
            "max_substeps" : 5
        }""")
        settings["convector_settings"].ValidateAndAssignDefaults(convector_settings)
        self.settings = settings

        if ( self.settings["model_part_name"].GetString() == "" ):
            KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","The value (string) of the parameter 'model_part_name' must not be empty.")

        self._fluid_model_part = Model[self.settings["model_part_name"].GetString()]
        self._write_to_log = settings["write_to_log_file"].GetBool()
        self._my_log_file = self.settings["log_file_name"].GetString()
        self._check_volume_loss_at_the_end = self.settings["check_volume_loss_at_the_end"].GetBool()
        self._correct_backwards = self.settings["correct_backwards"].GetBool()

        self._min_dt_factor = self.settings["min_dt_factor"].GetDouble()
        self._max_dt_factor = self.settings["max_dt_factor"].GetDouble()

        self._maximum_iterations = self.settings["maximum_iterations"].GetInt()
        self._is_printing_rank = ( self._fluid_model_part.GetCommunicator().MyPID() == 0 )
        self.mass_conservation_utility = KratosFluid.MassConservationUtility(self._fluid_model_part, self.settings["mass_correction_settings"])

        KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Construction finished.")

    def ExecuteFinalizeSolutionStepWriting(self):

        log_line_string = self.mass_conservation_utility.ComputeBalancedVolume()

        # writing first line in file
        if ( self._write_to_log and self._is_printing_rank ):
            with open(self._my_log_file, "a+") as logFile:
                logFile.write( log_line_string )

    # def ExecuteInitialize(self):
    #     self.forward_convection_process = self._set_levelset_convection_process()
    #     # first_lines_string=self.mass_conservation_utility.Initialize()
    #     if ( self._write_to_log and self._is_printing_rank ):
    #         with open(self._my_log_file, "w") as logFile:
    #              logFile.write( first_lines_string )
    #     KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Initialization finished (initial volumes calculated).")

    def  ExecuteBeforeSolutionLoop(self):
        self.forward_convection_process = self._set_levelset_convection_process()

        self.theoreticalVolume=self.mass_conservation_utility.CalculateWaterVolume()
        first_lines_string = self.mass_conservation_utility.CalculateInitialVolume()


        # if ( self._write_to_log and self._is_printing_rank ):
        #     with open(self._my_log_file, "w") as logFile:
        #          logFile.write( first_lines_string )

        KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Initialization finished (initial volumes calculated).")

    def Execute(self):
        self.ExecuteInitializeSolutionStep()
    
    def _DistanceGradientNorm(self,gradient):
        norm_distance_gradient=(gradient[0]**2+gradient[1]**2+gradient[2]**2)**(0.5)
        return norm_distance_gradient
    def _DotProductKratos(self, Vector1,Vector2):
        dot_product=Vector1[0]*Vector2[0]+Vector1[1]*Vector2[1]+Vector1[2]*Vector2[2]
        return dot_product
    def ExecuteInitializeSolutionStep(self):
        # Save current distance in case we do not improve
        KratosMultiphysics.VariableUtils().SaveScalarVar(
            KratosMultiphysics.DISTANCE,
            KratosFluid.AUX_DISTANCE,
            self._fluid_model_part.Nodes)
            
        # for node in self._fluid_model_part.Nodes:
        #     if node.Id in elementId:
        #         veloc=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
        #         print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        #         print(veloc)
        
        self.mass_conservation_utility.ComputeBalancedVolume()
        initial_volume_error = self.mass_conservation_utility.GetVolumeError()
        Flow_interface_corrected_old=self._CalculateFlowIntoAir()
        # print("2222222222222222222222222222222222222222222222____FLOW_____22222222222222222222222222222222222222222222222222")
        # print(Flow_interface_corrected_old)
        #Calculate de initial Volume error with the time step of the simulation
        dt_old = self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        # We change the initial error to old error in order to initialize mass conservation iteration
        print("initial_volume_error: ",initial_volume_error)
        old_volume_error = initial_volume_error
        
        # Volume error tolerance in order to correct mass
        user_error_volume_tolerance=1/100000
        tol_error = self.theoreticalVolume*user_error_volume_tolerance
        tol_flow=1e-05
        # tol_error=0.001
        print(tol_error)
        # mass conservation  correction 
        if abs(initial_volume_error) > tol_error and abs(Flow_interface_corrected_old)>tol_flow:
            iteration=0
            while iteration < self._maximum_iterations:
                print("we are correcting")
                flow_interface=self._CalculateFlowIntoAir() 
                volumetric_error=self.mass_conservation_utility.GetVolumeError()
                self._GetDistanceGradientProcess().Execute()
                global_corrector= self.mass_conservation_utility.GetInterfaceArea()
                print("AREA",global_corrector)
                for node in self._fluid_model_part.Nodes:
                    distance_gradient=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
                    distance_gradient_norm=self._DistanceGradientNorm(distance_gradient)
                    normalized_distance_gradient=distance_gradient*distance_gradient_norm
                    veolocity_nodal=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                    projected_velocity=self._DotProductKratos(veolocity_nodal,normalized_distance_gradient)
                    local_corrector=projected_velocity/flow_interface
                    
                    if abs(local_corrector)<abs(global_corrector):
                        corrector_distance=volumetric_error*abs(local_corrector)
                        
                    else: 
                        corrector_distance=volumetric_error/global_corrector
                    # corrector_distance=volumetric_error/global_corrector    
                    previous_distance=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                    new_distance= previous_distance+ corrector_distance
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,new_distance)
                
            
                self.mass_conservation_utility.ComputeBalancedVolume()
                new_error = self.mass_conservation_utility.GetVolumeError()

                if abs(new_error)< abs(initial_volume_error):
                    print("mass correction has been carried out succesfully")
                    print("PREVIUS_VOLUME", old_volume_error)
                    print("NEW_ERROR", new_error)
                else:
                    print("We are NOT improving mass conservation")
                    print("PREVIUS_VOLUME", old_volume_error)
                    print("NEW_ERROR", new_error)
                    return
                if abs(local_corrector)<abs(global_corrector):
                    print("WE ARE USING LOCAL CORRECTOR")
                else:
                    print("WE ARE USING GLOBAL CORRECTOR")

                if abs(new_error)>tol_error:
                    old_volume_error=new_error 
                    iteration+=1
                else:
                    print("MASS CONSERVATION CHECK IS COMPLETED")
                    return

                

                
    
        self.ExecuteFinalizeSolutionStepWriting()





    def _CreateDistanceGradientProcess(self):
        distance_gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(
                self._fluid_model_part,
                KratosMultiphysics.DISTANCE,
                KratosMultiphysics.DISTANCE_GRADIENT,
                KratosMultiphysics.NODAL_AREA)

        return distance_gradient_process



    def _GetDistanceGradientProcess(self):
        if not hasattr(self, '_distance_gradient_process'):
            self._distance_gradient_process = self._CreateDistanceGradientProcess()
        return self._distance_gradient_process




















        




        # # ADDING NEW THINGS 
        # # self.mass_conservation_utility.ReCheckTheMassConservation()
        # if  dt_max > 0:
        #     dt_low = 0.0
        #     dt_high = dt_max
        #     volume_error_low = initial_volume_error
        #     volume_error_high = new_volume_error
        # else:
        #     dt_low = dt_max
        #     dt_high = 0.0
        #     volume_error_low = new_volume_error
        #     volume_error_high = initial_volume_error
        # # self.iteration_is_on_min_side = abs(volume_error_high) > abs(volume_error_low)
        # # self.old_iteration_was_on_min_side = not self.iteration_is_on_min_side
        # # self.gamma = 1.0


        # if volume_error_low*volume_error_high > 0.0:
        #     if abs(volume_error_low) < abs(volume_error_high):
        #         self.mass_conservation_utility.RestoreDistanceValues(KratosFluid.AUX_DISTANCE)
        #     KratosMultiphysics.Logger.PrintInfo("Volume cannot be corrected, moving to global corrector")
        #     return

        # # iterations = 0

        # while iterations < self._maximum_iterations:
        #     # self.mass_conservation_utility.RestoreDistanceValues(KratosFluid.AUX_DISTANCE)
        #     dt = self._PredictNewDt(dt_low, dt_high, volume_error_low, volume_error_high)
        #     self._ConvectAuxiliaryDistance(dt)
        #     self.mass_conservation_utility.ReCheckTheMassConservation()
        #     new_volume_error = self.mass_conservation_utility.GetVolumeError()
        #     self.old_iteration_was_on_min_side = self.iteration_is_on_min_side
        #     if new_volume_error*volume_error_low > 0.0:
        #         volume_error_low = new_volume_error
        #         dt_low = dt  
        #     else:
        #         volume_error_high = new_volume_error
        #         dt_high = dt
                


                
        #     iterations += 1







        








       

    def _CalculateFlowIntoAir(self):
        print("entra")
        flow_orthogonal=self.mass_conservation_utility.OrthogonalFlowIntoAir()
        return flow_orthogonal

    # def _ComputeTimeStepForConvection(self,Flow_interface_corrected):
    #     """ REMARK: The pseudo time step dt has no meaning as a physical time step and is NOT accounted as such.
    #         Instead, it as purely used for an artificial convection that helps to conserve the mass.
    #         Accordingly, dt depends on the currently necessary mass convection.
    #     """
    #     dt = self.mass_conservation_utility.ComputeTimeStepForConvection(Flow_interface_corrected)

        # KratosMultiphysics.Logger.PrintInfo("Time step dt for convection in mass conservation: ", "{:e}".format(dt))
        # return dt


    def _ConvectAuxiliaryDistance(self, dt):
        """ The distance field is artificially convected to extrapolate the interface motion into the future.
            The time step considered for this artificial "forward convection" is dt.
            The previous distance field is saved in an auxiliary variable
        """
        i = 0
        phi_old_vect = KratosMultiphysics.Vector(self._fluid_model_part.NumberOfNodes())
        self.phi_prev_it_vect = KratosMultiphysics.Vector(self._fluid_model_part.NumberOfNodes())
        for node in self._fluid_model_part.Nodes:
            phi_now = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE, 0)
            phi_old = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE, 1)
            phi_old_vect[i] = phi_old
            self.phi_prev_it_vect[i] = phi_now
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1, phi_now)
            i += 1

        prev_dt = self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        revert_velocity = False
        if dt < 0.0:
            dt *= -1.0
            revert_velocity = True
        self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = dt
        if revert_velocity:
            self.mass_conservation_utility.RevertVelocityDirection()
        self.forward_convection_process.Execute()
        if revert_velocity:
            self.mass_conservation_utility.RevertVelocityDirection()
        self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = prev_dt

        for phi_old, node in zip(phi_old_vect, self._fluid_model_part.Nodes):
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1, phi_old)


    # def _ConvectAuxiliaryDistance(self, dt):
    #     """ The distance field is artificially convected to extrapolate the interface motion into the future.
    #         The time step considered for this artificial "forward convection" is dt.
    #         The previous distance field is saved in an auxiliary variable
    #     """
    #     prev_dt = self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
    #     revert_velocity = False
    #     if dt < 0.0:
    #         dt *= -1.0
    #         revert_velocity = True
    #     self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = dt
    #     if revert_velocity:
    #         self.mass_conservation_utility.RevertVelocityDirection()
    #     self.forward_convection_process.Execute()
    #     if revert_velocity:
    #         self.mass_conservation_utility.RevertVelocityDirection()
    #     self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = prev_dt


    def _PredictNewDt(self, dt_low, dt_high, volume_error_low, volume_error_high):
        return (volume_error_high*dt_low - volume_error_low*dt_high)/(volume_error_high - volume_error_low)

    def _set_levelset_convection_process(self):
        if self._fluid_model_part.IsDistributed():
            self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()
            self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["convector_settings"]["linear_solver_settings"])
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetConvectionProcess3D(
                self.EpetraCommunicator,
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.trilinos_linear_solver,
                self.settings["convector_settings"]["max_cfl"].GetDouble(),
                self.settings["convector_settings"]["cross_wind_stabilization_factor"].GetDouble(),
                self.settings["convector_settings"]["max_substeps"].GetInt())
        else:
            self.linear_solver = python_linear_solver_factory.ConstructSolver(self.settings["convector_settings"]["linear_solver_settings"])
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.linear_solver,
                self.settings["convector_settings"]["max_cfl"].GetDouble(),
                self.settings["convector_settings"]["cross_wind_stabilization_factor"].GetDouble(),
                self.settings["convector_settings"]["max_substeps"].GetInt())

        return level_set_convection_process

    def _set_new_waterlevel_surface(self,dt,elementId):
        prev_dt = self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        print("lastdt",prev_dt)
        print("newdt",dt)
        for node in self._fluid_model_part.Nodes:
            DISTANCE_NODE_n_1=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            DISTANCE_NODE_n=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE,1)
            Delta_Distance= DISTANCE_NODE_n_1-DISTANCE_NODE_n

            corrector_factor=(dt/prev_dt)*Delta_Distance
            new_distance=DISTANCE_NODE_n_1 + corrector_factor
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,new_distance)
            if node.Id in elementId:
                print("Identidad",node.Id)
                print(" ")
                print("Distancia paso actual",DISTANCE_NODE_n_1)
                print("Distancia paso anterior",DISTANCE_NODE_n)
                print("Incremento level_set",Delta_Distance)
                print("dt_simualtion",prev_dt)
                print("dt*",dt)
                print("corrector_factor",corrector_factor)
                print("La nueva distancia es phi n+1 corregida",new_distance)
