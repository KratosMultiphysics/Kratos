import KratosMultiphysics
from KratosMultiphysics.DEMApplication import *
import DEM_analysis_stage
import os
import shutil
import sys


class TimeStepTester(object):
    def __init__(self):
        #self.schemes_list = ["Forward_Euler", "Taylor_Scheme", "Symplectic_Euler", "Velocity_Verlet"]
        self.schemes_list = ["Beeman_Scheme","Symplectic_Euler", "Velocity_Verlet"]
        #self.schemes_list = ["Symplectic_Euler", "Velocity_Verlet","Cimne_Scheme","Gear_Scheme", "Beeman_Scheme"]
        #self.schemes_list = ["Symplectic_Euler", "Velocity_Verlet","Cimne_Scheme","Gear_Scheme"]
        self.stable_time_steps_list = []

        gnuplot_data = open("gnuplot_file.dem", 'w+')
        gnuplot_data.write('unset xrange' +'\n')
        gnuplot_data.write('unset yrange' +'\n')
        gnuplot_data.write('clear' +'\n')
        gnuplot_data.write('reset' +'\n')
        gnuplot_data.write('set style fill solid 1.0 noborder' +'\n')
        gnuplot_data.write("set terminal pngcairo size 700,250 enhanced font 'Verdana,8'" +'\n')
        gnuplot_data.write('set border 3 back ls 11' +'\n')
        gnuplot_data.write('set tics nomirror' +'\n')
        gnuplot_data.write('set grid back ls 12' +'\n')
        gnuplot_data.write('set boxwidth 0.8 absolute' +'\n')
        gnuplot_data.write('set style fill   solid 1.00 border' +'\n')
        gnuplot_data.write('set key inside right bottom vertical Right noreverse noenhanced autotitles columnhead nobox' +'\n')
        gnuplot_data.write('set key noinvert samplen 1 spacing 1 width 0 height 0 ' +'\n')
        gnuplot_data.write('set xtics border out scale 0.75 nomirror norotate  offset character 0, 0, 0 autojustify' +'\n')
        gnuplot_data.write('set xtics  norangelimit font ",8"' +'\n')
        gnuplot_data.write('set xtics   ()' +'\n')
        gnuplot_data.write('set ytics border out scale 0.75 nomirror norotate  offset character 0, 0, 0 #autojustify' +'\n')
        gnuplot_data.write('set ytics autofreq  norangelimit font ",8"' +'\n')
        gnuplot_data.write('\n')

        if os.path.exists('scheme_data'):
            shutil.rmtree('scheme_data')

        if not os.path.exists('scheme_data'):
            os.makedirs('scheme_data')


    def Run(self):

        for scheme in self.schemes_list:
            self.RunForACertainScheme(scheme)

        self.Finalize()

    def RunForACertainScheme(self, scheme):
        print("Computing stable time step for scheme: "+ scheme)
        tolerance = 0.5e-7
        dt = 1e-3    # changing the initial dt will change the number of iterations required to find a stable time step
        previous_dt = 0.0

        gnuplot_data = open("gnuplot_file.dem", 'a')
        gnuplot_data.write('\n'+'\n'+'\n'+'\n'+'\n'+'\n'+'\n')
        gnuplot_data.write('unset xrange' +'\n')
        gnuplot_data.write('unset yrange' +'\n')
        gnuplot_data.write('unset title' +'\n')
        gnuplot_data.write("set output '" +str(scheme) + "-time_vs_disp.png'" +'\n')
        gnuplot_data.write("set xrange [0:1]" +'\n')
        gnuplot_data.write("set title '" +str(scheme) + " Time|Displacement'" +'\n')
        gnuplot_data.write("set xlabel '{/Helvetica-Italic Time (s)}'" +'\n')
        gnuplot_data.write("set ylabel '{/Helvetica-Italic Displacement}'" +'\n')
        gnuplot_data.write('plot 	"a" , \\' +'\n')
        gnuplot_data.close()

        while dt > previous_dt + tolerance:
            try:
                print("current dt: " + str(dt))
                self.RunTestCaseWithCustomizedDtAndScheme(dt, scheme)

            except SystemExit:
                factor = min(0.2, 0.2*(dt-previous_dt))
                dt = factor * dt
                print("decreasing dt by " + str(factor))
                continue

            previous_dt = dt
            dt = dt * 1.2
            print("increasing dt by 1.2")

        self.stable_time_steps_list.append(previous_dt)

    @classmethod
    def RunTestCaseWithCustomizedDtAndScheme(self, dt, scheme):
        model = KratosMultiphysics.Model()
        CustomizedSolutionForTimeStepTesting(model, dt, scheme).Run()

    def Finalize(self):

        print("\n")
        print("#############################")
        print("List of tested schemes:")
        print(self.schemes_list)
        print("List of stable time steps:")
        print(self.stable_time_steps_list)
        print("#############################")
        print("\n")


class CustomizedSolutionForTimeStepTesting(DEM_analysis_stage.DEMAnalysisStage):

    def __init__(self, model, dt, scheme):
        self.customized_time_step = dt
        self.customized_scheme = scheme
        self.LoadParametersFile()
        self.previous_energy = 0.0

        import math
        def truncate(number, digits) -> float:
            stepper = 10.0 ** digits
            return math.trunc(stepper * number) / stepper
        dt_short = truncate(self.customized_time_step, 6)



        absolute_path_to_file = os.path.join('scheme_data', scheme + str(dt_short) + "_graph.grf")
        self.graph_export   = open(absolute_path_to_file, 'w')
        self.gnuplot_data = open("gnuplot_file.dem", 'a')
        # self.disp_ = []
        # self.energy_ = []
        # self.vel_ = []

        # self.disp_.append("        " +" '" + str(scheme) + str(dt_short) + "_graph.grf'" + "  " + "using ($1):($2) every 1 with lines lc 1 lw 2 lt 1  title '" + " " + str(scheme) + str(dt_short) + "_graph'" + ' , \\' +'\n')
        # self.energy_.append("        " +" '" + str(scheme) + str(dt_short) + "_graph.grf'" + "  " + "using ($1):($3) every 1 with lines lc 1 lw 2 lt 1  title '" + " " + str(scheme) + str(dt_short) + "_graph'" + ' , \\' +'\n')
        # self.vel_.append("        " +" '" + str(scheme) + str(dt_short) + "_graph.grf'" + "  " + "using ($1):($4) every 1 with lines lc 1 lw 2 lt 1  title '" + " " + str(scheme) + str(dt_short) + "_graph'" + ' , \\' +'\n')
        self.gnuplot_data.write("        " +" '" + str(scheme) + str(dt_short) + "_graph.grf'" + "  " + "using ($1):($2) every 1 with lines lc 1 lw 2 lt 1  title '" + " " + str(scheme) + str(dt_short) + "_graph'" + ' , \\' +'\n')
        super(CustomizedSolutionForTimeStepTesting, self).__init__(model, self.project_parameters)

    def LoadParametersFile(self):
        self.project_parameters = KratosMultiphysics.Parameters(
            """
            {
                "Dimension"                        : 3,
                "BoundingBoxOption"                : false,
                "BoundingBoxEnlargementFactor"     : 1.1,
                "AutomaticBoundingBoxOption"       : false,
                "BoundingBoxEnlargementFactor"     : 1.0,
                "BoundingBoxMaxX"                  : 1e3,
                "BoundingBoxMaxY"                  : 1e3,
                "BoundingBoxMaxZ"                  : 1e3,
                "BoundingBoxMinX"                  : -1e3,
                "BoundingBoxMinY"                  : -1e3,
                "BoundingBoxMinZ"                  : -1e3,
                "dem_inlet_option"                 : false,
                "GravityX"                         : 0.0,
                "GravityY"                         : 0.0,
                "GravityZ"                         : 0.0,
                "VelocityTrapOption"               : false,
                "RotationOption"                   : true,
                "CleanIndentationsOption"          : true,
                "RemoveBallsInEmbeddedOption"      : false,
                "solver_settings" :{
                    "strategy"                 : "sphere_strategy",
                    "RemoveBallsInitiallyTouchingWalls": false
                },

                "DeltaOption"                      : "Absolute",
                "SearchTolerance"                  : 0.0,
                "CoordinationNumber"               : 10,
                "AmplifiedSearchRadiusExtension"   : 1.10000e+00,
                "ModelDataInfo"                    : false,
                "VirtualMassCoefficient"           : 1.0,
                "RollingFrictionOption"            : false,
                "DontSearchUntilFailure"           : false,
                "ContactMeshOption"                : false,
                "OutputFileType"                   : "Binary",
                "Multifile"                        : "multiple_files",
                "TranslationalIntegrationScheme"   : "Forward_Euler",
                "RotationalIntegrationScheme"      : "Direct_Integration",
                "AutomaticTimestep"                : false,
                "DeltaTimeSafetyFactor"            : 1.0,
                "MaxTimeStep"                      : 1e-4,
                "FinalTime"                        : 1.0,
                "ControlTime"                      : 100,
                "NeighbourSearchFrequency"         : 1,
                "PeriodicDomainOption"             : false,
                "ElementType"                      : "SphericPartDEMElement3D",

                "GraphExportFreq"                  : 1e-5,
                "VelTrapGraphExportFreq"           : 1e-3,
                "OutputTimeStep"                   : 1e-4,
                "PostDisplacement"                 : true,
                "PostVelocity"                     : true,
                "PostElasticForces"                : false,
                "PostContactForces"                : false,
                "PostRigidElementForces"           : false,
                "PostTangentialElasticForces"      : false,
                "PostPressure"                     : false,
                "PostTotalForces"                  : true,
                "PostShearStress"                  : false,
                "PostNonDimensionalVolumeWear"     : false,
                "PostNodalArea"                    : false,
                "PostRHS"                          : false,
                "PostDampForces"                   : false,
                "PostAppliedForces"                : false,
                "PostRadius"                       : false,
                "PostGroupId"                      : false,
                "PostExportId"                     : false,
                "PostAngularVelocity"              : false,
                "PostParticleMoment"               : false,
                "PostEulerAngles"                  : false,
                "PostContactSigma"                 : false,
                "PostContactTau"                   : false,
                "PostLocalContactForce"            : false,
                "PostFailureCriterionState"        : false,
                "PostContactFailureId"             : false,
                "PostMeanContactArea"              : false,
                "PostStressStrainOption"           : false,
                "PostRollingResistanceMoment"      : false,
                "post_vtk_option"                  : false,
                "problem_name"                     : "TimeStepTests"
                }

            """
            )



    def SetDt(self):
        self.DEM_parameters["TranslationalIntegrationScheme"].SetString(self.customized_scheme)
        self.DEM_parameters["MaxTimeStep"].SetDouble(self.customized_time_step)
        default_input_parameters = self.GetDefaultInputParameters()
        self.DEM_parameters.ValidateAndAssignDefaults(default_input_parameters)
        self._GetSolver().dt = self.DEM_parameters["MaxTimeStep"].GetDouble()


    def ReadModelParts(self, max_node_Id=0, max_elem_Id=0, max_cond_Id=0):   # exists in DEM_analysis_stage
        properties = KratosMultiphysics.Properties(0)
        properties_walls = KratosMultiphysics.Properties(0)
        self.SetHardcodedProperties(properties, properties_walls)
        self.spheres_model_part.AddProperties(properties)
        self.rigid_face_model_part.AddProperties(properties_walls)

        DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
        DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, False)

        translational_scheme = ForwardEulerScheme()
        translational_scheme.SetTranslationalIntegrationSchemeInProperties(properties, True)
        rotational_scheme = ForwardEulerScheme()
        rotational_scheme.SetRotationalIntegrationSchemeInProperties(properties, True)

        element_name = "SphericParticle3D"
        PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        coordinates = KratosMultiphysics.Array3()
        coordinates[0] = 0.0
        coordinates[1] = -0.1
        coordinates[2] = 0.0
        radius = 0.1
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        coordinates[0] = 0.0
        coordinates[1] = 0.2
        coordinates[2] = 0.0
        radius = 0.1
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)


        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 10.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 5.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)


        self.rigid_face_model_part.CreateNewNode(11, -0.4, -0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(12, -0.4, -0.5, 0.5)

        self.rigid_face_model_part.CreateNewNode(13, 0.4, -0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(14, 0.4, -0.5, 0.5)


        self.rigid_face_model_part.CreateNewNode(15, 0.4, 0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(16, 0.4, 0.5, 0.5)

        self.rigid_face_model_part.CreateNewNode(17, -0.4, 0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(18, -0.4, 0.5, 0.5)

        condition_name = "RigidFace3D3N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 1, [11, 12, 13], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 2, [12, 13, 14], self.rigid_face_model_part.GetProperties()[0])

        self.rigid_face_model_part.CreateNewCondition(condition_name, 3, [13, 14, 15], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 4, [14, 15, 16], self.rigid_face_model_part.GetProperties()[0])

        self.rigid_face_model_part.CreateNewCondition(condition_name, 5, [15, 16, 17], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 6, [16, 17, 18], self.rigid_face_model_part.GetProperties()[0])

        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [17, 18, 11], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 8, [18, 11, 12], self.rigid_face_model_part.GetProperties()[0])

        self.initial_test_energy = self.ComputeEnergy()
        print("initial_energy: ", self.initial_test_energy)


    def ComputeEnergy(self):
        this_test_total_energy = 0.0

        for element in self.spheres_model_part.Elements:
            this_test_total_energy += element.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            this_test_total_energy += element.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            this_test_total_energy += element.Calculate(PARTICLE_ELASTIC_ENERGY, self.spheres_model_part.ProcessInfo)

        return this_test_total_energy

    def ComputeEnergyVariation(self):
        this_test_total_energy = 0.0
        energy_delta = 0.0

        for element in self.spheres_model_part.Elements:
            this_test_total_energy += element.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            this_test_total_energy += element.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            this_test_total_energy += element.Calculate(PARTICLE_ELASTIC_ENERGY, self.spheres_model_part.ProcessInfo)
        energy_delta = (this_test_total_energy - self.initial_test_energy)/self.initial_test_energy

        return energy_delta


    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        properties[PARTICLE_DENSITY] = 2650.0
        properties[KratosMultiphysics.YOUNG_MODULUS] = 7.0e6
        properties[KratosMultiphysics.POISSON_RATIO] = 0.30
        properties[FRICTION] = 0.0
        properties[PARTICLE_COHESION] = 0.0
        properties[COEFFICIENT_OF_RESTITUTION] = 1.0
        properties[KratosMultiphysics.PARTICLE_MATERIAL] = 1
        properties[ROLLING_FRICTION] = 0.0
        properties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEMContinuumConstitutiveLaw"
        properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"

        properties_walls[FRICTION] = 0.0
        properties_walls[WALL_COHESION] = 0.0
        properties_walls[COMPUTE_WEAR] = 0
        properties_walls[SEVERITY_OF_WEAR] = 0.001
        properties_walls[IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[BRINELL_HARDNESS] = 200.0
        properties_walls[KratosMultiphysics.YOUNG_MODULUS] = 7.0e10
        properties_walls[KratosMultiphysics.POISSON_RATIO] = 0.30


    def FinalizeTimeStep(self, time):
        super(CustomizedSolutionForTimeStepTesting, self).FinalizeTimeStep(time)
        self.energy_delta = self.ComputeEnergyVariation()
        if abs(self.energy_delta) > 0.5:
            print("ENERGY VARIATION OVER 50%")
            print("time step is:" + str(self.customized_time_step))
            import sys
            sys.exit()

        # self.current_test_energy = self.ComputeEnergy()
        # if self.current_test_energy/self.initial_test_energy > 1.4:
        #     print("GAINING ENERGY!!")
        #     print("time step is:" + str(self.customized_time_step))
        #     sys.exit()

        stime = self.spheres_model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.spheres_model_part.Nodes:
            if node.Id == 1:
                self.vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                self.disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                self.graph_export.write(str("%.8g"%stime).rjust(13) +"  "+str("%.6g"%self.disp).rjust(12) +"  "+str("%.6g"%self.vel).rjust(12) +"  "+str("%.6g"%self.energy_delta).rjust(12)+'\n')


        #if not self.step%200:
        #    print("Energy: "+str(current_test_energy))

        #elif self.initial_test_energy/current_test_energy > 1.5:
        #    print("LOSING ENERGY!!")
        #    print("time step is:" + str(self.customized_time_step))
        #    sys.exit()

    def Finalize(self):
        # total = self.disp_ + self.energy_ + self.vel_
        # for i in total:
        #     self.gnuplot_data.write(i)

        super(CustomizedSolutionForTimeStepTesting, self).Finalize()


    def PrintResultsForGid(self, time):
        pass


if __name__ == '__main__':
    TimeStepTester().Run()
