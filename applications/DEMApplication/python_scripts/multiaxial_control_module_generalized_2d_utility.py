import KratosMultiphysics
import KratosMultiphysics.DEMApplication as Dem


class MultiaxialControlModuleGeneralized2DUtility():
    def __init__(self, Model, settings):

        #NOTE: Negative target_stress means compression

        self.spheres_model_part = Model["SpheresPart"]
        self.fem_model_part = Model["RigidFacePart"] #rigid_walls_model_part

        self.parameters = settings["multiaxial_control_module_generalized_2d_utility"]

        self.cm_step = 0
        self.output_interval = self.parameters["Parameters"]["output_interval"].GetInt()

        self.cm_utility = Dem.MultiaxialControlModuleGeneralized2DUtilities(self.spheres_model_part,
                                                                            self.fem_model_part,
                                                                            self.parameters)

    def ExecuteInitialize(self):
        self.cm_utility.ExecuteInitialize()

        if self.output_interval == 0:
            return
        self.InitializePrintVariables()

    def ExecuteInitializeSolutionStep(self):
        self.cm_utility.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.cm_utility.ExecuteFinalizeSolutionStep()

        if self.output_interval == 0:
            return
        self.FillPrintVariables()

    def PrintResults(self):
        self.cm_step += 1

        if self.output_interval == 0:
            return
        if self.cm_step % self.output_interval != 0:
            return

        import matplotlib.pyplot as plt

        f = plt.figure()

        for name in self.actuator_names:
            plt.plot(self.times,
                    self.reactions[name],
                    color=self.colors[name],
                    linestyle=self.reactions['linestyle'],
                    linewidth=self.reactions['linewidth'],
                    label='reaction_' + name)
            # plt.plot(self.times,
            #         self.smoothed_reactions[name],
            #         color=self.colors[name],
            #         linestyle=self.smoothed_reactions['linestyle'],
            #         linewidth=self.smoothed_reactions['linewidth'],
            #         label='smoothed_reaction_' + name)
            # plt.plot(self.times,
            #         self.elastic_reactions[name],
            #         color=self.colors[name],
            #         linestyle=self.elastic_reactions['linestyle'],
            #         linewidth=self.elastic_reactions['linewidth'],
            #         label='elastic_reaction_' + name)
            # plt.plot(self.times,
            #         self.smoothed_elastic_reactions[name],
            #         color=self.colors[name],
            #         linestyle=self.smoothed_elastic_reactions['linestyle'],
            #         linewidth=self.smoothed_elastic_reactions['linewidth'],
            #         label='smoothed_elastic_reaction_' + name)
            plt.plot(self.times,
                    self.target_reactions[name],
                    color=self.colors[name],
                    linestyle=self.target_reactions['linestyle'],
                    linewidth=self.target_reactions['linewidth'],
                    label='target_reaction_' + name)

        # plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.legend()
        max_stress = 0
        for name in self.actuator_names:
            max_stress = max(max(abs(min(self.target_reactions[name])),max(self.target_reactions[name])), max_stress)
        plt.ylim(-0.1*max_stress,1.1*max_stress)

        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

        # naming the x axis
        plt.xlabel('Time (s)')
        # naming the y axis
        plt.ylabel('Stress (Psi)')
        # giving a title to my graph
        plt.title('Reaction vs target stresses')

        f.savefig("reactions.pdf", bbox_inches='tight')
        plt.close()

        f = plt.figure()

        for name in self.actuator_names:
            plt.plot(self.times,
                    self.velocities[name],
                    color=self.colors[name],
                    label='velocity_' + name)

        plt.legend()

        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

        # naming the x axis
        plt.xlabel('Time (s)')
        # naming the y axis
        plt.ylabel('Velocity (m/s)')
        # giving a title to my graph
        plt.title('Loading velocity')

        f.savefig("velocities.pdf", bbox_inches='tight')
        plt.close()

    def InitializePrintVariables(self):

        import matplotlib.pyplot as plt
        from matplotlib import colors as mcolors

        self.times = []
        self.actuator_names = []
        # Dictionaries containing one array per actuator
        self.colors = dict()
        self.reactions = dict()
        self.smoothed_reactions = dict()
        self.elastic_reactions = dict()
        self.smoothed_elastic_reactions = dict()
        self.target_reactions = dict()
        self.velocities = dict()
        self.dem_submodelparts = dict()
        self.fem_submodelparts = dict()

        self.colors_database = (c for c in mcolors.cnames.keys())
        self.colors_database = (c for c in ['k','r','b','g','c','m','y'])
        self.reactions['linestyle'] = 'solid'
        self.reactions['linewidth'] = 1
        self.smoothed_reactions['linestyle'] = 'solid'
        self.smoothed_reactions['linewidth'] = 3
        self.elastic_reactions['linestyle'] = 'dotted'
        self.elastic_reactions['linewidth'] = 1
        self.smoothed_elastic_reactions['linestyle'] = 'dotted'
        self.smoothed_elastic_reactions['linewidth'] = 3
        self.target_reactions['linestyle'] = 'dashed'
        self.target_reactions['linewidth'] = 3

        for actuator_parameters in self.parameters["list_of_actuators"]:
            name = actuator_parameters["Parameters"]["actuator_name"].GetString()
            self.actuator_names.append(name)
            self.reactions[name] = []
            self.colors[name] = next(self.colors_database)
            self.smoothed_reactions[name] = []
            self.elastic_reactions[name] = []
            self.smoothed_elastic_reactions[name] = []
            self.target_reactions[name] = []
            self.velocities[name] = []
            self.dem_submodelparts[name] = []
            self.fem_submodelparts[name] = []
            for boundary in actuator_parameters["list_of_dem_boundaries"]:
                mp = self.spheres_model_part.GetSubModelPart(boundary["model_part_name"].GetString())
                self.dem_submodelparts[name].append(mp)
            for boundary in actuator_parameters["list_of_fem_boundaries"]:
                mp = self.fem_model_part.GetSubModelPart(boundary["model_part_name"].GetString())
                self.fem_submodelparts[name].append(mp)

    def FillPrintVariables(self):
        time = self.spheres_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.times.append(time)

        # Pa to psi
        unit_trans = 0.000145038

        for name in self.actuator_names:
            if name == 'Z':
                for node in self.dem_submodelparts[name][0].Nodes:
                    self.reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Z))*unit_trans)
                    # self.smoothed_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_Z))*unit_trans)
                    # self.elastic_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_Z))*unit_trans)
                    # self.smoothed_elastic_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_Z))*unit_trans)
                    self.target_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Z))*unit_trans)
                    self.velocities[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Z)))
                    break

            elif name == 'X':
                for node in self.fem_submodelparts[name][0].Nodes:
                    self.reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_X))*unit_trans)
                    # self.smoothed_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_X))*unit_trans)
                    # self.elastic_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_X))*unit_trans)
                    # self.smoothed_elastic_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_X))*unit_trans)
                    self.target_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_X))*unit_trans)
                    self.velocities[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_X)))
                    break
            elif name == 'Y':
                for node in self.fem_submodelparts[name][0].Nodes:
                    self.reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Y))*unit_trans)
                    # self.smoothed_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_Y))*unit_trans)
                    # self.elastic_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_Y))*unit_trans)
                    # self.smoothed_elastic_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_Y))*unit_trans)
                    self.target_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Y))*unit_trans)
                    self.velocities[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Y)))
                    break
            elif name == 'Radial':
                for node in self.fem_submodelparts[name][0].Nodes:
                    value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_X))*unit_trans
                    value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Y))*unit_trans
                    norm = (value_x**2 + value_y**2)**0.5
                    self.reactions[name].append(norm)
                    # value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_X))*unit_trans
                    # value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_Y))*unit_trans
                    # norm = (value_x**2 + value_y**2)**0.5
                    # self.smoothed_reactions[name].append(norm)
                    # value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_X))*unit_trans
                    # value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_Y))*unit_trans
                    # norm = (value_x**2 + value_y**2)**0.5
                    # self.elastic_reactions[name].append(norm)
                    # value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_X))*unit_trans
                    # value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_Y))*unit_trans
                    # norm = (value_x**2 + value_y**2)**0.5
                    # self.smoothed_elastic_reactions[name].append(norm)
                    value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_X))*unit_trans
                    value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Y))*unit_trans
                    norm = (value_x**2 + value_y**2)**0.5
                    self.target_reactions[name].append(norm)
                    value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_X))
                    value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Y))
                    norm = (value_x**2 + value_y**2)**0.5
                    self.velocities[name].append(norm)
                    break