import KratosMultiphysics
import KratosMultiphysics.DEMApplication as Dem
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
from KratosMultiphysics.DEMApplication.multiaxial_control_module_generalized_2d_utility import MultiaxialControlModuleGeneralized2DUtility

class MultiaxialControlModuleFEMDEMGeneralized2DUtility(MultiaxialControlModuleGeneralized2DUtility):
    def __init__(self, Model, settings):

        #NOTE: Negative target_stress means compression

        self.spheres_model_part = Model["SpheresPart"]
        self.fem_model_part = Model["Structure"]

        self.parameters = settings["multiaxial_control_module_fem_dem_generalized_2d_utility"]

        self.cm_step = 0
        self.output_interval = self.parameters["Parameters"]["output_interval"].GetInt()

        self.cm_utility = DemFem.MultiaxialControlModuleFEMDEMGeneralized2DUtilities(self.spheres_model_part,
                                                                                self.fem_model_part,
                                                                                self.parameters)

    def ExecuteFinalizeSolutionStep(self):
        self.cm_utility.ExecuteFinalizeSolutionStep()

        if self.output_interval == 0:
            return
        self.FillPrintVariables()

        self.PrintResults()

    def FillPrintVariables(self):
        time = self.fem_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.times.append(time)

        # Pa to psi
        unit_trans = 0.000145038

        for name in self.actuator_names:
            if name == 'Z':
                for node in self.fem_submodelparts[name][0].Nodes:
                    self.reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Z))*unit_trans)
                    self.target_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Z))*unit_trans)
                    self.velocities[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Z))*unit_trans)
                    break
            elif name == 'X':
                for node in self.fem_submodelparts[name][0].Nodes:
                    self.reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_X))*unit_trans)
                    self.target_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_X))*unit_trans)
                    self.velocities[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_X))*unit_trans)
                    break
            elif name == 'Y':
                for node in self.fem_submodelparts[name][0].Nodes:
                    self.reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Y))*unit_trans)
                    self.target_reactions[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Y))*unit_trans)
                    self.velocities[name].append(abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Y))*unit_trans)
                    break
            elif name == 'Radial':
                for node in self.fem_submodelparts[name][0].Nodes:
                    value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_X))*unit_trans
                    value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Y))*unit_trans
                    norm = (value_x**2 + value_y**2)**0.5
                    self.reactions[name].append(norm)
                    value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_X))*unit_trans
                    value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Y))*unit_trans
                    norm = (value_x**2 + value_y**2)**0.5
                    self.target_reactions[name].append(norm)
                    value_x = abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_X))*unit_trans
                    value_y = abs(node.GetValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Y))*unit_trans
                    norm = (value_x**2 + value_y**2)**0.5
                    self.velocities[name].append(norm)
                    break

    def PrintResults(self):
        self.cm_step += 1

        if self.output_interval == 0:
            return
        if self.cm_step % self.output_interval != 0:
            return

        import matplotlib
        matplotlib.use('agg')
        import matplotlib.pyplot as plt

        f = plt.figure()

        for name in self.actuator_names:
            plt.plot(self.times,
                    self.reactions[name],
                    color=self.colors[name],
                    linestyle=self.reactions['linestyle'],
                    linewidth=self.reactions['linewidth'],
                    label='reaction_' + name)
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
        plt.ylabel('velocity (m/s)')
        # giving a title to my graph
        plt.title('Loading velocity')

        f.savefig("velocities.pdf", bbox_inches='tight')
        plt.close()