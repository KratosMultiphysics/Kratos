# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

class SaveVariablesUtility:

    def SaveMechanicalVariables(problem_name, ProjectParameters, main_model_part, time):

        OriginalMdpa = open(problem_name + ".mdpa" , 'r')
        OutputMdpa = open(problem_name + "_mechanical_" + str(time) + ".mdpa", 'w')

        for Line in OriginalMdpa:
            if Line.startswith("Begin SubModelPart"):
                OriginalMdpa.close()
                break
            OutputMdpa.write(Line)

        mechanical_loads_sub_model_part_list = ProjectParameters["solver_settings"]["mechanical_solver_settings"]["problem_domain_sub_model_part_list"]

        mechanical_parts = []

        if (ProjectParameters["problem_data"]["consider_construction"].GetBool()):
            for i in range(mechanical_loads_sub_model_part_list.size()):
                with open(ProjectParameters["construction_process"]["construction_input_file_name"].GetString(),'r') as construction_file:
                    for Line in construction_file:
                        Line = Line.strip('\n')# Remove the line-ending characters
                        if (str(mechanical_loads_sub_model_part_list[i].GetString()) == str("Parts_" + Line.split()[1])):
                            if int(Line.split()[0]) <= time:
                                mechanical_parts.append(main_model_part.GetSubModelPart(mechanical_loads_sub_model_part_list[i].GetString()))
                            break
        else:
            for i in range(mechanical_loads_sub_model_part_list.size()):
                mechanical_parts.append(main_model_part.GetSubModelPart(mechanical_loads_sub_model_part_list[i].GetString()))

        OutputMdpa.write('Begin NodalData DISPLACEMENT')
        for part in mechanical_parts:
            print(part)
            for node in part.Nodes:
                OutputMdpa.write('\n' + str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)))
        OutputMdpa.write('\nEnd NodalData\n')

        OutputMdpa.write('\nBegin NodalData NODAL_CAUCHY_STRESS_TENSOR')
        for part in mechanical_parts:
            for node in part.Nodes:
                OutputMdpa.write('\n' + str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosDam.NODAL_CAUCHY_STRESS_TENSOR)))
        OutputMdpa.write('\nEnd NodalData\n')

        OutputMdpa.close()

    def SaveFinalMechanicalVariables(problem_name, ProjectParameters, main_model_part):

        OriginalMdpa = open(problem_name + ".mdpa" , 'r')
        OutputMdpa = open(problem_name + "_mechanical.mdpa", 'w')

        for Line in OriginalMdpa:
            if Line.startswith("Begin SubModelPart"):
                OriginalMdpa.close()
                break
            OutputMdpa.write(Line)

        mechanical_loads_sub_model_part_list = ProjectParameters["solver_settings"]["mechanical_solver_settings"]["problem_domain_sub_model_part_list"]

        mechanical_parts = []
        for i in range(mechanical_loads_sub_model_part_list.size()):
            mechanical_parts.append(main_model_part.GetSubModelPart(mechanical_loads_sub_model_part_list[i].GetString()))

        OutputMdpa.write('Begin NodalData NODAL_CAUCHY_STRESS_TENSOR')
        for part in mechanical_parts:
            for node in part.Nodes:
                OutputMdpa.write('\n' + str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosDam.NODAL_CAUCHY_STRESS_TENSOR)))
        OutputMdpa.write('\nEnd NodalData\n')

        OutputMdpa.close()

    def SaveThermalVariables(problem_name, ProjectParameters, main_model_part, time):

        OriginalMdpa = open(problem_name + ".mdpa" , 'r')
        OutputMdpa = open(problem_name + "_thermal_" + str(time) + ".mdpa", 'w')

        for Line in OriginalMdpa:
            if Line.startswith("Begin SubModelPart"):
                OriginalMdpa.close()
                break
            OutputMdpa.write(Line)

        thermal_loads_sub_model_part_list = ProjectParameters["solver_settings"]["thermal_solver_settings"]["problem_domain_sub_model_part_list"]

        thermal_parts = []

        if (ProjectParameters["problem_data"]["consider_construction"].GetBool()):
            for i in range(thermal_loads_sub_model_part_list.size()):
                with open(ProjectParameters["construction_process"]["construction_input_file_name"].GetString(),'r') as construction_file:
                    for Line in construction_file:
                        Line = Line.strip('\n')# Remove the line-ending characters
                        if (str(thermal_loads_sub_model_part_list[i].GetString()) == str("Thermal_" + Line.split()[1])):
                            if int(Line.split()[0]) <= time:
                                thermal_parts.append(main_model_part.GetSubModelPart(thermal_loads_sub_model_part_list[i].GetString()))
                            break
        else:
            for i in range(thermal_loads_sub_model_part_list.size()):
                thermal_parts.append(main_model_part.GetSubModelPart(thermal_loads_sub_model_part_list[i].GetString()))

        OutputMdpa.write('Begin NodalData TEMPERATURE')
        for part in thermal_parts:
            for node in part.Nodes:
                OutputMdpa.write('\n' + str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)))
        OutputMdpa.write('\nEnd NodalData\n')

        OutputMdpa.write('\nBegin NodalData NODAL_REFERENCE_TEMPERATURE')
        for part in thermal_parts:
            for node in part.Nodes:
                OutputMdpa.write('\n' + str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosDam.NODAL_REFERENCE_TEMPERATURE)))
        OutputMdpa.write('\nEnd NodalData\n')

        OutputMdpa.close()

    def SaveFinalThermalVariables(problem_name, ProjectParameters, main_model_part):

        OriginalMdpa = open(problem_name + ".mdpa" , 'r')
        OutputMdpa = open(problem_name + "_thermal.mdpa", 'w')

        for Line in OriginalMdpa:
            if Line.startswith("Begin SubModelPart"):
                OriginalMdpa.close()
                break
            OutputMdpa.write(Line)

        thermal_loads_sub_model_part_list = ProjectParameters["solver_settings"]["thermal_solver_settings"]["problem_domain_sub_model_part_list"]

        thermal_parts = []
        for i in range(thermal_loads_sub_model_part_list.size()):
            thermal_parts.append(main_model_part.GetSubModelPart(thermal_loads_sub_model_part_list[i].GetString()))

        OutputMdpa.write('Begin NodalData TEMPERATURE')
        for part in thermal_parts:
            for node in part.Nodes:
                OutputMdpa.write('\n' + str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)))
        OutputMdpa.write('\nEnd NodalData\n')

        OutputMdpa.write('\nBegin NodalData NODAL_REFERENCE_TEMPERATURE')
        for part in thermal_parts:
            for node in part.Nodes:
                OutputMdpa.write('\n' + str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosDam.NODAL_REFERENCE_TEMPERATURE)))
        OutputMdpa.write('\nEnd NodalData\n')

        OutputMdpa.close()
