from xmc.distributedEnvironmentFramework import *

# Import packages
import numpy as np
import math


####################################################################################################
############################################# CLASSES ##############################################
####################################################################################################

OUTPUT_QUANTITIES = 1

class UnfolderManager(object):
    """
    Class used to organize a list of values into a list of sublists. Referring to f as future type, this class allows to pass from [f, f, f, ... , f] to [[f, f, f, ... , f], ... ,[f, ... , f]]. The length of the original list is "number", while the length of each sublist is group.
    A method managing multiple contributions is present as well.

    Attributes:
    - number: number of values of original list
    - group: desired length of sublists of output list.
    - groups: number of sublists of output list.

    Methods:
    - UnfoldNValues_Task: task method calling UnfoldNValues.
    - UnfoldNValues: method creating the list of sublists.
    - PostprocessContributionsPerInstance: task method summing together multiple contributions, if any. After summing all contributions, it calls UnfoldNValues to create the list of sublists.
    """

    def __init__(self, number, group):
        global OUTPUT_QUANTITIES
        self.groups = math.ceil(number/group)
        OUTPUT_QUANTITIES = self.groups
        self.number=number
        self.group=group

    def UnfoldNValues(self, number, group, values):
        """
        Method creating list of sublists.

        Inputs:
        - self: an instance of the class.
        - number: number of values of original list
        - group: desired length of sublists of output list.
        - values: original list of values.
        """

        partial_vals = []
        for val in range(1, number+1):
            partial_vals.append(values[val-1])
            if (val == number or val % group == 0):
                yield partial_vals
                partial_vals = []

    @ExaquteTask(target_direction=IN,returns='OUTPUT_QUANTITIES')
    def UnfoldNValues_Task(self, values):
        """
        Task method calling UnfoldNValues.

        Inputs:
        - self: an instance of the class.
        - values: original list of values.
        """

        list_unfolded = list(self.UnfoldNValues(self.number, self.group, values))
        if (self.groups == 1):
            list_unfolded = list_unfolded[0]
        return list_unfolded

    @ExaquteTask(aux_qoi_array_contributions={Type: COLLECTION_IN, Depth: 2},returns='OUTPUT_QUANTITIES')
    def PostprocessContributionsPerInstance(self,aux_qoi_array_contributions,number_qoi,number_combined_qoi):
        """
        Task method summing multiple contribution of a specific realization and calling UnfoldNValues.

        Inputs:
        - self: an instance of the class.
        - aux_qoi_array_contributions: original list of values with multiple contributions.
        - number_qoi: number of quantities of interest.
        - number_combined_qoi: number of combined quantities of interest.
        """

        aux_qoi_array = [[] for _ in range (0,number_qoi+number_combined_qoi)] # to store each qoi
        # append components to aux array
        for qoi_list in aux_qoi_array_contributions:
            for qoi_counter in range (0,number_qoi+number_combined_qoi):
                aux_qoi_array[qoi_counter].append(qoi_list[qoi_counter])
        assert(len(aux_qoi_array)==number_qoi+number_combined_qoi)
        qoi_list = []
        # expected value for steady state / time averaged qoi
        for qoi_counter in range (0,number_qoi):
            qoi_value = np.mean(aux_qoi_array[qoi_counter])
            qoi_list.append(qoi_value)
        # time power sums sum for time series qoi
        for combined_ps_counter in range (number_qoi,number_qoi+number_combined_qoi):
            S1 = 0 ; S2 = 0 ; S3 = 0 ; S4 = 0 ; S5 = 0 ; S6 = 0 ; S7 = 0 ; S8 = 0 ; S9 = 0 ; S10 = 0 ; contributions = 0
            for i in range (0,len(aux_qoi_array[combined_ps_counter])):
                S1 = S1 + aux_qoi_array[combined_ps_counter][i][0][0]
                S2 = S2 + aux_qoi_array[combined_ps_counter][i][1][0]
                if (len(aux_qoi_array[combined_ps_counter][i]) == 5):
                    S3 = S3 + aux_qoi_array[combined_ps_counter][i][2][0]
                    S4 = S4 + aux_qoi_array[combined_ps_counter][i][3][0]
                elif (len(aux_qoi_array[combined_ps_counter][i]) == 11):
                    S3 = S3 + aux_qoi_array[combined_ps_counter][i][2][0]
                    S4 = S4 + aux_qoi_array[combined_ps_counter][i][3][0]
                    S5 = S5 + aux_qoi_array[combined_ps_counter][i][4][0]
                    S6 = S6 + aux_qoi_array[combined_ps_counter][i][5][0]
                    S7 = S7 + aux_qoi_array[combined_ps_counter][i][6][0]
                    S8 = S8 + aux_qoi_array[combined_ps_counter][i][7][0]
                    S9 = S9 + aux_qoi_array[combined_ps_counter][i][8][0]
                    S10 = S10 + aux_qoi_array[combined_ps_counter][i][9][0]
                contributions = contributions + aux_qoi_array[combined_ps_counter][i][-1]
            if (len(aux_qoi_array[combined_ps_counter][-1]) == 11):
                qoi_list.append([[S1],[S2],[S3],[S4],[S5],[S6],[S7],[S8],[S9],[S10],contributions])
            elif (len(aux_qoi_array[combined_ps_counter][-1]) == 5):
                qoi_list.append([[S1],[S2],[S3],[S4],contributions])
            else:
                qoi_list.append([[S1],[S2],contributions])
        list_unfolded = list(self.UnfoldNValues(self.number,self.group,qoi_list))
        if (self.groups == 1):
            list_unfolded = list_unfolded[0]
        return list_unfolded
