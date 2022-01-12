# Import Python libraries
import numpy as np
import math

# Import distributed environment
from exaqute import *


####################################################################################################
############################################# CLASSES ##############################################
####################################################################################################

class UnfolderManager(object):
    """
    Class used to organize a list of values into a list of sublists. Referring to f as future type, this class allows to pass from [f, f, f, ... , f] to [[f, f, f, ... , f], ... ,[f, ... , f]]. The length of the original list is "number", while the length of each sublist is group.
    A method managing multiple contributions is present as well.

    Attributes:

    number: integer.
        Number of values of original list
    group: integer.
        Desired length of sublists of output list.
    groups: integer.
        Number of sublists of output list.

    Methods:

    UnfoldNValues_Task:
        Task method calling UnfoldNValues.
    UnfoldNValues:
        Method creating the list of sublists.
    PostprocessContributionsPerInstance:
        Task method summing together multiple contributions, if any. After summing all contributions, it calls UnfoldNValues to create the list of sublists.
    """

    def __init__(self, number, group):
        self.groups = math.ceil(number/group)
        self.number=number
        self.group=group

    def UnfoldNValues(self, number, group, values):
        """
        Method creating list of sublists.

        Inputs:

        number: integer.
            Number of values of original list
        group: integer.
            Desired length of sublists of output list.
        values: list.
            Original list of values.
        """

        partial_vals = []
        for val in range(1, number+1):
            partial_vals.append(values[val-1])
            if (val == number or val % group == 0):
                yield partial_vals
                partial_vals = []

    @task(keep=True, target_direction=IN,returns=1)
    def UnfoldNValues_Task(self, values):
        """
        Task method calling UnfoldNValues.

        Inputs:

        values: list.
            Original list of values.
        """

        list_unfolded = list(self.UnfoldNValues(self.number, self.group, values))
        if (self.groups == 1):
            list_unfolded = list_unfolded[0]
        return list_unfolded

    @task(keep=True, aux_qoi_array_contributions={Type: COLLECTION_IN, Depth: 2},returns=1)
    def PostprocessContributionsPerInstance_Task(self,aux_qoi_array_contributions,qoi_estimators):
        """
        Task method summing multiple contribution of a specific realization and calling UnfoldNValues.

        Inputs:

        aux_qoi_array_contributions: list.
            Original list of values with multiple contributions.
        qoi_estimators: list.
            List of strings. Each string is the corresponding moment estimator of each quantity of interest.
        """

        aux_qoi_array = [[] for _ in range (0,len(qoi_estimators))] # to store each qoi
        # loop over contributions and for each contribution over qoi
        # to append each contribution qoi in aux_qoi_array
        for qoi_list in aux_qoi_array_contributions:
            for qoi_counter in range (0,len(qoi_estimators)):
                aux_qoi_array[qoi_counter].append(qoi_list[qoi_counter])
        assert(len(aux_qoi_array)==len(qoi_estimators))
        qoi_list = []

        for estimator, qoi_counter in zip (qoi_estimators, range(len(qoi_estimators))):

            # expected value for moment estimators
            if estimator == "xmc.momentEstimator.MomentEstimator":
                qoi_value = np.mean(aux_qoi_array[qoi_counter])
                qoi_list.append(qoi_value)

            # combined power sums sum for combined moment estimators
            elif estimator == "xmc.momentEstimator.CombinedMomentEstimator":
                S1 = 0 ; S2 = 0 ; S3 = 0 ; S4 = 0 ; S5 = 0 ; S6 = 0 ; S7 = 0 ; S8 = 0 ; S9 = 0 ; S10 = 0 ; contributions = 0
                # loop over contributions of each CombinedMomentEstimator and sum power sums and contributions
                for i in range (0,len(aux_qoi_array[qoi_counter])):
                    S1 = S1 + aux_qoi_array[qoi_counter][i][0][0]
                    S2 = S2 + aux_qoi_array[qoi_counter][i][1][0]
                    if (len(aux_qoi_array[qoi_counter][i]) == 5):
                        S3 = S3 + aux_qoi_array[qoi_counter][i][2][0]
                        S4 = S4 + aux_qoi_array[qoi_counter][i][3][0]
                    elif (len(aux_qoi_array[qoi_counter][i]) == 11):
                        S3 = S3 + aux_qoi_array[qoi_counter][i][2][0]
                        S4 = S4 + aux_qoi_array[qoi_counter][i][3][0]
                        S5 = S5 + aux_qoi_array[qoi_counter][i][4][0]
                        S6 = S6 + aux_qoi_array[qoi_counter][i][5][0]
                        S7 = S7 + aux_qoi_array[qoi_counter][i][6][0]
                        S8 = S8 + aux_qoi_array[qoi_counter][i][7][0]
                        S9 = S9 + aux_qoi_array[qoi_counter][i][8][0]
                        S10 = S10 + aux_qoi_array[qoi_counter][i][9][0]
                    contributions = contributions + aux_qoi_array[qoi_counter][i][-1]
                # append results to qoi_list
                if (len(aux_qoi_array[qoi_counter][-1]) == 3): # order 1
                    qoi_list.append([[S1],[S2],contributions])
                elif (len(aux_qoi_array[qoi_counter][-1]) == 5): # order 2
                    qoi_list.append([[S1],[S2],[S3],[S4],contributions])
                elif (len(aux_qoi_array[qoi_counter][-1]) == 11): # order 5
                    qoi_list.append([[S1],[S2],[S3],[S4],[S5],[S6],[S7],[S8],[S9],[S10],contributions])
                else:
                    raise Exception("Ensemble average only supports CombinedMomentEstimator of order 1,2,5.")

            # expected value for multi moment estimators
            elif estimator == "xmc.momentEstimator.MultiMomentEstimator":
                # loop over scalar qoi of MultiMomentEstimator
                mme_members = []
                for member in range (0,len(aux_qoi_array[qoi_counter][0])):
                    # loop over contributions of each MultiMomentEstimator
                    mme_member = []
                    for i in range (0,len(aux_qoi_array[qoi_counter])):
                        mme_member.append(aux_qoi_array[qoi_counter][i][member])
                    mme_members.append(np.mean(mme_member))
                qoi_list.append(mme_members)

            # combined power sums sum for multi combined moment estimators
            elif estimator == "xmc.momentEstimator.MultiCombinedMomentEstimator":
                # loop over scalar qoi of CombinedMomentEstimator
                mcme_members = []
                for member in range (0,len(aux_qoi_array[qoi_counter][0])):
                    # initialize power sums and contributions of each member
                    S1 = 0 ; S2 = 0 ; S3 = 0 ; S4 = 0 ; S5 = 0 ; S6 = 0 ; S7 = 0 ; S8 = 0 ; S9 = 0 ; S10 = 0 ; contributions = 0
                    # loop over contributions of each MultiCombinedMomentEstimator and sum power sums and contributions
                    for i in range (0,len(aux_qoi_array[qoi_counter])):
                        S1 = S1 + aux_qoi_array[qoi_counter][i][member][0][0]
                        S2 = S2 + aux_qoi_array[qoi_counter][i][member][1][0]
                        if (len(aux_qoi_array[qoi_counter][i]) == 5):
                            S3 = S3 + aux_qoi_array[qoi_counter][i][member][2][0]
                            S4 = S4 + aux_qoi_array[qoi_counter][i][member][3][0]
                        elif (len(aux_qoi_array[qoi_counter][i]) == 11):
                            S3 = S3 + aux_qoi_array[qoi_counter][i][member][2][0]
                            S4 = S4 + aux_qoi_array[qoi_counter][i][member][3][0]
                            S5 = S5 + aux_qoi_array[qoi_counter][i][member][4][0]
                            S6 = S6 + aux_qoi_array[qoi_counter][i][member][5][0]
                            S7 = S7 + aux_qoi_array[qoi_counter][i][member][6][0]
                            S8 = S8 + aux_qoi_array[qoi_counter][i][member][7][0]
                            S9 = S9 + aux_qoi_array[qoi_counter][i][member][8][0]
                            S10 = S10 + aux_qoi_array[qoi_counter][i][member][9][0]
                        contributions = contributions + aux_qoi_array[qoi_counter][i][member][-1]
                    # append results to mcme_members
                    if (len(aux_qoi_array[qoi_counter][i][member]) == 3): # order 1
                        mcme_members.append([[S1],[S2],contributions])
                    elif (len(aux_qoi_array[qoi_counter][i][member]) == 5): # order 2
                        mcme_members.append([[S1],[S2],[S3],[S4],contributions])
                    elif (len(aux_qoi_array[qoi_counter][i][member]) == 11): # order 5
                        mcme_members.append([[S1],[S2],[S3],[S4],[S5],[S6],[S7],[S8],[S9],[S10],contributions])
                    else:
                        raise Exception("Ensemble average only supports MultiCombinedMomentEstimator of order 1,2,5.")
                qoi_list.append(mcme_members)

            else:
                err_msg =  "The moment estimator {} passed to the KratosSolverWrapper is not supported.\n".format(estimator)
                err_msg += "Available options are: \"MomentEstimator\", \"CombinedMomentEstimator\", \"MultiMomentEstimator\" and  \"MultiCombinedMomentEstimator.\""
                raise Exception(err_msg)

        list_unfolded = list(self.UnfoldNValues(self.number,self.group,qoi_list))
        if (self.groups == 1):
            list_unfolded = list_unfolded[0]
        return list_unfolded
