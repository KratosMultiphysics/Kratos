# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3
from pycompss.api.parameter import *

# import packages
from functools import wraps
import numpy as np

####################################################################################################
############################################# CLASSES ##############################################
####################################################################################################

OUTPUT_QUANTITIES = 1

class UnfolderManager(object):

    def __init__(self, number):
        global OUTPUT_QUANTITIES
        OUTPUT_QUANTITIES = self.number = number

    def UnfoldNValues(self, number, values):
      for val in range(0, number):
        yield values[val]

    @ExaquteTask(target_direction=IN,returns='OUTPUT_QUANTITIES')
    def UnfoldNValues_Task(self, values):
        return list(self.UnfoldNValues(self.number,values))


    ####################################################################################################
    ############################################ WRAPPERS ##############################################
    ####################################################################################################

def unfoldVales_Wrapper(number_values,values):
    if (number_values == 1):
        value_1 = Unfold1ValuesAux_Task(values)
        qoi_list = [value_1]
    elif (number_values == 2):
        value_1,value_2 = Unfold2ValuesAux_Task(values)
        qoi_list = [value_1,value_2]
    elif (number_values == 93):
        value_1,value_2,value_3,value_4,value_5,value_6,value_7,value_8,value_9,value_10,value_11,value_12,value_13,value_14,value_15,value_16,value_17,value_18,value_19,value_20,value_21,value_22,value_23,value_24,value_25,value_26,value_27,value_28,value_29,value_30,value_31,value_32,value_33,value_34,value_35,value_36,value_37,value_38,value_39,value_40,value_41,value_42,value_43,value_44,value_45,value_46,value_47,value_48,value_49,value_50,value_51,value_52,value_53,value_54,value_55,value_56,value_57,value_58,value_59,value_60,value_61,value_62,value_63,value_64,value_65,value_66,value_67,value_68,value_69,value_70,value_71,value_72,value_73,value_74,value_75,value_76,value_77,value_78,value_79,value_80,value_81,value_82,value_83,value_84,value_85,value_86,value_87,value_88,value_89,value_90,value_91,value_92,value_93 = Unfold93ValuesAux_Task(values)
        qoi_list = [value_1,value_2,value_3,value_4,value_5,value_6,value_7,value_8,value_9,value_10,value_11,value_12,value_13,value_14,value_15,value_16,value_17,value_18,value_19,value_20,value_21,value_22,value_23,value_24,value_25,value_26,value_27,value_28,value_29,value_30,value_31,value_32,value_33,value_34,value_35,value_36,value_37,value_38,value_39,value_40,value_41,value_42,value_43,value_44,value_45,value_46,value_47,value_48,value_49,value_50,value_51,value_52,value_53,value_54,value_55,value_56,value_57,value_58,value_59,value_60,value_61,value_62,value_63,value_64,value_65,value_66,value_67,value_68,value_69,value_70,value_71,value_72,value_73,value_74,value_75,value_76,value_77,value_78,value_79,value_80,value_81,value_82,value_83,value_84,value_85,value_86,value_87,value_88,value_89,value_90,value_91,value_92,value_93]
    else:
        raise Exception ("Number of QoI not supported. Add it in xmc/classDefs_solverWrapper/methodDefs_KratosSolverWrapper/utilities.py")
    return qoi_list

    ####################################################################################################
    ############################################## TASKS ###############################################
    ####################################################################################################

@ExaquteTask(returns=1)
def Unfold1ValuesAux_Task(values):
    return values[0]

@ExaquteTask(returns=2)
def Unfold2ValuesAux_Task(values):
    return values[0],values[1]

@ExaquteTask(returns=93)
def Unfold93ValuesAux_Task(values):
    return values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],values[10], \
        values[11],values[12],values[13],values[14],values[15],values[16],values[17],values[18],values[19],values[20], \
        values[21],values[22],values[23],values[24],values[25],values[26],values[27],values[28],values[29],values[30], \
        values[31],values[32],values[33],values[34],values[35],values[36],values[37],values[38],values[39],values[40], \
        values[41],values[42],values[43],values[44],values[45],values[46],values[47],values[48],values[49],values[50], \
        values[51],values[52],values[53],values[54],values[55],values[56],values[57],values[58],values[59],values[60], \
        values[61],values[62],values[63],values[64],values[65],values[66],values[67],values[68],values[69],values[70], \
        values[71],values[72],values[73],values[74],values[75],values[76],values[77],values[78],values[79],values[80], \
        values[81],values[82],values[83],values[84],values[85],values[86],values[87],values[88],values[89],values[90], \
        values[91],values[92]

@ExaquteTask(aux_qoi_array={Type: COLLECTION_IN, Depth: 2},returns=1)
def PostprocessContributionsPerInstance(aux_qoi_array,number_qoi,number_time_power_sums):
    assert(len(aux_qoi_array)==number_qoi+number_time_power_sums)
    qoi_list = []
    # expected value for steady state / time averaged qoi
    for qoi_counter in range (0,number_qoi):
        qoi_value = np.mean(aux_qoi_array[qoi_counter])
        qoi_list.append(qoi_value)
    # time power sums sum for time series qoi    
    for time_ps_counter in range (number_qoi,number_qoi+number_time_power_sums):
        S1 = 0 ; S2 = 0 ; S3 = 0 ; S4 = 0 ; contributions = 0
        for i in range (0,len(aux_qoi_array[time_ps_counter])):
            S1 = S1 + aux_qoi_array[time_ps_counter][i][0][0]
            S2 = S2 + aux_qoi_array[time_ps_counter][i][1][0]
            if (len(aux_qoi_array[time_ps_counter][i]) == 5):
                S3 = S3 + aux_qoi_array[time_ps_counter][i][2][0]
                S4 = S4 + aux_qoi_array[time_ps_counter][i][3][0]
            contributions = contributions + aux_qoi_array[time_ps_counter][i][-1]
        qoi_list.append([[S1],[S2],contributions])
    return qoi_list