# Import Kratos core and apps
import KratosMultiphysics as km
import KratosMultiphysics.ShapeOptimizationApplication as kso
import KratosMultiphysics.StructuralMechanicsApplication as kcsm
import KratosMultiphysics.MappingApplication as kma

# Additional imports
from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import time, os, shutil, csv, math
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.ShapeOptimizationApplication import custom_variable_utilities as cvu

# =================================================================================================================================================================
# Parameters
# =================================================================================================================================================================
# verification parameters
fd_start_level = 1
fd_end_level = 8
# list_of_move_nodes = [6362, 7680, 12169, 7744, 916, 698, 4547, 8135, 10613] # some selected nodes
# list_of_move_nodes = [65,70,75,82,94,108,119,132,144,157,170,185,196,209,224,239,253,267,281,295,309,324,337,351,366,383,386,387,391,393,399,403,411,413,421,436,448,449,468,485,487,510,528,533,561,577,598,624,634,675,681,722,728,756,766,796,813,841,855,882,894,920,930,958,971,1001,1015,1044,1057,1088,1098,1133,1139,1177,1180,1218,1224,1267,1270,1313,1314,1358,1359,1403,1404,1448,1450,1494,1495,1536,1537,1585,1588,1635,1638,1686,1697,1732,1749,1789,1804,1837,1859,1882,1910,1932,1969,1983,2023,2036,2086,2088,2144,2150,2189,2206,2243,2269,2290,2330,2345,2385,2397,2442,2453,2502,2505,2558,2561,2615,2616,2669,2673,2726,2732,2779,2783,2836,2840,2890,2897,2945,2948,3001,3004,3055,3056,3105,3110,3159,3165,3215,3225,3271,3276,3325,3335,3381,3393,3427,3444,3480,3501,3532,3555,3589,3611,3653,3681,3713,3747,3784,3823,3852,3896,3926,3971,3999,4054,4082,4145,4160,4232,4245,4320,4331,4402,4406,4445,4449,4452,4466,4473,4499,4500,4505,4519,4529,4547,4565,4586,4593,4610,4625,4651,4667,4670,4698,4727,4746,4747,4772,4804,4825,4827,4851,4894,4907,4908,4946,4985,4989,5003,5058,5071,5074,5112,5149,5157,5173,5234,5241,5254,5295,5324,5337,5355,5417,5420,5434,5490,5506,5525,5555,5597,5606,5623,5677,5680,5715,5738,5772,5807,5817,5857,5875,5923,5940,5949,6011,6024,6039,6082,6123,6127,6151,6212,6217,6233,6294,6306,6344,6362,6402,6436,6453,6485,6515,6575,6583,6601,6686,6690,6726,6802,6810,6868,6881,6907,6980,7003,7008,7061,7101,7136,7146,7172,7226,7255,7258,7303,7332,7361,7366,7390,7454,7465,7478,7504,7527,7558,7580,7583,7628,7647,7669,7680,7692,7722,7739,7767,7781,7783,7815,7829,7858,7870,7881,7894,7911,7938,7949,7979,7988,7991,8006,8016,8032,8041,8059,8066,8072,8089,8094,8107,8109,8119,8122,8129,8135,8140,8143,8146,8150,8152,8157,8160,8166,8170,8174,8177,8183,8185,8203,8230,8257,8284,8312,8341,8367,8393,8421,8449,8477,8505,8533,8560,8567,8569,8571,8573,8582,8587,8594,8600,8607,8612,8624,8633,8637,8647,8659,8669,8684,8702,8711,8724,8747,8751,8768,8791,8798,8831,8836,8864,8876,8887,8920,8925,8963,8970,8999,9009,9034,9063,9077,9110,9115,9152,9156,9190,9208,9231,9256,9273,9300,9322,9362,9366,9405,9414,9444,9463,9492,9510,9537,9557,9574,9602,9613,9644,9650,9690,9693,9738,9746,9782,9791,9830,9839,9873,9892,9915,9938,9958,9991,10009,10041,10059,10082,10100,10132,10145,10182,10199,10232,10249,10281,10299,10334,10349,10383,10395,10430,10445,10479,10492,10528,10541,10570,10589,10613,10637,10664,10684,10710,10729,10759,10775,10801,10817,10850,10864,10893,10905,10932,10953,10986,11001,11030,11041,11071,11084,11114,11130,11150,11167,11196,11213,11238,11249,11275,11301,11321,11335,11359,11371,11403,11410,11439,11448,11475,11485,11507,11522,11547,11562,11585,11596,11624,11632,11652,11670,11693,11701,11719,11729,11747,11757,11778,11792,11809,11818,11839,11847,11866,11872,11889,11900,11916,11925,11948,11952,11967,11974,11991,12002,12013,12018,12032,12039,12051,12053,12065,12073,12083,12086,12091,12095,12102,12106,12116,12123,12130,12131,12132] # entire wet interface
# list_of_move_nodes = [3,6,16,29,44,63,88,105,122,139,153,167,181,197,211,226,243,257,271,286,301,315,330,342,357,372,397,429,463,504,542,592,646,659,661,668,682,687,698,700,711,730,743,745,760,782,789,800,827,832,853,866,878,907,912,935,950,967,986,1003,1029,1040,1073,1081,1114,1123,1150,1166,1194,1209,1236,1257,1286,1303,1330,1349,1374,1399,1419,1442,1463,1487,1507,1532,1556,1582,1607,1628,1659,1680,1715,1730,1768,1786,1821,1834,1873,1879,1928,1931,1980,1988,2033,2049,2084,2110,2140,2170,2183,2231,2238,2280,2282,2334,2341,2377,2400,2428,2459,2469,2513,2520,2568,2573,2613,2624,2664,2681,2710,2738,2758,2789,2812,2849,2866,2902,2919,2955,2967,3007,3021,3060,3075,3112,3130,3161,3181,3214,3239,3269,3297,3323,3355,3377,3409,3422,3455,3478,3514,3530,3566,3583,3625,3630,3632,3637,3640,3646,3656,3665,3674,3690,3696,3701,3708,3729,3743,3763,3766,3772,3801,3829,3833,3836,3860,3895,3904,3914,3945,3976,3978,3991,4020,4048,4073,4076,4121,4128,4155,4192,4212,4244,4246,4290,4293,4332,4357,4378,4410,4421,4458,4470,4506,4539,4552,4594,4603,4632,4647,4666,4713,4719,4738,4768,4798,4809,4812,4873,4878,4892,4926,4958,4967,4975,5037,5038,5044,5095,5121,5126,5145,5200,5206,5209,5265,5274,5299,5325,5351,5387,5390,5427,5449,5486,5507,5516,5578,5580,5585,5638,5649,5670,5700,5721,5767,5768,5794,5833,5860,5867,5901,5928,5966,5969,5996,6027,6068,6076,6092,6117,6155,6177,6180,6205,6234,6261,6281,6292,6319,6339,6369,6394,6398,6419,6438,6462,6480,6504,6510,6528,6549,6559,6586,6613,6627,6640,6641,6656,6668,6682,6700,6716,6728,6735,6742,6753,6757,6765,6769,6777,6781,6792,6937,7058,7184,7307,7412,7528,7640,7729,7831,7933,8029,8116,8192,8222,8250,8276,8303,8332,8360,8387,8417,8442,8473,8502,8530,8559,8603,8639,8679,8721,8760,8808,8853,8894,8942,8991,9050,9096,9148,9207,9261,9309,9312,9317,9320,9324,9332,9345,9351,9368,9381,9385,9400,9424,9432,9437,9465,9486,9489,9521,9543,9546,9576,9593,9610,9637,9642,9679,9702,9721,9754,9768,9812,9817,9861,9876,9911,9928,9950,9989,10000,10048,10051,10095,10105,10141,10164,10194,10228,10248,10287,10304,10350,10360,10407,10412,10458,10474,10515,10539,10566,10596,10622,10656,10681,10723,10744,10788,10794,10854,10855,10906,10916,10974,10987,11032,11049,11099,11118,11158,11180,11226,11248,11295,11326,11361,11393,11434,11462,11500,11534,11578,11610,11654,11695,11731,11769,11813,11850,11897,11942,11999,12044,12118,12163,12207,12248,12296,12342,12397,12432,12474,12505,12547,12590,12643,12669,12702,12730,12763,12793,12838,12881,12910,12929,12957,12982,13017,13038,13080,13104,13137,13153,13176,13194,13218,13235,13264,13281,13304,13319,13343,13363,13386,13402,13425,13441,13465,13480,13502,13516,13537,13554,13576,13600,13618,13642,13662,13675,13692,13701,13715,13726,13742,13755,13773,13792,13812,13820,13830,13837,13847,13856,13864,13871,13882,13891,13899,13908,13916,13926,13939,13954,13957,13961,13966,13969,13972,13975,13978,13981,13984,13987,13990,13993,13996,13999,14002,14005,14008,14011,14014,14017,14020,14023,14026,14029,14032,14035,14039,14042,14046,14050,14054,14055,14056] # entire outer ring
list_of_move_nodes = [65,70,75,82,94,108,119,132,144,157,170,185,196,209,224,239,253,267,281,295,309,324,337,351,366,383,386,387,391,393,399,403,411,413,421,436,448,449,468,485,487,510,528,533,561,577,598,624,634,675,681,722,728,756,766,796,813,841,855,882,894,920,930,958,971,1001,1015,1044,1057,1088,1098,1133,1139,1177,1180,1218,1224,1267,1270,1313,1314,1358,1359,1403,1404,1448,1450,1494,1495,1536,1537,1585,1588,1635,1638,1686,1697,1732,1749,1789,1804,1837,1859,1882,1910,1932,1969,1983,2023,2036,2086,2088,2144,2150,2189,2206,2243,2269,2290,2330,2345,2385,2397,2442,2453,2502,2505,2558,2561,2615,2616,2669,2673,2726,2732,2779,2783,2836,2840,2890,2897,2945,2948,3001,3004,3055,3056,3105,3110,3159,3165,3215,3225,3271,3276,3325,3335,3381,3393,3427,3444,3480,3501,3532,3555,3589,3611,3653,3681,3713,3747,3784,3823,3852,3896,3926,3971,3999,4054,4082,4145,4160,4232,4245,4320,4331,4402,4406,4445,4449,4452,4466,4473,4499,4500,4505,4519,4529,4547,4565,4586,4593,4610,4625,4651,4667,4670,4698,4727,4746,4747,4772,4804,4825,4827,4851,4894,4907,4908,4946,4985,4989,5003,5058,5071,5074,5112,5149,5157,5173,5234,5241,5254,5295,5324,5337,5355,5417,5420,5434,5490,5506,5525,5555,5597,5606,5623,5677,5680,5715,5738,5772,5807,5817,5857,5875,5923,5940,5949,6011,6024,6039,6082,6123,6127,6151,6212,6217,6233,6294,6306,6344,6362,6402,6436,6453,6485,6515,6575,6583,6601,6686,6690,6726,6802,6810,6868,6881,6907,6980,7003,7008,7061,7101,7136,7146,7172,7226,7255,7258,7303,7332,7361,7366,7390,7454,7465,7478,7504,7527,7558,7580,7583,7628,7647,7669,7680,7692,7722,7739,7767,7781,7783,7815,7829,7858,7870,7881,7894,7911,7938,7949,7979,7988,7991,8006,8016,8032,8041,8059,8066,8072,8089,8094,8107,8109,8119,8122,8129,8135,8140,8143,8146,8150,8152,8157,8160,8166,8170,8174,8177,8183,8185,8203,8230,8257,8284,8312,8341,8367,8393,8421,8449,8477,8505,8533,8560,8567,8569,8571,8573,8582,8587,8594,8600,8607,8612,8624,8633,8637,8647,8659,8669,8684,8702,8711,8724,8747,8751,8768,8791,8798,8831,8836,8864,8876,8887,8920,8925,8963,8970,8999,9009,9034,9063,9077,9110,9115,9152,9156,9190,9208,9231,9256,9273,9300,9322,9362,9366,9405,9414,9444,9463,9492,9510,9537,9557,9574,9602,9613,9644,9650,9690,9693,9738,9746,9782,9791,9830,9839,9873,9892,9915,9938,9958,9991,10009,10041,10059,10082,10100,10132,10145,10182,10199,10232,10249,10281,10299,10334,10349,10383,10395,10430,10445,10479,10492,10528,10541,10570,10589,10613,10637,10664,10684,10710,10729,10759,10775,10801,10817,10850,10864,10893,10905,10932,10953,10986,11001,11030,11041,11071,11084,11114,11130,11150,11167,11196,11213,11238,11249,11275,11301,11321,11335,11359,11371,11403,11410,11439,11448,11475,11485,11507,11522,11547,11562,11585,11596,11624,11632,11652,11670,11693,11701,11719,11729,11747,11757,11778,11792,11809,11818,11839,11847,11866,11872,11889,11900,11916,11925,11948,11952,11967,11974,11991,12002,12013,12018,12032,12039,12051,12053,12065,12073,12083,12086,12091,12095,12102,12106,12116,12123,12130,12131,12132, 3,6,16,29,44,63,88,105,122,139,153,167,181,197,211,226,243,257,271,286,301,315,330,342,357,372,397,429,463,504,542,592,646,659,661,668,682,687,698,700,711,730,743,745,760,782,789,800,827,832,853,866,878,907,912,935,950,967,986,1003,1029,1040,1073,1081,1114,1123,1150,1166,1194,1209,1236,1257,1286,1303,1330,1349,1374,1399,1419,1442,1463,1487,1507,1532,1556,1582,1607,1628,1659,1680,1715,1730,1768,1786,1821,1834,1873,1879,1928,1931,1980,1988,2033,2049,2084,2110,2140,2170,2183,2231,2238,2280,2282,2334,2341,2377,2400,2428,2459,2469,2513,2520,2568,2573,2613,2624,2664,2681,2710,2738,2758,2789,2812,2849,2866,2902,2919,2955,2967,3007,3021,3060,3075,3112,3130,3161,3181,3214,3239,3269,3297,3323,3355,3377,3409,3422,3455,3478,3514,3530,3566,3583,3625,3630,3632,3637,3640,3646,3656,3665,3674,3690,3696,3701,3708,3729,3743,3763,3766,3772,3801,3829,3833,3836,3860,3895,3904,3914,3945,3976,3978,3991,4020,4048,4073,4076,4121,4128,4155,4192,4212,4244,4246,4290,4293,4332,4357,4378,4410,4421,4458,4470,4506,4539,4552,4594,4603,4632,4647,4666,4713,4719,4738,4768,4798,4809,4812,4873,4878,4892,4926,4958,4967,4975,5037,5038,5044,5095,5121,5126,5145,5200,5206,5209,5265,5274,5299,5325,5351,5387,5390,5427,5449,5486,5507,5516,5578,5580,5585,5638,5649,5670,5700,5721,5767,5768,5794,5833,5860,5867,5901,5928,5966,5969,5996,6027,6068,6076,6092,6117,6155,6177,6180,6205,6234,6261,6281,6292,6319,6339,6369,6394,6398,6419,6438,6462,6480,6504,6510,6528,6549,6559,6586,6613,6627,6640,6641,6656,6668,6682,6700,6716,6728,6735,6742,6753,6757,6765,6769,6777,6781,6792,6937,7058,7184,7307,7412,7528,7640,7729,7831,7933,8029,8116,8192,8222,8250,8276,8303,8332,8360,8387,8417,8442,8473,8502,8530,8559,8603,8639,8679,8721,8760,8808,8853,8894,8942,8991,9050,9096,9148,9207,9261,9309,9312,9317,9320,9324,9332,9345,9351,9368,9381,9385,9400,9424,9432,9437,9465,9486,9489,9521,9543,9546,9576,9593,9610,9637,9642,9679,9702,9721,9754,9768,9812,9817,9861,9876,9911,9928,9950,9989,10000,10048,10051,10095,10105,10141,10164,10194,10228,10248,10287,10304,10350,10360,10407,10412,10458,10474,10515,10539,10566,10596,10622,10656,10681,10723,10744,10788,10794,10854,10855,10906,10916,10974,10987,11032,11049,11099,11118,11158,11180,11226,11248,11295,11326,11361,11393,11434,11462,11500,11534,11578,11610,11654,11695,11731,11769,11813,11850,11897,11942,11999,12044,12118,12163,12207,12248,12296,12342,12397,12432,12474,12505,12547,12590,12643,12669,12702,12730,12763,12793,12838,12881,12910,12929,12957,12982,13017,13038,13080,13104,13137,13153,13176,13194,13218,13235,13264,13281,13304,13319,13343,13363,13386,13402,13425,13441,13465,13480,13502,13516,13537,13554,13576,13600,13618,13642,13662,13675,13692,13701,13715,13726,13742,13755,13773,13792,13812,13820,13830,13837,13847,13856,13864,13871,13882,13891,13899,13908,13916,13926,13939,13954,13957,13961,13966,13969,13972,13975,13978,13981,13984,13987,13990,13993,13996,13999,14002,14005,14008,14011,14014,14017,14020,14023,14026,14029,14032,14035,14039,14042,14046,14050,14054,14055,14056] # combined inner and outer ring

# analysis parameter
cfd_part_name = "ubend_2d"
cfd_interface_part_name = "design_surface"
cfd_wall_part_name = "WALL"
csm_interface_part_name = "wet_interface"

reference_pressure = 1#74571.25

load_mapper_settings = km.Parameters("""
{
    "mapper_type" : "nearest_element",
    "echo_level"  : 0
}""")
stress_response_settings = km.Parameters("""
{
        "response_type"                   : "adjoint_max_stress",
        "gradient_mode"                   : "semi_analytic",
        "step_size"                       : 1e-11,
        "critical_part_name"              : "stress_partition",
        "stress_type"                     : "VON_MISES_STRESS",
        "stress_treatment"                : "mean",
        "echo_level"                      : 1,
        "primal_settings"                 : "parameters_csm_analysis.json",
        "adjoint_settings"                : "auto",
        "primal_data_transfer_with_python": true,
        "sensitivity_settings" : {
            "sensitivity_model_part_name"     : "Parts_structure",
            "nodal_sensitivity_variables"     : ["SHAPE_SENSITIVITY"],
            "element_sensitivity_variables"   : [],
            "condition_sensitivity_variables" : [],
            "build_mode": "static"
        }
}""")

# Settings from SU2 config
SURFACE_FLOW_FILENAME = "surface_flow"

with open("parameters_fsi_optimization.json",'r') as parameter_file:
    parameters = km.Parameters(parameter_file.read())

# Change threading layer, otherwise the suprocess module used in the interface.py in SU2 will hang, (only)
# when Intel solvers are used as linear solvers for calculating the structure. This is a known issue of
# the Intel compiler 2018 and will be fixed in future releases. Compare:
# 1) https://github.com/numpy/numpy/issues/10060
# 2) https://software.intel.com/en-us/forums/intel-c-compiler/topic/758961
os.environ["MKL_THREADING_LAYER"] = "TBB"

# =================================================================================================================================================================
# Helper functions
# =================================================================================================================================================================
def OutputCFDAsVTK( output_mdpa, output_filename ):
    output_parameters = km.Parameters("""
    {
        "model_part_name": "XXX",
        "file_format": "ascii",
        "output_precision": 7,
        "output_control_type": "step",
        "output_frequency": 1.0,
        "output_sub_model_parts": true,
        "folder_name": "XXX",
        "save_output_files_in_folder": true,
        "nodal_solution_step_data_variables": ["SHAPE_UPDATE","NORMALIZED_SURFACE_NORMAL","POINT_LOAD"],
        "nodal_data_value_variables": [],
        "element_data_value_variables": [],
        "condition_data_value_variables": []
    }""")

    output_parameters["model_part_name"].SetString(output_mdpa.Name)
    output_parameters["folder_name"].SetString(output_filename)

    output_model = {}
    output_model[output_mdpa.Name] = output_mdpa

    vtk_output = VtkOutputProcess( output_model, output_parameters )
    vtk_output.ExecuteInitialize()
    vtk_output.ExecuteBeforeSolutionLoop()
    vtk_output.ExecuteInitializeSolutionStep()
    vtk_output.PrintOutput()
    vtk_output.ExecuteFinalizeSolutionStep()
    vtk_output.ExecuteFinalize()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
def OutputOPTAsVTK( output_mdpa, output_filename ):
    output_parameters = km.Parameters("""
    {
        "model_part_name": "XXX",
        "file_format": "ascii",
        "output_precision": 7,
        "output_control_type": "step",
        "output_frequency": 1.0,
        "output_sub_model_parts": true,
        "folder_name": "XXX",
        "save_output_files_in_folder": true,
        "nodal_solution_step_data_variables": ["CONTROL_POINT_UPDATE","SHAPE_UPDATE"],
        "nodal_data_value_variables": [],
        "element_data_value_variables": [],
        "condition_data_value_variables": []
    }""")

    output_parameters["model_part_name"].SetString(output_mdpa.Name)
    output_parameters["folder_name"].SetString(output_filename)

    output_model = {}
    output_model[output_mdpa.Name] = output_mdpa

    vtk_output = VtkOutputProcess( output_model, output_parameters )
    vtk_output.ExecuteInitialize()
    vtk_output.ExecuteBeforeSolutionLoop()
    vtk_output.ExecuteInitializeSolutionStep()
    vtk_output.PrintOutput()
    vtk_output.ExecuteFinalizeSolutionStep()
    vtk_output.ExecuteFinalize()

# =================================================================================================================================================================
# Define analyzer
# =================================================================================================================================================================
class CustomAnalyzer(AnalyzerBaseClass):
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __init__(self):
        self.cfd_interface_part = None
        self.cfd_part = None
        self.cfd_wall_part = None

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        # Initialize cfd data
        cfd_model = km.Model()
        self.cfd_part = cfd_model.CreateModelPart(cfd_part_name)
        self.cfd_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,2)
        self.cfd_part.AddNodalSolutionStepVariable(kcsm.POINT_LOAD)
        self.cfd_part.AddNodalSolutionStepVariable(kso.SHAPE_UPDATE)
        self.cfd_part.AddNodalSolutionStepVariable(km.NORMAL)
        self.cfd_part.AddNodalSolutionStepVariable(kso.NORMALIZED_SURFACE_NORMAL)

        model_part_io = km.ModelPartIO(cfd_part_name)
        model_part_io.ReadModelPart(self.cfd_part)

        self.cfd_interface_part = self.cfd_part.GetSubModelPart(cfd_interface_part_name)
        self.cfd_wall_part = self.cfd_part.GetSubModelPart(cfd_wall_part_name)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration):
        # obtain new csm mesh
        new_csm_mesh = {node.Id: [node.X, node.Y, node.Z] for node in current_design.Nodes}

        # Update fluid mesh (structure is controlled by the optimization algorithm) - Note that the mapper should be based on the previos design to match the forward map in the structure
        for node in current_design.Nodes:
            shape_update = node.GetSolutionStepValue(kso.SHAPE_UPDATE)
            node.X -= shape_update[0]
            node.Y -= shape_update[1]
            node.Z -= shape_update[2]
            node.X0 -= shape_update[0]
            node.Y0 -= shape_update[1]
            node.Z0 -= shape_update[2]

        vm_mapper = kso.MapperVertexMorphingMatrixFree(current_design, self.cfd_interface_part, parameters["optimization_settings"]["design_variables"]["filter"])
        vm_mapper.Map(kso.CONTROL_POINT_UPDATE, kso.SHAPE_UPDATE)

        for node in current_design.Nodes:
            shape_update = node.GetSolutionStepValue(kso.SHAPE_UPDATE)
            node.X += shape_update[0]
            node.Y += shape_update[1]
            node.Z += shape_update[2]
            node.X0 += shape_update[0]
            node.Y0 += shape_update[1]
            node.Z0 += shape_update[2]

        kso.MeshControllerUtilities(self.cfd_interface_part).UpdateMeshAccordingInputVariable(kso.SHAPE_UPDATE)
        kso.MeshControllerUtilities(self.cfd_interface_part).SetReferenceMeshToMesh()

        # Initialize SU2 interface
        from KratosMultiphysics.ShapeOptimizationApplication.interface_su2 import InterfaceSU2
        self.interface_su2 = InterfaceSU2(parameters["su2_interface_settings"])
        self.interface_su2.InitializeNewSU2Project()

        self.interface_su2.WriteNodesAsSU2MeshMotionFile(self.cfd_wall_part.GetNodes())

        # Run fluid
        update_mesh = True
        [cfd_response_value] = self.interface_su2.ComputeValues(["SURFACE_TOTAL_PRESSURE"], update_mesh, optimization_iteration)

        # Read and dimensionalize fluid results
        file_with_pressures = "DESIGNS/DSN_"+str(optimization_iteration).zfill(3)+"/DIRECT/"+SURFACE_FLOW_FILENAME+".csv"
        if self.interface_su2.su2_mesh_data["NDIME"] == 2:
            nodal_forces = self.interface_su2.ReadNodalValuesFromCSVFile(file_with_pressures,0,[5,6],1)
            nodal_forces = {key: [value[0],value[1],0.0] for key, value in nodal_forces.items()}
        else:
            nodal_forces = self.interface_su2.ReadNodalValuesFromCSVFile(file_with_pressures,0,[6,7,8],1)

        for node in self.cfd_interface_part.Nodes:
            force = nodal_forces[node.Id]
            dimensional_force = [value*reference_pressure for value in force]
            node.SetSolutionStepValue(kcsm.POINT_LOAD, dimensional_force)

        # Run CSM
        csm_response_value = self.__RunCSM(opt_itr=optimization_iteration, new_mesh=new_csm_mesh)

        # Some output
        OutputCFDAsVTK(self.cfd_interface_part, "cfd_interface_part")

        return cfd_response_value, csm_response_value

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __RunCSM(self, opt_itr, new_mesh):
        print("\n> Starting __RunCSM in analyzer...")
        start_time = time.time()

        analysis_model = km.Model()

        csm = structural_response_function_factory.CreateResponseFunction(stress_response_settings["response_type"].GetString(), stress_response_settings, analysis_model)

        csm_primal_part = csm.primal_model_part
        csm_primal_part.AddNodalSolutionStepVariable(kso.MESH_CHANGE)
        csm_primal_part.AddNodalSolutionStepVariable(kso.NORMALIZED_SURFACE_NORMAL)
        csm_primal_part.AddNodalSolutionStepVariable(km.NORMAL)

        # Initialize and set time to current optimization iteration
        csm.primal_analysis.project_parameters["problem_data"]["start_time"].SetDouble(opt_itr-1)
        csm.primal_analysis.project_parameters["problem_data"]["end_time"].SetDouble(opt_itr)

        csm.Initialize()

        csm_primal_part.ProcessInfo.SetValue(km.STEP, opt_itr-1)

        # apply mesh motion
        if new_mesh != {}:
            for node in csm_primal_part.Nodes:
                X_new = new_mesh[node.Id]
                mesh_change = [X_new[0]-node.X, X_new[1]-node.Y, X_new[2]-node.Z]
                node.SetSolutionStepValue(kso.MESH_CHANGE, mesh_change)
                node.X = X_new[0]
                node.Y = X_new[1]
                node.Z = X_new[2]
                node.X0 = X_new[0]
                node.Y0 = X_new[1]
                node.Z0 = X_new[2]

        # compute primal field

        # Map forces
        csm_interface_part = csm_primal_part.GetSubModelPart(csm_interface_part_name)
        csm_interface_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,2)
        kso.GeometryUtilities(csm_interface_part).ComputeUnitSurfaceNormals()
        csm_interface_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,3)

        cfd_to_csm_mapper = kma.MapperFactory.CreateMapper(csm_interface_part, self.cfd_interface_part, load_mapper_settings.Clone())
        cfd_to_csm_mapper.InverseMap(kcsm.POINT_LOAD, kcsm.POINT_LOAD, kma.Mapper.USE_TRANSPOSE)

        # Compute primals
        csm.InitializeSolutionStep()
        csm.CalculateValue()
        csm.FinalizeSolutionStep()
        csm.Finalize()

        print("> Finished __RunCSM in" ,round( time.time()-start_time, 3 ), " s.")

        return csm.GetValue()

# =================================================================================================================================================================
# Perform verification
# =================================================================================================================================================================
results_filename = "fd_results.txt"
file_to_write = open(results_filename,'w')
file_to_write.close()

for move_node_id in list_of_move_nodes:
    # Compute reference
    print("\n\n##########################################################################")
    print(">> Starting to compute reference of node "+str(move_node_id))
    print("##########################################################################")
    start_time = time.time()

    # Create case folder
    case_dir = "tmp_calculation_reference"
    if os.path.exists(case_dir):
        os.system("rm -rf "+case_dir+"/")
    os.system("mkdir "+case_dir)

    # Copy all necessary files to case folder
    shutil.copy("ubend_2d.mdpa", case_dir+"/")
    shutil.copy("CSM_model.mdpa", case_dir+"/")
    shutil.copy("OPT_model.mdpa", case_dir+"/")
    shutil.copy("parameters_csm_analysis.json", case_dir+"/")
    shutil.copy("parameters_material.json", case_dir+"/")
    shutil.copy("parameters_fsi_optimization.json", case_dir+"/")
    shutil.copy("restart_flow.dat", case_dir+"/")
    shutil.copy("ubend_coarse.cfg", case_dir+"/")
    shutil.copy("ubend_2d.su2", case_dir+"/")

    # Change to case directory
    top_dir = os.getcwd()
    os.chdir(case_dir)

    # Create analyzer
    analyzer = CustomAnalyzer()

    # Read optimization model part
    model = km.Model()
    opt_part = model.CreateModelPart("OPT_model")
    opt_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,3)
    opt_part.AddNodalSolutionStepVariable(kso.CONTROL_POINT_UPDATE)
    opt_part.AddNodalSolutionStepVariable(kso.SHAPE_UPDATE)

    model_part_io = km.ModelPartIO("OPT_model")
    model_part_io.ReadModelPart(opt_part)

    # Run analysis
    analyzer.InitializeBeforeOptimizationLoop()
    cfd_response_value_reference, csm_response_value_reference = analyzer.AnalyzeDesignAndReportToCommunicator(opt_part, 1)

    # Save values in file
    with open("cfd_value.txt",'w') as results_file:
        results_file.write(str(cfd_response_value_reference))
    with open("csm_value.txt",'w') as results_file:
        results_file.write(str(csm_response_value_reference))

    # Change back to top dir and save results
    os.chdir(top_dir)

    with open(results_filename,'a') as results_file:
        results_file.write("#####################################################################\n")
        results_file.write("node_id = "+str(move_node_id)+"\n")
        results_file.write("cfd_reference_value = "+str(cfd_response_value_reference)+"\n")
        results_file.write("csm_reference_value = "+str(csm_response_value_reference)+"\n")
        results_file.write("---\n")
        results_file.write("step_size,cfd_gradient_x,cfd_gradient_y,csm_gradient_x,csm_gradient_y,time\n")

    print("\n> Finished computing reference value in" ,round( time.time()-start_time, 3 ), " s.")
    print("> cfd reference response value =",cfd_response_value_reference)
    print("> csm reference response value =",csm_response_value_reference)
    print("\n")

    err("temporary break")

    # Compute pertubation
    for itr in range(fd_start_level,fd_end_level+1):

        # Compute pertubation
        current_delta = 1*10**(-itr)

        for dim_itr in range(2):

            print("\n\n##########################################################################")
            print(">> Starting finite difference level ",itr, " dimension ",dim_itr+1, "for node ", str(move_node_id))
            print("##########################################################################")

            # Create case folder
            case_dir = "tmp_calculation_"+str(itr)+"_"+str(dim_itr+1)
            if os.path.exists(case_dir):
                os.system("rm -rf "+case_dir+"/")
            os.system("mkdir "+case_dir)

            # Copy all necessary files to case folder
            shutil.copy("ubend_2d.mdpa", case_dir+"/")
            shutil.copy("CSM_model.mdpa", case_dir+"/")
            shutil.copy("OPT_model.mdpa", case_dir+"/")
            shutil.copy("parameters_csm_analysis.json", case_dir+"/")
            shutil.copy("parameters_material.json", case_dir+"/")
            shutil.copy("parameters_fsi_optimization.json", case_dir+"/")
            shutil.copy("restart_flow.dat", case_dir+"/")
            shutil.copy("ubend_coarse.cfg", case_dir+"/")
            shutil.copy("ubend_2d.su2", case_dir+"/")

            # Change to case directory
            top_dir = os.getcwd()
            os.chdir(case_dir)

            # Create analyzer
            analyzer = CustomAnalyzer()

            # Read optimization model part
            model = km.Model()
            opt_part = model.CreateModelPart("OPT_model")
            opt_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,3)
            opt_part.AddNodalSolutionStepVariable(kso.CONTROL_POINT_UPDATE)
            opt_part.AddNodalSolutionStepVariable(kso.SHAPE_UPDATE)

            model_part_io = km.ModelPartIO("OPT_model")
            model_part_io.ReadModelPart(opt_part)

            # Apply control point update
            if dim_itr == 0:
                opt_part.Nodes[move_node_id].SetSolutionStepValue(kso.CONTROL_POINT_UPDATE, [current_delta,0,0])
            elif dim_itr == 1:
                opt_part.Nodes[move_node_id].SetSolutionStepValue(kso.CONTROL_POINT_UPDATE, [0,current_delta,0])

            # Compute corresponding shape update
            vm_mapper = kso.MapperVertexMorphingMatrixFree(opt_part, opt_part, parameters["optimization_settings"]["design_variables"]["filter"])
            vm_mapper.Map(kso.CONTROL_POINT_UPDATE, kso.SHAPE_UPDATE)

            OutputOPTAsVTK(opt_part, "opt_part")

            # Update mesh
            kso.MeshControllerUtilities(opt_part).UpdateMeshAccordingInputVariable(kso.SHAPE_UPDATE)
            kso.MeshControllerUtilities(opt_part).SetReferenceMeshToMesh()

            # Run analysis
            analyzer.InitializeBeforeOptimizationLoop()
            cfd_response_value, csm_response_value = analyzer.AnalyzeDesignAndReportToCommunicator(opt_part, 1)

            # Evaluate fd gradients
            if dim_itr == 0:
                cfd_grad_x = (cfd_response_value - cfd_response_value_reference) / current_delta
                csm_grad_x = (csm_response_value - csm_response_value_reference) / current_delta
            elif dim_itr == 1:
                cfd_grad_y = (cfd_response_value - cfd_response_value_reference) / current_delta
                csm_grad_y = (csm_response_value - csm_response_value_reference) / current_delta

            print("\n> current cfd response value =",cfd_response_value)
            print("> current csm response value =",csm_response_value)
            print("\n")

            # Change back to top dir
            os.chdir(top_dir)

        # Write results
        with open(results_filename,'a') as results_file:
            results_file.write(str(current_delta)+","+str(cfd_grad_x)+","+str(cfd_grad_y)+","+str(csm_grad_x)+","+str(csm_grad_y)+","+str(time.ctime())+"\n")

    print("\n\n##########################################################################")
    print("\n> Finished verifiation process in" ,round( time.time()-start_time, 3 ), " s.")
    print("##########################################################################")
