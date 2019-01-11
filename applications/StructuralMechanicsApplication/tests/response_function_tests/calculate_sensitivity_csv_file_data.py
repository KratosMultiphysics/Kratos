import csv

response_value_1 = []
with open('response_values_original_geometry.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_1 = list(reader)

response_value_2 = []
with open('response_values_perturb_node1X_0.00000001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_2 = list(reader)

response_value_3 = []
with open('response_values_perturb_node1Y_0.00000001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_3 = list(reader)

response_value_4 = []
with open('response_values_perturb_node1Z_0.00000001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_4 = list(reader)

response_value_5 = []
with open('response_values_perturb_node2X_0.00000001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_5 = list(reader)

response_value_6 = []
with open('response_values_perturb_node2Y_0.00000001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_6 = list(reader)

response_value_7 = []
with open('response_values_perturb_node2Z_0.00000001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_7 = list(reader)

# Finite difference calculation
# delta = 0.00000001

# response_sensitivity_1_X = (response_value_2[49][1] - response_value_1[49][1]) / delta
# print("response_sensitivity_1_X : " , response_sensitivity_1_X)

# response_sensitivity_1_Y = (response_value_3[49][1] - response_value_1[49][1]) / delta
# print("response_sensitivity_1_Y : " , response_sensitivity_1_Y)

# response_sensitivity_1_Z = (response_value_4[49][1] - response_value_1[49][1]) / delta
# print("response_sensitivity_1_Z : " , response_sensitivity_1_Z)

# response_sensitivity_2_X = (response_value_5[49][1] - response_value_1[49][1]) / delta
# print("response_sensitivity_2_x : " , response_sensitivity_2_X)

# response_sensitivity_2_Y = (response_value_6[49][1] - response_value_1[49][1]) / delta
# print("response_sensitivity_2_Y : " , response_sensitivity_2_Y)

# response_sensitivity_2_Z = (response_value_7[49][1] - response_value_1[49][1]) / delta
# print("response_sensitivity_2_Z : " , response_sensitivity_2_Z)


response_value_8 = []
with open('response_values_perturb_single_truss.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_8 = list(reader)

response_value_9 = []
with open('response_values_perturb_single_truss_X.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_9 = list(reader)

delta = 0.00001

response_sensitivity_single_truss = (response_value_9[4][1] - response_value_8[4][1]) / delta
print("response_sensitivity_single truss : " , response_sensitivity_single_truss)