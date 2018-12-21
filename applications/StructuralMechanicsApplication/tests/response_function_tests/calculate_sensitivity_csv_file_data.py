import csv

response_value_1 = []
with open('response_values_original_geometry.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_1 = list(reader)

response_value_2 = []
with open('response_values_perturb_node1X_0.0001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_2 = list(reader)

response_value_3 = []
with open('response_values_perturb_node1Y_0.0001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_3 = list(reader)

response_value_4 = []
with open('response_values_perturb_node1Z_0.0001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_4 = list(reader)

response_value_5 = []
with open('response_values_perturb_node2X_0.0001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_5 = list(reader)

response_value_6 = []
with open('response_values_perturb_node2Y_0.0001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_6 = list(reader)

response_value_7 = []
with open('response_values_perturb_node2Z_0.0001.csv', "r") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    response_value_7 = list(reader)

# Finite difference calculation
delta = 0.0001

response_sensitivity_1_X = (response_value_2[9][1] - response_value_1[9][1]) / delta
print("response_sensitivity_1_X : " , response_sensitivity_1_X)

response_sensitivity_1_Y = (response_value_3[9][1] - response_value_1[9][1]) / delta
print("response_sensitivity_1_Y : " , response_sensitivity_1_Y)

response_sensitivity_1_Z = (response_value_4[9][1] - response_value_1[9][1]) / delta
print("response_sensitivity_1_Z : " , response_sensitivity_1_Z)

response_sensitivity_2_X = (response_value_5[9][1] - response_value_1[9][1]) / delta
print("response_sensitivity_2_x : " , response_sensitivity_2_X)

response_sensitivity_2_Y = (response_value_6[9][1] - response_value_1[9][1]) / delta
print("response_sensitivity_2_Y : " , response_sensitivity_2_Y)

response_sensitivity_2_Z = (response_value_7[9][1] - response_value_1[9][1]) / delta
print("response_sensitivity_2_Z : " , response_sensitivity_2_Z)