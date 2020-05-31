import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn')

# parent_dir = "/home/arming/dev/Examples/shape_optimization"
# case_name = "01_Strain_Energy_Minimization_3D_Hook"
histories = [
    "Optimization_Results",
    "Optimization_Results_05",
    "Optimization_Results_075",
    "Optimization_Results_old"
]


parent_dir = "/home/arming/dev/Kratos/applications/ShapeOptimizationApplication/tests"
case_name = "algorithm_trust_region_test"
histories = [
    "Optimization_Results"
]
column_names = [
    "c1: >",
    "c1_ref",
    "c2: =",
    "c2_ref"
    # "f"
]

for history in histories:
    path = f"{parent_dir}/{case_name}/{history}/optimization_log.csv"

    df = pd.read_csv(path, delimiter=",")
    df.columns = [x.strip() for x in df.columns]

    for column_name in column_names:
        plt.plot(df[column_name], label=column_name+history.replace("Optimization_Results", ""))

plt.legend()
plt.show()
