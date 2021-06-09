import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn')

history = "Optimization_Results"

plots = [
    {
        "name": "target",
        "column_names": [
            "f"
        ]
    },
    {
        "name": "Node 15: z <= -1.0",
        "column_names": [
            "c1: <=",
            "c1_ref",
        ]
    },
    {
        "name": "Node 40: z >= 1.0",
        "column_names": [
            "c2: >=",
            "c2_ref",
        ]
    },
    {
        "name": "Node 65: z = 1.0",
        "column_names": [
            "c3: =",
            "c3_ref",
        ]
    },
    # {
    #     "name": "step length",
    #     "column_names": [
    #         "step_size",
    #         "inf_norm_s",
    #         "inf_norm_c",
    #     ]
    # }
]

path = f"{history}/optimization_log.csv"

df = pd.read_csv(path, delimiter=",")
df.columns = [x.strip() for x in df.columns]


fig, axs = plt.subplots(1, len(plots), figsize=(18, 6))

for plot, ax in zip(plots, axs):
    for column_name in plot["column_names"]:
        ax.plot(df[column_name], label=column_name)
    ax.legend()
    ax.set_title(plot["name"])

plt.tight_layout()
plt.show()
