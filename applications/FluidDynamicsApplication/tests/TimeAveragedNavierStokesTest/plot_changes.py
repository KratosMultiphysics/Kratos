import matplotlib.pyplot as plt


datafile = open('cylinder/log.txt')
velocity_change = []
pressure_change = []
s_length = []
time = []

for line in datafile:
    if "TIME:  " in line:
        time.append(float(line[7:-1]))
    if "Averaging time length set to  " in line:
        s_length.append(float(line[30:36]))
    if "Change in velocity in percentage: " in line:
        velocity_change.append(float(line[34:-1]))
    if "Change in pressure in percentage: " in line:
        pressure_change.append(float(line[34:-1]))
        print(pressure_change[-1])

plt.figure(1)
plt.plot(time, velocity_change, label="velocity change")
plt.plot(time, pressure_change, label="pressure change")
plt.ylim([0, 0.1])
plt.xlabel("TIME")
plt.legend()
plt.show()

        

