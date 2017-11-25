import matplotlib.pyplot as plt

time = []
energy = []
with open('Adjoint_Energy.data', 'r') as finp:
    found_header = False
    for line in finp:
        temp = line.strip().split(',')
        if not found_header:
            header = temp
            found_header = True
            continue
    
        time.append(float(temp[0]))
        energy.append(float(temp[1]))
finp.close()
t = 0.0

import sys
residual_time = []
residual_error = []
list_errors = []
list_time = []
with open('log.kratos','r') as flog:
	for line in flog:
	    if line.find('CurrentTime =') != -1:
		t = float(line.strip().split()[-1])
	    if line[:6] == 'Error:':
		found_existing = False
		if len(list_time) > 0:
		    if list_time[-1]==t:
		        found_existing = True
		if not found_existing:
		    list_time.append(t)
		    list_errors.append(float(line.strip().split()[-1]))
		else:
		    list_errors[-1] = (float(line.strip().split()[-1]))

        
flog.close()

main_fig, ax = plt.subplots(nrows=2,ncols=1, figsize=[16,24])
fig = plt.subplot(2,1,1, title='Adjoint Energy for Singular Value Decomposition Method')
plt.semilogy(time, energy)
fig.set_ylabel('Adjoint energy')
fig.legend(loc=1)
fig.grid(True)

fig = plt.subplot(2,1,2, title='Approximation Error for Adjoint Fields', sharex = fig)
plt.semilogy(list_time, list_errors)
fig.set_ylabel('Residual')
fig.legend(loc=1)
fig.grid(True)
fig.set_xlabel('Time [s]')


plt.savefig('adjoint_energy.pdf', bbox_inches='tight')
plt.show()
            
            
