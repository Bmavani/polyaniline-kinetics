# Polyaniline Polymerization Kinetics Code Sample
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from initial_solverfinal import initial_solverfinal

# Constants
ANI0_list = [0.02, 0.035, 0.05, 0.1, 0.2]
k_base_1  = [0.0008, 0.0010, 0.0011, 0.008, 0.01]
k_base_2  = [5.2,    4.6,    4.7,    5.2,   2.4 ]
n_cases   = len(ANI0_list)
t_span    = (0, 180)
t_eval    = np.linspace(0, 180, 90)
labels    = ['0.02M', '0.035M', '0.05M', '0.1M', '0.2M']

AN   = np.zeros((len(t_eval), n_cases))
X    = np.zeros((len(t_eval), n_cases))
APS  = np.zeros((len(t_eval), n_cases))
APSR = np.zeros((len(t_eval), n_cases))   
Mn   = np.zeros((len(t_eval), n_cases))

# Main loop
i = 0
for j in range(n_cases):
    ANI0 = ANI0_list[j]
    APS0 = ANI0_list[j]
    i += 1
    if i - 1 == j:              
        k1        = k_base_1[i - 1]
        k_prime_2 = k_base_2[i - 1]

        sol = solve_ivp(
            initial_solverfinal,
            t_span,
            [ANI0],             
            t_eval=t_eval,
            method='RK45',      
            args=(k1, k_prime_2, APS0, ANI0),   
            dense_output=False,
        )

        ANI2 = sol.y[0]         

        AN[:, i - 1]   = ANI2 / ANI0
        X[:, i - 1]    = (ANI0 - ANI2) / ANI0
        APS[:, i - 1]  = APS0 - ((ANI0 - ANI2) / 1.25)
        APSR[:, i - 1] = APS[:, i - 1] / APS0   
        Mn[:, i - 1]   = ((ANI0 - ANI2) * 93.13 * 180) / APS0

# Figures 
plt.figure(1)
plt.plot(t_eval, X)
plt.xlabel('time (min)')
plt.ylabel('Conversion (X)')
plt.title('Initial concentration effect')
plt.legend(labels)

plt.figure(2)
plt.plot(t_eval, AN)
plt.xlabel('time (min)')
plt.ylabel('[NA]/[NA0]')
plt.title('Initial concentration effect')
plt.legend(labels)

plt.figure(3)
plt.plot(t_eval, APS)
plt.xlabel('time (min)')
plt.ylabel('APS concentration (M or mol/L)')
plt.title('APS concentration vs time')
plt.legend(labels)

plt.figure(4)
plt.plot(X, Mn)
plt.xlabel('Conversion (X)')
plt.ylabel('Number average molecular weight (cumulative) (g/mol)')
plt.title('Chain length')
plt.legend(labels)

plt.figure(5)
plt.plot(t_eval, APS)
plt.xlabel('time (min)')
plt.ylabel('Initiator concentration (M)')
plt.legend(labels)

plt.figure(6)
plt.plot(ANI0_list, Mn[-1, :], '-o')
plt.xlabel('Initial ANI concentration (M)')
plt.ylabel('Number average molecular weight (cumulative) (g/mol)')

plt.figure(7)
plt.plot(X, APS)
plt.xlabel('Conversion (x)')
plt.ylabel('Initiator concentration (M)')
plt.legend(labels)

plt.figure(8)
plt.plot(t_eval, APSR)  
plt.xlabel('time (min)')
plt.ylabel('Initiator ratio')

plt.show()
