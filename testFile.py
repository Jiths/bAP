'''
Created on Feb 13, 2013

@author: raphaelholca
'''
import sys
import numpy as np
import matplotlib.pyplot as plt

def main(argv=None):
    
    stepDur = 6
    timeDur = 300
    dt = 0.01
    theta = 20.
    NMDA = True
    if NMDA:
        tau_rise = 2.
        tau_fast = 10.
        a_fast = 0.527
        tau_slow = 45.
        a_slow = 0.473
        tau_f = 25. #original: 50.
        tau_s = 55. #original: 200.
    else:
        tau_rise = 0.58
        tau_fast = 7.6
        a_fast = 0.903
        tau_slow = 25.69
        a_slow = 0.097
        tau_f = 5.
        tau_s = 25. #original: 50.
    t1=20
    t2=100
    time = np.linspace(0,timeDur,timeDur/dt+1)
#    step = np.append(np.zeros(t1/dt), np.append(np.ones(stepDur/dt),np.zeros(timeDur/dt-stepDur/dt+1-t1/dt)))
    step = np.append(np.zeros(t1/dt), np.append(np.ones(stepDur/dt), np.append(np.zeros((t2-stepDur-t1)/dt), np.append(np.ones(stepDur/dt),np.zeros((timeDur-t2-stepDur)/dt+1)))))
    S = np.zeros(timeDur/dt+1)
    Ca_hard = np.zeros(timeDur/dt+1)
    Ca_easy = np.zeros(timeDur/dt+1)
    I_NMDA = np.zeros(timeDur/dt+1)
    S_rise = 0.
    S_fast = 0.
    S_slow = 0.
    for t in range(int(timeDur/dt)):
        S_rise += (-theta*(1.-S_fast-S_slow)*step[t]-S_rise/tau_rise)*dt
        S_fast += (theta*(a_fast-S_fast)*step[t]-S_fast/tau_fast)*dt
        S_slow += (theta*(a_slow-S_slow)*step[t]-S_slow/tau_slow)*dt
        S[t] = S_rise+S_fast+S_slow 
        Ca_hard[t] = Ca_hard[t-1] + S[t]*25*0.01-0.083*(Ca_hard[t-1]-0.07)-(0.083/6)*pow(Ca_hard[t-1],2)
        
        if t1/dt<=t: I_NMDA[t]  = 2.2e4*0.5*0.002*(0.5*np.exp(-(t*dt-t1)/tau_f)+0.5*np.exp(-(t*dt-t1)/tau_s))
        if t2/dt<=t: I_NMDA[t] += 2.2e4*0.5*0.002*(0.5*np.exp(-(t*dt-t2)/tau_f)+0.5*np.exp(-(t*dt-t2)/tau_s))
        
        Ca_easy[t] = Ca_easy[t-1] + (0.1*I_NMDA[t]-Ca_easy[t-1])*dt
    
    plt.plot(time,S*25)
    plt.plot(time,Ca_hard, 'r')
    plt.plot(time,Ca_easy, '--r')
#    plt.plot(time,step)
    plt.plot(time,I_NMDA)
    plt.show()
    
    
if __name__ == "__main__":
    sys.exit(main())
