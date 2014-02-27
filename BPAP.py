'''
Created on Feb 14, 2013

@author: raphaelholca
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

def main(argv=None):
    dt = 0.1
    duration = 200
    time = np.linspace(0,duration,duration/dt+1)
    
    G_L = 0.003
    E_L = -54.387#-17.0
    G_Na = 1.2
    E_Na = 50.0#55.0 
    G_K = 0.36#0.2
    E_K = -77.0#-72.0
    C_m = 0.1
    I_e = 0.0
    m = np.zeros(np.size(time))
    h = np.zeros(np.size(time))
    n = np.zeros(np.size(time))
    I = np.zeros(np.size(time))
    V = np.zeros(np.size(time))
    V[-1] = -65.0
    m[-1] = 0.053
    h[-1] = 0.596
    n[-1] = 0.3184
    
    for t in range(np.size(time)):
        
        if t*dt>=50: I_e= 1.0
        if t*dt>=52: I_e= 0.0
        
        m[t] = m[t-1] + M(V[t-1], m[t-1])*dt
        h[t] = h[t-1] + H(V[t-1], h[t-1])*dt
        n[t] = n[t-1] + N(V[t-1], n[t-1])*dt
        I[t] = G_L*(V[t-1]-E_L) + G_Na*pow(m[t],3)*h[t]*(V[t-1]-E_Na) + G_K*pow(n[t],4)*(V[t-1]-E_K)
        V[t] = V[t-1] + ((-I[t]+I_e)/C_m)*dt
        
    plt.figure()
    plt.plot(time,V)
    plt.xlim(45,60)
    plt.ylim(-100, 100)
    plt.figure()
    plt.plot(time, m, 'm')
    plt.plot(time, h, 'r')
    plt.plot(time, n, 'b')
    plt.xlim(45,60)
    plt.ylim(0, 1)
    plt.show()  

def M(V, m):
    alpha = 0.1*(V+40.0)/(1-np.exp(-0.1*(V+40.0)))#0.38*(V+29.7)/(1-np.exp(-0.1*(V+29.7)))
    beta = 4.0*np.exp(-0.0556*(V+65.0))#15.2*np.exp(-0.0556*(V+54.7))
    return alpha*(1-m) - beta*m

def H(V, h):
    alpha = 0.07*np.exp(-0.05*(V+65.0))#0.266*np.exp(-0.05*(V+48.0))
    beta = 1.0/(1.0+np.exp(-0.1*(V+35.0)))#3.8/(1.0+np.exp(-0.1*(V+18)))
    return alpha*(1-h) - beta*h

def N(V, n):
    alpha = 0.01*(V+55.0)/(1.0-np.exp(-0.1*(V+55.0)))#0.02*(V+45.7)/(1-np.exp(-0.1*(V+45.7)))
    beta = 0.125*np.exp(-0.0125*(V+65.0))#0.25*np.exp(-0.0125*(V+55.7))
    return alpha*(1-n) - beta*n

if __name__ == "__main__":
    sys.exit(main())

