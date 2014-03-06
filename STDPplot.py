'''
Created on Jan 31, 2013

@author: raphaelholca
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import pickle

from CompSim import CompSim

def main(argv=None):
        
        tic = time.time()
        
        global sim
        print "plotting STDP curve...\r\r" 
        sim = CompSim(duration=100, dt=0.1)
        sim.setParam()
        sim.setParamSTDPplot()
        sim.createSynapses()
        sim.createPopulations()
#        isi = np.linspace(-25,25,101) #interval between the pre- and post-synaptic spikes (post-pre)
        isi = np.array([-5])
        deltaW = np.zeros(np.shape(isi))
        P = np.zeros(np.shape(isi))
        D = np.zeros(np.shape(isi))
        for i in range(np.size(isi)):
            sim.resetValueSTDPplot()
            deltaW[i] = sim.runSimulation(isi=isi[i])
            P[i] = np.sum(sim.STDPtracker[3,:])
            D[i] = np.sum(sim.STDPtracker[4,:])
            
        if True: #displays Vm, etc.
            plt.figure()
            plt.plot(sim.time_array, sim.STDPtracker[0,:], 'k') #soma Vm
            plt.plot(sim.time_array, sim.STDPtracker[1,:], 'b') #spine Vm
            plt.plot(sim.time_array, sim.STDPtracker[2,:], 'k') #spine calcium
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[3,:]*50-20, 'b') #P
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[4,:]*100, 'r') #D
            plt.plot(sim.time_array, sim.STDPtracker[5,:]*1000-40, 'k') #delta g
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[6,:], 'r') #FS Vm
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[7,:], 'r') #inhibibitory conductance
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[8,:], 'b') #excitatory conductance
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[9,:]*10, 'b') #OP NMDAR due to Glu
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[10,:]*10, 'r') #OP NMDAR due to bAP
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[11,:]*10, 'k') #OP VGCC
#            plt.Plot_Print(sim.time_array, sim.STDPtracker[12,:], 'b') #B
#            plt.ylim(-1,10)

        print "delta W =", sim.STDPtracker[5,-1]
                     
        if False: #displays the STDP curve
            plt.figure()
            plt.plot(isi, deltaW, 'ok')
#            plt.Plot_Print(isi, P/2000, 'ob-')
#            plt.Plot_Print(isi, D/2000, 'or-')
            plt.xlabel('interspike interval (post-pre)')
            plt.ylabel('change in synaptic weight')
            plt.axhline(0)
            plt.axvline(0)
        
            ltp=np.sum(deltaW[deltaW>=0])
            ltd=-np.sum(deltaW[deltaW<0])
            print "LTP/LTD", ltp/ltd, "\r\r"
            
        #save variables to files
        if False:
            pickleFile = open("calcium_low", 'w')
            pickle.dump(sim.STDPtracker[2,:], pickleFile)
            pickleFile.close()

#            pickleFile = open("B", 'w')
#            pickle.dump(sim.STDPtracker[12,:], pickleFile)
#            pickleFile.close()
#            pickleFile = open("bAP_NMDA", 'w')
#            pickle.dump(sim.STDPtracker[10,:], pickleFile)
#            pickleFile.close()
#            pickleFile = open("calcium", 'w')
#            pickle.dump(sim.STDPtracker[2,:], pickleFile)
#            pickleFile.close()
#            pickleFile = open("D", 'w')
#            pickle.dump(sim.STDPtracker[4,:], pickleFile)
#            pickleFile.close()
#            pickleFile = open("OP_NMDA", 'w')
#            pickle.dump(sim.STDPtracker[9,:], pickleFile)
#            pickleFile.close()
#            pickleFile = open("OP_VGCC", 'w')
#            pickle.dump(sim.STDPtracker[11,:], pickleFile)
#            pickleFile.close()
#            pickleFile = open("time", 'w')
#            pickle.dump(sim.time_array, pickleFile)
#            pickleFile.close()
#            pickleFile = open("Vm_spine", 'w')
#            pickle.dump(sim.STDPtracker[1,:], pickleFile)
#            pickleFile.close()
#            pickleFile = open("Vm_soma", 'w')
#            pickle.dump(sim.STDPtracker[0,:], pickleFile)
#            pickleFile.close()
        
        print "STDP curve plotting time:", time.time()-tic
        plt.show()
        
if __name__ == "__main__":
    sys.exit(main())