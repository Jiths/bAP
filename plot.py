'''
Created on Feb 28, 2014

@author: raphaelholca
'''

import pylab as plt
import matplotlib as mpl
import numpy as np

#plots to display        
RS1_Vm      = []
RS2_Vm      = []
RS3_Vm      = []
FS1_Vm      = []
calcium     = []
raster      = [1]
Vm_TC       = []
TC_RS_init  = []
TC_RS_final = []
TC_FS_init  = []
TC_FS_final = []
RS_RS_init  = []
RS_RS_final = []
RS_FS_init  = []
RS_FS_final = []
ODC         = []
RS1_weight  = []
RS2_weight  = []
RS3_weight  = []
FS1_weight  = []
weightChange= []
extraPlot   = []
        
if RS1_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
    plt.figure()
    spikes = np.argwhere(self.TC.spikeTimes[0,:])
    plt.plot(spikes*self.dt,98*np.ones(len(spikes)), 'vr')
    spikes = np.argwhere(self.TC.spikeTimes[3,:])
    plt.plot(spikes*self.dt,97*np.ones(len(spikes)), 'vb')
    plt.plot(self.time_array, self.valueTracker[7,:], 'r') #synaptic Vm
    plt.plot(self.time_array, self.RS.Vm[self.RS1,:], 'k') #soma Vm
    plt.plot(self.time_array, self.valueTracker[8,:], 'b') #synpatic lateral excitation
    plt.plot(self.time_array, self.valueTracker[9,:], 'r') #synpatic lateral inhibition
    plt.plot(self.time_array, self.valueTracker[10,:]*1000+80, 'b') #potentiation variable
    plt.plot(self.time_array, -self.valueTracker[11,:]*1000+80, 'r') #depression variable
    plt.ylim(-100,100)
    plt.title('RS '+ str(self.RS1))
    plt.savefig('../output/' + 'RS1_Vm.png')

if RS2_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
    plt.figure()
    spikes = np.argwhere(self.TC.spikeTimes[0,:])
    plt.plot(spikes*self.dt,98*np.ones(len(spikes)), 'vr')
    spikes = np.argwhere(self.TC.spikeTimes[3,:])
    plt.plot(spikes*self.dt,97*np.ones(len(spikes)), 'vb')
    plt.plot(self.time_array, self.valueTracker[12,:], 'r') #synaptic Vm
    plt.plot(self.time_array, self.RS.Vm[self.RS2,:], 'k') #soma Vm
    plt.plot(self.time_array, self.valueTracker[13,:], 'b') #synpatic lateral excitation
    plt.plot(self.time_array, self.valueTracker[14,:], 'r') #synpatic lateral inhibition
    plt.plot(self.time_array, self.valueTracker[15,:]*1000+80, 'b') #potentiation variable
    plt.plot(self.time_array, -self.valueTracker[16,:]*1000+80, 'r') #depression variable
    plt.ylim(-100,100)
    plt.title('RS ' + str(self.RS2))
    plt.savefig('../output/' + 'RS2_Vm.png')
    
if RS3_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
    plt.figure()
    spikes = np.argwhere(self.TC.spikeTimes[0,:])
    plt.plot(spikes*self.dt,98*np.ones(len(spikes)), 'vr')
    spikes = np.argwhere(self.TC.spikeTimes[3,:])
    plt.plot(spikes*self.dt,97*np.ones(len(spikes)), 'vb')
    plt.plot(self.time_array, self.valueTracker[17,:], 'r') #synaptic Vm
    plt.plot(self.time_array, self.RS.Vm[self.RS3,:], 'k') #soma Vm
    plt.plot(self.time_array, self.valueTracker[18,:], 'b') #synpatic lateral excitation
    plt.plot(self.time_array, self.valueTracker[19,:], 'r') #synpatic lateral inhibition
    plt.plot(self.time_array, self.valueTracker[20,:]*1000+80, 'b') #potentiation variable
    plt.plot(self.time_array, -self.valueTracker[21,:]*1000+80, 'r') #depression variable
    plt.ylim(-100,100)
    plt.title('RS ' + str(self.RS3))
    plt.savefig('../output/' + 'RS3_Vm.png')

if FS1_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
    plt.figure()
    plt.plot(self.time_array, self.valueTracker[22,:], 'b') #synaptic Vm
    plt.plot(self.time_array, self.FS.Vm[self.FS1,:], 'k') #soma Vm
    plt.plot(self.time_array, self.valueTracker[23,:], 'b') #synpatic lateral excitation
    plt.plot(self.time_array, self.valueTracker[24,:], 'r') #synpatic lateral inhibition
    plt.plot(self.time_array, self.valueTracker[25,:]*1000+80, 'b') #potentiation variable
    plt.plot(self.time_array, -self.valueTracker[26,:]*1000+80, 'r') #depression variable
    plt.ylim(-100,100)
    plt.title('FS ' + str(self.FS1))
    plt.savefig('../output/' + 'FS1_Vm.png')

if calcium:#displays evolution of the calcium detectors
    plt.figure()
    plt.plot(self.time_array, self.calciumTracker[0,:], 'k')
    plt.plot(self.time_array, self.CaDetectorsTracker[0,:]*10+2, 'b') #P
    plt.plot(self.time_array, self.CaDetectorsTracker[1,:]*10+2, 'r') #D
    plt.plot(self.time_array, self.CaDetectorsTracker[2,:], 'c') #B
    plt.plot(self.time_array, self.CaDetectorsTracker[3,:], 'm') #V
    plt.ylim(-0.1,3.3)
    plt.savefig('../output/' + 'calcium.png')

if raster:#display spike raster plot
    plt.figure()
    neuronCounter = 0
    popCounter = 0
    color = ['|r', '|b', '|k']
    for pops in self.population:
        for neurons in range(np.size(pops.spikeTimes,0)):
            spikes = np.argwhere(pops.spikeTimes[neurons,:])
            plt.plot(spikes*self.dt,neuronCounter*np.ones(len(spikes)), color[popCounter])
            plt.hold(True)
            neuronCounter+=1
        popCounter+=1
    plt.xlim(0,self.trialDuration)
    plt.ylim(-1,neuronCounter)
    plt.xlabel('time(ms)')
    plt.ylabel('neuron')
    plt.savefig('../output/' + 'raster.png')

if Vm_TC:#displays Vm for all TC cells
    plt.figure()
    for TC in range(np.size(self.TC.Vm,0)):
        plt.subplot(int(np.floor(np.sqrt(np.size(self.TC.Vm,0)))),int(np.ceil(np.sqrt(np.size(self.TC.Vm,0)))),TC+1)
        plt.plot(self.time_array, self.TC.Vm[TC,:])
        plt.ylim(-90,50)
    plt.savefig('../output/' + 'Vm_TC') 

if TC_RS_init: #display initial TC->RS weights            
    plt.figure()
    rootSize = np.sqrt(np.size(self.TC_RS.g[:,0]))
    for RS in range(self.RS.size):
        plt.subplot(int(np.ceil(np.sqrt(self.RS.size))),int(np.ceil(np.sqrt(self.RS.size))), RS+1)
        plt.imshow(np.reshape(self.TC_RS_initWeight[:,RS],[rootSize, rootSize]), interpolation='nearest', vmin=0, vmax=self.TC_RS.g_max, cmap='bwr')
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)
    plt.suptitle('Initial TC->RS Weights')
    plt.savefig('../output/' + 'TC_RS_init')

if TC_RS_final: #display final TC->RS weights
    plt.figure()
        rootSize = np.sqrt(np.size(self.TC_RS.g[:,0]))
        for RS in range(self.RS.size):
            plt.subplot(int(np.ceil(np.sqrt(self.RS.size))),int(np.ceil(np.sqrt(self.RS.size))), RS+1)
            plt.imshow(np.reshape(self.TC_RS.g[:,RS],[rootSize, rootSize]), interpolation='nearest', vmin=0, vmax=self.TC_RS.g_max, cmap='bwr')
            plt.gca().axes.get_xaxis().set_visible(False)
            plt.gca().axes.get_yaxis().set_visible(False)
            plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        plt.colorbar(cax=cax)
        plt.suptitle('Final TC->RS Weights')
        plt.savefig('../output/' + 'TC-RS_final.png')
    
if TC_FS_final: #display initial TC->FS weights
    plt.figure()
    rootSize = np.sqrt(np.size(self.TC_FS.g[:,0]))
    for FS in range(self.FS.size):
        plt.subplot(int(np.ceil(np.sqrt(self.FS.size))),int(np.ceil(np.sqrt(self.FS.size))), FS+1)
        plt.imshow(np.reshape(self.TC_FS_initWeight[:,FS],[rootSize, rootSize]), interpolation='nearest', vmin=0, vmax=self.TC_FS.g_max, cmap='bwr')
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)
    plt.suptitle('Initial TC->FS Weights')
    plt.savefig('../output/' + 'TC_FS_final')

if TC_FS_final: #display final TC->FS weights
    plt.figure()
    rootSize = np.sqrt(np.size(self.TC_FS.g[:,0]))
    for FS in range(self.FS.size):
        plt.subplot(int(np.ceil(np.sqrt(self.FS.size))),int(np.ceil(np.sqrt(self.FS.size))), FS+1)
        plt.imshow(np.reshape(self.TC_FS.g[:,FS],[rootSize, rootSize]), interpolation='nearest', vmin=0, vmax=self.TC_FS.g_max, cmap='bwr')
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)
    plt.suptitle('Final TC->FS Weights')
    plt.savefig('../output/' + 'TC_FS_final')
    
if RS_RS_init:
    plt.figure()
    for projecting in range(np.size(self.RS_RS.g,0)):
        plt.subplot(int(np.sqrt(np.size(self.RS_RS.g,0))),int(np.sqrt(np.size(self.RS_RS.g,0))),projecting+1)
        plt.imshow(np.reshape(self.RS_RS_initWeight[projecting,:],[np.sqrt(self.RS_RS.size[1]),np.sqrt(self.RS_RS.size[1])]), interpolation='nearest', vmin=0, vmax=self.RS_RS.g_max)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)
    plt.suptitle('Initial RS->RS Weights')
    plt.savefig('../output/' + 'RS_RS_init.png')

if RS_RS_final:
    plt.figure()
    for projecting in range(np.size(self.RS_RS.g,0)):
        plt.subplot(int(np.sqrt(np.size(self.RS_RS.g,0))),int(np.sqrt(np.size(self.RS_RS.g,0))),projecting+1)
        plt.imshow(np.reshape(self.RS_RS.g[projecting,:],[np.sqrt(self.RS_RS.size[1]),np.sqrt(self.RS_RS.size[1])]), interpolation='nearest', vmin=0, vmax=self.RS_RS.g_max)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)
    plt.suptitle('Final RS->RS Weights')
    plt.savefig('../output/' + 'RS_RS_final.png')
   
if RS_FS_init:
    plt.figure()
    for projecting in range(np.size(self.RS_RS.g,0)):
        plt.subplot(int(np.sqrt(np.size(self.RS_FS.g,0))),int(np.sqrt(np.size(self.RS_RS.g,0))),projecting+1)
        plt.imshow(np.reshape(self.RS_RS_initWeight[projecting,:],[np.sqrt(self.RS_RS.size[1]),np.sqrt(self.RS_RS.size[1])]), interpolation='nearest', vmin=0, vmax=self.RS_RS.g_max)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)
    plt.suptitle('Initial RS->FS Weights')
    plt.savefig('../output/' + 'RS_FS_init.png')
 
if RS_FS_final:
    plt.figure()
    for projecting in range(np.size(self.RS_FS.g,0)):
        plt.subplot(int(np.sqrt(np.size(self.RS_FS.g,0))),int(np.sqrt(np.size(self.RS_FS.g,0))),projecting+1)
        plt.imshow(np.reshape(self.RS_FS.g[projecting,:],[np.sqrt(self.RS_FS.size[1]),np.sqrt(self.RS_FS.size[1])]), interpolation='nearest', vmin=0, vmax=self.RS_FS.g_max)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)
    plt.suptitle('Final RS->FS Weights')
    plt.savefig('../output/' + 'RS_FS_final.png')
    
if ODC: #dipslays an ocular dominance map
    plt.figure(figsize=(7,7))
        rootSize = np.sqrt(self.RS.size)
        ODC_mat = np.zeros(self.RS.size)
        alpha_mat = np.zeros(self.RS.size)
        #create white color map with transparency gradient
        cmap_trans = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',['black','black'],256) 
        cmap_trans._init()
        alphas = np.linspace(1.0, 0, cmap_trans.N+3)
        cmap_trans._lut[:,-1] = alphas
        
        for RS in range(self.RS.size):
            prefPattern = [np.sum(self.TC_RS.g[[0,4,8,12],RS]),np.sum(self.TC_RS.g[[1,5,9,13],RS]),np.sum(self.TC_RS.g[[2,6,10,14],RS]),np.sum(self.TC_RS.g[[3,7,11,15],RS])]
            ODC_mat[RS] = np.argmax(prefPattern)
            alpha_mat[RS] = np.max(prefPattern)-(np.sum(prefPattern)-np.max(prefPattern))/self.numPattern
#            ODC_mat[RS] = np.mean(self.TC_RS.g[[0,5,10,15],RS]) - np.mean(self.TC_RS.g[[3,6,9,12],RS])   
        plt.imshow(np.reshape(ODC_mat, [rootSize, rootSize]),interpolation='nearest', cmap='Spectral', vmin=0,vmax=3) #color
        plt.imshow(np.reshape(alpha_mat, [rootSize, rootSize]),interpolation='nearest', cmap=cmap_trans, vmin=-0.25,vmax=1.5) #transparency
        plt.title('Ocular Dominance at ' + np.str(t) + 'ms')
        plt.savefig('../output/' + 'OcularDominance_final.png')

if RS1_weight: #display weight evolution of RS neuron 1
    plt.figure()
    plt.plot(self.time_array, self.weightTracker[0,:], 'b')
#            plt.plot(self.time_array, self.weightTracker[1,:], 'r') 
#            plt.plot(self.time_array, self.weightTracker[2,:], 'k') 
    plt.plot(self.time_array, self.weightTracker[3,:], 'r') 
    plt.ylim(0,self.TC_RS.g_max)
    spikes = np.argwhere(self.TC.spikeTimes[0,:])
    plt.plot(spikes*self.dt,(self.TC_RS.g_max-0.01)*np.ones(len(spikes)), 'vr')
    spikes = np.argwhere(self.TC.spikeTimes[3,:])
    plt.plot(spikes*self.dt,(self.TC_RS.g_max-0.02)*np.ones(len(spikes)), 'vb')
    plt.legend("\\")
    plt.title("RS " + np.str(self.RS1))
    plt.savefig('../output/' + 'RS_1.png')

if RS2_weight: #display weight evolution of RS neuron 2
    plt.figure()
    plt.plot(self.time_array, self.weightTracker[4,:], 'b')
#            plt.plot(self.time_array, self.weightTracker[5,:], 'r') 
#            plt.plot(self.time_array, self.weightTracker[6,:], 'k') 
    plt.plot(self.time_array, self.weightTracker[7,:], 'r') 
    spikes = np.argwhere(self.TC.spikeTimes[0,:])
    plt.plot(spikes*self.dt,(self.TC_RS.g_max-0.01)*np.ones(len(spikes)), 'vr')
    spikes = np.argwhere(self.TC.spikeTimes[3,:])
    plt.plot(spikes*self.dt,(self.TC_RS.g_max-0.02)*np.ones(len(spikes)), 'vb')
    plt.ylim(0,self.TC_RS.g_max)
    plt.legend('\\')
    plt.title("RS " + np.str(self.RS2))
    plt.savefig('../output/' + 'RS_2.png')
    
if RS3_weight: #display weight evolution of RS neuron 3
    plt.figure()
    plt.plot(self.time_array, self.weightTracker[8,:],  'b')
#            plt.plot(self.time_array, self.weightTracker[9,:],  'r') 
#            plt.plot(self.time_array, self.weightTracker[10,:], 'k') 
    plt.plot(self.time_array, self.weightTracker[11,:], 'r') 
    spikes = np.argwhere(self.TC.spikeTimes[0,:])
    plt.plot(spikes*self.dt,(self.TC_RS.g_max-0.01)*np.ones(len(spikes)), 'vr')
    spikes = np.argwhere(self.TC.spikeTimes[3,:])
    plt.plot(spikes*self.dt,(self.TC_RS.g_max-0.02)*np.ones(len(spikes)), 'vb')
    plt.ylim(0,self.TC_RS.g_max)
    plt.legend('\\')
    plt.title("RS " + np.str(self.RS3))
    plt.savefig('../output/' + 'RS_3.png')
    
if FS1_weight: #display weight evolution of FS neuron 1
    plt.figure()
    plt.plot(self.time_array, self.weightTracker[12,:],  'b')
#            plt.plot(self.time_array, self.weightTracker[13,:],  'r') 
#            plt.plot(self.time_array, self.weightTracker[14,:], 'k') 
    plt.plot(self.time_array, self.weightTracker[15,:], 'r') 
    plt.ylim(0,self.TC_RS.g_max)
    plt.legend('\\')
    plt.title("FS " + np.str(self.FS1))
    plt.savefig('../output/' + 'FS_1.png')
    
if weightChange: #displays the percent of total weight change every at every 100ms
    plt.figure()
    x=np.linspace(0,self.trialDuration,int((self.trialDuration/500)/self.dt)+1)
    plt.plot(x, self.percentChange, '-ob')
    plt.title('percent weight change')
    plt.savefig('../output/' + 'percent_weight_change.png')
    
if extraPlot: #displays an extra, custom plot
    plt.figure()
    plt.plot(self.time_array, self.valueTracker[0], 'k')
#            plt.plot(self.time_array, self.valueTracker[1]+10, 'b')
    plt.plot(self.time_array, self.valueTracker[2], 'b')
#            plt.plot(self.time_array, self.valueTracker[3]+30, 'b')
#            plt.plot(self.time_array, self.valueTracker[4]+40, 'b')
#            plt.plot(self.time_array, self.valueTracker[5]+50, 'b')
#            plt.plot(self.time_array, self.valueTracker[6]+60, 'b')
    plt.ylim(-100,50)
    plt.savefig('../output/' + 'extraPlot.png')