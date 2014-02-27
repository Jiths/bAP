'''
Created on Nov 13, 2012

@author: raphaelholca
'''

import scipy.io
import pylab as plt
import matplotlib as mpl
import time
import numpy as np
from scipy import fftpack
from Population import Population
from Synapses import Synapses
import pickle


class CompSim():
    '''
    This is the main simulation function.
    '''
    def __init__(self, duration=200., dt=0.2):
        
        print 'simulation started:', time.ctime()

        self.tic_prep = time.time()
        
        """ general simulation parameters """
        self.trialDuration = duration #time to simulate (msec)
        self.dt = dt #length of each simulation time step (msec)
        self.time_array = np.linspace(0,self.trialDuration,self.trialDuration/self.dt+1) #time array
        self.naturalImage = [] #whether to use a natural image as an input - leave empty to use artificial stimulus
        self.imageSize = 4 #pixel size of the SIDE of the input image (total size = imageSize^2)
        self.popSize = 5 #side of the size of the population of neurons (total size = popSize^2)
        self.numPattern = 3 #number of patterns embedded in the input
        self.FS_RS = [1] #turns off lateral inhibition to excitatory cells - leave empty for OFF, put 1 for ON
        self.RS_RS = [1] #turns off lateral excitation to excitatory cells - leave empty for OFF, put 1 for ON
        self.RS_FS = [1] #turns off lateral excitation to inhibitory cells - leave empty for OFF, put 1 for ON
        self.FS_FS = [1] #turns off lateral inhibition to inhibitory cells - leave empty for OFF, put 1 for ON
        self.gaussian = [] #whether 'intracortical' connections are distributed according to a Gaussian distribution. If not, follows a square wave
#        np.random.seed(83251) #seeds the random number generator
        
        #neurons whose stats should be displayed
        self.RS1=0
        self.RS2=5
        self.RS3=3
        self.FS1=0
        
        #plots to display        
        self.RS1_Vm      = []
        self.RS2_Vm      = []
        self.RS3_Vm      = []
        self.FS1_Vm      = []
        self.calcium     = []
        self.raster      = []
        self.Vm_TC       = []
        self.TC_RS_init  = []
        self.TC_RS_evol  = [1]
        self.TC_RS_final = []
        self.TC_FS_init  = []
        self.TC_FS_final = []
        self.RS_RS_init  = []
        self.RS_RS_final = []
        self.RS_FS_init  = []
        self.RS_FS_final = []
        self.ODC         = []
        self.ODC_evol    = [1]
        self.RS1_weight  = []
        self.RS2_weight  = []
        self.RS3_weight  = []
        self.FS1_weight  = []
        self.weightChange= []
        self.extraPlot   = []
      
    def setParam(self):
        
        """ value trackers and others """
        
        self.weightTracker = np.zeros([20,self.trialDuration/self.dt+1]) #tracks value of weights over time
        self.OPTracker=np.zeros([3,self.trialDuration/self.dt+1]) #array space to record OP
        self.calciumTracker = np.zeros([3,self.trialDuration/self.dt+1]) #tracks value of the calcium transient over time
        self.valueTracker = np.zeros([30,self.trialDuration/self.dt+1]) #tracks other values...
        self.CaDetectorsTracker = np.zeros([5,self.trialDuration/self.dt+1]) #tracks value of the calcium detector variables
        self.STDPtracker = np.zeros([15,self.trialDuration/self.dt+1]) #tracks value for STDP plots
        self.percentChange = np.zeros(int(self.trialDuration/self.dt/500)+1) #tracks how much weight changed during the last time step
        self.STDPplot = False #boolean used to plot the STDP curve - set to True by setParamSTDP()
        self.stimPresentation = 500.0 #duration of presentation of a particular image patch 
        
        """ general neuron parameters """ 
        self.E_leak = -78. #rest membrane potential (mV)
        self.E_e = +0. #excitatory synaptic reversal potential (mV)
        self.E_i = -80. #inhibitory synaptic reversal potential (mV)
        self.E_Ca = +130 #Ca reversal potential (mV) 
        self.g_NMDA = 0.005 #conductance of the NMDA channel
        self.g_leak = 1.0 #conductance of the leakage current
        self.Vth = -50 #membrane threshold (mV)
        self.Vspike = +40. #max voltage of an action potential (mV)
        self.GrefMax = 1.0/self.dt #refractory period max conductance (constant)
        self.Vreset = self.E_leak #reset voltage after an action potential (mV)
        self.Crest = 0.07 #resting level of [Ca2+] (uM) 
        self.tau_ampaF = 5. #time constant of the fast AMPA current
        self.tau_ampaS = 10 #original: 50. #time constant of the slow AMPA current
        self.tau_nmdaF = 7.#25. #original: 50. #time constant of the fast NMDA current
        self.tau_nmdaS = 12.#55. #original: 200. #time constant of the slow NMDA current
        self.tau_bpapF = 1.#3. #time constant of the fast BPAP current
        self.tau_bpapS = 25. #time constant of the slow BPAP current
        self.tau_Ca = 30.#90. #time constant of the calcium concentration
        self.Vmax_BPAP = 63#67. #maximum depolarization due to a BPAP
        self.lr_LTP = 0.004 #learning rate of the synapse potentiation
        self.lr_LTD = 0.004 #learning rate of the synapse depression
        
        """ TC-specific parameters """
        self.TC_size = pow(self.imageSize,2) #number of TC neurons
        self.TC_Rm = 1. #specific membrane resistance (MOhm)
        self.TC_Cm = 30. #specific membrane capacitance (uF)
        self.TC_Ie = 40. + np.zeros(self.TC_size) #injected current
        self.TC_tau_Gref = 10.0#15. #time constant for refractory period current
        
        """ RS-specific parameters """
        self.RS_size = pow(self.popSize,2) #number of RS neurons
        self.RS_Rm = 3.90 #specific membrane resistance (MOhm)
        self.RS_Cm = 0.714 #specific membrane capacitance (uF)
        self.RS_Ie = 0. + np.zeros(self.RS_size) #injected current
        self.RS_tau_Gref = 5.0#15. #time constant for refractory period current
        
        """ FS-specific parameters """
        self.FS_size = pow(self.popSize,2) #number of FS neurons
        self.FS_Rm = 0.789 #specific membrane resistance (MOhm)
        self.FS_Cm = 1.141 #specific membrane capacitance (uF)
        self.FS_Ie = 0. + np.zeros(self.FS_size) #injected current
        self.FS_tau_Gref = 4.0 #time constant for refractory period current
        
        """ synapse parameters """
        
        """ TC->RS synapses parameters """
        self.TC_RS_size = [self.TC_size, self.RS_size]
        self.TC_RS_G = 0.13 #maximal unitary conductance of the thalamocortical synapses on RS cells
        self.TC_RS_weightRand = 0.03 #variability in weight initialization
        self.TC_RS_rise = 5 #rise time constant of the excitatory synaptic conductance
        self.TC_RS_decay = 0.5 #decay time constant of the excitatory synaptic conductance
        self.TC_RS_E = self.E_e #synaptic reversial potential
        self.TC_RS_G_max = 0.35 #maximal synaptic conductance - weights cannot grow stronger than this value
        
        """ TC->FS synapses parameters """
        self.TC_FS_size = [self.TC_size, self.FS_size]
        self.TC_FS_G = 0.0#0.1 #maximal unitary conductance of the thalamocortical synapses on FS cells
        self.TC_FS_weightRand = 0.0#0.02 #variability in weight initialization
        self.TC_FS_rise = 10.0 #rise time constant of the excitatory synaptic conductance
        self.TC_FS_decay = 0.5 #decay time constant of the excitatory synaptic conductance
        self.TC_FS_E = self.E_e #synaptic reversial potential
        self.TC_FS_G_max = 0.0#0.2 #maximal synaptic conductance - weights cannot grow stronger than this value
        
        """ RS->FS synapses parameters """
        self.RS_FS_size = [self.RS_size, self.FS_size]
        self.RS_FS_G = 1.0 #maximal unitary conductance, when OP=1
        if not self.RS_FS: self.RS_FS_G = 0.0
        self.RS_FS_weightRand = 0.0 #MUST BE 0 #variability in weight initialization
        self.RS_FS_rise = 10.0 #rise time constant of the excitatory synaptic conductance
        self.RS_FS_decay = 2.2 #decay time constant of the excitatory synaptic conductance
        self.RS_FS_E = self.E_e #synaptic reversial potential
        self.RS_FS_G_max = 1.0 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.RS_FS_extent = 1.0#3.0 #extent of lateral connections
        
        """ FS->RS synapses parameters """
        self.FS_RS_size = [self.FS_size, self.RS_size]
        self.FS_RS_G = 0.006 #maximal unitary conductance of the inhibitory synapses on RS cells
        if not self.FS_RS: self.FS_RS_G = 0. #turns inhibition off
        self.FS_RS_weightRand = 0 #MUST BE 0 #variability in weight initialization
        self.FS_RS_rise = 5.5 #rise time constant of the inhibitory synaptic conductance
        self.FS_RS_decay = 0.5 #decay time constant of the inhibitory synaptic conductance
        self.FS_RS_E = self.E_i #synaptic reversial potential
        self.FS_RS_G_max = 0.013#0.13 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.FS_RS_extent = 5 #extent of lateral connections
        
        """ RS->RS synapses parameters """
        self.RS_RS_size = [self.RS_size, self.RS_size]
        self.RS_RS_G = 0.06#0.12 #initial maximal unitary conductance
        if not self.RS_RS: self.RS_RS_G = 0.
        self.RS_RS_weightRand = 0. #MUST BE 0 #variability in weight initialization
        self.RS_RS_rise = 5 #rise time constant of the excitatory synaptic conductance
        self.RS_RS_decay = 0.5 #decay time constant of the excitatory synaptic conductance
        self.RS_RS_E = self.E_e #synaptic reversial potential
        self.RS_RS_G_max = 0.12 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.RS_RS_extent = 3#1.0 #extent of lateral connections
        
        """ FS->FS synapses parameters """
        self.FS_FS_size = [self.FS_size, self.FS_size]
        self.FS_FS_G = 0.0#0.2 #maximal unitary conductance of the inhibitory synapses on RS cells
        if not self.FS_FS: self.FS_FS_G = 0. #turns inhibition off
        self.FS_FS_weightRand = 0 #MUST BE 0 #variability in weight initialization
        self.FS_FS_rise = 5.5 #rise time constant of the inhibitory synaptic conductance
        self.FS_FS_decay = 0.5 #decay time constant of the inhibitory synaptic conductance
        self.FS_FS_E = self.E_i #synaptic reversial potential
        self.FS_FS_G_max = 0.2 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.FS_FS_extent = 1.#0.8 #extent of lateral connections
        
        if not self.STDPplot: print "duration:", int(self.trialDuration), "ms\nTC:", self.TC_size, "\nRS:", self.RS_size, "\nFS:", self.FS_size
    
    def setParamSTDPplot(self):
        """ sets the  parameter values for the STDP plot """
               
        self.RS_RS_G = 0.0
        self.RS_FS_G = 0.0
        self.FS_RS_G = 0.03
        self.RS_RS_G_max = 0.0
        self.RS_FS_G_max = 0.0
        self.FS_RS_G_max = 0.03
        self.STDPplot_initialG_TC_RS = self.TC_RS_G
        self.STDPplot_initialG_TC_FS = self.TC_FS_G
        self.STDPplot_initialG_FS_RS = self.FS_RS_G
        self.TC_size = 1
        self.RS_size = 1
        self.FS_size = 1
        self.TC_RS_size = [1,1]
        self.TC_FS_size = [1,1]
        self.RS_FS_size = [1,1]
        self.RS_RS_size = [1,1]
        self.FS_RS_size = [1,1]
        self.FS_FS_size = [1,1]
        self.TC_RS_weightRand = 0.
        self.TC_FS_weightRand = 0.
        self.RS_FS_weightRand = 0.
        self.FS_RS_weightRand = 0.
        self.TC_Ie = np.zeros(self.FS_size)
        self.RS_Ie = np.zeros(self.FS_size)
        self.FS_Ie = np.zeros(self.FS_size)
        self.STDPplot = True
        
    def resetValueSTDPplot(self):
        """ resets variables between each iteration of the STDP plot function """        
        
        for neuron in [self.TC, self.RS, self.FS]: 
            neuron.Vm = np.zeros([1,self.trialDuration/self.dt+1])
            neuron.Vm[:,-1] = self.E_leak
            neuron.lastCellSpike = np.ones(1)*-1e99
            neuron.spikeTimes = np.zeros([1,self.trialDuration/self.dt+1])
            neuron.OP = np.zeros([1,1])
            neuron.OP_NMDA = np.zeros([1,1])
            neuron.BPAP = np.zeros(1)
        
        for synapses in [self.TC_RS, self.TC_FS, self.FS_RS]:
            synapses.OP = np.zeros([1,1])
            synapses.OP_NMDA = np.zeros([1,1])
            synapses.lastSynsSpike = np.ones([1,1])*-1e99
            synapses.risingBool[:,:] = False
            synapses.I_NMDA = np.zeros([1,1]) #Ca2+ current through NMDAr
            synapses.I_BPAP = np.zeros([1,1]) #current due to the back-propagating action-potential
            synapses.I_Ca = np.zeros([1,1]) #Ca2+ current through NMDAr
            synapses.I_VGCC = np.zeros([1,1]) #Ca2+ current through voltage-gated calcium channels
            synapses.I = np.zeros([1,1]) #current through AMPAr and GABAr
            synapses.Vm = np.ones([1,1])*self.E_leak #local membrane potential at the synapse
            synapses.Mg = np.zeros([1,1]) #extent of the Mg blockade of the NMDAr
            synapses.lastSynsSpike = np.ones([1,1])*-1e99 #time of last pre-synaptic spike
            synapses.calcium = np.zeros([1,1]) #internal "calcium" concentration, used to compute LTD       
            synapses.P = np.zeros([1,1])
            synapses.V = np.zeros([1,1])
            synapses.B = np.zeros([1,1])
            synapses.D = np.zeros([1,1])
            synapses.W = np.zeros([1,1])
        self.TC_RS.g = np.ones([1,1])*self.STDPplot_initialG_TC_RS
        self.TC_FS.g = np.ones([1,1])*self.STDPplot_initialG_TC_FS
        self.FS_RS.g = np.ones([1,1])*self.STDPplot_initialG_FS_RS
        
    def createSynapses(self):
        """ create synapses """
        
        self.TC_RS = Synapses(self, 'TC', 'RS', {}, {}, self.TC_RS_size, self.TC_RS_G, self.TC_RS_weightRand, self.TC_RS_G_max, self.TC_RS_rise, self.TC_RS_decay, self.TC_RS_E)
        self.TC_FS = Synapses(self, 'TC', 'FS', {}, {}, self.TC_FS_size, self.TC_FS_G, self.TC_FS_weightRand, self.TC_FS_G_max, self.TC_FS_rise, self.TC_FS_decay, self.TC_FS_E)
        self.FS_RS = Synapses(self, 'FS', 'RS', {}, {}, self.FS_RS_size, self.FS_RS_G, self.FS_RS_weightRand, self.FS_RS_G_max, self.FS_RS_rise, self.FS_RS_decay, self.FS_RS_E, self.FS_RS_extent)
        self.FS_FS = Synapses(self, 'FS', 'FS', {}, {}, self.FS_FS_size, self.FS_FS_G, self.FS_FS_weightRand, self.FS_FS_G_max, self.FS_FS_rise, self.FS_FS_decay, self.FS_FS_E, self.FS_FS_extent)
        self.RS_FS = Synapses(self, 'RS', 'FS', {}, {}, self.RS_FS_size, self.RS_FS_G, self.RS_FS_weightRand, self.RS_FS_G_max, self.RS_FS_rise, self.RS_FS_decay, self.RS_FS_E, self.RS_FS_extent)
        self.RS_RS = Synapses(self, 'RS', 'RS', {}, {}, self.RS_RS_size, self.RS_RS_G, self.RS_RS_weightRand, self.RS_RS_G_max, self.RS_RS_rise, self.RS_RS_decay, self.RS_RS_E, self.RS_RS_extent)
        self.synapses = {self.TC_RS, self.TC_FS, self.FS_RS, self.FS_FS, self.RS_FS, self.RS_RS}
        
    def createPopulations(self):
        """ create neuron populations and connections between neurons """
        #all neuron matrix are projectsFrom x projectsTo 
        
        """ neural populations """
        self.TC = Population(self, 'TC', self.TC_size, self.TC_Rm, self.TC_Cm, self.TC_Ie, self.TC_tau_Gref)
        self.RS = Population(self, 'RS', self.RS_size, self.RS_Rm, self.RS_Cm, self.RS_Ie, self.RS_tau_Gref)
        self.FS = Population(self, 'FS', self.FS_size, self.FS_Rm, self.FS_Cm, self.FS_Ie, self.FS_tau_Gref)
        
        """ receiving synapses """
        self.TC.receivesFrom = {} #synapse population TC is post-synaptic to
        self.RS.receivesFrom = {self.TC_RS, self.FS_RS, self.RS_RS} #synapse population RS is post-synaptic to
        self.FS.receivesFrom = {self.TC_FS, self.RS_FS} #synapse population FS is post-synaptic to
        
        """ projecting synapses """
        self.TC.projectsTo = {self.TC_RS, self.TC_FS} #synapse population TC is pre-synaptic to
        self.RS.projectsTo = {self.RS_FS, self.RS_RS} #synapse population RS is pre-synaptic to
        self.FS.projectsTo = {self.FS_RS} #synapse population FS is pre-synaptic to
        
        """ pointer to pre- and post-synaptic population """
        self.TC_RS.pre = self.TC
        self.TC_RS.post = self.RS
        self.TC_FS.pre = self.TC
        self.TC_FS.post = self.FS
        self.FS_RS.pre = self.FS
        self.FS_RS.post = self.RS
        self.RS_FS.pre = self.RS
        self.RS_FS.post = self.FS
        self.RS_RS.pre = self.RS
        self.RS_RS.post = self.RS
        
        self.population = {self.TC,self.RS,self.FS}
    
    def createConnectionPattern(self):
        """ assign connection weights to connections between neurons """
        #connections can either follow a Guassian decay with distance or a 'square' decay (connections equal to a given weight up to a distance and then are set to zero)
        
        for syns in self.synapses:
            if syns.extent!=None:

                weights = np.zeros([np.sqrt(syns.size[1]),np.sqrt(syns.size[1])]) #initialize a square matrix the size of the receiving population 
                
                if self.gaussian: #create a 3D Gaussian function, centered at the center of the matrix. The width of the function depends on syns.extent
                    for y in range(np.size(weights,0)):
                        for x in range(np.size(weights,1)):
                            weights[y,x]=np.exp(-pow((x-np.int(np.floor(np.size(weights,0)/2)))/syns.extent,2))*np.exp(-pow((y-np.int(np.floor(np.size(weights,0)/2)))/syns.extent,2))
                else: #creates a square wave form, with ones at the center of the matrix and zeros outside. This replaces the guassian decay
                    low=int(np.sqrt(syns.size[1])/2)-int(syns.extent/2)
                    high=int(np.sqrt(syns.size[1])/2)+int(syns.extent/2)+1
                    weights[low:high,low:high]=np.ones([syns.extent,syns.extent])
                
                #now roll the matrix in the correct position for each projecting neuron, and then reshape to 1D array
                if syns.preName=='FS' and syns.postName=='RS': #FS->RS
                    FSroot=int(np.sqrt(syns.size[0]))
                    RSroot=int(np.sqrt(self.RS.size))
                    for projecting in range(syns.size[0]):
                        rolledWeights = np.roll(weights, -int(np.floor(RSroot/2))+int(np.floor(RSroot/FSroot))*(1+np.int(projecting/FSroot))-1, 0)
                        rolledWeights = np.roll(rolledWeights, -int(np.floor(RSroot/2))+int(np.floor(RSroot/FSroot))*(1+np.int(np.mod(projecting,FSroot)))-1, 1)
                        syns.g_max_matrix[projecting,:] = syns.g_max*np.reshape(rolledWeights,np.size(syns.g[projecting,:]))
                        syns.g[projecting,:] *= np.reshape(rolledWeights,np.size(syns.g[projecting,:]))
                elif syns.preName=='RS' and syns.postName=='FS': #RS->FS
                    syns.g_max_matrix=np.zeros(syns.size)
                    syns.g=np.zeros(syns.size)
                    np.fill_diagonal(syns.g, syns.g_max)
                    np.fill_diagonal(syns.g_max_matrix, syns.g_max)
                    if self.gaussian:
                        FSroot=int(np.sqrt(self.FS.size))
                        RSroot=int(np.sqrt(self.RS.size))
                        FSpositionX = np.zeros(self.FS.size)
                        FSpositionY = np.zeros(self.FS.size)
                        for i in range(self.FS.size): #finds the position of the FS cells
                            FSpositionX[i]=int(np.floor(RSroot/FSroot))*(1+int(i/FSroot))-1
                            FSpositionY[i]=int(np.floor(RSroot/FSroot))*(1+np.int(np.mod(i,FSroot)))-1
                        for proj in range(self.RS.size): #finds the euclidian distance between each FS and RS cell. Uses this distance to find the weight in the 3D Gaussian (wraps around)
                            for receiv in range(self.FS.size):
                                x=FSpositionX[receiv]-int(proj/RSroot) - int((FSpositionX[receiv]-int(proj/RSroot))/int(RSroot/2))*np.mod(FSpositionX[receiv]-int(proj/RSroot),int(RSroot/2))*2
                                y=FSpositionY[receiv]-int(np.mod(proj,RSroot)) - int((FSpositionY[receiv]-int(np.mod(proj,RSroot)))/int(RSroot/2))*np.mod(FSpositionY[receiv]-int(np.mod(proj,RSroot)),int(RSroot/2))*2 
                                syns.g_max_matrix[proj,receiv] = syns.g_max*np.exp(-pow(x/syns.extent,2))*np.exp(-pow(y/syns.extent,2)) 
                                syns.g[proj,receiv] *= np.exp(-pow(x/syns.extent,2))*np.exp(-pow(y/syns.extent,2))
                else: #FS->FS and RS->RS
                    for projecting in range(np.size(syns.g,0)):
                        root=np.size(weights,0)
                        rolledWeights = np.roll(weights, -int(np.floor(root/2))+np.int(projecting/root), 0)
                        rolledWeights = np.roll(rolledWeights, -int(np.floor(root/2))+np.int(np.mod(projecting,root)), 1)
                        syns.g_max_matrix[projecting,:] = syns.g_max*np.reshape(rolledWeights,np.size(syns.g[projecting,:]))
                        syns.g[projecting,:] *= np.reshape(rolledWeights,np.size(syns.g[projecting,:]))
#                        if syns.preName=='RS' and np.all(projecting!=self.plast):syns.g[projecting,:]=np.zeros(np.size(syns.g[projecting,:])) #only a few central neurons have lateral connections
#                        if syns.preName=='RS':syns.g[projecting,:]=np.zeros(np.size(syns.g[projecting,:]))
                
                #set diagonal of self-projecting weights to zero
                if True and syns.preName==syns.postName:
                    np.fill_diagonal(syns.g, 0)
                    np.fill_diagonal(syns.g_max_matrix, 0)
                       
                #plot weights
                if False:
                    print "plotting weights of", syns.preName, "->", syns.postName
                    plt.figure()
                    for projecting in range(np.size(syns.g,0)):
                        plt.subplot(int(np.sqrt(np.size(syns.g,0))),int(np.sqrt(np.size(syns.g,0))),projecting+1)
                        plt.imshow(np.reshape(syns.g[projecting,:],[np.size(weights,0),np.size(weights,0)]), interpolation='nearest', vmin=0, vmax=syns.g_max)
                        plt.gca().axes.get_xaxis().set_visible(False)
                        plt.gca().axes.get_yaxis().set_visible(False)
                    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    plt.colorbar(cax=cax)
                    plt.suptitle(syns.preName+' -> '+syns.postName)    
            
        plt.show()

    def createInput(self):
        """ create an image input """   
        
        if not self.naturalImage:
            #creates inputs: either four vertical lines or two diagonal lines
            isi = np.random.uniform(low=29, high=50, size=self.TC.size)
            for i in range(self.TC.size):
                if self.numPattern == 2:
                    #if np.mod(i,np.sqrt(self.TC.size)+1)==0: isi[i] = 41 # /
                    #if np.mod(i,np.sqrt(self.TC.size)-1)==0 and i!=0 and i!=self.TC.size-1: isi[i] = 41 # \
                    if np.mod(i,4)==0:isi[i] = 41
                    if np.mod(i,4)==1:isi[i] = 41
                elif self.numPattern == 3:
                    if np.mod(i,4)==0:isi[i] = 41
                    if np.mod(i,4)==1:isi[i] = 41
                    if np.mod(i,4)==2:isi[i] = 41
                elif self.numPattern == 4:
                    if np.mod(i,4)==0:isi[i] = 41
                    if np.mod(i,4)==1:isi[i] = 41
                    if np.mod(i,4)==2:isi[i] = 41
                    if np.mod(i,4)==3:isi[i] = 41

            #fill in the spike time array based on the interspike interval
            isi = np.reshape(np.floor(isi/self.dt), [self.TC.size,1])
            self.TC.spikeTimes = np.ones([self.TC.size,self.trialDuration/self.dt+1])*np.linspace(0,self.trialDuration/self.dt,self.trialDuration/self.dt+1) 
            spikeTimes_Mask = np.mod(self.TC.spikeTimes,isi)==0
            self.TC.spikeTimes[spikeTimes_Mask]=1 #put a spike in the spikeTime array based on the interspike interval
            self.TC.spikeTimes[~spikeTimes_Mask]=0
            self.TC.spikeTimes[:,0]=0          
            
            for i in range(self.TC.size): # roll the firing times of each pattern to create an offset in the patterns
                if self.numPattern == 2:
                    #if np.mod(i,np.sqrt(self.TC.size)+1)==0: self.TC.spikeTimes[i] = np.roll(self.TC.spikeTimes[i],int((41/2)/self.dt))
                    if np.mod(i,4)==1: self.TC.spikeTimes[i] = np.roll(self.TC.spikeTimes[i],int((41/2)/self.dt))
                if self.numPattern == 3:
                    if np.mod(i,4)==1: self.TC.spikeTimes[i] = np.roll(self.TC.spikeTimes[i],int((41/3)/self.dt))
                    if np.mod(i,4)==2: self.TC.spikeTimes[i] = np.roll(self.TC.spikeTimes[i],int((41*2/3)/self.dt))
                if self.numPattern == 4:
                    if np.mod(i,4)==1: self.TC.spikeTimes[i] = np.roll(self.TC.spikeTimes[i],int((41/4)/self.dt))
                    if np.mod(i,4)==2: self.TC.spikeTimes[i] = np.roll(self.TC.spikeTimes[i],int((41/2)/self.dt))
                    if np.mod(i,4)==3: self.TC.spikeTimes[i] = np.roll(self.TC.spikeTimes[i],int((41*3/4)/self.dt))
                    
        
        if self.naturalImage: #constructs input to TC from a natual image

            highISI = 65.0 #maximal ISI
            lowISI = 64.0 #minimal ISI
            timeLow = 0 #time interval during which to keep the input (ignore spike times before this point)
            timeHigh = timeLow+self.stimPresentation
            self.allImages = scipy.io.loadmat('IMAGES.mat')['IMAGES'] #load the image set
            
            if np.mod(self.trialDuration,self.stimPresentation)!=0: print "!!! training duration is not a multiple of stimulus presentation time !!!"
            
            for i in range(int(self.trialDuration/self.stimPresentation)):
                #compute an ISI for each TC based on the intensity of its input pixel
#                imageIndex = 5
#                imageX = 200
#                imageY = 200
                imageIndex = np.floor(np.random.uniform(0,np.size(self.allImages,2))) #select a random image from the image set
                imageX = np.floor(np.random.uniform(0,np.size(self.allImages,0)-self.imageSize)) #X index of the top left corner of an image patch
                imageY = np.floor(np.random.uniform(0,np.size(self.allImages,0)-self.imageSize)) #Y index of the top left corner of an image patch
                self.inputM = self.allImages[imageY:imageY+self.imageSize, imageX:imageX+self.imageSize, imageIndex] #select image patch
                self.inputM = np.clip(self.inputM,0,1) #creates an ON channel
                isi = np.reshape(np.abs(self.inputM-1.0)*(highISI-lowISI)+lowISI,[self.TC.size,1]) #compute interspike interval from pixel intensity
                isi[isi==highISI] = self.trialDuration/self.dt #make sure that the neurons that receive zero input never spike
                isi = np.floor(isi/self.dt) #convert time to time steps, and round to use as index
                
                #fill in the spike time array based on the interspike interval; on consider times between timeLow and timeHigh
                tmpSpikeTimes = np.ones([self.TC.size,timeHigh/self.dt])*np.linspace(0,timeHigh/self.dt-1,timeHigh/self.dt) 
                spikeTimes_Mask = np.mod(tmpSpikeTimes,isi)==0 
                tmpSpikeTimes[spikeTimes_Mask]=1 #put a spike in the spikeTime array based on the interspike interval
                tmpSpikeTimes[~spikeTimes_Mask]=0
                self.TC.spikeTimes[:,i*self.stimPresentation/self.dt:(i+1)*self.stimPresentation/self.dt] = tmpSpikeTimes[:,timeLow/self.dt:timeHigh/self.dt] #only consider times between timeLow and timeHigh
                
                if False: #plot image 
                    plt.figure()
                    plt.imshow(self.allImages[:,:,imageIndex], cmap = 'gray', interpolation='nearest')
                    plt.title('whole image')
            if True: #plot image patch
                plt.figure()
                plt.imshow(self.inputM, cmap = 'gray', interpolation='nearest', vmin=-1, vmax=1)
                plt.title('image patch')
#                    plt.show()
                
            if False: #plot image statistics
                plt.hist(self.allImages.reshape(-1), bins=100)
                F1 = fftpack.fft2(self.allImages[:,:,5]) #performs fourier transform
                F2 = fftpack.fftshift( F1 ) #shifts the FT so that low frequencies are at the center
                psd2D = abs(F2)**2 #calculate the 2D power spectral density
                plt.figure()
                plt.imshow(np.log10(psd2D))
                plt.figure()
        
    def runSimulation(self, isi=None):
        """ simulate the propagation of activity in the neural network """
        
        """ The simulation implements a multi-compartment neuron model (soma and spine) to explicitly simulate back-propagating action potentials (BPAPs) in spines.
            At the spines, the depolarization associated with BPAPs activates voltage-gated calcium channels (VGCCs) and removes the Mg2+ blocker from
            NMDA receptors. Pre-synaptic activity increases the open probabilty (OP) of ion channels at the spines, changing the conductance of the membrane.
            The interaction between pre-synaptic acivity and the BPAP results in calcium currents. Changes in intracellular calcium concentration at the 
            spines are read by four calcium detectors (P, D, B and V) controlling synaptic plasticity. Plasticity is reflected in modifications of the 
            excitatory conductance of individual spines  """
        
        #copy initial synaptic weights
        self.TC_RS_initWeight = np.copy(self.TC_RS.g[:,:])
        self.TC_FS_initWeight = np.copy(self.TC_FS.g[:,:])
        self.RS_RS_initWeight = np.copy(self.RS_RS.g[:,:])
        self.RS_RS_initWeight = np.copy(self.RS_FS.g[:,:])
        self.FS_RS_initWeight = np.copy(self.FS_RS.g[:,:])
        
        if not self.STDPplot: print "prep time:", int((time.time()-self.tic_prep)/60), 'min,',  int(np.mod((time.time()-self.tic_prep),60)), 'sec'
        tic_run=time.time()
        
        for t in range(np.size(self.time_array)): #loop through all time steps of the trial
            
            #used to compute weight change at each time step
            TC_RSmemory = np.copy(self.TC_RS.g[:,:])
            changeBin = 0.0
            
            if self.STDPplot: #triggers spikes in the TC and RS; used to plot the STDP curve
                if np.mod(t*self.dt,1000)==40: self.TC.lastCellSpike=np.ones(np.shape(self.TC.lastCellSpike))*t #triggers a pre-synaptic spike at t=40ms
                if np.mod(t*self.dt,1000)==40-isi: self.RS.Vm[0,t-1] = self.Vth + 10. #triggers a post-synaptic spike at t=40ms-isi
#                if np.mod(t*self.dt,1000)==40-isi-2: self.FS.Vm[0,t-1] = self.Vth + 10. #triggers a spike in a FS neuron at t=40ms-isi-2ms
            
            for pops in {self.RS, self.FS}: #loop through the different populations and compute membrane currents based on conductance change from last time step
                #1- reset neurons that spiked at last time step
                pastSpikes = pops.Vm[:,t-1]>=self.Vspike-1 #find the neurons that spiked at last time step
                pops.Vm[pastSpikes,t]=self.Vreset #reset Vm to Vreset after a spike
                
                #2- compute currents & associated voltage changes for the neurons that didn't spike during last time step
                I_leak = self.g_leak*(self.E_leak-pops.Vm[~pastSpikes,t-1]) #leakage current
                I_inj = pops.Ie[~pastSpikes] #injected current
                I_ref = pops.Rm*pops.Gref[~pastSpikes]*(self.E_i-pops.Vm[~pastSpikes,t-1]) #refractory period current
                I_syn = 0 #total synaptic currents (initialize)
                for syns in pops.receivesFrom: #compute all the synaptic currents received by a population of neurons
                    syns.Mg = 1./(1.+np.exp(-0.092*(syns.Vm+13))*0.56) #compute the voltage-dependent blockade level of NMDAr by Mg
                    syns.I = syns.g*syns.pre.OP*(syns.E-syns.Vm) #AMPA and GABA current
                    I_syn += sum(syns.I,0)[~pastSpikes]  #sum of AMPA and GABA synaptic currents
                    syns.V_BPAP = self.Vmax_BPAP*syns.post.BPAP #BPAP depolarization
                    syns.VmPSP += ((self.E_leak-syns.VmPSP) + (pops.Rm/pops.Rm)*syns.I*1.1)/pops.Tau_m*self.dt#syns.VmPSP += ((self.E_leak-syns.VmPSP) + pops.Rm*syns.I*0.28)/pops.Tau_m*self.dt #compute the synaptically local EPSP
                    syns.Vm = syns.VmPSP + syns.V_BPAP #superposition of the synaptic PSP and bAP
                    syns.I_NMDA = self.g_NMDA*syns.pre.OP_NMDA*syns.Mg*(self.E_Ca-syns.Vm) #NMDA calcium current
                    syns.I_VGCC = 1./(1+np.exp((syns.Vm-25.)/7.))-1./(1.+np.exp((syns.Vm+15.)/7.)) #VGCC calcium current 
                    syns.calcium += 5*(syns.I_NMDA+0.1*syns.I_VGCC)-(syns.calcium-self.Crest)/self.tau_Ca #update intracellular calcium concentration                     
                pops.Vm[~pastSpikes, t] = pops.Vm[~pastSpikes,t-1] + (I_leak + I_inj + I_ref + I_syn)/pops.Tau_m*self.dt #increment Vm at the cell soma
                
                #3- make the neurons that reached threshold during the ongoing time step spike
                currentSpike = pops.Vm[:, t] >= self.Vth #find the neurons that reached threshold during the ongoing time step
                pops.Vm[currentSpike,t] = self.Vspike #make those neurons spike
                pops.spikeTimes[currentSpike,t] = 1 #record all spike times for raster plot
                pops.lastCellSpike[currentSpike] = t #record last spike time

            self.TC.lastCellSpike[self.TC.spikeTimes[:,t]==1] = t #make the TC neuron spike based on their pre-determined firing rate
                
            for pops in self.population: #loop through the populations and compute changes in conductance due to activity during the ongoing time step
                #4- compute changes in refractory period conductance due to post-synaptic spikes
                pops.Gref[(t-pops.lastCellSpike)<1./self.dt] = self.GrefMax
                pops.Gref[(t-pops.lastCellSpike)>=1./self.dt] -= (pops.Gref[(t-pops.lastCellSpike)>=1./self.dt]/pops.tau_Gref)*self.dt
                #5- compute the changes in conductance (OP) due to pre-synaptic spikes
                deltaT = -(t-pops.lastCellSpike)*self.dt #time since last pre-synaptic spike
                pops.OP = np.reshape(0.5*np.exp(deltaT/self.tau_ampaF)+0.5*np.exp(deltaT/self.tau_ampaS),[pops.size,1]) #compute AMPA OP
                pops.OP_NMDA = np.reshape(0.85*np.exp(deltaT/self.tau_nmdaF)+0.15*np.exp(deltaT/self.tau_nmdaS),[pops.size,1]) #compute NMDA OP
                if pops.name == 'RS': #compute the difference between inhibitory and excitatory conductance; used to suppress BPAP
                    pops.g_excit = np.sum(95*self.RS_RS.g*self.RS.OP,0) #total excitatory conductance
                    pops.g_inhib = np.sum(375*self.FS_RS.g*self.FS.OP,0) #total inhibitory conductance 
#                    I_tot = 1.0/(1.0+np.exp((pops.g_inhib-pops.g_excit-8)/2)) #sigmoid relationship between g_inhib-g_excit and BPAP amplitude decrease 
                    I_tot = np.clip(-0.072*(pops.g_inhib-pops.g_excit)+1.0,0,1) #linear relationship between g_inhib-g_excit and BPAP amplitude decrease, clipped to max=1
                else: I_tot = np.ones(pops.size)
                mask = np.logical_and(-deltaT>=1.0,-deltaT<=2.0) #triggers a BPAP in the spines from 1ms to 2ms after a spike at the soma
                pops.BPAP[mask] = I_tot[mask] #the amplitude of the BPAP is decreased proportionally to EPSC and IPSC.
                pops.BPAP[~mask] = 0. #BPAPs are simulated as square waves; set to zero 2ms after spike initiation at the soma 

            #7- compute changes in synaptic efficacy (plasticity) - only TC->RS synapses exhibit plasticity
            for syns in self.synapses:
                if (syns.preName == 'TC' and syns.postName == 'RS'):
                    #compute changes in the different calcium detector (D, B, P, V)
                    syns.D += (2.4/(1.+np.exp((syns.B-2.5)/-0.20))-syns.D*15)*self.dt/275. #D: depression detector; depends on the peak of B
                    syns.B += (25./(1.+np.exp((syns.calcium-0.65)/-0.01))-syns.B*3-4.*syns.B*syns.V*8)*self.dt/40. #B: reads of [Ca] timecourse; highest for continuous [Ca]>0.65
                    syns.P += ((7.*pow(syns.calcium/4.,4)/(1.+pow(syns.calcium/4.,4)))-26.*syns.P)*self.dt/500. #P: potentiation; increases propotionally to [Ca]
                    syns.V += (4./(1.+np.exp((syns.calcium-1.82)/-0.05))-syns.V*2.4)*self.dt/10. #V: veto; prevents depression for high [Ca]
                    syns.g = np.clip(syns.g + (self.lr_LTP*syns.P-self.lr_LTD*syns.D)*self.dt, 0, syns.g_max_matrix) #update synaptic weights based on P and D
               
            #value trackers
            changeBin += np.sum(np.abs(TC_RSmemory-self.TC_RS.g))/np.sum(self.TC_RS.g)
            if np.mod(t,500)==0: 
                self.percentChange[int(t/500)]=changeBin/500
                changeBin=0.0  
            self.STDPtracker[0,t] = self.RS.Vm[0,t]
            self.STDPtracker[1,t] = self.TC_RS.Vm[0,0]
            self.STDPtracker[2,t] = self.TC_RS.calcium[0,0]
            self.STDPtracker[3,t] = self.TC_RS.P[0,0]
            self.STDPtracker[4,t] = self.TC_RS.D[0,0]
            self.STDPtracker[5,t] = self.TC_RS.g[0,0]-self.TC_RS_initWeight[0,0]
            self.STDPtracker[6,t] = self.FS.Vm[0,t]
            self.STDPtracker[7,t] = self.RS.g_inhib[0]
            self.STDPtracker[8,t] = self.RS.g_excit[0]
            self.STDPtracker[9,t] = self.TC.OP_NMDA[0,0]
            self.STDPtracker[10,t] = self.TC_RS.Mg[0,0]
            self.STDPtracker[11,t] = self.TC_RS.I_VGCC[0,0]
            self.STDPtracker[12,t] = self.TC_RS.B[0,0]
            
            self.calciumTracker[0,t] = self.TC_RS.calcium[0,0]
            self.CaDetectorsTracker[0,t] = self.TC_RS.P[0,0]
            self.CaDetectorsTracker[1,t] = self.TC_RS.D[0,0]
            self.CaDetectorsTracker[2,t] = self.TC_RS.B[0,0]
            self.CaDetectorsTracker[3,t] = self.TC_RS.V[0,0]
            
            if not self.STDPplot: 
                self.trackValue(t)
                if self.ODC_evol and np.mod(t*self.dt,1000)==0: self.showODC(t*self.dt)
                if self.TC_RS_evol and np.mod(t*self.dt,1000)==0:self.showTC_RS(t*self.dt)
        
        # returns the change in synaptic conductance
        if self.STDPplot: return self.TC_RS.g[0,0]-self.TC_RS_initWeight[0,0]
        print "run time:", int((time.time()-tic_run)/60), 'min,',  int(np.mod((time.time()-tic_run),60)), 'sec'
     
    def trackValue(self, t):
        """ track different values """
        
        self.OPTracker[0,t]=self.TC.OP_NMDA[0]
        
        self.valueTracker[0,t] = self.RS.Vm[0,t]                   
        self.valueTracker[1,t] = self.TC_RS.Vm[0,0]
        self.valueTracker[2,t] = self.TC_RS.Vm[1,0]
        self.valueTracker[3,t] = self.TC_RS.Vm[2,0]
        self.valueTracker[4,t] = self.TC_RS.Vm[3,0]
        self.valueTracker[5,t] = self.TC_RS.Vm[4,0]
        self.valueTracker[6,t] = self.TC_RS.Vm[5,0]
        
        self.weightTracker[0,t]  = np.mean(self.TC_RS.g[[0,5,10,15],self.RS1])#RS neuron 1 \
        self.weightTracker[1,t]  = np.mean(self.TC_RS.g[[1,5, 9,13],self.RS1])
        self.weightTracker[2,t]  = np.mean(self.TC_RS.g[[2,6,10,14],self.RS1])
        self.weightTracker[3,t]  = np.mean(self.TC_RS.g[[3,6, 9,12],self.RS1])
        self.weightTracker[4,t]  = np.mean(self.TC_RS.g[[0,5,10,15],self.RS2])#RS neuron 2 \
        self.weightTracker[5,t]  = np.mean(self.TC_RS.g[[1,5, 9,13],self.RS2])
        self.weightTracker[6,t]  = np.mean(self.TC_RS.g[[2,6,10,14],self.RS2])
        self.weightTracker[7,t]  = np.mean(self.TC_RS.g[[3,6, 9,12],self.RS2])
        self.weightTracker[8,t]  = np.mean(self.TC_RS.g[[0,5,10,15],self.RS3])#RS neuron 3 \
        self.weightTracker[9,t]  = np.mean(self.TC_RS.g[[1,5, 9,13],self.RS3])
        self.weightTracker[10,t] = np.mean(self.TC_RS.g[[2,6,10,14],self.RS3])
        self.weightTracker[11,t] = np.mean(self.TC_RS.g[[3,6, 9,12],self.RS3])
        self.weightTracker[12,t] = np.mean(self.TC_FS.g[[0,5,10,15],self.FS1])#FS neuron 1 \
        self.weightTracker[13,t] = np.mean(self.TC_FS.g[[1,5, 9,13],self.FS1])
        self.weightTracker[14,t] = np.mean(self.TC_FS.g[[2,6,10,14],self.FS1])
        self.weightTracker[15,t] = np.mean(self.TC_FS.g[[3,6, 9,12],self.FS1])
             
        self.valueTracker[7,t]  = self.TC_RS.Vm[0,self.RS1]#RS neuron 1
        self.valueTracker[8,t]  = self.RS.g_excit[self.RS1]
        self.valueTracker[9,t]  = self.RS.g_inhib[self.RS1]
        self.valueTracker[10,t] = self.TC_RS.P[0,self.RS1]
        self.valueTracker[11,t] = self.TC_RS.D[0,self.RS1]
        self.valueTracker[12,t] = self.TC_RS.Vm[0,self.RS2]#RS neuron 2
        self.valueTracker[13,t] = self.RS.g_excit[self.RS2]
        self.valueTracker[14,t] = self.RS.g_inhib[self.RS2]
        self.valueTracker[15,t] = self.TC_RS.P[0,self.RS2]
        self.valueTracker[16,t] = self.TC_RS.D[0,self.RS2]
        self.valueTracker[17,t] = self.TC_RS.Vm[0,self.RS3]#RS neuron 3
        self.valueTracker[18,t] = self.RS.g_excit[self.RS3]
        self.valueTracker[19,t] = self.RS.g_inhib[self.RS3]
        self.valueTracker[20,t] = self.TC_RS.P[0,self.RS3]
        self.valueTracker[21,t] = self.TC_RS.D[0,self.RS3]
        self.valueTracker[22,t] =  self.TC_FS.Vm[0,self.FS1]#FS neuron 1
        self.valueTracker[23,t] =  self.FS.g_excit[self.FS1]
        self.valueTracker[24,t] =  self.FS.g_inhib[self.FS1]
        self.valueTracker[25,t] = self.TC_FS.P[0,self.FS1]
        self.valueTracker[26,t] = self.TC_FS.D[0,self.FS1]
    
    def showResults(self, t=[]):
        """ display plots """
        
        tic_plot=time.time()
                
        if self.RS1_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
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
            plt.savefig('../output/' + 'RS1_Vm')
        
        if self.RS2_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
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
            plt.savefig('../output/' + 'RS2_Vm')
            
        if self.RS3_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
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
            plt.savefig('../output/' + 'RS3_Vm')
        
        if self.FS1_Vm:#displays evolution of membrane and synaptic potentials of neuron RS1
            plt.figure()
            plt.plot(self.time_array, self.valueTracker[22,:], 'b') #synaptic Vm
            plt.plot(self.time_array, self.FS.Vm[self.FS1,:], 'k') #soma Vm
            plt.plot(self.time_array, self.valueTracker[23,:], 'b') #synpatic lateral excitation
            plt.plot(self.time_array, self.valueTracker[24,:], 'r') #synpatic lateral inhibition
            plt.plot(self.time_array, self.valueTracker[25,:]*1000+80, 'b') #potentiation variable
            plt.plot(self.time_array, -self.valueTracker[26,:]*1000+80, 'r') #depression variable
            plt.ylim(-100,100)
            plt.title('FS ' + str(self.FS1))
            plt.savefig('../output/' + 'FS1_Vm')
        
        if self.calcium:#displays evolution of the calcium detectors
            plt.figure()
            plt.plot(self.time_array, self.calciumTracker[0,:], 'k')
            plt.plot(self.time_array, self.CaDetectorsTracker[0,:]*10+2, 'b') #P
            plt.plot(self.time_array, self.CaDetectorsTracker[1,:]*10+2, 'r') #D
            plt.plot(self.time_array, self.CaDetectorsTracker[2,:], 'c') #B
            plt.plot(self.time_array, self.CaDetectorsTracker[3,:], 'm') #V
            plt.ylim(-0.1,3.3)
            plt.savefig('../output/' + 'calcium')
        
        if self.raster:#display spike raster plot
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
            plt.savefig('../output/' + 'raster')

        if self.Vm_TC:#displays Vm for all TC cells
            plt.figure()
            for TC in range(np.size(self.TC.Vm,0)):
                plt.subplot(int(np.floor(np.sqrt(np.size(self.TC.Vm,0)))),int(np.ceil(np.sqrt(np.size(self.TC.Vm,0)))),TC+1)
                plt.plot(self.time_array, self.TC.Vm[TC,:])
                plt.ylim(-90,50)
            plt.savefig('../output/' + 'Vm_TC') 
        
        if self.TC_RS_init: #display initial TC->RS weights            
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
        
        if self.TC_RS_final: #display final TC->RS weights
            self.showTC_RS()
            
        if self.TC_FS_final: #display initial TC->FS weights
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
        
        if self.TC_FS_final: #display final TC->FS weights
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
            
        if self.RS_RS_init:
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
            plt.savefig('../output/' + 'RS_RS_init')
        
        if self.RS_RS_final:
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
            plt.savefig('../output/' + 'RS_RS_final')
           
        if self.RS_FS_init:
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
            plt.savefig('../output/' + 'RS_FS_init')
         
        if self.RS_FS_final:
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
            plt.savefig('../output/' + 'RS_FS_final')
            
        if self.ODC: #dipslays an ocular dominance map
            self.showODC()
        
        if self.RS1_weight: #display weight evolution of RS neuron 1
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
            plt.savefig('../output/' + 'RS_1')
        
        if self.RS2_weight: #display weight evolution of RS neuron 2
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
            plt.savefig('../output/' + 'RS_2')
            
        if self.RS3_weight: #display weight evolution of RS neuron 3
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
            plt.savefig('../output/' + 'RS_3')
            
        if self.FS1_weight: #display weight evolution of FS neuron 1
            plt.figure()
            plt.plot(self.time_array, self.weightTracker[12,:],  'b')
#            plt.plot(self.time_array, self.weightTracker[13,:],  'r') 
#            plt.plot(self.time_array, self.weightTracker[14,:], 'k') 
            plt.plot(self.time_array, self.weightTracker[15,:], 'r') 
            plt.ylim(0,self.TC_RS.g_max)
            plt.legend('\\')
            plt.title("FS " + np.str(self.FS1))
            plt.savefig('../output/' + 'FS_1')
            
        if self.weightChange: #displays the percent of total weight change every at every 100ms
            plt.figure()
            x=np.linspace(0,self.trialDuration,int((self.trialDuration/500)/self.dt)+1)
            plt.plot(x, self.percentChange, '-ob')
            plt.title('percent weight change')
            plt.savefig('../output/' + 'percent_weight_change')
            
        if self.extraPlot: #displays an extra, custom plot
            plt.figure()
            plt.plot(self.time_array, self.valueTracker[0], 'k')
#            plt.plot(self.time_array, self.valueTracker[1]+10, 'b')
            plt.plot(self.time_array, self.valueTracker[2], 'b')
#            plt.plot(self.time_array, self.valueTracker[3]+30, 'b')
#            plt.plot(self.time_array, self.valueTracker[4]+40, 'b')
#            plt.plot(self.time_array, self.valueTracker[5]+50, 'b')
#            plt.plot(self.time_array, self.valueTracker[6]+60, 'b')
            plt.ylim(-100,50)
            plt.savefig('../output/' + 'extraPlot')
            
        if not self.STDPplot: print "plot time:", int((time.time()-tic_plot)/60), 'min,',  int(np.mod((time.time()-tic_plot),60)), 'sec'

    def showODC(self, t=[]):
        plt.figure(figsize=(7,7))
        rootSize = np.sqrt(self.RS.size)
        ODC_mat = np.zeros(self.RS.size)
        alpha_mat = np.zeros(self.RS.size)
        #creat white color map with transparency gradient
        cmap_trans = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',['black','black'],256) 
        cmap_trans._init()
        alphas = np.linspace(1.0, 0, cmap_trans.N+3)
        cmap_trans._lut[:,-1] = alphas
        
        for RS in range(self.RS.size):
            prefPattern = [np.sum(self.TC_RS.g[[0,4,8,12],RS]),np.sum(self.TC_RS.g[[1,5,9,13],RS]),np.sum(self.TC_RS.g[[2,6,10,14],RS]),np.sum(self.TC_RS.g[[3,7,11,15],RS])]
            ODC_mat[RS] = np.argmax(prefPattern)
            alpha_mat[RS] = np.max(prefPattern)-(np.sum(prefPattern)-np.max(prefPattern))/self.numPattern
#            ODC_mat[RS] = np.mean(self.TC_RS.g[[0,5,10,15],RS]) - np.mean(self.TC_RS.g[[3,6,9,12],RS])   
        plt.imshow(np.reshape(ODC_mat, [rootSize, rootSize]),interpolation='nearest', cmap='Spectral') #color
        plt.imshow(np.reshape(alpha_mat, [rootSize, rootSize]),interpolation='nearest', cmap=cmap_trans, vmin=-0.25,vmax=1.5) #transparency
        plt.title('Ocular Dominance at ' + np.str(t) + 'ms')
        if t == []: plt.savefig('../output/' + 'OcularDominance_final')
        else: plt.savefig('../output/' + 'OcularDominance_' + np.str(int(t))) 
        
        #print np.reshape(-alpha_mat, [rootSize, rootSize])
        
        for RS in range(self.RS.size):
            prefPattern = [np.sum(self.TC_RS.g[[0,4,8,12],RS]),np.sum(self.TC_RS.g[[1,5,9,13],RS]),np.sum(self.TC_RS.g[[2,6,10,14],RS]),np.sum(self.TC_RS.g[[3,7,11,15],RS])]
            ODC_mat[RS] = np.argmax(prefPattern)  
        plt.imshow(np.reshape(ODC_mat, [rootSize, rootSize]),interpolation='nearest', cmap='Spectral') #color
        plt.title('no transparency')
        plt.savefig('../output/' + 'OcularDominance_NT_' + np.str(int(t)))
#         file = open('DATA/12-2_OD_multi3-2/' + 'TC_RS_g_' + np.str(int(t)), 'w')
#         pickle.dump(self.TC_RS.g, file)
#         file.close()
    
    def showTC_RS(self, t=[]):
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
        if t == []: plt.suptitle('Final TC->RS Weights')
        else: plt.suptitle('TC->RS Weights at ' + np.str(t) + 'ms')
        print t
        if t == []: plt.savefig('../output/' + 'TC-RS_final')
        else: plt.savefig('../output/' + 'TC-RS_' + np.str(int(t)))
    