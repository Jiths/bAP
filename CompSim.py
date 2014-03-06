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
from SynTracker import SynTracker
import cPickle as pickle
import sys


class CompSim():
    '''
    This is the main simulation function.
    '''
    def __init__(self, duration=1000., dt=0.2):

        self.tic_prep = time.time()
        
        """ general simulation parameters """
        self.trialDuration = duration #time to simulate (msec)
        self.dt = dt #length of each simulation time step (msec)
        self.time_array = np.linspace(0,self.trialDuration,self.trialDuration/self.dt+1) #time array
        self.naturalImage = [] #whether to use a natural image as an input - leave empty to use artificial stimulus
        self.imageSize = 4 #pixel size of the SIDE of the input image (total size = imageSize^2)
        self.popSize = 5 #side of the size of the population of neurons (total size = popSize^2)
        self.numPattern = 2 #number of patterns embedded in the input
        self.FS_RS = [1] #turns off lateral inhibition to excitatory cells - leave empty for OFF, put 1 for ON
        self.RS_RS = [1] #turns off lateral excitation to excitatory cells - leave empty for OFF, put 1 for ON
        self.RS_FS = [1] #turns off lateral excitation to inhibitory cells - leave empty for OFF, put 1 for ON
        self.FS_FS = [1] #turns off lateral inhibition to inhibitory cells - leave empty for OFF, put 1 for ON
        self.gaussian = [] #whether 'intracortical' connections are distributed according to a Gaussian distribution. If not, follows a square wave
        np.random.seed(83251) #seeds the random number generator
        self.neuronsToTrack = ['0', '1', '3', '4', '5', '8'] #Indices of RS neurons whose synaptic variables should be tracked
        self.RF_sampleTime = 50. #time interval at which to sample (take a snapshot) of the TC_RS weight (RF) -- in ms
        self.plotRFbool = False #plot all RFs in .png file at the end of the simulation
      
    def setParam(self):
        
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
        self.TC_RS_weightRand = 0.03 #variability in weight initialization - weight initialization
        self.TC_RS_rise = 5 #rise time constant of the excitatory synaptic conductance
        self.TC_RS_decay = 0.5 #decay time constant of the excitatory synaptic conductance
        self.TC_RS_E = self.E_e #synaptic reversial potential
        self.TC_RS_G_max = 0.35 #maximal synaptic conductance - weights cannot grow stronger than this value
        
        """ TC->FS synapses parameters """#not in use
        self.TC_FS_size = [self.TC_size, self.FS_size]
        self.TC_FS_G = 0.0#0.1 #maximal unitary conductance of the thalamocortical synapses on FS cells
        self.TC_FS_weightRand = 0.0#0.02 #variability in weight initialization
        self.TC_FS_rise = 10.0 #rise time constant of the excitatory synaptic conductance
        self.TC_FS_decay = 0.5 #decay time constant of the excitatory synaptic conductance
        self.TC_FS_E = self.E_e #synaptic reversial potential
        self.TC_FS_G_max = 0.0#0.2 #maximal synaptic conductance - weights cannot grow stronger than this value
        
        """ RS->FS synapses parameters """
        self.RS_FS_size = [self.RS_size, self.FS_size]
        self.RS_FS_G = 1.0 #maximal unitary conductance, when OP=1 - weight initialization
        if not self.RS_FS: self.RS_FS_G = 0.0
        self.RS_FS_weightRand = 0.0 #MUST BE 0 #variability in weight initialization
        self.RS_FS_rise = 10.0 #rise time constant of the excitatory synaptic conductance
        self.RS_FS_decay = 2.2 #decay time constant of the excitatory synaptic conductance
        self.RS_FS_E = self.E_e #synaptic reversial potential
        self.RS_FS_G_max = 1.0 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.RS_FS_extent = 1.0#3.0 #extent of lateral connections
        
        """ FS->RS synapses parameters """
        self.FS_RS_size = [self.FS_size, self.RS_size]
        self.FS_RS_G = 0.006 #maximal unitary conductance of the inhibitory synapses on RS cells - weight initialization
        if not self.FS_RS: self.FS_RS_G = 0. #turns inhibition off
        self.FS_RS_weightRand = 0 #MUST BE 0 #variability in weight initialization
        self.FS_RS_rise = 5.5 #rise time constant of the inhibitory synaptic conductance
        self.FS_RS_decay = 0.5 #decay time constant of the inhibitory synaptic conductance
        self.FS_RS_E = self.E_i #synaptic reversial potential
        self.FS_RS_G_max = 0.013#0.13 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.FS_RS_extent = 5 #extent of lateral connections
        
        """ RS->RS synapses parameters """
        self.RS_RS_size = [self.RS_size, self.RS_size]
        self.RS_RS_G = 0.06#0.12 #initial maximal unitary conductance - weight initialization
        if not self.RS_RS: self.RS_RS_G = 0.
        self.RS_RS_weightRand = 0. #MUST BE 0 #variability in weight initialization
        self.RS_RS_rise = 5 #rise time constant of the excitatory synaptic conductance
        self.RS_RS_decay = 0.5 #decay time constant of the excitatory synaptic conductance
        self.RS_RS_E = self.E_e #synaptic reversial potential
        self.RS_RS_G_max = 0.12 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.RS_RS_extent = 3#1.0 #extent of lateral connections
        
        """ FS->FS synapses parameters """ #not in use
        self.FS_FS_size = [self.FS_size, self.FS_size]
        self.FS_FS_G = 0.0#0.2 #maximal unitary conductance of the inhibitory synapses on RS cells - weight initialization
        if not self.FS_FS: self.FS_FS_G = 0. #turns inhibition off
        self.FS_FS_weightRand = 0 #MUST BE 0 #variability in weight initialization
        self.FS_FS_rise = 5.5 #rise time constant of the inhibitory synaptic conductance
        self.FS_FS_decay = 0.5 #decay time constant of the inhibitory synaptic conductance
        self.FS_FS_E = self.E_i #synaptic reversial potential
        self.FS_FS_G_max = 0.2 #maximal synaptic conductance - weights cannot grow stronger than this value
        self.FS_FS_extent = 1.#0.8 #extent of lateral connections
        
        """ value trackers and others """
        
        self.weightTracker = np.zeros([20,self.trialDuration/self.dt+1]) #tracks value of weights over time
        self.OPTracker=np.zeros([3,self.trialDuration/self.dt+1]) #array space to record OP
        self.calciumTracker = np.zeros([3,self.trialDuration/self.dt+1]) #tracks value of the calcium transient over time
        self.valueTracker = np.zeros([30,self.trialDuration/self.dt+1]) #tracks other values...
        self.CaDetectorsTracker = np.zeros([5,self.trialDuration/self.dt+1]) #tracks value of the calcium detector variables
        self.STDPtracker = np.zeros([15,self.trialDuration/self.dt+1]) #tracks value for STDP plots
        self.percentChange = np.zeros(int(self.trialDuration/self.dt/500)+1) #tracks how much weight changed during the last time step
        self.STDPplot = False #boolean used to Plot_Print the STDP curve - set to True by setParamSTDP()
        self.stimPresentation = 500.0 #duration of presentation of a particular image patch 
        
        if not self.STDPplot: print "duration:", int(self.trialDuration), "ms\nTC:", self.TC_size, "\nRS:", self.RS_size, "\nFS:", self.FS_size
    
    def setParamSTDPplot(self):
        """ sets the  parameter values for the STDP Plot_Print """
               
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
        """ resets variables between each iteration of the STDP Plot_Print function """        
        
        for neuron in [self.TC, self.RS, self.FS]: 
            neuron.Vm = np.zeros([1,self.trialDuration/self.dt+1])
            neuron.Vm[:,-1] = self.E_leak
            neuron.lastCellSpike = np.ones(1)*-1e99
            neuron.spikeTimes = np.zeros([1,self.trialDuration/self.dt+1])
            neuron.OP = np.zeros([1,1])
            neuron.OP_NMDA = np.zeros([1,1])
            neuron.BPAP = np.zeros(1)
        
        for synapses in [self.TC_RS, self.TC_FS, self.FS_RS]:
            synapses.OP = np.zeros([1,1,self.trialDuration/self.dt+1])
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
        self.TC_RS.g = np.ones([1,1])*self.STDPplot_initialG_TC_RS
        self.TC_FS.g = np.ones([1,1])*self.STDPplot_initialG_TC_FS
        self.FS_RS.g = np.ones([1,1])*self.STDPplot_initialG_FS_RS
        
    def createSynapses(self):
        """ create synapses """
        
        self.TC_RS = Synapses(self, 'TC', 'RS', {}, {}, self.TC_RS_size, self.TC_RS_G, self.TC_RS_weightRand, self.TC_RS_G_max, self.TC_RS_rise, self.TC_RS_decay, self.TC_RS_E)
        self.FS_RS = Synapses(self, 'FS', 'RS', {}, {}, self.FS_RS_size, self.FS_RS_G, self.FS_RS_weightRand, self.FS_RS_G_max, self.FS_RS_rise, self.FS_RS_decay, self.FS_RS_E, self.FS_RS_extent)
        self.RS_FS = Synapses(self, 'RS', 'FS', {}, {}, self.RS_FS_size, self.RS_FS_G, self.RS_FS_weightRand, self.RS_FS_G_max, self.RS_FS_rise, self.RS_FS_decay, self.RS_FS_E, self.RS_FS_extent)
        self.RS_RS = Synapses(self, 'RS', 'RS', {}, {}, self.RS_RS_size, self.RS_RS_G, self.RS_RS_weightRand, self.RS_RS_G_max, self.RS_RS_rise, self.RS_RS_decay, self.RS_RS_E, self.RS_RS_extent)
        self.synapses = {self.TC_RS, self.FS_RS, self.RS_FS, self.RS_RS} #some synapse population are not in use
        
    def createPopulations(self):
        """ create neuron populations and connections between neurons """
        #all neuron matrix are projectsFrom x projectsTo 
        
        """ neural populations """
        self.TC = Population(self, 'TC', self.TC_size, self.TC_Rm, self.TC_Cm, self.TC_Ie, self.TC_tau_Gref)
        self.RS = Population(self, 'RS', self.RS_size, self.RS_Rm, self.RS_Cm, self.RS_Ie, self.RS_tau_Gref)
        self.FS = Population(self, 'FS', self.FS_size, self.FS_Rm, self.FS_Cm, self.FS_Ie, self.FS_tau_Gref)
        
        """ receiving synapses """
        self.TC.receivesFrom = [] #synapse population TC is post-synaptic to
        self.RS.receivesFrom = [self.TC_RS, self.FS_RS, self.RS_RS] #synapse population RS is post-synaptic to
        self.FS.receivesFrom = [self.RS_FS] #synapse population FS is post-synaptic to
        
        """ projecting synapses """
        self.TC.projectsTo = [self.TC_RS] #synapse population TC is pre-synaptic to
        self.RS.projectsTo = [self.RS_FS, self.RS_RS] #synapse population RS is pre-synaptic to
        self.FS.projectsTo = [self.FS_RS] #synapse population FS is pre-synaptic to
        
        """ pointer to pre- and post-synaptic population """
        self.TC_RS.pre = self.TC
        self.TC_RS.post = self.RS
        self.FS_RS.pre = self.FS
        self.FS_RS.post = self.RS
        self.RS_FS.pre = self.RS
        self.RS_FS.post = self.FS
        self.RS_RS.pre = self.RS
        self.RS_RS.post = self.RS
        
        self.population = {self.TC, self.FS, self.RS}
    
    def createSynTracker(self):
        """ create objects to track the synaptic variables for a few selected neurons """
        
        self.allSynTracker = dict([])
        for n in self.neuronsToTrack:
            self.allSynTracker[str(n)] = SynTracker(n, [self.TC_size, np.size(self.time_array)])
            
        #creates an array to store synaptic weights over time
        self.TC_RS_gTracker = np.zeros([self.TC_size*self.RS_size, np.size(self.time_array)/250+1])
        self.timeStamp = np.zeros(np.size(self.time_array)/250+1)
    
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
                       
                #Plot_Print weights
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
                
                if False: #Plot_Print image 
                    plt.figure()
                    plt.imshow(self.allImages[:,:,imageIndex], cmap = 'gray', interpolation='nearest')
                    plt.title('whole image')
            if True: #Plot_Print image patch
                plt.figure()
                plt.imshow(self.inputM, cmap = 'gray', interpolation='nearest', vmin=-1, vmax=1)
                plt.title('image patch')
#                    plt.show()
                
            if False: #Plot_Print image statistics
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
        
        print 'Simulation:'
        #copy initial synaptic weights for STDP Plot_Print
        self.TC_RS_initWeight = np.copy(self.TC_RS.g[:,:])

        tic_run=time.time()
        for t in range(np.size(self.time_array)): #loop through all time steps of the trial
            
            if self.STDPplot: #triggers spikes in the TC and RS; used to Plot_Print the STDP curve
                if np.mod(t*self.dt,1000)==40: self.TC.lastCellSpike=np.ones(np.shape(self.TC.lastCellSpike))*t #triggers a pre-synaptic spike at t=40ms
                if np.mod(t*self.dt,1000)==40-isi: self.RS.Vm[0,t-1] = self.Vth + 10. #triggers a post-synaptic spike at t=40ms-isi
#                if np.mod(t*self.dt,1000)==40-isi-2: self.FS.Vm[0,t-1] = self.Vth + 10. #triggers a spike in a FS neuron at t=40ms-isi-2ms
            
            for pops in [self.RS, self.FS]: #loop through the different populations and compute membrane currents based on conductance change from last time step
                #1- reset neurons that spiked at last time step
                pastSpikes = pops.Vm[:,t-1]>=self.Vspike-1 #find the neurons that spiked at last time step
                pops.Vm[pastSpikes,t]=self.Vreset #reset Vm to Vreset after a spike
                
                #2- compute currents & associated voltage changes for the neurons that didn't spike during last time step
                I_leak = self.g_leak*(self.E_leak-pops.Vm[~pastSpikes,t-1]) #leakage current
                I_inj = pops.Ie[~pastSpikes] #injected current
                I_ref = pops.Rm*pops.Gref[:,t-1][~pastSpikes]*(self.E_i-pops.Vm[~pastSpikes,t-1]) #refractory period current
                I_syn = 0 #total synaptic currents (initialize)
                for syns in pops.receivesFrom: #compute all the synaptic currents received by a population of neurons
                    syns.Mg = 1./(1.+np.exp(-0.092*(syns.Vm+13))*0.56) #compute the voltage-dependent blockade level of NMDAr by Mg
                    syns.I = syns.g*syns.pre.OP[:,t-1][:,np.newaxis]*(syns.E-syns.Vm) #AMPA and GABA current
                    I_syn += sum(syns.I,0)[~pastSpikes]  #sum of AMPA and GABA synaptic currents
                    syns.V_BPAP = self.Vmax_BPAP*syns.post.BPAP[:,t-1] #BPAP depolarization
                    syns.VmPSP += ((self.E_leak-syns.VmPSP) + (pops.Rm/pops.Rm)*syns.I*1.1)/pops.Tau_m*self.dt#syns.VmPSP += ((self.E_leak-syns.VmPSP) + pops.Rm*syns.I*0.28)/pops.Tau_m*self.dt #compute the synaptically local EPSP
                    syns.Vm = syns.VmPSP + syns.V_BPAP #superposition of the synaptic PSP and bAP
                    syns.I_NMDA = self.g_NMDA*syns.pre.OP_NMDA[:,t-1][:,np.newaxis]*syns.Mg*(self.E_Ca-syns.Vm) #NMDA calcium current
                    syns.I_VGCC = 1./(1+np.exp((syns.Vm-25.)/7.))-1./(1.+np.exp((syns.Vm+15.)/7.)) #VGCC calcium current 
                    syns.calcium += 5*(syns.I_NMDA+0.1*syns.I_VGCC)-(syns.calcium-self.Crest)/self.tau_Ca #update intracellular calcium concentration                     
                pops.Vm[~pastSpikes, t] = pops.Vm[~pastSpikes,t-1] + (I_leak + I_inj + I_ref + I_syn)/pops.Tau_m*self.dt #increment Vm at the cell soma
                
                #3- make the neurons that reached threshold during the ongoing time step spike
                if t*self.dt>6000 and pops.name == 'FS' and (t-self.TC.lastCellSpike[0])*self.dt<10.:thrashThis=0 #added to remove inhibition after a spike in the TC neuron 0 (input pattern 1)
                else:
                        currentSpike = pops.Vm[:, t] >= self.Vth #find the neurons that reached threshold during the ongoing time step
                        pops.Vm[currentSpike,t] = self.Vspike #make those neurons spike
                        pops.spikeTimes[currentSpike,t] = 1 #record all spike times for raster Plot_Print
                        pops.lastCellSpike[currentSpike] = t #record last spike time

            self.TC.lastCellSpike[self.TC.spikeTimes[:,t]==1] = t #make the TC neuron spike based on their pre-determined firing rate
                
            for pops in self.population: #loop through the populations and compute changes in conductance due to activity during the ongoing time step
                #4- compute changes in refractory period conductance due to post-synaptic spikes
                pops.Gref[:,t][(t-pops.lastCellSpike)<1./self.dt] = self.GrefMax
                pops.Gref[:,t][(t-pops.lastCellSpike)>=1./self.dt] -= (pops.Gref[:,t][(t-pops.lastCellSpike)>=1./self.dt]/pops.tau_Gref)*self.dt
                #5- compute the changes in conductance (OP) due to pre-synaptic spikes
                deltaT = -(t-pops.lastCellSpike)*self.dt #time since last pre-synaptic spike
                pops.OP[:,t] = (0.5*np.exp(deltaT/self.tau_ampaF)+0.5*np.exp(deltaT/self.tau_ampaS)) #compute AMPA OP
                pops.OP_NMDA[:,t] = (0.85*np.exp(deltaT/self.tau_nmdaF)+0.15*np.exp(deltaT/self.tau_nmdaS)) #compute NMDA OP
                
            for pops in self.population:
                if pops.name == 'RS': #compute the difference between inhibitory and excitatory conductance; used to suppress BPAP
                    pops.g_excit[:,t] = np.sum(95*self.RS_RS.g*self.RS.OP[:,t],0) #total excitatory conductance
                    pops.g_inhib[:,t] = np.sum(375*self.FS_RS.g*self.FS.OP[:,t],0) #total inhibitory conductance 
                    I_tot = np.clip(-0.072*(pops.g_inhib[:,t]-pops.g_excit[:,t])+1.0,0,1) #linear relationship between g_inhib-g_excit and BPAP amplitude decrease, clipped to max=1
                else: I_tot = np.ones(pops.size)
                mask = np.logical_and(-deltaT>=1.2,-deltaT<=2.2) #triggers a BPAP in the spines from 1.2ms to 2.2ms after a spike at the soma
                pops.BPAP[:,t][mask] = I_tot[mask] #the amplitude of the BPAP is decreased proportionally to EPSC and IPSC.
                pops.BPAP[:,t][~mask] = 0. #BPAPs are simulated as square waves; set to zero 2.2ms after spike initiation at the soma 

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
            if not self.STDPplot: self.trackValue(t)
            
            if np.mod(t*self.dt,int(self.trialDuration/50))==0: 
                i=int(t*self.dt/int(self.trialDuration/50))
                timeLeft = -1
                if i==0: timeStart = time.time()
                else: timeLeft = (50-i)*(time.time()-timeStart)/i
                sys.stdout.write('\r')
                sys.stdout.write("[%-50s] %d%% done in: %d min %d sec" % ('='*i, 2*i, int((timeLeft)/60), int(np.mod((timeLeft),60)) ))
                sys.stdout.flush()
        
        # returns the change in synaptic conductance for STDP Plot_Print
        if self.STDPplot: return self.TC_RS.g[0,0]-self.TC_RS_initWeight[0,0]
        print "\nrun time:", int((time.time()-tic_run)/60), 'min,',  int(np.mod((time.time()-tic_run),60)), 'sec'
        if self.plotRFbool: self.plotRF() 

    def pickleValue(self):
        print 'Pickle:'
        #save neuron variables to file
        tic_pickle = time.time()
        pFile = open('../output/neurons', 'w')
        pickle.dump({'TC':self.TC, 'FS':self.FS, 'RS':self.RS}, pFile, protocol=2)
        pFile.close()
        
        #save synapse parameters to file
        pFile = open('../output/synParam', 'w')
        pickle.dump({'TC_RS':self.TC_RS, 'FS_RS':self.FS_RS, 'RS_FS':self.RS_FS, 'RS_RS':self.RS_RS}, pFile, protocol=2)
        pFile.close()
        
        #save synapse variables to file
        pFile = open('../output/synTracker', 'w')
        pickle.dump(self.allSynTracker, pFile, protocol=2)
        pFile.close()
        
        #save TC_RS weights to file
        pFile = open('../output/weights', 'w')
        pickle.dump({'w':self.TC_RS_gTracker, 'time':self.timeStamp}, pFile, protocol=2)
        pFile.close()
        
        #save simulation paramters to file
        pFile = open('../output/genParam', 'w')
        pickle.dump({'numPattern':self.numPattern, 'dt':self.dt, 'trialDuration':self.trialDuration, 'timeArray':self.time_array,
                     'RF_sampleTime':self.RF_sampleTime, 'neuronsToTrack':self.neuronsToTrack}, pFile, protocol=2)
        pFile.close()
        print "pickle time:", int((time.time()-tic_pickle)/60), 'min,',  int(np.mod((time.time()-tic_pickle),60)), 'sec'
     
    def trackValue(self, t):
        """ track synaptic variables """
        
        for n in self.allSynTracker.keys():
            self.allSynTracker[n].I_NMDA[:,t]  = self.TC_RS.I_NMDA  [:,int(n)]  
            self.allSynTracker[n].I_VGCC[:,t]  = self.TC_RS.I_VGCC  [:,int(n)]
            self.allSynTracker[n].I[:,t]       = self.TC_RS.I       [:,int(n)]
            self.allSynTracker[n].Vm[:,t]      = self.TC_RS.Vm      [:,int(n)]
            self.allSynTracker[n].Mg[:,t]      = self.TC_RS.Mg      [:,int(n)]
            self.allSynTracker[n].g[:,t]       = self.TC_RS.g       [:,int(n)]
            self.allSynTracker[n].calcium[:,t] = self.TC_RS.calcium [:,int(n)]
            self.allSynTracker[n].P[:,t]       = self.TC_RS.P       [:,int(n)]
            self.allSynTracker[n].V[:,t]       = self.TC_RS.V       [:,int(n)]
            self.allSynTracker[n].B[:,t]       = self.TC_RS.B       [:,int(n)]
            self.allSynTracker[n].D[:,t]       = self.TC_RS.D       [:,int(n)]
        
        step_sample = self.RF_sampleTime/self.dt
        if np.mod(t,step_sample)==0: 
            self.TC_RS_gTracker[:,t/step_sample] = np.reshape(np.copy(self.TC_RS.g),-1,1)
            self.timeStamp[t/step_sample] = t*self.dt

    def plotRF(self):
        if True: return
        tic_plot = time.time()
        print 'Plot:'
        for t,i in zip(self.timeStamp,range(np.size(self.timeStamp))):
#             if np.mod(t*self.dt,int(self.trialDuration/50))==0: 
#                 j=int(t*self.dt/int(self.trialDuration/50))
#                 timeLeft = -1
#                 if j==0: timeStart = time.time()
#                 else: timeLeft = (50-j)*(time.time()-timeStart)/j
#                 sys.stdout.write('\r')
#                 sys.stdout.write("[%-50s] %d%% done in: %d min %d sec" % ('='*j, 2*j, int((timeLeft)/60), int(np.mod((timeLeft),60)) ))
#                 sys.stdout.flush()
                 
            if np.mod(t,200) == 0:
                #plot color-coded map
                squareW = np.reshape(self.TC_RS_gTracker[:,i], (self.TC.size, self.RS.size))    
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
                    prefPattern = [np.sum(squareW[[0,4,8,12],RS]),np.sum(squareW[[1,5,9,13],RS]),np.sum(squareW[[2,6,10,14],RS]),np.sum(squareW[[3,7,11,15],RS])]
                    ODC_mat[RS] = np.argmax(prefPattern)
                    alpha_mat[RS] = np.max(prefPattern)-(np.sum(prefPattern)-np.max(prefPattern))/self.numPattern
                #            ODC_mat[RS] = np.mean(self.TC_RS.g[[0,5,10,15],RS]) - np.mean(self.TC_RS.g[[3,6,9,12],RS])   
                plt.imshow(np.reshape(ODC_mat, [rootSize, rootSize]),interpolation='nearest', cmap='Spectral', vmin=0,vmax=3) #color
                plt.imshow(np.reshape(alpha_mat, [rootSize, rootSize]),interpolation='nearest', cmap=cmap_trans, vmin=-0.25,vmax=1.5) #transparency
                plt.title('Ocular Dominance at ' + np.str(t) + 'ms')
                plt.savefig('../output/' + 'codedMap/' + np.str(int(t)) + '.png')
                 
                #plot detailed map
                plt.figure()
                rootSize = np.sqrt(self.TC.size)
                for RS in range(self.RS.size):
                    plt.subplot(int(np.ceil(np.sqrt(self.RS.size))),int(np.ceil(np.sqrt(self.RS.size))), RS+1)
                    plt.imshow(np.reshape(squareW[:,RS],[rootSize, rootSize]), interpolation='nearest', vmin=0, vmax=self.TC_RS.g_max, cmap='bwr')
                    plt.gca().axes.get_xaxis().set_visible(False)
                    plt.gca().axes.get_yaxis().set_visible(False)
                    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                plt.colorbar(cax=cax)
                plt.suptitle('TC->RS Weights at ' + np.str(t) + 'ms')
                plt.savefig('../output/' + 'detailedMap/' + np.str(int(t)) + '.png')
        print "plot time:", int((time.time()-tic_plot)/60), 'min,',  int(np.mod((time.time()-tic_plot),60)), 'sec'