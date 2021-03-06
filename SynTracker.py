'''
Created on Mar 3, 2014

@author: raphaelholca
'''

from numpy import *

class SynTracker(object):
    
    def __init__(self, neuronID, size_time):
        ''' 
        Object Constructor; tracker of synapses variables 
        '''
        
        self.neuronID = neuronID
        self.synapseNum = size_time[0]
        self.time = size_time[1]
        
        #synapse variables
#         self.V_BPAP = zeros(1) #depolarization due to the back-propagating action-potential
        self.Vm = zeros(size_time) #local membrane potential at the synapse
        self.g = zeros(size_time) #synaptic weights 
        self.calcium = zeros(size_time) #internal "calcium" concentration, used to compute LTD
        self.Mg = zeros(size_time) #extent of the Mg blockade of the NMDAr
        self.I = zeros(size_time) #current through AMPAr and GABAr
        self.I_NMDA = zeros(size_time) #Ca2+ current through NMDAr
        self.I_VGCC = zeros(size_time) #Ca2+ current through voltage-gated calcium channels
             
        
        #calcium detectors
        self.P=zeros(size_time)
        self.V=zeros(size_time)
        self.B=zeros(size_time)
        self.D=zeros(size_time)
