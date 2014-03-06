'''
Created on Nov 13, 2012

@author: raphaelholca
'''

from numpy import *

class Population(object):

    def __init__(self, sim, name, size, Rm, Cm, Ie, tau_Gref, receivesFrom={}, projectsTo={}):
        ''' 
        Object Constructor 
        '''
        
        size_time = [size,sim.trialDuration/sim.dt+1]
        
        #neuron parameters
        self.name = name
        self.size = size
        self.Rm = Rm
        self.Cm = Cm
        self.Tau_m = Cm*Rm
        self.Ie = Ie
        self.tau_Gref = tau_Gref
        
        #Synapse projections
        self.receivesFrom = receivesFrom
        self.projectsTo = projectsTo
        
        #neuron variables
        self.Vm = zeros(size_time)
        self.spikeTimes = zeros(size_time)
        self.OP = zeros(size_time)
        self.OP_NMDA = zeros(size_time)
        self.BPAP = zeros(size_time)
        self.g_excit = zeros(size_time)
        self.g_inhib = zeros(size_time)
        self.Gref = zeros(size_time)
        
        self.Vm[:,-1] = sim.E_leak
        self.lastCellSpike = ones(size)*-1e99