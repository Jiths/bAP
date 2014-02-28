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
        size_t = [size,sim.trialDuration/sim.dt+1]
        #Neuron parameters
        self.name = name
        self.size = size
        self.Rm = Rm
        self.Cm = Cm
        self.Tau_m = Cm*Rm
        self.Ie = Ie
        self.I_inhib = zeros(size)
        self.Vm = zeros(size_t)
        self.Vm[:,-1] = sim.E_leak
        self.Gref = zeros(size)
        self.lastCellSpike = ones(size)*-1e99
        self.spikeTimes = zeros([size,sim.trialDuration/sim.dt+1])
        self.OP = zeros([size,1])
        self.OP_NMDA = zeros([size,1])
        self.BPAP = zeros(size)
        self.g_excit = zeros(size)
        self.g_inhib = zeros(size)
        self.tau_Gref = tau_Gref
        
        #Synapse projections
        self.receivesFrom = receivesFrom
        self.projectsTo = projectsTo
        
        