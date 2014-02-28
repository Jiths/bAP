'''
Created on Nov 13, 2012

@author: raphaelholca
'''

import sys
import cProfile

#Force matplotlib to not use any Xwindows backend; comment out to plot on screen
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from CompSim import CompSim
#@profile  
def main(argv=None):
        
        global sim
        sim = CompSim()
        sim.setParam()
        sim.createSynapses()
        sim.createPopulations()
        sim.createConnectionPattern()
        sim.createInput()
        sim.runSimulation()
#         cProfile.run('sim.runSimulation()')
#         sim.showResults()
        plt.show()
    
if __name__ == "__main__":
    sys.exit(main())
