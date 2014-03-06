'''
Created on Feb 28, 2014

@author: raphaelholca
'''

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cPickle as pickle
from CompSim import CompSim

sim = CompSim()

def main(argv=None):
    #load data
    global neurons, synapses, synParam, weights, genParam
    
    pFile = open('../output/neurons', 'r')
    neurons = pickle.load(pFile)
    pFile.close()
    
    pFile = open('../output/synTracker', 'r')
    synapses = pickle.load(pFile)
    pFile.close()
    
    pFile = open('../output/synParam', 'r')
    synParam = pickle.load(pFile)
    pFile.close()
    
    pFile = open('../output/weights', 'r')
    weights = pickle.load(pFile)
    pFile.close()
    
    pFile = open('../output/genParam', 'r')
    genParam = pickle.load(pFile)
    pFile.close()
    
    toPlot()
    plot()

def toPlot():
    global traces, ODM
    global g_excit, g_inhib, OP, calcium, BPAP
    global time, ylim
    
    #plots to display
    traces      = [1]
    ODM         = [1]
     
    #traces to display        
    g_excit     = [1]
    g_inhib     = [1]
    OP          = [1]
    BPAP        = [1]
    calcium     = [1]
    
    #plot param
    time        = np.array([0.,500.])/genParam['dt']
    ylim        = [-10,60]

def plot():
    if traces: 
        plt.figure()
        if g_excit:  plt.plot(genParam['timeArray'],  neurons['RS'].g_excit[0,:]  , 'r')
        if g_inhib:  plt.plot(genParam['timeArray'],  neurons['RS'].g_inhib[0,:]  , 'b')
        if OP:       plt.plot(genParam['timeArray'],  neurons['RS'].OP[0,:]*30    , 'r.')
        if OP:       plt.plot(genParam['timeArray'],  neurons['FS'].OP[0,:]*30    , 'b.')
        if BPAP:     plt.plot(genParam['timeArray'],  neurons['RS'].BPAP[0,:]*50  , 'k')
        
        plt.ylim(ylim)
        plt.savefig('../output/' + 'g_record_delay.png')
     
    #plot stimulus preference map
    if ODM:
        for t,i in zip(weights['time'],range(np.size(weights['time']))): 
            if np.mod(t,1000) == 0:
                #color-coded map
                squareW = np.reshape(weights['w'][:,i], (neurons['TC'].size, neurons['RS'].size))    
                plt.figure(figsize=(7,7))
                rootSize = np.sqrt(neurons['RS'].size)
                ODC_mat = np.zeros(neurons['RS'].size)
                alpha_mat = np.zeros(neurons['RS'].size)
                #create white color map with transparency gradient
                cmap_trans = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',['black','black'],256) 
                cmap_trans._init()
                alphas = np.linspace(1.0, 0, cmap_trans.N+3)
                cmap_trans._lut[:,-1] = alphas
                 
                for RS in range(neurons['RS'].size):
                    prefPattern = [np.sum(squareW[[0,4,8,12],RS]),np.sum(squareW[[1,5,9,13],RS]),np.sum(squareW[[2,6,10,14],RS]),np.sum(squareW[[3,7,11,15],RS])]
                    ODC_mat[RS] = np.argmax(prefPattern)
                    alpha_mat[RS] = np.max(prefPattern)-(np.sum(prefPattern)-np.max(prefPattern))/genParam['numPattern']
                #            ODC_mat[RS] = np.mean(self.TC_RS.g[[0,5,10,15],RS]) - np.mean(self.TC_RS.g[[3,6,9,12],RS])   
                plt.imshow(np.reshape(ODC_mat, [rootSize, rootSize]),interpolation='nearest', cmap='Spectral', vmin=0,vmax=3) #color
                plt.imshow(np.reshape(alpha_mat, [rootSize, rootSize]),interpolation='nearest', cmap=cmap_trans, vmin=-0.25,vmax=1.5) #transparency
                plt.title('Ocular Dominance at ' + np.str(t) + 'ms')
                plt.savefig('../output/' + 'OcularDominance_' + np.str(int(t)) + '.png')
                 
                #full stimulus preference
                plt.figure()
                rootSize = np.sqrt(neurons['TC'].size)
                for RS in range(neurons['RS'].size):
                    plt.subplot(int(np.ceil(np.sqrt(neurons['RS'].size))),int(np.ceil(np.sqrt(neurons['RS'].size))), RS+1)
                    plt.imshow(np.reshape(squareW[:,RS],[rootSize, rootSize]), interpolation='nearest', vmin=0, vmax=synParam['TC_RS'].g_max, cmap='bwr')
                    plt.gca().axes.get_xaxis().set_visible(False)
                    plt.gca().axes.get_yaxis().set_visible(False)
                    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                plt.colorbar(cax=cax)
                plt.suptitle('TC->RS Weights at ' + np.str(t) + 'ms')
                plt.savefig('../output/' + 'TC-RS_' + np.str(int(t)) + '.png')

if __name__ == "__main__":
    sys.exit(main())