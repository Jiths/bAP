'''
Created on Mar 4, 2014

@author: raphaelholca
'''

from Tkinter import *
import Tkinter as Tk
import ttk
import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import cPickle as pickle

from Plot_Traces import Plot_Traces
from Plot_RFs import Plot_RFs

class Plot_Root(Tk.Tk):
    
    def __init__(self):
        Tk.Tk.__init__(self)
        self.title("ROOT")
        self.wm_attributes("-topmost",1)
        
        self.loadData()
        self.createRootGUI()
        self.allGUI = []

    def loadData(self):
        #load data
        
        pFile = open('../output/neurons', 'r')
        self.neurons = pickle.load(pFile)
        pFile.close()
        
        pFile = open('../output/synTracker', 'r')
        self.synapses = pickle.load(pFile)
        pFile.close()
        
        pFile = open('../output/synParam', 'r')
        self.synParam = pickle.load(pFile)
        pFile.close()
        
        pFile = open('../output/weights', 'r')
        self.weights = pickle.load(pFile)
        pFile.close()
        
        pFile = open('../output/genParam', 'r')
        self.genParam = pickle.load(pFile)
        pFile.close()
    
    def createRootGUI(self):
        
        #create mainFrame widget
        mainFrame = Tk.Frame(self, borderwidth=5, relief='sunken')
        mainFrame.grid(column=0, row=0, sticky=(N, W, E, S), padx=(5,5), pady=(5,5))
        
        #create click buttons
        button = Tk.Button(mainFrame, text='new trace GUI',  command=self.launchTracesGUI)
        button.grid(column=0, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))
        
        button = Tk.Button(mainFrame, text='new RF GUI',  command=self.launchRFsGUI)
        button.grid(column=1, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))
        
        button = Tk.Button(mainFrame, text='refresh root',  command=self.loadData)
        button.grid(column=2, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))
        
        button = Tk.Button(mainFrame, text='close all GUI',  command=self.closeAll)
        button.grid(column=3, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))

    
    def launchTracesGUI(self):
        #creates new Plot_Traces object
        traces = Plot_Traces(self)
        self.allGUI.append(traces)
        traces.mainloop()
        
    def launchRFsGUI(self):
        #creates new Plot_Traces object
        RFs = Plot_RFs(self)
        self.allGUI.append(RFs)
        RFs.mainloop()
        
    def closeAll(self):
        #close all open GUIs
        for gui in self.allGUI: gui.closeAll()
        self.destroy()