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
import os

from Plot_Traces import Plot_Traces
from Plot_RFs import Plot_RFs

class Plot_Root(Tk.Tk):
    
    def __init__(self):
        Tk.Tk.__init__(self)
        self.title("ROOT")
        self.wm_attributes("-topmost",1)
        
        self.createRootGUI()
        self.allGUI = []
        self.neurons = {}
        self.synapses = {}
        self.synParam = {}
        self.weights = {}
        self.genParam = {}
        self.alreadyLoaded = []

    def chooseData(self, command):
        #creates dialog window to choose the data to view
        
        path = '../output/'
        folders = []
        for x in os.walk(path+'.').next()[1]: 
            folders.append(x)
        
        self.dataDialog = Tk.Toplevel()
        self.dataDialog.title('choose data to view')
        
        mainFrame = Tk.Frame(self.dataDialog, borderwidth=5, relief='sunken')
        mainFrame.grid(column=0, row=0, sticky=(N, W, E, S), padx=(5,5), pady=(5,5))
        
        rowCount = 0
        for f in folders:
            button = Tk.Button(mainFrame, text=f, command=lambda f_call=f:self.loadData(path, f_call, command))
            button.grid(column=0, row=rowCount, sticky=(W, E), padx=(5,5), pady=(5,5))
            rowCount += 1
         
    def loadData(self, path, folder, command):
        #loads the requested data

        self.dataDialog.destroy()
        
        if folder not in self.alreadyLoaded: #makes sure to not uselessly load same data twice 
        
            pFile = open(path + folder + '/neurons', 'r')
            self.neurons[folder] = pickle.load(pFile)
            pFile.close()
               
            pFile = open(path + folder + '/synTracker', 'r')
            self.synapses[folder] = pickle.load(pFile)
            pFile.close()
               
            pFile = open(path + folder + '/synParam', 'r')
            self.synParam[folder] = pickle.load(pFile)
            pFile.close()
               
            pFile = open(path + folder + '/weights', 'r')
            self.weights[folder] = pickle.load(pFile)
            pFile.close()
               
            pFile = open(path + folder + '/genParam', 'r')
            self.genParam[folder] = pickle.load(pFile)
            pFile.close()
            
            self.alreadyLoaded.append(folder)
           
        if command:
            command(folder)
    
    def createRootGUI(self):
        
        #create mainFrame widget
        mainFrame = Tk.Frame(self, borderwidth=5, relief='sunken')
        mainFrame.grid(column=0, row=0, sticky=(N, W, E, S), padx=(5,5), pady=(5,5))
        
        #create click buttons
        button = Tk.Button(mainFrame, text='new trace GUI',  command=lambda:self.chooseData(self.launchTracesGUI))
        button.grid(column=0, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))
        
        button = Tk.Button(mainFrame, text='new RF GUI',  command=lambda:self.chooseData(self.launchRFsGUI))
        button.grid(column=1, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))
        
        button = Tk.Button(mainFrame, text='refresh root',  command=lambda:self.chooseData(False))
        button.grid(column=2, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))
        
        button = Tk.Button(mainFrame, text='close all GUI',  command=self.closeAll)
        button.grid(column=3, row=0, sticky=(W, E), padx=(5,5), pady=(5,5))

    
    def launchTracesGUI(self, folder):
        #creates new Plot_Traces object
        traces = Plot_Traces(self, folder)
        self.allGUI.append(traces)
        traces.mainloop()
        
    def launchRFsGUI(self, folder):
        #creates new Plot_Traces object
        RFs = Plot_RFs(self, folder)
        self.allGUI.append(RFs)
        RFs.mainloop()
        
    def closeAll(self):
        #close all open GUIs
        for gui in self.allGUI: gui.closeGUI(allGUI=True)
        self.destroy()