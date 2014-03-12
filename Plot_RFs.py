'''
Created on Mar 4, 2014

@author: raphaelholca
'''

from Tkinter import *
import Tkinter as Tk
import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib as mpl

class Plot_RFs(Tk.Tk):
    
    def __init__(self, root, folder):
        Tk.Tk.__init__(self)
        self.root = root
        self.folder = folder
        self.title('RFs - ' + folder)
        self.protocol('WM_DELETE_WINDOW', self.closeGUI)
        
        self.createRFsGUI()
    
    def createRFsGUI(self):
        
        #declare button variables
        self.plotContent = Tk.StringVar()
        self.plotContent.set('color-coded')
        
        self.RFtime = Tk.StringVar()
        self.RFtime.set('0')
        
        #create plotframe widget
        mainFrame = Tk.Frame(self)
        mainFrame.grid(column=0, row=0)
        # mainFrame.columnconfigure(0, weight=1) #window resizing
        # mainFrame.rowconfigure(0, weight=1) #window resizing
        self.plotFrame = Tk.Frame(mainFrame, borderwidth=5, relief='sunken', width=500, height=500)
        self.plotFrame.grid(column=0, row=0, columnspan=3, rowspan=2, sticky=(N, W, E), padx=(10,10), pady=(10,10))
        buttonFrame = Tk.Frame(mainFrame, borderwidth=5, relief='sunken', width=500, height=100)
        buttonFrame.grid(column=0, row=3, columnspan=3, rowspan=1, sticky=(W, E, S), padx=(10,10), pady=(10,10))
        
        #create main figure
        self.mainFig = Figure()
        
        #create main plot canvas
        self.canvas = FigureCanvasTkAgg(self.mainFig, master=self.plotFrame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        
        #create radio buttons
        self.color_coded = Tk.Radiobutton(buttonFrame, text='color-coded',  command=self.codedMap, variable=self.plotContent, value='color-coded', width=12)
        self.color_coded.grid(column=0, row=0, sticky=(W, E), padx=(0,20))
        self.color_coded.select()
          
        self.detailed    = Tk.Radiobutton(buttonFrame, text='detailed',     command=self.detailedMap, variable=self.plotContent, value='detailed', width=12)
        self.detailed.grid(column=1, row=0, sticky=(W, E), padx=(0,20))
        self.detailed.deselect()
          
        self.individual  = Tk.Radiobutton(buttonFrame, text='individual',   command=self.individualRF, variable=self.plotContent, value='individual', width=12)
        self.individual.grid(column=2, row=0, sticky=(W, E), padx=(0,20))
        self.individual.deselect()
        
        #create time entry field
        self.timeInput = Tk.Entry(buttonFrame, textvariable=self.RFtime, width=12)
        self.timeInput.grid(column=0, row=1, sticky=(W))
        self.timeInput.insert(0,'0')
        self.timeInput.bind('<Return>', self.roundTime)
        
        #create scale widget for time selection
        self.timeScale = Tk.Scale(buttonFrame, orient=HORIZONTAL, from_=0, to=self.root.genParam[self.folder]['trialDuration'], resolution=self.root.genParam[self.folder]['RF_sampleTime'], 
                                  command=self.timeScale_C, length=500, showvalue=0)
        self.timeScale.grid(column=1, row=1, columnspan=5, sticky=(E,W))
        
        self.allCommand = {'color-coded':self.color_coded, 'detailed':self.detailed, 'individual':self.individual}
        self.codedMap()
        self.timeInput.focus_set()
        
    def timeScale_C(self, event):
        self.timeInput.delete(0, END)
        self.timeInput.insert(0,str(self.timeScale.get()))
        self.RFtime.set(str(self.timeScale.get()))
        self.allCommand[self.plotContent.get()].invoke()
    
    def roundTime(self, event):
        ''' executed on <Return>; rounds the requested time and executes command associated with radio button'''
        #find the closest saved weight snapshot to requested time
        self.RFtime.set(self.timeInput.get())
        m = np.mod(float(self.RFtime.get()),self.root.genParam[self.folder]['RF_sampleTime'])
        round_t = np.floor(float(self.RFtime.get())/self.root.genParam[self.folder]['RF_sampleTime'])*self.root.genParam[self.folder]['RF_sampleTime'] + np.round(m*2,-2)/2
        
        if round_t/self.root.genParam[self.folder]['RF_sampleTime']>=np.size(self.root.weights[self.folder]['w'],2): #checks for out of bound requested time
            round_t = (np.size(self.root.weights[self.folder]['w'],2)-1)*self.root.genParam[self.folder]['RF_sampleTime']
            self.mainFig.text(0.5,0.5,'time out-of-bound',bbox=dict(facecolor='red', alpha=0.8), horizontalalignment='center', verticalalignment='center')
            self.canvas.show()
        
        self.RFtime.set(str(int(round_t)))
        self.timeScale.set(int(round_t))
        self.timeInput.delete(0, END)
        self.timeInput.insert(0, str(int(round_t)))
        
        #execute the command associated with the current check radiobutton
        if event: self.allCommand[self.plotContent.get()].invoke()
       
    def codedMap(self):
        #plot color-coded receptive fields
        self.plotContent.set('color-coded')
        self.mainFig.clf() #clear existing figure
        self.plt = self.mainFig.add_subplot(111) #creates subplot
        
        #retrieve time
        self.roundTime(False)
        t = int(self.RFtime.get())
        i = int(t/self.root.genParam[self.folder]['RF_sampleTime'])
        
        squareW = self.root.weights[self.folder]['w'][:,:,i]    
        rootSize = np.sqrt(self.root.neurons[self.folder]['RS'].size)
        ODC_mat = np.zeros(self.root.neurons[self.folder]['RS'].size)
        alpha_mat = np.zeros(self.root.neurons[self.folder]['RS'].size)
        #create white color map with transparency gradient
        cmap_trans = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',['black','black'],256) 
        cmap_trans._init()
        alphas = np.linspace(1.0, 0, cmap_trans.N+3)
        cmap_trans._lut[:,-1] = alphas
          
        for RS in range(self.root.neurons[self.folder]['RS'].size):
            prefPattern = [np.sum(squareW[[0,4,8,12],RS]),np.sum(squareW[[1,5,9,13],RS]),np.sum(squareW[[2,6,10,14],RS]),np.sum(squareW[[3,7,11,15],RS])]
            ODC_mat[RS] = np.argmax(prefPattern)
            alpha_mat[RS] = np.max(prefPattern)-(np.sum(prefPattern)-np.max(prefPattern))/self.root.genParam[self.folder]['numPattern']
        #            ODC_mat[RS] = np.mean(self.TC_RS.g[[0,5,10,15],RS]) - np.mean(self.TC_RS.g[[3,6,9,12],RS])   
        self.plt.imshow(np.reshape(ODC_mat, [rootSize, rootSize]),interpolation='nearest', cmap='Spectral', vmin=0,vmax=3) #color
        self.plt.imshow(np.reshape(alpha_mat, [rootSize, rootSize]),interpolation='nearest', cmap=cmap_trans, vmin=-0.25,vmax=1.5) #transparency
        self.plt.set_xticks([])
        self.plt.set_yticks([])
        self.canvas.show()
        
    def detailedMap(self):
        #plot detailed receptive fields
        self.plotContent.set('detailed')
        self.mainFig.clf() #clear existing figure
        
        #retrieve time
        self.roundTime(False)
        t = int(self.RFtime.get())
        i = t/self.root.genParam[self.folder]['RF_sampleTime']
        
        rootSize = np.sqrt(self.root.neurons[self.folder]['TC'].size)
        squareW = self.root.weights[self.folder]['w'][:,:,i]
        for RS in range(self.root.neurons[self.folder]['RS'].size):
            subplt = self.mainFig.add_subplot(int(np.ceil(np.sqrt(self.root.neurons[self.folder]['RS'].size))),int(np.ceil(np.sqrt(self.root.neurons[self.folder]['RS'].size))), RS+1)
            subplt.imshow(np.reshape(squareW[:,RS],[rootSize, rootSize]), interpolation='nearest', vmin=0, vmax=self.root.synParam[self.folder]['TC_RS'].g_max, cmap='bwr')
            subplt.set_xticks([])
            subplt.set_yticks([])
#             self.plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#         cax = self.plt.axes([0.85, 0.1, 0.075, 0.8])
#         self.plt.colorbar(cax=cax)
    
        self.canvas.show()
        
    def individualRF(self):
        print '!!NOT IMPLEMENTED YET!!'
    
    def closeGUI(self, allGUI=False):
        if not allGUI: self.root.allGUI.remove(self)
        self.destroy()
        