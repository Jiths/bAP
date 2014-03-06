'''
Created on Mar 4, 2014

@author: raphaelholca
'''

from Tkinter import *
import Tkinter as Tk
import ttk
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import cPickle as pickle

class Plot_Traces(Tk.Tk):
    
    def __init__(self, root):
        Tk.Tk.__init__(self)
        self.title("trace GUI")
        
        self.root = root
        self.createTracesGUI()
    
    def createTracesGUI(self):
        
        #declare button variables
        #control variables for neuron
        self.neuronNumber  = Tk.StringVar()
        self.Vm_neuron     = Tk.BooleanVar()
        self.spikeTimes    = Tk.BooleanVar()
        self.OP            = Tk.BooleanVar()
        self.OP_NMDA       = Tk.BooleanVar()
        self.g_excit       = Tk.BooleanVar()
        self.g_inhib       = Tk.BooleanVar() 
        self.BPAP          = Tk.BooleanVar()
        self.G_ref         = Tk.BooleanVar()
        self.neuronAvail   = True
        
        #control variables for synapses
        self.synapseNumber = Tk.StringVar()
        self.Vm_syn        = Tk.BooleanVar()
        self.g             = Tk.BooleanVar()
        self.calcium       = Tk.BooleanVar()
        self.Mg            = Tk.BooleanVar()
        self.I_syn         = Tk.BooleanVar()
        self.I_NMDA        = Tk.BooleanVar()
        self.I_VGCC        = Tk.BooleanVar()
        self.P             = Tk.BooleanVar()
        self.B             = Tk.BooleanVar()
        self.D             = Tk.BooleanVar()
        self.V             = Tk.BooleanVar()
        
        #control variables for axes
        self.y_min         = Tk.StringVar()
        self.y_max         = Tk.StringVar()
        self.x_min         = Tk.StringVar()
        self.x_max         = Tk.StringVar()
        self.step_y        = 10
        self.step_x        = 100
        self.offValue      = False #False = use provided values
        self.legend        = False #True = display label info
        
        #create plotframe widget
        mainFrame = Tk.Frame(self)
        mainFrame.grid(column=0, row=0)
        plotFrame = Tk.Frame(mainFrame, borderwidth=5, relief='sunken', width=500, height=500)
        plotFrame.grid(column=0, row=0, columnspan=3, rowspan=2, sticky=(N, W, E), padx=(10,10), pady=(10,10))
        buttonFrame = Tk.Frame(mainFrame, borderwidth=5, relief='sunken', width=500, height=125)
        buttonFrame.grid(column=0, row=3, columnspan=7, rowspan=5, sticky=(N, S, E, W), padx=(10,10), pady=(10,10))
        
        #create main figure
        self.mainFig = Figure()
        
        #create main plot canvas
        self.canvas = FigureCanvasTkAgg(self.mainFig, master=plotFrame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        
        #create control for neuron variables
        Tk.Label(buttonFrame, text='neuron: ').grid(column=0, row=0, sticky=(W))
        self.neuronInput = Tk.Entry(buttonFrame, textvariable=self.neuronNumber, width=11)
        self.neuronInput.grid(column=1, row=0, sticky=(W))
        self.neuronNumber.set('0')
        self.neuronInput.insert(0,'0')
        self.neuronInput.bind('<Return>', self.renewCanvas)
        
        check = Tk.Checkbutton(buttonFrame, text='neuron Vm',    command=lambda:self.checked(self.Vm_neuron),  variable=self.Vm_neuron,  onvalue=True, offvalue=False)
        check.grid(column=0, row=1, sticky=(W), padx=(0,10))
        self.Vm_neuron.set(False)
        
        check = Tk.Checkbutton(buttonFrame, text='spike times',  command=lambda:self.checked(self.spikeTimes), variable=self.spikeTimes,  onvalue=True, offvalue=False)
        check.grid(column=1, row=1, sticky=(W), padx=(0,10))
        self.spikeTimes.set(False)
        
        check = Tk.Checkbutton(buttonFrame, text='OP',           command=lambda:self.checked(self.OP),         variable=self.OP,          onvalue=True, offvalue=False)
        check.grid(column=0, row=2, sticky=(W))
        self.OP.set(False)
        
        check = Tk.Checkbutton(buttonFrame, text='OP NMDA',      command=lambda:self.checked(self.OP_NMDA),    variable=self.OP_NMDA,     onvalue=True, offvalue=False)
        check.grid(column=1, row=2, sticky=(W))
        self.OP_NMDA.set(False)
        
        check = Tk.Checkbutton(buttonFrame, text='g excit',      command=lambda:self.checked(self.g_excit),    variable=self.g_excit,     onvalue=True, offvalue=False)
        check.grid(column=0, row=3, sticky=(W))
        self.g_excit.set(False)
          
        check = Tk.Checkbutton(buttonFrame, text='g inhib',      command=lambda:self.checked(self.g_inhib),    variable=self.g_inhib,     onvalue=True, offvalue=False)
        check.grid(column=1, row=3, sticky=(W))
        self.g_inhib.set(False)
          
        check = Tk.Checkbutton(buttonFrame, text='BPAP',         command=lambda:self.checked(self.BPAP),       variable=self.BPAP,        onvalue=True, offvalue=False)
        check.grid(column=0, row=4, sticky=(W))
        self.BPAP.set(False)
        
        check = Tk.Checkbutton(buttonFrame, text='G ref',        command=lambda:self.checked(self.G_ref),      variable=self.G_ref,       onvalue=True, offvalue=False)
        check.grid(column=1, row=4, sticky=(W))
        self.G_ref.set(False)
        
        s = ttk.Separator(buttonFrame, orient=VERTICAL)
        s.grid(column=2, row=0,rowspan=5, sticky="snew", pady=(2,2), padx=(2,3))
        
        #create controls for synapse variables
        Tk.Label(buttonFrame, text='synapse: ').grid(column=3, row=0, sticky=(W))
        self.synapseInput = Tk.Entry(buttonFrame, textvariable=self.synapseNumber, width=14)
        self.synapseInput.grid(column=4, row=0, columnspan=2, sticky=(W))
        self.synapseNumber.set('0')
        self.synapseInput.insert(0,'0')
        self.synapseInput.bind('<Return>', self.renewCanvas)
        
        self.CB_Vm_syn = Tk.Checkbutton(buttonFrame, text='syn Vm',    command=lambda:self.checked(self.Vm_syn),     variable=self.Vm_syn,      onvalue=True, offvalue=False)
        self.CB_Vm_syn.grid(column=3, row=1, sticky=(W), padx=(0,10))
        self.Vm_syn.set(False)
        
        self.CB_g = Tk.Checkbutton(buttonFrame, text='g syn',               command=lambda:self.checked(self.g),          variable=self.g,           onvalue=True, offvalue=False)
        self.CB_g.grid(column=3, row=2, sticky=(W))
        self.g.set(False)
        
        self.CB_calcium = Tk.Checkbutton(buttonFrame, text='calcium',   command=lambda:self.checked(self.calcium),    variable=self.calcium,     onvalue=True, offvalue=False)
        self.CB_calcium.grid(column=3, row=3, sticky=(W))
        self.calcium.set(False)
        
        self.CB_Mg = Tk.Checkbutton(buttonFrame, text='Mg',             command=lambda:self.checked(self.Mg),         variable=self.Mg,          onvalue=True, offvalue=False)
        self.CB_Mg.grid(column=3, row=4, sticky=(W))  
        self.Mg.set(False)
        
        self.CB_I_syn = Tk.Checkbutton(buttonFrame, text='I syn',      command=lambda:self.checked(self.I_syn),      variable=self.I_syn,       onvalue=True, offvalue=False)
        self.CB_I_syn.grid(column=4, row=1, sticky=(W))
        self.I_syn.set(False)
          
        self.CB_I_NMDA = Tk.Checkbutton(buttonFrame, text='I NMDA',     command=lambda:self.checked(self.I_NMDA),     variable=self.I_NMDA,      onvalue=True, offvalue=False)
        self.CB_I_NMDA.grid(column=4, row=2, sticky=(W), padx=(0,10))
        self.I_NMDA.set(False)
          
        self.CB_I_VGCC = Tk.Checkbutton(buttonFrame, text='I VGCC',     command=lambda:self.checked(self.I_VGCC),     variable=self.I_VGCC,      onvalue=True, offvalue=False)
        self.CB_I_VGCC.grid(column=4, row=3, sticky=(W))
        self.I_VGCC.set(False)
        
        self.legendButton = Tk.Button(buttonFrame, text='info ON',     command=self.legendClicked, width=6)
        self.legendButton.grid(column=4, row=4, sticky=(W))

        self.CB_P = Tk.Checkbutton(buttonFrame, text='P',               command=lambda:self.checked(self.P),          variable=self.P,           onvalue=True, offvalue=False)
        self.CB_P.grid(column=5, row=1, sticky=(W), padx=(0,10))
        self.P.set(False)
        
        self.CB_B = Tk.Checkbutton(buttonFrame, text='B',               command=lambda:self.checked(self.B),          variable=self.B,           onvalue=True, offvalue=False)
        self.CB_B.grid(column=5, row=2, sticky=(W), padx=(0,10))
        self.B.set(False)
        
        self.CB_D = Tk.Checkbutton(buttonFrame, text='D',               command=lambda:self.checked(self.D),          variable=self.D,           onvalue=True, offvalue=False)
        self.CB_D.grid(column=5, row=3, sticky=(W), padx=(0,10))
        self.D.set(False)
        
        self.CB_V = Tk.Checkbutton(buttonFrame, text='V',               command=lambda:self.checked(self.V),          variable=self.V,           onvalue=True, offvalue=False)
        self.CB_V.grid(column=5, row=4, sticky=(W), padx=(0,10))
        self.V.set(False)
        
        s = ttk.Separator(buttonFrame, orient=VERTICAL)
        s.grid(column=6, row=0,rowspan=5, sticky="snew", pady=(2,2), padx=(2,3))
        
        #create controls for time variables
        self.offButton = Tk.Button(buttonFrame,              command=self.disable, width=3)
        self.offButton.grid(column=7, row=0, rowspan=2, sticky=(N,S))
        
        self.y_upButton = Tk.Button(buttonFrame, text='^',               command=lambda:self.move('y_up'), width=3)
        self.y_upButton.grid(column=8, row=0)
        
        self.y_downButton = Tk.Button(buttonFrame, text='v',             command=lambda:self.move('y_down'), width=3)
        self.y_downButton.grid(column=8, row=1)
        
        self.ymaxInput = Tk.Entry(buttonFrame, textvariable=self.y_max, width=8)
        self.ymaxInput.grid(column=9, row=0, columnspan=2, sticky=(E,W))
        self.y_max.set('80')
        self.ymaxInput.insert(0,self.y_max.get())
        self.ymaxInput.bind('<Return>', self.moveReturn)

        self.yminInput = Tk.Entry(buttonFrame, textvariable=self.y_min, width=8)
        self.yminInput.grid(column=9, row=1, columnspan=2, sticky=(E,W))
        self.y_min.set('-100')
        self.yminInput.insert(0,self.y_min.get())
        self.yminInput.bind('<Return>', self.moveReturn)

        s = ttk.Separator(buttonFrame, orient=HORIZONTAL)
        s.grid(column=7, row=2,columnspan=5, sticky="snew", pady=(2,2), padx=(2,2))
        
        self.xminInput = Tk.Entry(buttonFrame, textvariable=self.x_min, width=12)
        self.xminInput.grid(column=7, row=3, columnspan=2, sticky=(E,W))
        self.x_min.set('0')
        self.xminInput.insert(0,self.x_min.get())
        self.xminInput.bind('<Return>', self.moveReturn)

        self.xmaxInput = Tk.Entry(buttonFrame, textvariable=self.x_max, width=12)
        self.xmaxInput.grid(column=9, row=3, columnspan=2, sticky=(E,W))
        self.x_max.set('200')
        self.xmaxInput.insert(0,self.x_max.get())
        self.xmaxInput.bind('<Return>', self.moveReturn)
        
        button = Tk.Button(buttonFrame, text='<',                command=lambda:self.move('x_left'), width=6)
        button.grid(column=7, row=4, columnspan=2, sticky=(E,W))
        
        button = Tk.Button(buttonFrame, text='>',                command=lambda:self.move('x_right'), width=6)
        button.grid(column=9, row=4, columnspan=2, sticky=(E,W))
        
        self.offValue = not self.offValue
        self.disable()
        self.renewCanvas(False)
        
    def legendClicked(self):
        self.legend = not self.legend
        if self.legend: self.legendButton.config(text='info OFF')
        else: self.legendButton.config(text='info ON')
        self.renewCanvas(False)
    
    def checked(self, var):
        var.set(not var.get())
        self.renewCanvas(False)
        
    def disable(self):
        self.offValue = not self.offValue
        if self.offValue:
            self.offButton.config(text='ON')
            self.y_upButton.config(state=DISABLED)
            self.y_downButton.config(state=DISABLED)
            self.yminInput.config(state=DISABLED)
            self.ymaxInput.config(state=DISABLED)
        else:
            self.offButton.config(text='OFF')
            self.y_upButton.config(state=NORMAL)
            self.y_downButton.config(state=NORMAL)
            self.yminInput.config(state=NORMAL)
            self.ymaxInput.config(state=NORMAL)
            
            self.yminInput.delete(0, END)
            self.yminInput.insert(0,self.y_min.get())
            
            self.ymaxInput.delete(0, END)
            self.ymaxInput.insert(0,self.y_max.get())
        self.renewCanvas(False)
                    
    def moveReturn(self, event):
        self.y_min.set(self.yminInput.get())
        self.y_max.set(self.ymaxInput.get())
        self.x_min.set(self.xminInput.get())
        self.x_max.set(self.xmaxInput.get())
        self.renewCanvas(False)
        
    def move(self, moveDir):
        if moveDir == 'y_up':
            self.y_min.set(str(int(self.y_min.get())+self.step_y))
            self.y_max.set(str(int(self.y_max.get())+self.step_y))
        if moveDir == 'y_down':
            self.y_min.set(str(int(self.y_min.get())-self.step_y))
            self.y_max.set(str(int(self.y_max.get())-self.step_y))
        if moveDir == 'x_left':
            self.x_min.set(str(int(self.x_min.get())-self.step_x))
            self.x_max.set(str(int(self.x_max.get())-self.step_x))
        if moveDir == 'x_right':
            self.x_min.set(str(int(self.x_min.get())+self.step_x))
            self.x_max.set(str(int(self.x_max.get())+self.step_x))
        
        self.yminInput.delete(0, END)
        self.yminInput.insert(0,self.y_min.get())
        
        self.ymaxInput.delete(0, END)
        self.ymaxInput.insert(0,self.y_max.get())
        
        self.xminInput.delete(0, END)
        self.xminInput.insert(0,self.x_min.get())
        
        self.xmaxInput.delete(0, END)
        self.xmaxInput.insert(0,self.x_max.get())
        
        self.renewCanvas(False)
                  
    def renewCanvas(self, event):
        #clear previous figure
        self.mainFig.clf()
        self.plt = self.mainFig.add_subplot(111)
        
        #get input variables
        self.neuronNumber.set(self.neuronInput.get())
        self.synapseNumber.set(self.synapseInput.get())
        
        n, s = self.checkAvail(int(self.neuronNumber.get()), int(self.synapseNumber.get()))
        
        #plot variables
        p1=p2=p3=p4=p5=p6=p7=p8=p9=p10=p11=p12=p13=p14=p15=p16=p17=p18=p19=[]
        spikes = np.argwhere(self.root.neurons['RS'].spikeTimes[n,:])
        #neuron variables:
        if self.Vm_neuron.get():  p1, = self.plt.plot(self.root.genParam['timeArray'],  self.root.neurons['RS'].Vm[n,:],            'k')
        if self.spikeTimes.get(): p2, = self.plt.plot(spikes*self.root.genParam['dt'],  10*np.ones(len(spikes)),                    'kd', markersize=4)
        if self.OP.get():         p3, = self.plt.plot(self.root.genParam['timeArray'],  self.root.neurons['RS'].OP[n,:]*40,         'r.')
        if self.OP_NMDA.get():    p4, = self.plt.plot(self.root.genParam['timeArray'],  self.root.neurons['RS'].OP_NMDA[n,:]*40,    'r.')
        if self.g_excit.get():    p5, = self.plt.plot(self.root.genParam['timeArray'],  self.root.neurons['RS'].g_excit[n,:],       'r')
        if self.g_inhib.get():    p6, = self.plt.plot(self.root.genParam['timeArray'],  self.root.neurons['RS'].g_inhib[n,:],       'b')
        if self.BPAP.get():       p7, = self.plt.plot(self.root.genParam['timeArray'],  self.root.neurons['RS'].BPAP[n,:]*40,       'k')
        if self.G_ref.get():      p8, = self.plt.plot(self.root.genParam['timeArray'],  self.root.neurons['RS'].Gref[n,:]*8,        'k')
        
        #synapse variables
        if self.neuronAvail:
            if self.Vm_syn.get():     p9,  = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].Vm[s,:],            'k')
            if self.g.get():          p10, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].g[s,:]*100,         'k')
            if self.calcium.get():    p11, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].calcium[s,:]*8,     'k')
            if self.Mg.get():         p12, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].Mg[s,:]*50,         'k')
            if self.I_syn.get():      p13, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].I[s,:]*5,           'k')
            if self.I_NMDA.get():     p14, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].I_NMDA[s,:]*50,     'k')
            if self.I_VGCC.get():     p15, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].I_VGCC[s,:]*50,     'k')
            if self.P.get():          p16, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].P[s,:]*1000,        'k')
            if self.B.get():          p17, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].B[s,:]*10,          'k')
            if self.D.get():          p18, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].D[s,:]*500,         'k')
            if self.V.get():          p19, = self.plt.plot(self.root.genParam['timeArray'],  self.root.synapses[str(n)].V[s,:]*30,          'k')
        
        #plot parameters
        self.plt.set_xlabel('time (ms)')
        self.plt.set_xlim(int(self.x_min.get()),int(self.x_max.get()))
        if not self.offValue: self.plt.set_ylim(int(self.y_min.get()),int(self.y_max.get()))
        #make plot legend
        legendPlots = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19])
        legendNames = np.array(['Vm', 'spikes', 'OP', 'OP NDMA', 'g excit', 'g inhib', 'BPAP', 'G ref', 'Vm syn', 'g', 'calcium', 'Mg',
                                'I syn', 'I NMDA', 'I VGCC', 'P', 'V', 'B', 'D'])
        legendMask = np.array(map(bool,[self.Vm_neuron.get(), self.spikeTimes.get(), self.OP.get(), self.OP_NMDA.get(), self.g_excit.get(), self.g_inhib.get(), 
                      self.BPAP.get(), self.G_ref.get(), self.Vm_syn.get(), self.g.get(), self.calcium.get(), self.Mg.get(), self.I_syn.get(), 
                      self.I_NMDA.get(), self.I_VGCC.get(), self.P.get(), self.V.get(), self.B.get(), self.D.get()]))
        if self.legend: self.mainFig.legend(legendPlots[legendMask], legendNames[legendMask], mode="expand", ncol=6)
        
        self.canvas.show()

    def checkAvail(self, n, s):
        #check for out-of-bound neuron
        if n>=np.size(self.root.neurons['RS'].g_excit,0):
            n = np.size(self.root.neurons['RS'].g_excit,0)-1
            self.neuronNumber.set(str(n))
            self.neuronInput.delete(0, END)
            self.neuronInput.config(foreground='red')
            self.neuronInput.insert(0,str(n))
            self.mainFig.text(0.5,0.5,'neuron out-of-bound',bbox=dict(facecolor='red', alpha=0.8), horizontalalignment='center', verticalalignment='center')
            self.canvas.show()
        else: self.neuronInput.config(foreground='black')
        
        #check for out-of-bound neuron
        if s>=self.root.neurons['TC'].size:
            s = self.root.neurons['TC'].size-1
            self.synapseNumber.set(str(s))
            self.synapseInput.delete(0, END)
            self.synapseInput.config(foreground='red')
            self.synapseInput.insert(0,str(s))
            self.mainFig.text(0.5,0.5,'synapse out-of-bound',bbox=dict(facecolor='red', alpha=0.8), horizontalalignment='center', verticalalignment='center')
            self.canvas.show()
        else: self.synapseInput.config(foreground='black')
        
        #check for neuron with unavailable synapse information
        if self.neuronNumber.get() not in self.root.genParam['neuronsToTrack']:
            tt = 'synapse information not available; use neurons '
            for k in self.root.genParam['neuronsToTrack']: tt += k + ' '
            self.mainFig.text(0.02,0.95,tt,bbox=dict(facecolor='red', alpha=0.35))
            self.canvas.show()
            self.neuronAvail = False
            self.CB_Vm_syn.config(state=DISABLED)
            self.CB_g.config(state=DISABLED)
            self.CB_calcium.config(state=DISABLED)
            self.CB_Mg.config(state=DISABLED)
            self.CB_I_syn.config(state=DISABLED)
            self.CB_I_NMDA.config(state=DISABLED)
            self.CB_I_VGCC.config(state=DISABLED)
            self.CB_P.config(state=DISABLED)
            self.CB_B.config(state=DISABLED)
            self.CB_D.config(state=DISABLED)
            self.CB_V.config(state=DISABLED)
            self.synapseInput.config(state=DISABLED)
        else: 
            self.neuronAvail = True
            self.CB_Vm_syn.config(state=NORMAL)
            self.CB_g.config(state=NORMAL)
            self.CB_calcium.config(state=NORMAL)
            self.CB_Mg.config(state=NORMAL)
            self.CB_I_syn.config(state=NORMAL)
            self.CB_I_NMDA.config(state=NORMAL)
            self.CB_I_VGCC.config(state=NORMAL)
            self.CB_P.config(state=NORMAL)
            self.CB_B.config(state=NORMAL)
            self.CB_D.config(state=NORMAL)
            self.CB_V.config(state=NORMAL)
            self.synapseInput.config(state=NORMAL)
            
        return n, s

    def closeAll(self):
        self.destroy()   
