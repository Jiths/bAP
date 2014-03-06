'''
Created on Mar 4, 2014

@author: raphaelholca
'''

from Tkinter import *
import Tkinter as ttk
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np

def calculate(*args): #defines the calculation procedure; called with 'return' key of 'calculate' button
    try:
        value = float(feet.get()) #gets value from entry widget
        meters.set((0.3048 * value * 10000.0 + 0.5)/10000.0) #place calculate value into the label widget
    except ValueError:
        pass
    
root = Tk()
root.title("Feet to Meters")

mainframe = ttk.Frame(root) #create frame widget
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))#, padx=(100,100), pady=(100,100))
mainframe.columnconfigure(0, weight=1) #window resizing
mainframe.rowconfigure(0, weight=1) #window resizing

feet = StringVar() #this is a global variable; any time its value changes it will be automatically updated globally
meters = StringVar()

# mainplot = FigureCanvasTkAgg(f, master=root)
# mainplot.get_tk_widget().pack(side=ttk.TOP, fill=ttk.BOTH, expand=1)
# mainplot.show()

# mainplot = ttk.Canvas(mainframe)
# mainplot.grid(column=0, row=0)

feet_entry = ttk.Entry(mainframe, width=7, textvariable=feet) #create entry widget
feet_entry.grid(column=2, row=1, sticky=(W, E)) #place it on the screen

ttk.Label(mainframe, textvariable=meters).grid(column=2, row=2, sticky=(W, E)) #create label widget and place on the screen
ttk.Button(mainframe, text="Calculate", command=calculate).grid(column=3, row=3, sticky=W) #create button widget and place on the screen

ttk.Label(mainframe, text="feet").grid(column=3, row=1, sticky=W) #create label widget and place on the screen
ttk.Label(mainframe, text="is equivalent to").grid(column=1, row=2, sticky=E) #create label widget and place on the screen
ttk.Label(mainframe, text="meters").grid(column=3, row=2, sticky=W) #create label widget and place on the screen

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5) #makes sure the widget are not too close to each other
feet_entry.focus() #puts focus (cursor) on the feet_entry widget
root.bind('<Return>', calculate) #behavior of 'return' key

root.mainloop() #makes tk enter its event loop



        