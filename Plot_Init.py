'''
Created on Mar 4, 2014

@author: raphaelholca
'''

import sys
import matplotlib
matplotlib.use('GTK')
from Plot_Root import Plot_Root
from Plot_Traces import Plot_Traces

def main(argv=None):
    root = Plot_Root()
    root.mainloop()
#     traces = Plot_Traces()
#     traces.mainloop()
    
if __name__ == "__main__":
    sys.exit(main())