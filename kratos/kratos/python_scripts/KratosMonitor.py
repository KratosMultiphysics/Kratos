import Tkinter 


class KratosMonitor:

    def DeleteRoot(self):
        self.root.destroy()
        print "Monitor Closed"

    def __init__(self):
        self.root = Tkinter.Tk()
        self.root.title("Kratos")
#        self.root.protocol("WM_DELETE_WINDOW", self.DeleteRoot())

        self.frame = Tkinter.Frame(self.root)
        self.frame.grid()

        self.Variables = dict()
        self.Values = dict()

    def __del__(self):
        self.DeleteRoot()

    def Set(self, VariableName, Value):
        if(self.root):
            if(VariableName in self.Variables):
                self.Values[VariableName].config(text=Value)
            else:
                new_label = Tkinter.Label(self.frame, text=VariableName, anchor=Tkinter.W, justify=Tkinter.LEFT)
                new_label.grid(row=len(self.Variables), column=0, padx=8, sticky=Tkinter.W)
                new_label.config(anchor=Tkinter.W, justify=Tkinter.LEFT)

                self.Variables[VariableName] = new_label;
            
                new_label = Tkinter.Label(self.frame, text=Value, anchor=Tkinter.W, justify=Tkinter.LEFT)
                new_label.grid(row=len(self.Values), column=1, padx=8, sticky=Tkinter.W)
                self.Values[VariableName] = new_label;


 

kratos_monitor = KratosMonitor()


#kratos_monitor.root.mainloop()


#######################This part is for testing

#kratos_monitor.Set("Project", "test")

#from time import gmtime, strftime, mktime, sleep
#start_time = gmtime()
#kratos_monitor.Set("Started at", strftime("%a, %d %b %Y %H:%M:%S ", start_time))
#sleep(4)
