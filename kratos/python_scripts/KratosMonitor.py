import Tkinter 

class KratosMonitor:
    def __init__(self, master):

        self.frame = Tkinter.Frame(master)
        self.frame.grid()

        self.Variables = dict()
        self.Values = dict()


    def SetVariable(self, VariableName, Value):
        if(VariableName in self.Variables):
            self.Values[VariableName].config(text=Value)
        else:
            new_label = Tkinter.Label(self.frame, text=VariableName, anchor=Tkinter.W, justify=Tkinter.LEFT)
            new_label.grid(row=len(self.Variables), column=0)
            self.Variables[VariableName] = new_label;
            
            new_label = Tkinter.Label(self.frame, text=Value, anchor=Tkinter.W, justify=Tkinter.LEFT)
            new_label.grid(row=len(self.Values), column=1)
            self.Values[VariableName] = new_label;

 
root = Tkinter.Tk()
monitor = KratosMonitor(root)

monitor.SetVariable("Time", 34)
monitor.SetVariable("Delt Time", 0.1)


#root.mainloop()

