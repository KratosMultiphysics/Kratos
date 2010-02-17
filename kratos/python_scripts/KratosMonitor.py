import Tkinter 

class KratosMonitor:
    def __init__(self, master):

        self.frame = Tkinter.Frame(master)
        self.frame.grid()

        self.Variables = dict()
        self.Values = dict()


    def Set(self, VariableName, Value):
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

 
root = Tkinter.Tk()
kratos_monitor = KratosMonitor(root)


#root.mainloop()

