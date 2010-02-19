import Tkinter 


class KratosMonitor:

    def DeleteRoot(self):
        self.toplevel.update()
        self.toplevel.destroy();
        self.root.destroy();
        self.toplevel = 0;
        print "Monitor Closed"

    def __init__(self):
        self.root = Tkinter.Tk()
        self.root.withdraw();
        self.toplevel = Tkinter.Toplevel( self.root);
        self.toplevel.title("Kratos")
        self.toplevel.protocol("WM_DELETE_WINDOW", self.DeleteRoot)
        self.frame = Tkinter.Frame( self.toplevel)
        self.frame.grid();
        self.toplevel.bind( "<Control-c>", exit);

        self.Variables = dict()
        self.Values = dict()

    def __del__(self):
        self.DeleteRoot()

    def Set(self, VariableName, Value):
        if(self.toplevel != 0):
            bg_color="white"
            if(VariableName in self.Variables):
                self.Values[VariableName].config(text=Value)
            else:
                new_label = Tkinter.Label(self.toplevel, text=VariableName, anchor=Tkinter.W, justify=Tkinter.LEFT)
                new_label.grid(row=len(self.Variables), column=0, padx=8, pady=1, sticky=Tkinter.N+Tkinter.W+Tkinter.E)
##                if ( len(self.Variables) % 2):
##                    new_label.config(anchor=Tkinter.W, justify=Tkinter.LEFT)
##                else:
##                    new_label.config(anchor=Tkinter.W, justify=Tkinter.LEFT, background=bg_color)

                self.Variables[VariableName] = new_label;
            
                new_label = Tkinter.Label(self.toplevel, text=Value, anchor=Tkinter.W, justify=Tkinter.LEFT)
                new_label.grid(row=len(self.Values), column=1, padx=8, pady=1, sticky=Tkinter.N+Tkinter.W+Tkinter.E)
##                if ( len(self.Values) % 2):
##                    new_label.config(anchor=Tkinter.W, justify=Tkinter.LEFT)
##                else:
                new_label.config(anchor=Tkinter.W, justify=Tkinter.LEFT, background=bg_color)
                self.Values[VariableName] = new_label;
            self.root.update()
        else:
            print "KratosMonitor : ", VariableName, " = ", Value


 

kratos_monitor = KratosMonitor()


# kratos_monitor.root.mainloop()
# kratos_monitor.root.update( )

#######################This part is for testing

#kratos_monitor.Set("Project", "test")

#from time import gmtime, strftime, mktime, sleep
#start_time = gmtime()
#kratos_monitor.Set("Started at", strftime("%a, %d %b %Y %H:%M:%S ", start_time))
#sleep(4)

#kratos_monitor.root.mainloop()
