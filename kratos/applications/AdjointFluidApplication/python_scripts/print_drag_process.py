import KratosMultiphysics
#import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PrintDragProcess(Model, settings["Parameters"])

class PrintDragProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        
        self.drag_model_part = Model[settings["drag_model_part_name"].GetString()]
        self.drag_file_name = settings["drag_file_name"].GetString()

    def ExecuteInitialize(self):
        self.drag_file = open(self.drag_file_name, 'w')
        self.drag_file.write("#time         RX                     RY                     RZ\n")

    def ExecuteFinalizeSolutionStep(self):
        time = self.drag_model_part.ProcessInfo[KratosMultiphysics.TIME]
        dx = 0.0
        dy = 0.0
        dz = 0.0
        for node in self.drag_model_part.Nodes:
            reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)
            dx -= reaction[0]
            dy -= reaction[1]
            dz -= reaction[2]
        self.drag_file.write('{0:12.5e} {1:22.15e} {2:22.15e} {3:22.15e}\n'.format(time,dx,dy,dz))
        self.drag_file.flush()

    def ExecuteFinalize(self):
        self.drag_file.close()
