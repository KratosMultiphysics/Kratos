#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()

import matplotlib
import collections

from numpy import *
from pylab import *


class GraphPlotUtility:
    #######################################################################
    def __init__(self,model_part,problem_path,plot_active,plot_frequency):

        self.model_part   = model_part
        self.problem_path = problem_path
        
        self.plot_active = False
        if(plot_active == "True"):
            self.plot_active = True
            
        self.plot_frequency = plot_frequency

        # graph variables
        self.mesh_id   = 0
        self.x_var     = "TIME"
        self.y_var     = "FORCE_CONTACT_NORMAL"
        self.plot_step = 0

        # graph data limits
        self.x_limit = 100
        self.y_limit = 100
        self.x_dir = 0
        self.y_dir = 0

        # graph data containers
        self.X   = []
        self.Y   = []
        self.Y_x = [] 
        self.Y_y = [] 
        

    #######################################################################
    def Initialize(self,initial_plot_step):
        clf()
        self.X   = []
        self.Y   = []
        self.Y_x = []
        self.Y_y = []

        #set initial plot step
        self.plot_step   = initial_plot_step

        #print file headers
        if(initial_plot_step == 0):
            self.graph_path  = self.problem_path + "/ContactForce_vs_Displacement.cvs"
            graphfile  = open(self.graph_path,"w")
            linevalue  = "Displacement TotalForce ForceX ForceY\n"
            graphfile.write(linevalue)
            graphfile.close()


    #######################################################################
    def SetPlotVariables(self,x_variable,y_variable,mesh_id):
        # graph variables
        self.x_var   = x_variable
        self.y_var   = y_variable
        self.mesh_id  = mesh_id

        # graph name
        #self.graph_name = str(self.y_var)+"_vs_"+str(self.x_var)
        self.graph_name = str(self.y_var)+"_vs_DISPLACEMENT"
 
    #######################################################################   
    def SetStepResult(self,velocity):
        
        if(self.plot_active == True):

            if(self.x_var == "TIME"):
                X_value  = self.model_part.ProcessInfo[TIME]

                velocity = sqrt( velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]) 

                if(velocity != 0):
                    X_value  = X_value * velocity;

                Y_value = []
                Y_value_x = 0
                Y_value_y = 0

                initial_mesh  = 0
                number_meshes = 2 #numberofmeshes
            
                if(self.mesh_id > 0):
                    initial_mesh  = 1
                    number_meshes = 3
                
                for mesh in range(int(initial_mesh),int(number_meshes)):
                    #print "Mesh", mesh
                    if(mesh == self.mesh_id):
                        for node in self.model_part.GetNodes(mesh):
                            value   = node.GetSolutionStepValue(FORCE_CONTACT_NORMAL);
                            Y_value = Y_value + value
                            Y_value_x = Y_value_x + value[0]
                            Y_value_y = Y_value_y + value[1]
                            
                        self.X.append(X_value)
                        totalforce = sqrt( Y_value[0]*Y_value[0] + Y_value[1]*Y_value[1] + Y_value[2]*Y_value[2]) 
                        self.Y.append(totalforce)
                        self.Y_x.append(Y_value_x)
                        self.Y_y.append(Y_value_y)
                
                        forcepath = self.problem_path + "/ContactForce_vs_Time.cvs"
                        forcefile  = open(forcepath,"a")
                        linevalue = str(X_value) + " " + str(totalforce) + " " + str(Y_value_x) + " " + str(Y_value_y)+"\n"
                        forcefile.write(linevalue)
                        forcefile.close()
                
            
            #print " X ", self.X

    #######################################################################   
    def Plot(self,write_id):
        
        if(self.plot_active == True and self.plot_step == self.plot_frequency ):
            clf()
            plot(self.X,self.Y,'g-o',self.X,self.Y_x,'b-s', self.X,self.Y_y,'r-^')
            grid(True)
            title('TITLE: '+str(self.graph_name))
            
            #write_id = self.model_part.ProcessInfo[WRITE_ID]
            
            xlabel(str(self.x_var))
            ylabel(str(self.y_var))
            
            force_extra = 10
            time_extra  = 10 * self.model_part.ProcessInfo[DELTA_TIME]

            x_limits = self.SearchLimits(self.X,time_extra)
            y_limits = self.SearchLimits(self.Y,force_extra)
            y_xlimits = self.SearchLimits(self.Y_x,force_extra)
            y_ylimits = self.SearchLimits(self.Y_y,force_extra)
            
            y_max = y_limits.max
            y_min = y_limits.min
            if(y_max<y_xlimits.max):
                y_max = y_xlimits.max
            if(y_max<y_ylimits.max):
                y_max = y_ylimits.max

            if(y_min>y_xlimits.min):
                y_min = y_xlimits.min
            if(y_min>y_ylimits.min):
                y_min = y_ylimits.min
                         
            xlim(x_limits.min,x_limits.max)
            ylim(y_min,y_max)

            figure_name = os.path.join(self.problem_path,self.graph_name + "_" + str(write_id))
            #figure_name = str(self.problempath)+"/"+str(self.graph_name)+"_"+str(write_id)
            savefig(figure_name)

            self.plot_step = 0;
        else:
            self.plot_step = self.plot_step + 1;

    #######################################################################   
    def SearchLimits(self,X,extra):
        
        limits = collections.namedtuple('Limits', ['min', 'max'])
        limits.min = X[0]
        limits.max = X[0]
        for x in X:
          if(limits.min>x):
             limits.min=x
          if(limits.max<x):
             limits.max=x

        limits.min = limits.min - extra;
        limits.max = limits.max + extra;

        return limits
            
            
            
