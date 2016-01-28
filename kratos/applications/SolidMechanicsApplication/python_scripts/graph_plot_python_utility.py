from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()

#import matplotlib
import collections

from numpy import *
from pylab import *


class GraphPlotUtility:
    #

    def __init__(self, model_part, problem_path):

        self.model_part = model_part

        self.problem_path = problem_path

        # graph variables
        self.mesh_id = 0
        self.x_var = "DISPLACEMENT"
        self.y_var = "REACTION"

        # graph data limits
        self.x_limit = 100
        self.y_limit = 100
        self.x_dir = 0
        self.y_dir = 0

        # graph data containers
        self.Time = []

        self.X = []
        self.X_x = []
        self.X_y = []
        self.X_z = []

        self.Y = []
        self.Y_x = []
        self.Y_y = []
        self.Y_z = []

    #
    def Initialize(self, x_variable, y_variable, mesh_id):

        clf()

        self.SetPlotVariables(x_variable, y_variable, mesh_id)

        figure_path = os.path.join(self.problem_path, str(self.x_var) + "_vs_" + str(self.y_var) + ".post.csv")

        if(os.path.exists(figure_path) == False):
            # print file headers
            figure_file = open(figure_path, "w")
            line_header  = "Time," + str(self.x_var) + "," + str(self.y_var) + "," +str(self.x_var)+"_X,"+str(self.x_var)+"_Y,"+str(self.x_var)+"_Z,"+str(self.y_var)+"_X,"+str(self.y_var)+"_Y,"+str(self.y_var)+"_Z"+"\n"
            figure_file.write(line_header)
            line_header  = "0,0,0,0,0,0,0,0,0\n"
            figure_file.write(line_header)
            figure_file.close()

        time = self.GetStepTime()
        
        X_value = self.GetStepVariable(self.x_var)
        Y_value = self.GetStepVariable(self.y_var)

        X_value_norm = self.Get3DArrayModulus(X_value)
        Y_value_norm = self.Get3DArrayModulus(Y_value)

        self.Time.append(time)
       
        self.X.append(X_value_norm)
        self.X_x.append(X_value[0])
        self.X_y.append(X_value[1])
        self.X_z.append(X_value[2])

        self.Y.append(Y_value_norm)
        self.Y_x.append(Y_value[0])
        self.Y_y.append(Y_value[1])
        self.Y_z.append(Y_value[2])

    #
    def SetPlotVariables(self, x_variable, y_variable, mesh_id):
        # graph variables
        self.x_var = x_variable
        self.y_var = y_variable
        self.mesh_id = mesh_id

        # graph name
        self.graph_name = str(self.y_var) + "_vs_" + str(self.x_var)
        self.plot_name = str(self.y_var) + "_vs_TIME"

        # initialize variables
        self.InitializePlotVariables()

    #
    def InitializePlotVariables(self):

        self.Time = []

        self.X = []
        self.X_x = []
        self.X_y = []
        self.X_z = []

        self.Y = []
        self.Y_x = []
        self.Y_y = []
        self.Y_z = []

    #
    def GetStepTime(self):

        return self.model_part.ProcessInfo[TIME]

    #
    def GetStepDeltaTime(self):

        return self.model_part.ProcessInfo[DELTA_TIME]

    #
    def GetStepVariable(self, variable):

        variable_value = []
        for node in self.model_part.GetNodes(self.mesh_id):
            kratos_variable = globals()[variable]
            nodal_value = node.GetSolutionStepValue(kratos_variable);
            variable_value = variable_value + nodal_value
 
        return variable_value

    #
    def Get3DArrayModulus(self, variable):

        modulus = 0
        for var in variable:
            modulus = modulus + var * var

        return sqrt(modulus)

    #
    def SetStepResult(self):

        time = self.GetStepTime()

        X_value = self.GetStepVariable(self.x_var)
        Y_value = self.GetStepVariable(self.y_var)

        X_value_norm = self.Get3DArrayModulus(X_value)
        Y_value_norm = self.Get3DArrayModulus(Y_value)

        self.Time.append(time)

        self.X.append(X_value_norm)
        self.X_x.append(X_value[0])
        self.X_y.append(X_value[1])
        self.X_z.append(X_value[2])

        self.Y.append(Y_value_norm)
        self.Y_x.append(Y_value[0])
        self.Y_y.append(Y_value[1])
        self.Y_z.append(Y_value[2])

        figure_path = os.path.join(self.problem_path, str(self.x_var) + "_vs_" + str(self.y_var) + ".post.csv")
        figure_file = open(figure_path, "a")
        line_value = str(time) + "," + str(X_value_norm) + "," + str(Y_value_norm) + "," + str(X_value[0]) + "," + str(X_value[1]) + "," + str(X_value[2]) + "," + str(Y_value[0]) + "," + str(Y_value[1]) + "," + str(Y_value[2]) + "\n"
        figure_file.write(line_value)
        figure_file.close()

    #
    def Plot(self, write_id):

    #
    def PlotMP(self, write_id): #matplot needed
        clf()
        plot(self.Time, self.Y, 'g-o')
        # plot(self.Time,self.Y,'g-o',self.Time,self.Y_x,'b-s', self.Time,self.Y_y,'r-^')
        grid(True)
        title('TITLE: ' + str(self.plot_name))

        # xlabel(str(self.x_var))
        xlabel("TIME")
        ylabel(str(self.y_var))

        force_extra = 10
        time_extra = 10 * self.GetStepDeltaTime()

        x_limits = self.SearchLimits(self.Time, time_extra)
        # x_limits = self.SearchLimits(self.X,time_extra)
        y_limits = self.SearchLimits(self.Y, force_extra)
        y_xlimits = self.SearchLimits(self.Y_x, force_extra)
        y_ylimits = self.SearchLimits(self.Y_y, force_extra)

        y_max = y_limits.max
        y_min = y_limits.min
        if(y_max < y_xlimits.max):
            y_max = y_xlimits.max
        if(y_max < y_ylimits.max):
            y_max = y_ylimits.max

        if(y_min > y_xlimits.min):
            y_min = y_xlimits.min
        if(y_min > y_ylimits.min):
            y_min = y_ylimits.min

        xlim(x_limits.min, x_limits.max)
        ylim(y_min, y_max)

        figure_name = os.path.join(self.problem_path, self.plot_name + "_" + str(write_id))+ ".graph.png"

        savefig(figure_name)

    #
    def SearchLimits(self, X, extra):

        limits = collections.namedtuple('Limits', ['min', 'max'])
        limits.min = X[0]
        limits.max = X[0]
        for x in X:
            if(limits.min > x):
                limits.min = x
            if(limits.max < x):
                limits.max = x

        limits.min = limits.min - extra;
        limits.max = limits.max + extra;

        return limits
