import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlow
import math
import sys
import pprint
import numpy as np
import scipy
from scipy.optimize import linprog

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InitializeGeometryProcess(Model, settings["Parameters"])
def in_hull(points, x):
    n_points = len(points)
    n_dim = len(x)
    c = np.zeros(n_points)
    A = np.r_[points.T,np.ones((1,n_points))]
    b = np.r_[x, np.ones(1)]
    lp = linprog(c, A_eq=A, b_eq=b)
    return lp.success
def inside(points, point):
    x=point[0]
    y=point[1]
    counter=0
    for i in range(0,len(points)-1):
        x_geometry1=points[i,0]
        y_geometry1=points[i,1]
        x_geometry2=points[i+1,0]
        y_geometry2=points[i+1,1]        
        if x_geometry1>=x or x_geometry2>=x:
            check_y_sign=(y_geometry1-y)*(y_geometry2-y)
            if check_y_sign<=0:
                counter += 1
    if counter % 2 == 0:
        is_inside=False 
    else:
        is_inside=True  

    return is_inside
def rotate(origin, point_array, angleInDegrees):   
	ox, oy = origin
	px=point_array[:,0]
	py=point_array[:,1]
	angle=math.radians(angleInDegrees)
	qx = ox+math.cos(angle)*(px - ox)-math.sin(angle)*(py - oy)
	qy = oy+math.sin(angle)*(px - ox)+math.cos(angle)*(py - oy)
	point_array[:,0]=qx
	point_array[:,1]=qy
	return point_array

## All the processes python should be derived from "Process"
class InitializeGeometryProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        print("Initialize Geometry process")
        default_parameters = KratosMultiphysics.Parameters( """
            {   
                "model_part_name": "insert_model_part",
                "initial_case": "circle",
                "geometry_parameter": 1.0,
                "initial_point": [0,0]
            }  """ );

        settings.ValidateAndAssignDefaults(default_parameters)
        self.model_part = Model.GetModelPart(settings["model_part_name"].GetString())
        self.initial_case=settings["initial_case"].GetString()
        self.geometry_parameter = settings["geometry_parameter"].GetDouble();        
        self.initial_point = KratosMultiphysics.Vector(2)
        self.initial_point[0] = settings["initial_point"][0].GetDouble()
        self.initial_point[1] = settings["initial_point"][1].GetDouble()

    def Execute(self):
        print("Executing Initialize Geometry")
        # points=np.loadtxt('n0012.dat')
        # levelset=CompressiblePotentialFlow.InitializeLevelSetProcess(self.model_part,points)
        # levelset.Execute()
        if self.initial_case == "0012":
            multiplier = 1
            angle=self.geometry_parameter
            points = np.loadtxt('n0012.dat')*multiplier
            points = rotate([0, 0], points, angle)
            points[:, 0] = points[:, 0]+self.initial_point[0]
            points[:, 1] = points[:, 1]+self.initial_point[1]
            for node in self.model_part.Nodes:
                point = [node.X, node.Y]
                min_dist=1000
                for i in range(0,len(points)):
                    dist=math.sqrt((node.X-points[i,0])**2+(node.Y-points[i,1])**2)
                    if dist<min_dist:
                        min_dist=dist
                if abs(point[0])<2 and abs(point[1])<2: 
                    if not in_hull(points, point):  # positive distance for fluid
                        node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,min_dist)
                    else:  # negative distance for solid
                        node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,-min_dist)
                else:  # positive distance for fluid
                        node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,min_dist)
        elif self.initial_case == "0010":
            T=0.1
            a0=0.2969
            a1=-0.126
            a2=-0.3516
            a3=0.2843
            a4=-0.1036        
            for node in self.model_part.Nodes:
                x=node.X 
                y=node.Y
                if x>=0:
                    if y>=0:
                        in_airfoil=T/0.2*(math.sqrt(x)*a0+a1*x+a2*x**2+a3*x**3+a4*x**4)
                        node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,y-in_airfoil)
                    else:
                        in_airfoil=-T/0.2*(math.sqrt(x)*a0+a1*x+a2*x**2+a3*x**3+a4*x**4)
                        node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,in_airfoil-y)
                else:
                    node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,1)
        elif self.initial_case == "circle":
            radius=self.geometry_parameter
            for node in self.model_part.Nodes:
                in_circle=(node.X-self.initial_point[0])**2+(node.Y-self.initial_point[1])**2
                if abs(in_circle)==0.0:
                    in_circle=-1e-5+1
                node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,in_circle-radius**2)
        elif self.initial_case == "ellipse":
            angle=math.radians(self.geometry_parameter)
            a=1
            b=0.2*a
            for node in self.model_part.Nodes:
                diffx=node.X-self.initial_point[0]
                diffy=node.Y-self.initial_point[1]
                in_ellipse=((diffx*math.cos(angle)+diffy*math.sin(angle))/a)**2+((diffx*math.sin(angle)+diffy*math.cos(angle))/b)**2
                if abs(in_ellipse-1)==0.0:
                    in_ellipse=-1e-5+1
                node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,in_ellipse-1)
        else:
            raise Exception("Initial geometry case not added")

        self.ApplyFlags()
        print("Level Set geometry initialized")
        
        
    def ExecuteInitialize(self):
        self.Execute()
        

    def ApplyFlags(self):

        for element in self.model_part.Elements:
            IsPositive=False
            IsNegative=False
            for node in element.GetNodes():
                distance=node.GetSolutionStepValue(KratosMultiphysics.NODAL_H)
                if distance>0:
                    IsPositive=True
                else:
                    IsNegative=True           
            if IsPositive and IsNegative:
                element.Set(KratosMultiphysics.BOUNDARY,True)
            elif IsPositive:
                element.Set(KratosMultiphysics.FLUID,True)
            else:
                element.Set(KratosMultiphysics.FLUID,False)
                element.Set(KratosMultiphysics.ACTIVE,False)
        # for element in self.model_part.Elements:
        #     IsPositive=False
        #     IsNegative=False
        #     for node in element.GetNodes():
        #         distance=node.GetSolutionStepValue(KratosMultiphysics.NODAL_H)
        #         if distance>0:
        #             IsPositive=True
        #         else:
        #             IsNegative=True     

        #     if not IsPositive and IsNegative:
        #         element.Set(KratosMultiphysics.FLUID,False)
        #         element.Set(KratosMultiphysics.ACTIVE,False)
                # for node in element.GetNodes():
                #     if node.IsNot(KratosMultiphysics.VISITED):
                #         node.Set(KratosMultiphysics.VISITED,True)
                        # node.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0)
                        # node.Fix(KratosMultiphysics.POSITIVE_FACE_PRESSURE)

        
