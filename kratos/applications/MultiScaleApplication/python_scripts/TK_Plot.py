from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *
import os
import glob
CheckForPreviousImport()


class GNUPlot_Folder:

	def __init__(self, Name):
		self.Name = Name
		
	def Initialize(self, Model):
		if not os.path.exists(self.Name):
			os.mkdir(self.Name)
		
	def OnBeforeSolutionStep(self, Model):
		pass
		
	def OnSolutionStepCompleted(self, Model):
		pass
		
	def Finalize(self, Model):
		pass


class GNUPlot_Nodal:

	def __init__(self,
				 BaseName=None,
				 Name=None,
				 VarX=None,
				 VarY=None,
				 NodeX=None,
				 NodeY=None,
				 ):
		self.BaseName = BaseName
		self.Name = Name
		self.VarX = VarX
		self.VarY = VarY
		self.NodeX = NodeX
		self.NodeY = NodeY
		self.ofile = None
		
	def Initialize(self, Model):
		if not os.path.exists(self.BaseName):
			os.mkdir(self.BaseName)
		filename = self.BaseName + os.sep + self.Name + "_" + str(self.NodeX) + "_" + str(self.NodeY) + ".grf"
		self.ofile = open(filename, "w+")
		xx = Model.Nodes[self.NodeX].GetSolutionStepValue(self.VarX)
		yy = Model.Nodes[self.NodeY].GetSolutionStepValue(self.VarY)
		self.ofile.write(str(xx) + "   " + str(yy) + '\n')
		self.ofile.flush()
	
	def OnBeforeSolutionStep(self, Model):
		pass
	
	def OnSolutionStepCompleted(self, Model):
		xx = Model.Nodes[self.NodeX].GetSolutionStepValue(self.VarX)
		yy = Model.Nodes[self.NodeY].GetSolutionStepValue(self.VarY)
		self.ofile.write(str(xx) + "   " + str(yy) + '\n')
		self.ofile.flush()
	
	def Finalize(self, Model):
		self.ofile.close()


class GNUPlot_Nodal_Multiple:

	def __init__(self,
				 BaseName=None,
				 Name=None,
				 VarX=None,
				 VarY=None,
				 NodesX=None,
				 NodesY=None,
				 ):
		self.BaseName = BaseName
		self.Name = Name
		self.VarX = VarX
		self.VarY = VarY
		self.NodesX = NodesX
		self.NodesY = NodesY
		self.ofile = None
		
	def Initialize(self, Model):
		if not os.path.exists(self.BaseName):
			os.mkdir(self.BaseName)
		filename = self.BaseName + os.sep + self.Name + ".grf"
		self.ofile = open(filename, "w+")
		sum_x = 0.0
		sum_y = 0.0
		for i in self.NodesX:
			sum_x += Model.Nodes[i].GetSolutionStepValue(self.VarX)
		for i in self.NodesY:
			sum_y += Model.Nodes[i].GetSolutionStepValue(self.VarY)
		self.ofile.write(str(sum_x) + "   " + str(sum_y) + '\n')
		self.ofile.flush()
	
	def OnBeforeSolutionStep(self, Model):
		pass
	
	def OnSolutionStepCompleted(self, Model):
		sum_x = 0.0
		sum_y = 0.0
		for i in self.NodesX:
			sum_x += Model.Nodes[i].GetSolutionStepValue(self.VarX)
		for i in self.NodesY:
			sum_y += Model.Nodes[i].GetSolutionStepValue(self.VarY)
		self.ofile.write(str(sum_x) + "   " + str(sum_y) + '\n')
		self.ofile.flush()
	
	def Finalize(self, Model):
		self.ofile.close()


class GNUPlot_Elemental:

	def __init__(self,
				 BaseName=None,
				 Name=None,
				 VarX=None,
				 VarY=None,
				 ComponentX=None,
				 ComponentY=None,
				 ElementID=None,
				 GpID=None,
				 ):
		self.BaseName = BaseName
		self.Name = Name
		self.VarX = VarX
		self.VarY = VarY
		self.ElementID = ElementID
		self.GpID = GpID
		self.ofile = None
		self.ComponentX = ComponentX
		self.ComponentY = ComponentY
		
	def Initialize(self, Model):
		if not os.path.exists(self.BaseName):
			os.mkdir(self.BaseName)
		filename = self.BaseName + os.sep + self.Name + "_" + str(self.ElementID) + "_" + str(self.GpID) + ".grf"
		self.ofile = open(filename, "w+")
		xx = 0.0
		yy = 0.0
		self.ofile.write(str(xx) + "   " + str(yy) + '\n')
		self.ofile.flush()
	
	def OnBeforeSolutionStep(self, Model):
		pass
	
	def OnSolutionStepCompleted(self, Model):
		allx = Model.Elements[self.ElementID].GetValuesOnIntegrationPoints(self.VarX, Model.ProcessInfo)
		ally = Model.Elements[self.ElementID].GetValuesOnIntegrationPoints(self.VarY, Model.ProcessInfo)
		mx = allx[self.GpID]
		my = ally[self.GpID]
		if self.ComponentX is None:
			xx = mx
		else:
			xx = mx[self.ComponentX]
		if self.ComponentY is None:
			yy = my
		else:
			yy = my[self.ComponentY]
		self.ofile.write(str(xx) + "   " + str(yy) + '\n')
		self.ofile.flush()
	
	def Finalize(self, Model):
		self.ofile.close()


class GNUPlot_YieldSurf2D:

	def __init__(self,
				 BaseName=None,
				 Name=None,
				 ElementID=None,
				 GpID=None,
				 ):
		self.BaseName = BaseName
		self.Name = Name
		self.ElementID = ElementID
		self.GpID = GpID
		
	def Initialize(self, Model):
		if not os.path.exists(self.BaseName):
			os.mkdir(self.BaseName)
	
	def OnBeforeSolutionStep(self, Model):
		pass
	
	def OnSolutionStepCompleted(self, Model):
		time = Model.ProcessInfo[TIME]
		stime = str(time).replace('.','_')
		filename = self.BaseName + os.sep + self.Name + "_" + stime + ".grf"
		ofile = open(filename, "w+")
		allx = Model.Elements[self.ElementID].GetValuesOnIntegrationPoints(YIELD_SURFACE_DATA_2D_X, Model.ProcessInfo)
		ally = Model.Elements[self.ElementID].GetValuesOnIntegrationPoints(YIELD_SURFACE_DATA_2D_Y, Model.ProcessInfo)
		xx = allx[self.GpID]
		yy = ally[self.GpID]
		for i in range(len(xx)):
			x = str(xx[i])
			y = str(yy[i])
			ofile.write(x + "   " + y + '\n')
		ofile.close()
		
	def Finalize(self, Model):
		pass
