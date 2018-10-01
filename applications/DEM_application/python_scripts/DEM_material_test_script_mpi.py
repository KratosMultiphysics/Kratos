from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from functools import reduce

import math
import datetime
import shutil

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.mpi import *

import DEM_material_test_script

class MaterialTest(DEM_material_test_script.MaterialTest):

  def __init__(self, DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, RigidFace_model_part):
      super(MaterialTest,self).__init__(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, RigidFace_model_part)

  def Initialize(self):
      super(MaterialTest,self).Initialize()

  def Flush(self,a):
      pass

  def PrepareTests(self):
      ##Fixing horizontally top and bot
      if(self.parameters.TestType != "BTS"):

        for node in self.TOP:

          node.SetSolutionStepValue(VELOCITY_X, 0.0);
          node.SetSolutionStepValue(VELOCITY_Z, 0.0);
          node.Fix(VELOCITY_X);
          node.Fix(VELOCITY_Z);

        for node in self.BOT:

          node.SetSolutionStepValue(VELOCITY_X, 0.0);
          node.SetSolutionStepValue(VELOCITY_Z, 0.0);
          node.Fix(VELOCITY_X);
          node.Fix(VELOCITY_Z);

      if(self.parameters.TestType == "BTS"):

        self.bts_export = open(self.parameters.problem_name + "_bts" + ".grf", 'w');
        self.BtsSkinDetermination()

      else:

        self.graph_export = open(self.parameters.problem_name +"_graph.grf", 'w')
        self.graph_export_1 = open(self.parameters.problem_name +"_graph_top.grf", 'w')
        self.graph_export_2 = open(self.parameters.problem_name +"_graph_bot.grf", 'w')

        if( self.parameters.TestType =="Hydrostatic"):
          self.graph_export_volumetric = open(self.parameters.problem_name+"_graph_VOL.grf",'w')

        print ('Initial Height of the Model: ' + str(self.height)+'\n')

        #if(self.parameters.PredefinedSkinOption == "ON" ):
        #  print ("ERROR: in Concrete Test Option the Skin is automatically predefined. Switch the Predefined Skin Option OFF")

        (xtop_area,xbot_area,xlat_area,xtopcorner_area,xbotcorner_area,y_top_total,weight_top, y_bot_total, weight_bot) = self.CylinderSkinDetermination()

        xtop_area_gath        = mpi.allgather(mpi.world, xtop_area)
        xbot_area_gath        = mpi.allgather(mpi.world, xbot_area)
        xlat_area_gath        = mpi.allgather(mpi.world, xlat_area)
        xtopcorner_area_gath  = mpi.allgather(mpi.world, xtopcorner_area)
        xbotcorner_area_gath  = mpi.allgather(mpi.world, xbotcorner_area)

        xtop_area = reduce(lambda x, y: x + y, xtop_area_gath)
        xbot_area = reduce(lambda x, y: x + y, xbot_area_gath)
        xlat_area = reduce(lambda x, y: x + y, xlat_area_gath)
        xtopcorner_area = reduce(lambda x, y: x + y, xtopcorner_area_gath)
        xbotcorner_area = reduce(lambda x, y: x + y, xbotcorner_area_gath)

        weight_top_gath = mpi.allgather(mpi.world, weight_top)
        weight_bot_gath = mpi.allgather(mpi.world, weight_bot)
        y_top_total_gath = mpi.allgather(mpi.world, y_top_total)
        y_bot_total_gath = mpi.allgather(mpi.world, y_bot_total)

        weight_top = reduce(lambda x, y: x + y, weight_top_gath)
        weight_bot = reduce(lambda x, y: x + y, weight_bot_gath)
        y_top_total = reduce(lambda x, y: x + y, y_top_total_gath)
        y_bot_total = reduce(lambda x, y: x + y, y_bot_total_gath)

        initial_height_top = y_top_total/weight_top
        initial_height_bot = y_bot_total/weight_bot

        inner_initial_height = initial_height_top - initial_height_bot

        specimen_length = self.parameters.SpecimenLength
        extended_length = specimen_length + (specimen_length - inner_initial_height)

        self.length_correction_factor = specimen_length/extended_length

  def PrepareDataForGraph(self):

    prepare_check = [0,0,0,0]
    prepare_check_gath = [0,0,0,0]
    self.total_check = 0

    for smp in self.RigidFace_model_part.SubModelParts:
        if smp[TOP]:
            self.top_mesh_nodes = smp.Nodes
            prepare_check[0] = 1
        if smp[BOTTOM]:
            self.bot_mesh_nodes = smp.Nodes
            prepare_check[1] = 1

    for smp in self.spheres_model_part.SubModelParts:
        if smp[TOP]:
            self.top_mesh_nodes = smp.Nodes
            prepare_check[2] = -1

        if smp[BOTTOM]:
            self.bot_mesh_nodes = smp.Nodes
            prepare_check[3] = -1

    prepare_check_gath[0] = mpi.gather(mpi.world,prepare_check[0],0)
    prepare_check_gath[1] = mpi.gather(mpi.world,prepare_check[1],0)
    prepare_check_gath[2] = mpi.gather(mpi.world,prepare_check[2],0)
    prepare_check_gath[3] = mpi.gather(mpi.world,prepare_check[3],0)

    if(mpi.rank == 0 ):
      prepare_check[0] = reduce(lambda x,y: max(x,y), prepare_check_gath[0])
      prepare_check[1] = reduce(lambda x,y: max(x,y), prepare_check_gath[1])
      prepare_check[2] = reduce(lambda x,y: min(x,y), prepare_check_gath[2])
      prepare_check[3] = reduce(lambda x,y: min(x,y), prepare_check_gath[3])

      for it in range(len(prepare_check)):

        self.total_check += prepare_check[it]

      if(math.fabs(self.total_check)!=2):

        print(" ERROR in the definition of TOP BOT groups. Both groups are required to be defined, they have to be either on FEM groups or in DEM groups")

  def MeasureForcesAndPressure(self):

    dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)

    #if(mpi.rank == 0 ):
    self.strain += -100*self.length_correction_factor*1.0*self.parameters.LoadingVelocityTop*dt/self.parameters.SpecimenLength

    if( self.parameters.TestType =="BTS"):

      total_force_bts = 0.0

      for node in self.top_mesh_nodes:

        force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
        total_force_bts += force_node_y

      if(mpi.rank == 0 ):
        total_force_bts_gather = mpi.gather(mpi.world, total_force_bts, 0)
        total_force_bts = reduce(lambda x, y: x + y, total_force_bts_gather)

        self.total_stress_bts = 2.0*total_force_bts/(3.14159*self.parameters.SpecimenLength*self.parameters.SpecimenDiameter*1e6)
        self.strain_bts += -100*2*self.parameters.LoadingVelocityTop*dt/self.parameters.SpecimenDiameter

    else:

      if( self.parameters.TestType =="Hydrostatic"):

        radial_strain = -100*self.MeasureRadialStrain()
        if(mpi.rank == 0):
          self.volumetric_strain = self.strain + 2.0*radial_strain

      total_force_top = 0.0
      total_force_bot = 0.0

      for node in self.top_mesh_nodes:

        force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]

        total_force_top += force_node_y

      total_force_top_gath = mpi.allgather(mpi.world, total_force_top)
      total_force_top = reduce(lambda x, y: x + y, total_force_top_gath)

      self.total_stress_top = total_force_top/(self.parameters.MeasuringSurface*1000000)

      for node in self.bot_mesh_nodes:

        force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]

        total_force_bot += force_node_y

      total_force_bot_gath = mpi.allgather(mpi.world, total_force_bot)
      total_force_bot = reduce(lambda x, y: x + y, total_force_bot_gath)

      self.total_stress_bot = total_force_bot/(self.parameters.MeasuringSurface*1000000)
      self.total_stress_mean = 0.5*(self.total_stress_bot + self.total_stress_top)

      if( ( (self.parameters.TestType == "Triaxial") or (self.parameters.TestType == "Hydrostatic") ) and (self.parameters.ConfinementPressure != 0.0) ):

          self.Pressure = min(self.total_stress_mean*1e6, self.parameters.ConfinementPressure * 1e6)

          if( self.parameters.TestType == "Hydrostatic"):

              self.Pressure = self.total_stress_mean*1e6

          self.ApplyLateralPressure(self.Pressure, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)

  def PrintGraph(self,step):
      if(mpi.rank == 0 ):
          super(MaterialTest,self).PrintGraph(step)

  def PrintChart(self):
      if(mpi.rank == 0 ):
          super(MaterialTest,self).PrintChart()

  def FinalizeGraphs(self):
      if(mpi.rank == 0):
          super(MaterialTest,self).FinalizeGraphs()

  def MeasureRadialStrain(self):

    mean_radial_strain = 0.0
    radial_strain = 0.0
    weight = 0.0

    for node in self.XLAT:

      r = node.GetSolutionStepValue(RADIUS)
      x = node.X
      z = node.Z

      x0 = node.X0
      z0 = node.Z0

      dist_initial = math.sqrt(x0 * x0 + z0 * z0)
      dist_now = math.sqrt(x * x + z * z)

      node_radial_strain = (dist_now - dist_initial) / dist_initial

      mean_radial_strain += node_radial_strain*r*r

      weight += r*r

    mean_radial_strain_gath = mpi.allgather(mpi.world,mean_radial_strain)
    weight_gath = mpi.allgather(mpi.world,weight)

    mean_radial_strain = reduce(lambda x, y: x + y, mean_radial_strain_gath)
    weight = reduce(lambda x, y: x + y, weight_gath)

    if(mpi.rank == 0 ):
      if(weight == 0.0):
        print ("Error in MeasureRadialStrain. Lateral skin particles not well defined")

    radial_strain = mean_radial_strain/weight

    return radial_strain
