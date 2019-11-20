from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from functools import reduce

import math
import datetime
import shutil

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import KratosMultiphysics.DEMApplication.DEM_material_test_script as DEM_material_test_script

class SPTest(DEM_material_test_script.MaterialTest):

  def __init__(self, DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part):
      super(SPTest,self).__init__(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

  def Initialize(self):
      super(SPTest,self).Initialize()
      self.PrepareTestSP()

  def PrepareTestSP(self):
    if self.test_type == "SandP":

        ####### Correction Coefs  TODO 0.25* for cylinder section EXXON
        self.alpha_top = 3.141592*self.diameter*self.diameter*0.25/(self.xtop_area + 0.70710678*self.xtopcorner_area)
        self.alpha_bot = 3.141592*self.diameter*self.diameter*0.25/(self.xbot_area + 0.70710678*self.xbotcorner_area)
        self.alpha_lat = 3.141592*self.diameter*self.height/(self.xlat_area + 0.70710678*self.xtopcorner_area + 0.70710678*self.xbotcorner_area)

  def MeasureForcesAndPressure(self):
    super(SPTest,self).MeasureForcesAndPressure()

    total_force_top = 0.0
    z_tensor_value = 0.0

    for node in self.top_mesh_nodes:

      force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]

      total_force_top += force_node_y

    self.total_stress_top = total_force_top/(self.MeasuringSurface)


    if ( (self.test_type == "SandP") and (self.ConfinementPressure != 0.0) ):

      # self.Pressure = min(self.total_stress, self.ConfinementPressure * 1e6)
      # self.aux = AuxiliaryUtilities() from material test init
      self.aux.ComputeAverageZStressFor2D(self.spheres_model_part, z_tensor_value)

      # aixo estara immposat desde stage
      # self.spheres_model_part.ProcessInfo.SetValue(IMPOSED_Z_STRAIN_VALUE, -1e-6 * self.time)

      self.customapply(z_tensor_value, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)






  def customapply(self, Pressure, XLAT, XBOT, XTOP, XBOTCORNER, XTOPCORNER, alpha_top, alpha_bot, alpha_lat):

      for node in XLAT:

          r = node.GetSolutionStepValue(RADIUS)
          x = node.X
          y = node.Y
          z = node.Z

          values = Array3()
          vect = Array3()

          cross_section = 3.141592 * r * r

          # vector normal al centre:
          vect_moduli = math.sqrt(x * x + z * z)

          if(vect_moduli > 0.0):
              vect[0] = -x / vect_moduli
              vect[1] = 0
              vect[2] = -z / vect_moduli

          values[0] = cross_section * alpha_lat * Pressure * vect[0]
          values[1] = 0.0
          values[2] = cross_section * alpha_lat * Pressure * vect[2]

          node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)


  def Flush(self,a):
      pass
