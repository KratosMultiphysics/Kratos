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
    self.PrepareTests()
    self.DefineBoundaries()


  def PrepareTests(self):

    ##Fixing vertical
    if self.test_type == "SandP":
        for node in self.TOP:
            node.SetSolutionStepValue(VELOCITY_Z, 0.0)
            node.Fix(VELOCITY_Z)

    if self.test_type == "SandP":
        absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        self.graph_export   = open(absolute_path_to_file, 'w')

        self.Procedures.KratosPrintInfo('Initial Height of the Model: ' + str(self.height)+'\n')

        (self.xlat_area) = self.CircularSkinDetermination()


  def CircularSkinDetermination(self):

        # Cylinder dimensions
        d = self.diameter
        eps = 2.0

        perimeter = 3.141592 * d
        xlat_area = 0.0

        for element in self.spheres_model_part.Elements:

            element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 0)

            node = element.GetNode(0)
            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y

            cross_section = 2.0 * r

            if ((x * x + y * y) >= ((d / 2 - eps * r) * (d / 2 - eps * r))):

                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)
                self.LAT.append(node)
                xlat_area = xlat_area + cross_section

        if(len(self.LAT)==0):
            self.Procedures.KratosPrintWarning("ERROR! in Circular Skin Determination - NO LATERAL PARTICLES" + "\n")

        return (xlat_area)


  def MeasureForcesAndPressure(self):
    super(SPTest,self).MeasureForcesAndPressure()
    average_zstress_value = 0.0

    if self.test_type == "SandP":
      average_zstress_value = self.aux.ComputeAverageZStressFor2D(self.spheres_model_part)
      print (average_zstress_value)
      self.ApplyLateralStress(average_zstress_value, self.LAT)



  def ApplyLateralStress(self, average_zstress_value, LAT):

      for node in LAT:

          r = node.GetSolutionStepValue(RADIUS)
          x = node.X
          y = node.Y

          values = Array3()
          vect = Array3()

          cross_section = 2.0 * r

          # vector normal al centre:
          vect_moduli = math.sqrt(x * x + y * y)

          if(vect_moduli > 0.0):
              vect[0] = -x / vect_moduli
              vect[1] = -y / vect_moduli

          values[0] = cross_section * average_zstress_value * vect[0]
          values[1] = cross_section * average_zstress_value * vect[1]

          node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)


  def Flush(self,a):
      pass
