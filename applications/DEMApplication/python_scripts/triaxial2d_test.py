import math
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as Dem

class Triaxial2D(Dem.DEM_material_test_script.MaterialTest):

    def __init__(self, DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part):
        super().__init__(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

    def Initialize(self):
        self.PrepareTests()
        self.PrepareTestTriaxial2D()

    def PrepareTestTriaxial2D(self):
        self.alpha_lat = self.perimeter/(self.xlat_area)

    def PrepareTests(self):
        ##Fixing vertical
        if self.test_type == "Triaxial2D":
            for element in self.spheres_model_part.Elements:
                node = element.GetNode(0)
                node.SetSolutionStepValue(Kratos.VELOCITY_Z, 0.0)
                node.Fix(Kratos.VELOCITY_Z)

        absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        self.graph_export   = open(absolute_path_to_file, 'w')
        self.xlat_area = self.CircularSkinDetermination()

    def CircularSkinDetermination(self):

        # Cylinder dimensions
        d = self.diameter
        eps = 2.0

        self.perimeter = 3.141592 * d
        self.xlat_area = 0.0

        for element in self.spheres_model_part.Elements:
            element.GetNode(0).SetSolutionStepValue(Dem.SKIN_SPHERE, 0)
            node = element.GetNode(0)
            r = node.GetSolutionStepValue(Kratos.RADIUS)
            x = node.X
            y = node.Y

            cross_section = 2.0 * r

            if ((x * x + y * y) >= ((d / 2 - eps * r) * (d / 2 - eps * r))):
                element.GetNode(0).SetSolutionStepValue(Dem.SKIN_SPHERE, 1)
                self.LAT.append(node)
                self.xlat_area = self.xlat_area + cross_section

        if len(self.LAT)==0:
            self.Procedures.KratosPrintWarning("ERROR! in Circular Skin Determination - NO LATERAL PARTICLES" + "\n")

        return self.xlat_area

    def ApplyLateralStress(self, average_zstress_value, LAT, alpha_lat):

        for node in LAT:
            r = node.GetSolutionStepValue(Kratos.RADIUS)
            x = node.X
            y = node.Y

            values = Kratos.Array3()
            vect = Kratos.Array3()

            cross_section = 2.0 * r

            vect_moduli = math.sqrt(x * x + y * y)

            if vect_moduli > 0.0:
                vect[0] = x / vect_moduli
                vect[1] = y / vect_moduli

            values[0] = cross_section * average_zstress_value * vect[0] * alpha_lat
            values[1] = cross_section * average_zstress_value * vect[1] * alpha_lat
            node.SetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE, values)

    def MeasureForcesAndPressure(self):
        super().MeasureForcesAndPressure()
        average_zstress_value = 0.0

        if self.test_type == "Triaxial2D":
            average_zstress_value = self.aux.ComputeAverageZStressFor2D(self.spheres_model_part)
            self.ApplyLateralStress(average_zstress_value, self.LAT, self.alpha_lat)

    def PrintGraph(self, time):
        pass
        #if self.graph_counter == self.graph_frequency:
        #    self.graph_counter = 0
        #    self.graph_export.write(str("%.6g"%self.strain).rjust(13) + "  " +
        #                            str("%.6g"%(self.total_stress_mean*1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12)+'\n')
        #    self.Flush(self.graph_export)

        #self.graph_counter += 1

    def Flush(self,a):
        pass
