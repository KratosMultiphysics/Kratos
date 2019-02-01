from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import RVEAnalysis


class TestPatchTestShells(KratosUnittest.TestCase):

    def test_isotropic_rve(self):
        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        boundary_mp_name = "Structure.SurfacePressure3D_Pressure_on_surfaces_Auto2"
        averaging_mp_name = "Structure.computing_domain"

        model = KratosMultiphysics.Model()
        simulation = RVEAnalysis.RVEAnalysis(model,parameters,boundary_mp_name,averaging_mp_name)
        simulation.Run()

        Cestimated = model["Structure.computing_domain"].GetValue(StructuralMechanicsApplication.EIGENVECTOR_MATRIX)

        Canalytic = KratosMultiphysics.Matrix(6,6)
        Canalytic.fill(0.0)
        E = 1e6
        nu = 0.3
        l = E*nu/((1+nu)*(1-2*nu))
        G = E/(2.0*(1.0+nu))
        Canalytic[0,0] = l+2*G
        Canalytic[0,1] = l
        Canalytic[0,2] = l

        Canalytic[1,0] = l
        Canalytic[1,1] = l+2*G
        Canalytic[1,2] = l

        Canalytic[2,0] = l
        Canalytic[2,1] = l
        Canalytic[2,2] = l+2*G

        Canalytic[3,3] = G
        Canalytic[4,4] = G
        Canalytic[5,5] = G

        for i in range(0,Cestimated.Size1()):
            for j in range(0,Cestimated.Size2()):
                print(i,j,Cestimated[i,j],Canalytic[i,j])
                self.assertAlmostEqual(abs(Cestimated[i,j] - Canalytic[i,j])/(l+2*G),0.0,5)

        print(Canalytic)
        print(Cestimated)




if __name__ == '__main__':
    KratosUnittest.main()
