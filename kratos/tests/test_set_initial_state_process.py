import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.set_initial_state_process as set_initial_state_process


class TestSetInitialStateProcess(KratosUnittest.TestCase):
    def test_set_initial_state_process(self):
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        # Set up model part
        model = KM.Model()
        modelpart = model.CreateModelPart("Main")
        node1 = modelpart.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = modelpart.CreateNewNode(2, 1.0, 0.0, 0.0)
        node3 = modelpart.CreateNewNode(3, 0.0, 1.0, 0.0)
        node4 = modelpart.CreateNewNode(4, 0.0, 0.0, 1.0)
        modelpart.AddProperties(KM.Properties(1))
        elem = modelpart.CreateNewElement(
            "Element3D4N", 1, [1, 2, 3, 4], modelpart.GetProperties()[1]
        )

        # Set up constitutive law
        cl = CLA.SmallStrainIsotropicDamage3DLaw()
        props = modelpart.Properties[0]
        props.SetValue(KM.CONSTITUTIVE_LAW, cl)
        props.SetValue(KM.YOUNG_MODULUS, 6)
        props.SetValue(KM.POISSON_RATIO, 0.3)
        props.SetValue(CLA.HARDENING_CURVE, 1)
        props.SetValue(CLA.STRESS_LIMITS, [1.5, 2.3, 3.0])
        props.SetValue(CLA.HARDENING_PARAMETERS, [0.6, 0.4, 0.0])
        geom = KM.Tetrahedra3D4(node1, node2, node3, node4)
        cl_options = KM.Flags()
        cl_options.Set(KM.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
        cl_options.Set(KM.ConstitutiveLaw.COMPUTE_STRESS, True)
        cl_options.Set(KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
        stress_vector = KM.Vector(cl.GetStrainSize())
        strain_vector = KM.Vector(cl.GetStrainSize())
        constitutive_matrix = KM.Matrix(cl.GetStrainSize(), cl.GetStrainSize())
        cl_params = KM.ConstitutiveLawParameters()
        cl_params.SetOptions(cl_options)
        cl_params.SetStrainVector(strain_vector)
        cl_params.SetStressVector(stress_vector)
        cl_params.SetConstitutiveMatrix(constitutive_matrix)
        cl_params.SetProcessInfo(modelpart.ProcessInfo)
        cl_params.SetMaterialProperties(props)
        cl_params.SetElementGeometry(geom)
        F = KM.Matrix(3, 3, 0.0)
        F[0, 0] = F[1, 1] = F[2, 2] = 1
        cl_params.SetDeformationGradientF(F)
        cl_params.SetDeterminantF(1)

        # Initialize
        N = KM.Vector(4)
        elem.Initialize(modelpart.ProcessInfo)
        cl.InitializeMaterial(props, geom, N)
        imposed_strain = [-0.0759, 0.7483, 0.1879, 0.5391, 0.0063, -0.3292]
        imposed_strain_vector = KM.Vector(imposed_strain)
        parameters = KM.Parameters(
            """{
            "Parameters" : {
                "model_part_name": "Main",
                "imposed_strain_multiplier": "1.0 * t",
                "imposed_strain": """
            + "{}".format(imposed_strain)
            + """
            }
        }"""
        )
        process = set_initial_state_process.Factory(parameters, model)

        # Some iterations
        time = 0.0
        step = 0
        while time < 0.99:
            time += 0.1
            print(time)
            step += 1
            modelpart.ProcessInfo[KM.STEP] += 1
            modelpart.CloneTimeStep(time)

            zero_strain = KM.Vector([0, 0, 0, 0, 0, 0])
            cl_params.SetStrainVector(zero_strain)
            process.ExecuteInitializeSolutionStep()
            cl.InitializeMaterialResponseCauchy(cl_params)

            cl_params.SetStrainVector(zero_strain)
            cl.CalculateMaterialResponseCauchy(cl_params)

            cl_params.SetStrainVector(zero_strain)
            cl.FinalizeMaterialResponseCauchy(cl_params)

        strain = cl_params.GetStrainVector()
        stress = cl_params.GetStressVector()
        reference_stress_vector = KM.Vector(
            [1.28659, 3.14916, 1.88274, 0.60914, 0.00712, -0.37197]
        )
        self.assertVectorAlmostEqual(strain, imposed_strain_vector, prec=5)
        self.assertVectorAlmostEqual(stress, reference_stress_vector, prec=5)


if __name__ == "__main__":
    KratosUnittest.main()
