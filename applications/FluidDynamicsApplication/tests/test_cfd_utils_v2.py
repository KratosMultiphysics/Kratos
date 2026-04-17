import numpy as np
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication.cfd_utils as cfd_utils_module

cfd_utils_module.configure("OpenMP", "float64")

class TestCFDUtilsV2(UnitTest.TestCase):
    def ElemData(self, v, connectivities):
        return np.take(v, connectivities, axis=0)

    def setUp(self):
        #read mdpa
        self.dim = 3
        current_model = KM.Model()

        main_model_part = current_model.CreateModelPart("Main")
        self.model_part = main_model_part.CreateSubModelPart("submodelpart")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)

        model_part_io = KM.ModelPartIO("cfd_utils")
        model_part_io.ReadModelPart(self.model_part)

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, node.X)
            node.SetSolutionStepValue(KM.VELOCITY, [node.X**2, node.Y**2, node.Z**2])

        v_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.VELOCITY,data_shape=[self.dim],step_index=0)
        p_adaptor = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,KM.PRESSURE,0)

        v_adaptor.CollectData()
        self.v = v_adaptor.data.reshape((len(self.model_part.Nodes),self.dim))

        p_adaptor.CollectData()
        self.p = p_adaptor.data

        connectivity_adaptor = KM.TensorAdaptors.ConnectivityIdsTensorAdaptor(self.model_part.Elements)
        connectivity_adaptor.CollectData()
        self.connectivity = connectivity_adaptor.data
        self.connectivity -= 1 #have indices to start in

        # Preallocation of local array of shape functions
        self.N = np.ones(self.dim+1)/(self.dim+1)
        self.Ngauss = np.array([self.N])

        # Obtain the shape function derivatives
        geometry_adaptor_DN = KM.TensorAdaptors.GeometriesTensorAdaptor(
            self.model_part.Elements, #choose gometries if model_part contains geometries
            KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.ShapeFunctionDerivatives,
            KM.GeometryData.IntegrationMethod.GI_GAUSS_1)
        geometry_adaptor_DN.CollectData()
        self.DN = np.squeeze(geometry_adaptor_DN.data).copy() #this has shape nel*1*nnodes_in_el*dim - the copy is important as we need to own the data

        is_float32 = (np.dtype(cfd_utils_module.PRECISION).itemsize == 4)
        self.rtol = 1e-5 if is_float32 else 1e-12
        self.atol = 1e-5 if is_float32 else 1e-12

        self.velem = self.ElemData(self.v, self.connectivity)
        self.pelem = self.ElemData(self.p, self.connectivity)

        self.cfd_utils = cfd_utils_module.CFDUtils()

        # Stabilization parameters and prerequisites
        self.rho = 1.225
        self.elemental_volumes = np.array([elem.GetGeometry().DomainSize() for elem in self.model_part.Elements])

        Mscalar = np.zeros(len(self.model_part.Nodes), dtype=np.float64)
        Mel = np.einsum("e,i->ei", self.elemental_volumes, self.N)
        self.cfd_utils.AssembleVector(self.connectivity, Mel, Mscalar)
        M_expanded = np.tile(Mscalar[:, np.newaxis], (1, self.dim))
        self.Minv = 1.0 / M_expanded.ravel()

        det_J_volume_factor = 2.0 if self.dim == 2 else 6.0
        self.w_int_order = det_J_volume_factor * self.cfd_utils.GetGaussIntegrationWeights(self.dim, 1)

    def test_ComputeElementalDivergence(self):
        div_v = self.cfd_utils.ComputeElementalDivergence(self.DN, self.velem)
        ref_div_v = np.array([3.7, 4.100000000000001, 3.6, 3.45, 3.5, 3.6])
        np.testing.assert_allclose(div_v, ref_div_v, rtol=self.rtol, atol=self.atol)

    def test_ComputeElementwiseNodalDivergence(self):
        nodal_div_v = self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, self.velem)
        ref_nodal_div_v = np.array([[0.925 , 0.925 , 0.925 , 0.925 ],
                                    [1.025 , 1.025 , 1.025 , 1.025 ],
                                    [0.9   , 0.9   , 0.9   , 0.9   ],
                                    [0.8625, 0.8625, 0.8625, 0.8625],
                                    [0.875 , 0.875 , 0.875 , 0.875 ],
                                    [0.9   , 0.9   , 0.9   , 0.9   ]])
        np.testing.assert_allclose(nodal_div_v, ref_nodal_div_v, rtol=self.rtol, atol=self.atol)

    def test_ComputeNDN(self):
        n_dn = self.cfd_utils.ComputeNDN(self.N, self.DN, self.pelem)
        ref_n_dn = np.array([[[0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0]],
                             [[0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0]],
                             [[0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0]],
                             [[0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0]],
                             [[0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0]],
                             [[0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0],
                              [0.25, 0.0, 0.0]]])
        np.testing.assert_allclose(n_dn, ref_n_dn, rtol=self.rtol, atol=1e-14)

    def test_ComputeDNN(self):
        dn_n = self.cfd_utils.ComputeDNN(self.N, self.DN, self.pelem)
        ref_dn_n = np.array([[[ 0.591666666666667,  0.070902203856749, -0.537878787878788],
                              [ 0.               , -0.953512396694215,  0.806818181818182],
                              [ 0.               ,  0.806818181818182,  0.               ],
                              [-0.591666666666667,  0.075792011019284, -0.268939393939394]],
                             [[-0.608333333333333, -0.187179487179487, -0.052194280078895],
                              [ 0.608333333333333, -0.51474358974359 , -0.055793885601578],
                              [ 0.               ,  0.701923076923077, -0.593934911242604],
                              [ 0.               ,  0.               ,  0.701923076923077]],
                             [[-0.027173913043478, -0.195652173913044,  0.1875           ],
                              [-0.239130434782609,  0.228260869565217,  0.               ],
                              [ 0.25             ,  0.               ,  0.               ],
                              [ 0.016304347826087, -0.032608695652174, -0.1875           ]],
                             [[-0.009548611111111, -0.21875          , -0.009114583333333],
                              [ 0.030381944444444,  0.21875          , -0.209635416666667],
                              [-0.270833333333333,  0.               ,  0.21875          ],
                              [ 0.25             ,  0.               ,  0.               ]],
                             [[-0.109821428571429,  0.071550324675325, -0.366071428571429],
                              [-0.402678571428571, -0.04825487012987 ,  0.366071428571429],
                              [ 0.5125           , -0.489204545454545,  0.               ],
                              [ 0.               ,  0.465909090909091,  0.               ]],
                             [[-0.529356060606061,  0.447916666666667,  0.048550407925408],
                              [ 0.488636363636364,  0.               , -0.394667832167832],
                              [-0.               ,  0.               ,  0.413461538461538],
                              [ 0.040719696969697, -0.447916666666667, -0.067344114219114]]])
        np.testing.assert_allclose(dn_n, ref_dn_n, rtol=self.rtol, atol=self.atol)

    def test_ComputeElementalConvectiveOperator(self):
        grad_u = self.cfd_utils.ComputeElementalGradient(self.DN, self.velem)
        # Vector case
        conv_op = self.cfd_utils.ComputeElementalConvectiveOperator(self.velem, grad_u)
        ref_conv_op = np.array([[[3.375            , 0.               , 0.               ],
                                 [0.95             , 0.               , 1.331            ],
                                 [1.015818181818182, 1.331            , 2.145            ],
                                 [0.               , 0.               , 0.               ]],
                                [[0.               , 0.               , 0.               ],
                                 [3.375            , 0.               , 0.               ],
                                 [1.243000000000001, 2.197000000000001, 0.               ],
                                 [1.113961538461539, 1.287000000000001, 2.197000000000001]],
                                [[0.               , 0.               , 2.743999999999999],
                                 [0.               , 1.520875         , 1.74             ],
                                 [1.157625         , 1.33375          , 2.218021739130435],
                                 [0.               , 0.               , 0.               ]],
                                [[0.               , 0.               , 0.               ],
                                 [0.               , 1.728            , 0.               ],
                                 [0.               , 1.517999999999999, 1.728            ],
                                 [1.157625         , 1.320927083333333, 2.1645           ]],
                                [[0.               , 0.               , 0.               ],
                                 [0.               , 0.               , 2.743999999999999],
                                 [1.               , 0.               , 1.364            ],
                                 [1.16025          , 1.331            , 2.240325000000001]],
                                [[0.               , 1.728            , 0.               ],
                                 [1.331            , 2.171            , 0.               ],
                                 [1.1445           , 1.277977272727273, 2.197000000000001],
                                 [0.               , 0.               , 0.               ]]])
        np.testing.assert_allclose(conv_op, ref_conv_op, rtol=self.rtol, atol=self.atol)

        # Scalar case
        grad_p = self.cfd_utils.ComputeElementalGradient(self.DN, self.pelem)
        conv_op_scalar = self.cfd_utils.ComputeElementalConvectiveOperator(self.velem, grad_p)
        ref_conv_op_scalar = np.array([[2.25  , 1.    , 1.1025, 0.    ],
                                       [0.    , 2.25  , 1.21  , 1.1025],
                                       [0.    , 0.    , 1.1025, 0.    ],
                                       [0.    , 0.    , 0.    , 1.1025],
                                       [0.    , 0.    , 1.    , 1.1025],
                                       [0.    , 1.21  , 1.1025, 0.    ]])
        np.testing.assert_allclose(conv_op_scalar, ref_conv_op_scalar, rtol=self.rtol, atol=self.atol)

    def test_ComputeElementalGradient(self):
        grad_u = self.cfd_utils.ComputeElementalGradient(self.DN, self.velem)
        ref_grad_u = np.array([[[ 1.5              ,  0.107644628099173, -0.454545454545454],
                                [ 0.               ,  1.1              ,  0.               ],
                                [ 0.               ,  0.236363636363636,  1.1              ]],
                               [[ 1.5              , -0.338461538461538, -0.07707100591716 ],
                                [ 0.               ,  1.3              , -0.169230769230769],
                                [ 0.               ,  0.               ,  1.3              ]],
                               [[ 1.05             ,  0.               ,  0.               ],
                                [-0.052380952380952,  1.15             ,  0.               ],
                                [ 0.094824016563147, -0.208695652173913,  1.4              ]],
                               [[ 1.05             ,  0.               ,  0.               ],
                                [-0.045436507936508,  1.2              , -0.047916666666667],
                                [ 0.123809523809524,  0.               ,  1.2              ]],
                               [[ 1.               ,  0.047727272727273,  0.               ],
                                [ 0.               ,  1.1              ,  0.               ],
                                [-0.33             ,  0.196818181818182,  1.4              ]],
                               [[ 1.1              ,  0.               , -0.040384615384615],
                                [ 0.118181818181818,  1.2              , -0.18006993006993 ],
                                [ 0.               ,  0.               ,  1.3              ]]])
        np.testing.assert_allclose(grad_u, ref_grad_u, rtol=self.rtol, atol=self.atol)

    def test_InterpolateValue(self):
        # Vector case
        a_gauss = self.cfd_utils.InterpolateValue(self.Ngauss, self.velem)
        ref_a_gauss = np.array([[[1.088125, 0.3025  , 0.725   ]],
                               [[1.140625, 0.725   , 0.4225  ]],
                               [[0.275625, 0.633125, 1.2725  ]],
                               [[0.275625, 0.993125, 0.7825  ]],
                               [[0.525625, 0.3025  , 1.215   ]],
                               [[0.578125, 1.085   , 0.4225  ]]])
        np.testing.assert_allclose(a_gauss, ref_a_gauss, rtol=self.rtol, atol=self.atol)

        # Scalar case
        a_gauss_scalar = self.cfd_utils.InterpolateValue(self.Ngauss, self.pelem)
        ref_a_gauss_scalar = np.array([[0.8875],
                                       [0.9125],
                                       [0.2625],
                                       [0.2625],
                                       [0.5125],
                                       [0.5375]])
        np.testing.assert_allclose(a_gauss_scalar, ref_a_gauss_scalar, rtol=self.rtol, atol=self.atol)

    def test_ComputeConvectiveContribution(self):
        grad_u = self.cfd_utils.ComputeElementalGradient(self.DN, self.velem)
        conv_op_elem = self.cfd_utils.ComputeElementalConvectiveOperator(self.velem, grad_u)
        conv_contrib = self.cfd_utils.ComputeConvectiveContribution(conv_op_elem)
        ref_conv_contrib = np.array([
            [[0.4357909090909091, 0.06655000000000003, 0.17380000000000004],
            [0.31454090909090915, 0.06655000000000003, 0.24035000000000006],
            [0.31783181818181827, 0.13310000000000005, 0.28105],
            [0.26704090909090916, 0.06655000000000003, 0.17380000000000004]],
            [[0.28659807692307704, 0.17420000000000008, 0.10985000000000003],
            [0.4553480769230771, 0.17420000000000008, 0.10985000000000003],
            [0.3487480769230771, 0.28405000000000014, 0.10985000000000003],
            [0.342296153846154, 0.23855000000000012, 0.21970000000000006]],
            [[0.05788125000000001, 0.14273124999999998, 0.4723010869565217],
            [0.05788125000000001, 0.21877499999999997, 0.4221010869565217],
            [0.11576250000000002, 0.20941875, 0.4460021739130434],
            [0.05788125000000001, 0.14273124999999998, 0.3351010869565217]],
            [[0.05788125000000001, 0.22834635416666663, 0.19462500000000002],
            [0.05788125000000001, 0.31474635416666663, 0.19462500000000002],
            [0.05788125000000001, 0.30424635416666657, 0.281025],
            [0.11576250000000002, 0.2943927083333333, 0.30285000000000006]],
            [[0.10801250000000001, 0.06655000000000003, 0.31741625],
            [0.10801250000000001, 0.06655000000000003, 0.45461625],
            [0.15801250000000003, 0.06655000000000003, 0.38561625000000005],
            [0.16602500000000003, 0.13310000000000005, 0.42943250000000005]],
            [[0.12377500000000002, 0.3452488636363637, 0.10985000000000003],
            [0.19032500000000005, 0.3673988636363637, 0.10985000000000003],
            [0.18100000000000005, 0.3227477272727273, 0.21970000000000006],
            [0.12377500000000002, 0.25884886363636367, 0.10985000000000003]]])
        np.testing.assert_allclose(conv_contrib, ref_conv_contrib, rtol=self.rtol, atol=self.atol)

    def test_ApplyLaplacian(self):
        laplacian_v = self.cfd_utils.ApplyLaplacian(self.DN, self.velem)
        ref_laplacian_v = np.array([[[ 1.284081802244837,  0.087878787878788, -0.64778362133734 ],
        [-0.528874393825558, -1.181818181818182,  0.746055597295267],
        [ 0.09785875281743 ,  1.               ,  0.214876033057851],
        [-0.85306616123671 ,  0.093939393939394, -0.313148009015778]],
       [[-0.92616359137752 , -0.256986800182066, -0.074358974358974],
        [ 1.195639450065941, -0.722985889849795, -0.07948717948718 ],
        [-0.210190469521375,  1.110150204824761, -0.846153846153846],
        [-0.059285389167046, -0.130177514792899,  1.               ]],
       [[-0.108695652173913, -0.851720398304249,  1.145733403632404],
        [-0.956521739130435,  1.047717637779749, -0.267856607041052],
        [ 1.               , -0.049886621315193,  0.090308587202997],
        [ 0.065217391304348, -0.146110618160308, -0.968185383794349]],
       [[-0.038194444444444, -0.996683443825061, -0.046170319979844],
        [ 0.121527777777778,  1.033007927322164, -0.944003527336861],
        [-1.083333333333333,  0.006948381204333,  0.872260015117158],
        [ 1.               , -0.043272864701436,  0.117913832199547]],
       [[-0.207622491145218,  0.153571428571428, -0.901807851239669],
        [-0.790208087367178, -0.103571428571428,  1.240754132231405],
        [ 0.954442148760331, -1.05             , -0.517871900826446],
        [ 0.043388429752066,  1.               ,  0.178925619834711]],
       [[-1.086981127846513,  0.86734375764096 ,  0.117424242424242],
        [ 1.029653039268424,  0.239657195950902, -0.954545454545454],
        [-0.031065088757397, -0.138515330823023,  1.               ],
        [ 0.088393177335485, -0.96848562276884 , -0.162878787878788]]])
        np.testing.assert_allclose(laplacian_v, ref_laplacian_v, rtol=self.rtol, atol=self.atol)

        laplacian_p = self.cfd_utils.ApplyLaplacian(self.DN, self.pelem)
        ref_laplacian_p = np.array([[ 0.666666666666667,  0.               , -0.               ,
        -0.666666666666667],
       [-0.666666666666667,  0.666666666666667,  0.               ,
         0.               ],
       [-0.10351966873706 , -0.910973084886128,  0.952380952380952,
         0.062111801242236],
       [-0.036375661375661,  0.115740740740741, -1.031746031746032,
         0.952380952380952],
       [-0.214285714285714, -0.785714285714286,  1.               ,
         0.               ],
       [-0.984848484848485,  0.909090909090909,  0.               ,
         0.075757575757576]])
        np.testing.assert_allclose(laplacian_p, ref_laplacian_p, rtol=self.rtol, atol=self.atol)

    def test_ComputeMomentumStabilization(self):
        grad_v = self.cfd_utils.ComputeElementalGradient(self.DN, self.velem)
        conv_op_elem = self.cfd_utils.ComputeElementalConvectiveOperator(self.velem, grad_v)
        convective = self.rho * self.cfd_utils.ComputeConvectiveContribution(conv_op_elem)
        convective *= self.elemental_volumes[:, None, None]

        pi_conv = np.zeros((len(self.model_part.Nodes), self.dim), dtype=np.float64)
        self.cfd_utils.AssembleVector(self.connectivity, convective, pi_conv)
        pi_conv = pi_conv.ravel() * self.Minv
        pi_conv = pi_conv.reshape((len(self.model_part.Nodes), self.dim))

        conv_proj_el = self.ElemData(pi_conv, self.connectivity)
        stab_mom = self.cfd_utils.ComputeMomentumStabilization(self.DN, self.velem, conv_op_elem, conv_proj_el)
        # import sys
        # np.set_printoptions(threshold=sys.maxsize, precision=15, suppress=True)
        # print(repr(stab_mom))
        ref_stab_mom = np.array([
            [[ 0.241165647869191,  0.013774781008715,  0.053069369820089],
            [ 0.152320472769996,  0.052033571185278,  0.116770422009721],
            [-0.361769818005304, -0.068541295676325, -0.195464890524385],
            [-0.031716302633884,  0.002732943482331,  0.025625098694574]],
            [[-0.06936240955999 , -0.002253043857513, -0.006435899786187],
            [ 0.233364763152239,  0.04029544781584 ,  0.038219495703796],
            [ 0.149381562346675,  0.129839887570243,  0.019406715274716],
            [-0.313383915938924, -0.167882291528569, -0.051190311192325]],
            [[ 0.06691774558499 ,  0.122004000695001,  0.318149650257302],
            [ 0.040522598627701,  0.115284967677564,  0.188780905835649],
            [-0.069997635438296, -0.18355769822628 , -0.392630777866901],
            [-0.037442708774394, -0.053731270146285, -0.11429977822605 ]],
            [[ 0.001337742928752,  0.018645861286902,  0.005392977118481],
            [ 0.010599387980355,  0.14020064542423 ,  0.066903993108358],
            [ 0.060923988285783,  0.163016533181729,  0.159348388524561],
            [-0.07286111919489 , -0.321863039892861, -0.231645358751401]],
            [[-0.021454293024849, -0.029755708870864, -0.064658344330006],
            [ 0.157916862668863,  0.101270206536038,  0.410427995290483],
            [-0.002295271980404, -0.009136344763324,  0.003895796518437],
            [-0.13416729766361 , -0.062378152901849, -0.349665447478914]],
            [[ 0.181247002800016,  0.352071254148047,  0.116037736825715],
            [-0.005736846984176,  0.007733049874852, -0.038971725483889],
            [-0.133030346337386, -0.301131355106468, -0.055339139920372],
            [-0.042479809478454, -0.058672948916431, -0.021726871421454]]])
        np.testing.assert_allclose(stab_mom, ref_stab_mom, rtol=self.rtol, atol=self.atol)

    def test_ComputeDivDivStabilization(self):
        aux_scalar = self.cfd_utils.ComputeElementwiseNodalDivergence(self.N, self.DN, self.velem)
        aux_scalar *= self.elemental_volumes[:, np.newaxis]
        pi_div = np.zeros((len(self.model_part.Nodes)), dtype=np.float64)
        self.cfd_utils.AssembleVector(self.connectivity, aux_scalar, pi_div)
        pi_div *= self.Minv[::self.dim]

        div_proj_el = self.ElemData(pi_div, self.connectivity)
        stab_div = self.cfd_utils.ComputeDivDivStabilization(self.N, self.DN, self.velem, div_proj_el)

        ref_stab_div = np.array([[[-0.023159903208659, -0.002775360301864,  0.021054457462417],
        [ 0.               ,  0.037323810956103, -0.031581686193626],
        [ 0.               , -0.031581686193626,  0.               ],
        [ 0.023159903208659, -0.002966764460613,  0.010527228731209]],
       [[-0.195178694074141, -0.060054982792043, -0.016746100970858],
        [ 0.195178694074141, -0.165151202678119, -0.01790100448609 ],
        [ 0.               ,  0.225206185470163, -0.190559080013215],
        [ 0.               ,  0.               ,  0.225206185470163]],
       [[ 0.002048885567121,  0.014751976083268, -0.014137310413131],
        [ 0.018030192990661, -0.017210638763812,  0.               ],
        [-0.018849747217509,  0.               ,  0.               ],
        [-0.001229331340272,  0.002458662680545,  0.014137310413131]],
       [[ 0.00597087706625 ,  0.136787365517738,  0.005699473563239],
        [-0.018998245210797, -0.136787365517738,  0.131087891954499],
        [ 0.169355785879104,  0.               , -0.136787365517738],
        [-0.156328417734558,  0.               ,  0.               ]],
       [[ 0.029902785629794, -0.01948211791032 ,  0.099675952099314],
        [ 0.109643547309245,  0.013139102776728, -0.099675952099314],
        [-0.139546332939039,  0.133203317805446,  0.               ],
        [ 0.               , -0.126860302671854,  0.               ]],
       [[ 0.10477905558717 , -0.088659200881451, -0.009609913382255],
        [-0.096719128234311,  0.               ,  0.078119295881558],
        [ 0.               ,  0.               , -0.081839262352109],
        [-0.008059927352859,  0.088659200881451,  0.013329879852806]]])
        np.testing.assert_allclose(stab_div, ref_stab_div, rtol=self.rtol, atol=self.atol)

    def test_ComputePressureStabilizationProjectionTerm(self):
        pi_press_el = self.cfd_utils.ComputeNDN(self.N, self.DN, self.pelem)
        pi_press_el *= self.elemental_volumes[:, np.newaxis, np.newaxis]
        pi_press = np.zeros((len(self.model_part.Nodes), self.dim), dtype=np.float64)
        self.cfd_utils.AssembleVector(self.connectivity, pi_press_el, pi_press)
        pi_press = pi_press.ravel() * self.Minv
        pi_press = pi_press.reshape((len(self.model_part.Nodes), self.dim))

        pres_proj_el = self.ElemData(pi_press, self.connectivity)

        stab_press = self.cfd_utils.ComputePressureStabilizationProjectionTerm(self.N, self.DN, pres_proj_el)
        ref_stab_press = np.array([[ 0.666666666666667,  0.               ,  0.               ,
        -0.666666666666667],
       [-0.666666666666667,  0.666666666666667, -0.               ,
         0.               ],
       [-0.10351966873706 , -0.910973084886128,  0.952380952380952,
         0.062111801242236],
       [-0.036375661375661,  0.115740740740741, -1.031746031746032,
         0.952380952380952],
       [-0.214285714285714, -0.785714285714286,  1.               ,
         0.               ],
       [-0.984848484848485,  0.909090909090909,  0.               ,
         0.075757575757576]])
        np.testing.assert_allclose(stab_press, ref_stab_press, rtol=self.rtol, atol=self.atol)

if __name__ == '__main__':
    UnitTest.main()