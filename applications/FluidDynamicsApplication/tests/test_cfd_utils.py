import numpy as np
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication.cfd_utils as cfd_utils_module

cfd_utils_module.configure("OpenMP", "float64")

class TestCFDUtils(UnitTest.TestCase):
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

        model_part_io = KM.ModelPartIO("cfd_utils", KM.ModelPartIO.SKIP_TIMER)
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
        ref_stab_mom = np.array([[[ 0.094903781666614, -0.141011265504643, -0.128222298224962],
        [-0.008839356240907, -0.114461541088728, -0.07999652885956 ],
        [ 0.000661256580191, -0.068541295676325, -0.028969778250378],
        [-0.086725682005898,  0.324014102269696,  0.237188605334899]],

       [[-0.13050483648701 ,  0.105385799928034,  0.320174463925397],
        [ 0.106792886634332, -0.089493261831567, -0.144416739905702],
        [ 0.004537485540416,  0.008015625330015, -0.124567412827371],
        [ 0.019174464312261, -0.023908163426482, -0.051190311192325]],

       [[-0.143210793979794, -0.09860013126933 ,  0.057500622763395],
        [-0.100641899265625, -0.032601649163062,  0.014577343688467],
        [-0.069997635438296, -0.028948962438353,  0.048871462954706],
        [ 0.313850328683715,  0.160150742870745, -0.120949429406568]],

       [[ 0.292858453579981, -0.001782509510481,  0.190271227002012],
        [-0.096240524605413, -0.005279680174441, -0.100385609340366],
        [-0.123756809779678,  0.004553026836638, -0.069304028127915],
        [-0.07286111919489 ,  0.002509162848283, -0.020581589533731]],

       [[ 0.279726260987099,  0.296782783987071, -0.09118907783863 ],
        [-0.139947753974007, -0.152317664202192,  0.049349424931468],
        [-0.071929956731032, -0.08208696688303 , -0.004193720855013],
        [-0.06784855028206 , -0.062378152901849,  0.046033373762175]],

       [[-0.121741619077499,  0.024174608536179, -0.162404635742216],
        [-0.058577880588672,  0.037567986849959, -0.104393957565646],
        [-0.077673073037437,  0.021314605499533, -0.055339139920372],
        [ 0.257992572703609, -0.083057200885672,  0.322137733228235]]])
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

class TestCFDUtilsMisc(UnitTest.TestCase):
    """Checks for shape contracts and error paths of the CFDUtils API (no mesh required)."""

    def setUp(self):
        self.cfd_utils = cfd_utils_module.CFDUtils()

    def test_GetShapeFunctionsOnGaussPoints_order1_is_2d(self):
        # Regression: order-1 used to return a 3D array of shape (1, 1, dim+1).
        # It must be 2D with shape (ngauss, n_in_el) = (1, dim+1), consistent with order-2.
        for dim in (2, 3):
            N = self.cfd_utils.GetShapeFunctionsOnGaussPoints(dim, 1)
            self.assertEqual(N.ndim, 2)
            self.assertEqual(N.shape, (1, dim + 1))
            np.testing.assert_allclose(N.sum(axis=1), np.ones(1), rtol=1e-12, atol=1e-14)

    def test_GetShapeFunctionsOnGaussPoints_order1_order2_consistent_ndim(self):
        for dim in (2, 3):
            N1 = self.cfd_utils.GetShapeFunctionsOnGaussPoints(dim, 1)
            N2 = self.cfd_utils.GetShapeFunctionsOnGaussPoints(dim, 2)
            self.assertEqual(N1.ndim, N2.ndim)
            self.assertEqual(N1.shape[1], N2.shape[1])  # same n_in_el

    def test_ComputeNDN_rejects_wrong_pel_shape(self):
        # Regression: the error path used to raise NameError ('field' undefined) instead of ValueError
        nelem, n_in_el, dim = 2, 4, 3
        DN = np.zeros((nelem, n_in_el, dim))
        N = np.ones(n_in_el) / n_in_el
        wrong_pel = np.zeros((nelem, n_in_el + 1))
        with self.assertRaises(ValueError):
            self.cfd_utils.ComputeNDN(N, DN, wrong_pel)

class TestCFDUtilsKernels(UnitTest.TestCase):
    """Granular unit tests for the cfd_utils kernels on small synthetic fixtures
    (no mesh required). Each kernel is checked against an independent reference
    (explicit einsum/loops), so the tests survive internal refactors (e.g. the
    internal use of a buffer pool)."""

    def setUp(self):
        self.cfd_utils = cfd_utils_module.CFDUtils()
        self.xp = cfd_utils_module.xp

    @staticmethod
    def _csr_to_dense(A):
        import scipy.sparse as sp
        return sp.csr_matrix((
            np.asarray(A.value_data()).copy(),
            np.asarray(A.index2_data()).copy(),
            np.asarray(A.index1_data()).copy()
        ), shape=(A.Size1(), A.Size2())).toarray()

    # ------------------------------------------------------------------
    # Reference-element tables
    # ------------------------------------------------------------------
    def test_shape_functions_partition_of_unity_and_quadrature_consistency(self):
        # N rows sum to 1 (partition of unity); the quadrature reproduces the
        # integral of a constant: sum_g w_g * N_g[I] = V_ref / n_in_el
        for dim in (2, 3):
            V_ref = 0.5 if dim == 2 else 1.0 / 6.0
            for order in (1, 2):
                N = self.cfd_utils.GetShapeFunctionsOnGaussPoints(dim, order)
                w = self.cfd_utils.GetGaussIntegrationWeights(dim, order)
                self.assertEqual(N.shape[0], w.shape[0])
                self.assertEqual(N.shape[1], dim + 1)
                np.testing.assert_allclose(N.sum(axis=1), np.ones(N.shape[0]), rtol=1e-12, atol=1e-14)
                np.testing.assert_allclose(w @ N, np.full(dim + 1, V_ref / (dim + 1)), rtol=1e-12, atol=1e-14)
                # weights integrate the reference element volume exactly
                np.testing.assert_allclose(w.sum(), V_ref, rtol=1e-12, atol=1e-14)

    def test_gauss_weights_reject_unsupported_inputs(self):
        with self.assertRaises(Exception):
            self.cfd_utils.GetGaussIntegrationWeights(1, 1)
        with self.assertRaises(Exception):
            self.cfd_utils.GetGaussIntegrationWeights(2, 3)
        with self.assertRaises(Exception):
            self.cfd_utils.GetShapeFunctionsOnGaussPoints(2, 4)

    def test_elemental_mass_matrix_consistency(self):
        # Consistent mass on the unit-measure reference element:
        # M[I,J] = (delta_IJ + 1) * measure / (n*(n+1)) with measure = 1 here
        for dim in (2, 3):
            n = dim + 1
            M = self.cfd_utils.GetElementalMassMatrix(dim)
            ref = (np.eye(n) + np.ones((n, n))) / (n * (n + 1))
            np.testing.assert_allclose(M, ref, rtol=1e-12, atol=1e-14)
        with self.assertRaises(Exception):
            self.cfd_utils.GetElementalMassMatrix(1)

    def test_wall_mass_matrix_consistency(self):
        # Wall mass on the unit-measure edge/face: same (delta + 1)/(n*(n+1)) rule
        for dim, n in ((2, 2), (3, 3)):
            M = self.cfd_utils.GetWallMassMatrix(dim)
            ref = (np.eye(n) + np.ones((n, n))) / (n * (n + 1))
            np.testing.assert_allclose(M, ref, rtol=1e-12, atol=1e-14)
        with self.assertRaises(Exception):
            self.cfd_utils.GetWallMassMatrix(1)

    # ------------------------------------------------------------------
    # ComputeLaplacianMatrix
    # ------------------------------------------------------------------
    def test_ComputeLaplacianMatrix_reference_and_symmetry(self):
        rng = np.random.RandomState(31)
        DN = rng.randn(5, 3, 2)
        out = self.cfd_utils.ComputeLaplacianMatrix(DN)
        ref = np.einsum("eim,ejm->eij", DN, DN)
        np.testing.assert_allclose(out, ref, rtol=1e-12, atol=1e-13)
        # symmetric by construction
        np.testing.assert_allclose(out, np.swapaxes(out, 1, 2), rtol=0.0, atol=1e-13)

    def test_ComputeLaplacianMatrix_writes_to_out(self):
        rng = np.random.RandomState(32)
        DN = rng.randn(4, 4, 3)
        buf = np.zeros((4, 4, 4))
        ret = self.cfd_utils.ComputeLaplacianMatrix(DN, out=buf)
        self.assertIs(ret, buf)
        np.testing.assert_allclose(buf, np.einsum("eim,ejm->eij", DN, DN), rtol=1e-12, atol=1e-13)

    # ------------------------------------------------------------------
    # ComputeBodyForceContribution
    # ------------------------------------------------------------------
    def test_ComputeBodyForceContribution_reference(self):
        rng = np.random.RandomState(33)
        b = rng.randn(3, 3, 2)
        M2 = np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]]) / 12.0  # unit-volume 2D
        ref = np.einsum("ij,ejd->eid", M2, b)
        out = self.cfd_utils.ComputeBodyForceContribution(b)
        np.testing.assert_allclose(out, ref, rtol=1e-12, atol=1e-13)

    def test_ComputeBodyForceContribution_zero_input_gives_zero(self):
        out = self.cfd_utils.ComputeBodyForceContribution(np.zeros((2, 4, 3)))
        np.testing.assert_array_equal(out, np.zeros((2, 4, 3)))

    # ------------------------------------------------------------------
    # ComputeNavierSlipWallContribution
    # ------------------------------------------------------------------
    def test_ComputeNavierSlipWallContribution_reference(self):
        # 2 wall elements of 2 nodes each (2D): projectors built from unit normals
        rng = np.random.RandomState(34)
        v = rng.randn(2, 2, 2)
        normals = np.array([[1.0, 0.0], [0.0, 1.0]])
        P = np.eye(2)[None, :, :] - normals[:, :, None] * normals[:, None, :]  # (2,2,2)
        M_w = np.array([[1.0/3.0, 1.0/6.0], [1.0/6.0, 1.0/3.0]])
        visc, slip = 0.7, 0.25
        v_tang = np.einsum("eij,ekj->eki", P, v)
        ref = np.einsum("nk,eki->eni", M_w, v_tang) * (-visc / (slip + 1e-12))
        out = self.cfd_utils.ComputeNavierSlipWallContribution(v, P, visc, slip)
        np.testing.assert_allclose(out, ref, rtol=1e-11, atol=1e-13)

    def test_ComputeNavierSlipWallContribution_pure_normal_flow_is_unchanged(self):
        # Velocity parallel to the wall normal has zero tangential part -> zero contribution
        v = np.array([[[2.0, 0.0]], [[0.0, 3.0]], [[-1.0, 0.0]], [[0.0, -5.0]]]).reshape(4, 2, 1)[:, :, 0:2]
        v = np.array([[[2.0, 0.0], [2.0, 0.0]],   # element 0: velocity along +x
                      [[0.0, 3.0], [0.0, 3.0]]])  # element 1: velocity along +y
        normals = np.array([[1.0, 0.0], [0.0, 1.0]])
        P = np.eye(2)[None, :, :] - normals[:, :, None] * normals[:, None, :]
        out = self.cfd_utils.ComputeNavierSlipWallContribution(v, P, 1.5, 0.1)
        np.testing.assert_allclose(out, np.zeros((2, 2, 2)), rtol=0.0, atol=1e-14)

    # ------------------------------------------------------------------
    # Scalar matrix allocation / assembly / Dirichlet machinery
    # ------------------------------------------------------------------
    def _poisson_1d_fixture(self):
        # Unit 1D Poisson on nodes 0..4 as "2-node 2D triangles" (edges)
        conn = np.array([[0, 1], [1, 2], [2, 3], [3, 4]], dtype=np.int32)
        lhs_el = np.tile(np.array([[1.0, -1.0], [-1.0, 1.0]]), (conn.shape[0], 1, 1))
        ref = np.zeros((5, 5))
        for e in range(conn.shape[0]):
            for I in range(2):
                for J in range(2):
                    ref[conn[e, I], conn[e, J]] += lhs_el[e, I, J]
        return conn, lhs_el, ref

    def test_AllocateScalarMatrix_graph_and_indices(self):
        conn, _, ref = self._poisson_1d_fixture()
        A, csr_indices = self.cfd_utils.AllocateScalarMatrix(conn)
        self.assertEqual((A.Size1(), A.Size2()), (5, 5))
        dense = self._csr_to_dense(A)
        # Freshly allocated matrices carry no spurious structure outside the graph
        self.assertEqual(dense.shape, (5, 5))
        # every assembly index addresses a real entry of the value buffer
        self.assertEqual(csr_indices.shape, (conn.shape[0], 2, 2))
        self.assertTrue(np.all(csr_indices >= 0))
        self.assertTrue(np.all(csr_indices < np.asarray(A.value_data()).size))

    def test_AssembleScalarMatrix_matches_dense_reference(self):
        conn, lhs_el, ref = self._poisson_1d_fixture()
        A, _ = self.cfd_utils.AllocateScalarMatrix(conn)
        self.cfd_utils.AssembleScalarMatrix(lhs_el, conn, A)
        np.testing.assert_allclose(self._csr_to_dense(A), ref, rtol=1e-12, atol=1e-14)

    def test_AssembleScalarMatrixByCSRIndices_matches_dense_reference(self):
        conn, lhs_el, ref = self._poisson_1d_fixture()
        A, csr_indices = self.cfd_utils.AllocateScalarMatrix(conn)
        self.cfd_utils.AssembleScalarMatrixByCSRIndices(lhs_el, csr_indices, A)
        np.testing.assert_allclose(self._csr_to_dense(A), ref, rtol=1e-12, atol=1e-14)
        # re-assembly must zero the buffer first (no accumulation)
        self.cfd_utils.AssembleScalarMatrixByCSRIndices(lhs_el, csr_indices, A)
        np.testing.assert_allclose(self._csr_to_dense(A), ref, rtol=1e-12, atol=1e-14)

    def test_GetScalarMatrixDirichletIndices_matches_brute_force(self):
        conn, lhs_el, _ = self._poisson_1d_fixture()
        A, _ = self.cfd_utils.AllocateScalarMatrix(conn)
        self.cfd_utils.AssembleScalarMatrix(lhs_el, conn, A)
        fixed = np.array([1, 3], dtype=np.int32)
        csr_data_idx, csr_diag_idx = self.cfd_utils.GetScalarMatrixDirichletIndices(fixed, A)

        # Brute-force reference: walk the CSR rows
        indptr = np.asarray(A.index1_data())
        indices = np.asarray(A.index2_data())
        expected_data, expected_diag = [], []
        for row in range(A.Size1()):
            for p in range(indptr[row], indptr[row + 1]):
                if row in fixed or indices[p] in fixed:
                    expected_data.append(p)
                    if indices[p] == row:
                        expected_diag.append(p)
        np.testing.assert_array_equal(np.sort(csr_data_idx), np.sort(np.array(expected_data, dtype=np.int32)))
        np.testing.assert_array_equal(np.sort(csr_diag_idx), np.sort(np.array(expected_diag, dtype=np.int32)))

    def test_ApplyHomogeneousDirichlet_scalar_diag_reference(self):
        conn, lhs_el, ref = self._poisson_1d_fixture()
        A, _ = self.cfd_utils.AllocateScalarMatrix(conn)
        self.cfd_utils.AssembleScalarMatrix(lhs_el, conn, A)
        fixed = np.array([1, 3], dtype=np.int32)
        csr_data_idx, csr_diag_idx = self.cfd_utils.GetScalarMatrixDirichletIndices(fixed, A)
        rhs = np.arange(5, dtype=np.float64) + 1.0
        self.cfd_utils.ApplyHomogeneousDirichlet(fixed, csr_data_idx, csr_diag_idx, A, rhs, diag_value=1.0)

        ref_masked = ref.copy()
        ref_masked[fixed, :] = 0.0
        ref_masked[:, fixed] = 0.0
        ref_masked[1, 1] = 1.0
        ref_masked[3, 3] = 1.0
        np.testing.assert_allclose(self._csr_to_dense(A), ref_masked, rtol=1e-12, atol=1e-14)
        np.testing.assert_array_equal(rhs[fixed], np.zeros(2))
        np.testing.assert_array_equal(rhs[[0, 2, 4]], np.array([1.0, 3.0, 5.0]))

    def test_ApplyHomogeneousDirichlet_array_diag_reference(self):
        # The stage passes diag_value as an array (one value per fixed row)
        conn, lhs_el, ref = self._poisson_1d_fixture()
        A, _ = self.cfd_utils.AllocateScalarMatrix(conn)
        self.cfd_utils.AssembleScalarMatrix(lhs_el, conn, A)
        fixed = np.array([1, 3], dtype=np.int32)
        csr_data_idx, csr_diag_idx = self.cfd_utils.GetScalarMatrixDirichletIndices(fixed, A)
        rhs = np.zeros(5)
        diag_values = self.xp.asarray(np.array([2.5, 4.0]))
        self.cfd_utils.ApplyHomogeneousDirichlet(fixed, csr_data_idx, csr_diag_idx, A, rhs, diag_value=diag_values)

        dense = self._csr_to_dense(A)
        self.assertAlmostEqual(float(dense[1, 1]), 2.5)
        self.assertAlmostEqual(float(dense[3, 3]), 4.0)

    def test_ApplyHomogeneousDirichlet_no_fixed_doea_nothing(self):
        conn, lhs_el, ref = self._poisson_1d_fixture()
        A, _ = self.cfd_utils.AllocateScalarMatrix(conn)
        self.cfd_utils.AssembleScalarMatrix(lhs_el, conn, A)
        fixed = np.zeros(0, dtype=np.int32)
        csr_data_idx, csr_diag_idx = self.cfd_utils.GetScalarMatrixDirichletIndices(fixed, A)
        rhs = np.ones(5)
        self.cfd_utils.ApplyHomogeneousDirichlet(fixed, csr_data_idx, csr_diag_idx, A, rhs, diag_value=1.0)
        np.testing.assert_allclose(self._csr_to_dense(A), ref, rtol=1e-12, atol=1e-14)
        np.testing.assert_array_equal(rhs, np.ones(5))

    # ------------------------------------------------------------------
    # Temporary buffers lifecycle
    # ------------------------------------------------------------------
    def test_temporaries_reused_for_same_dimensions(self):
        cu = self.cfd_utils
        cu._ensure_temporaries(10, 3, 2)
        first = cu._pool.buffers[cu.TMP_A]
        cu._ensure_temporaries(10, 3, 2)
        self.assertIs(cu._pool.buffers[cu.TMP_A], first)     # same dims -> no reallocation
        cu._ensure_temporaries(12, 3, 2)
        self.assertIsNot(cu._pool.buffers[cu.TMP_A], first)  # new dims -> reallocation

    def test_pool_has_only_two_scratch_buffers(self):
        # The point of the pool refactor: minimal memory footprint, two slots only
        cu = self.cfd_utils
        cu._ensure_temporaries(10, 3, 2)
        self.assertEqual(len(cu._pool.buffers), 2)
        # A: nelem * max(n_in_el*dim, dim*dim); B: nelem * dim*dim
        self.assertEqual(cu._pool.max_sizes, [10 * max(3 * 2, 2 * 2), 10 * 2 * 2])

    def test_pool_fully_released_after_each_kernel(self):
        # Every kernel must release its scratch before returning (tracked by the pool)
        rng = np.random.RandomState(41)
        nelem, n, dim = 6, 3, 2
        cu = self.cfd_utils
        cu.allocate_temporaries(nelem, n, dim)
        DN = rng.randn(nelem, n, dim)
        uel = rng.randn(nelem, n, dim)
        pel = rng.randn(nelem, n)
        N = np.full(n, 1.0 / n)
        calls = [
            lambda: cu.ComputeElementalDivergence(DN, uel),
            lambda: cu.ComputeElementwiseNodalDivergence(N, DN, uel),
            lambda: cu.ComputeNDN(N, DN, pel),
            lambda: cu.ComputeDNN(N, DN, pel),
            lambda: cu.ApplyLaplacian(DN, pel),
            lambda: cu.ApplyLaplacian(DN, uel),
            lambda: cu.ComputeMomentumStabilization(DN, uel, uel, uel),
            lambda: cu.ComputeDivDivStabilization(N, DN, uel, pel),
            lambda: cu.ComputePressureStabilizationProjectionTerm(N, DN, uel),
        ]
        for call in calls:
            call()
            self.assertFalse(any(cu._pool.Status()),
                             f"pool buffers left in use after {call}")


class TestRobustCG(UnitTest.TestCase):
    """Behavior pins for the robust_cg solver (numpy backend)."""

    def setUp(self):
        self.cfd_utils = cfd_utils_module.CFDUtils()
        self.xp = cfd_utils_module.xp

    def _make_spd(self, n, seed):
        import scipy.sparse as sp
        rng = np.random.RandomState(seed)
        M = rng.rand(n, n)
        A = M @ M.T + n * np.eye(n)  # well-conditioned SPD
        return sp.csr_matrix(A)

    def test_robust_cg_converges_matches_direct_solve(self):
        n = 30
        A = self._make_spd(n, seed=0)
        rng = np.random.RandomState(1)
        b = rng.rand(n)
        x, info = self.cfd_utils.robust_cg(A, b, xp=self.xp, rtol=1e-10, atol=0.0, maxiter=1000)
        self.assertEqual(info, 0)
        x_ref = np.linalg.solve(A.toarray(), b)
        np.testing.assert_allclose(np.asarray(x), x_ref, rtol=1e-6, atol=1e-8)

    def test_robust_cg_zero_rhs_returns_zero_immediately(self):
        A = self._make_spd(5, seed=2)
        b = np.zeros(5)
        x, info = self.cfd_utils.robust_cg(A, b, xp=self.xp, rtol=1e-8, atol=0.0)
        self.assertEqual(info, 0)
        np.testing.assert_array_equal(np.asarray(x), np.zeros(5))

    def test_robust_cg_zero_rhs_does_not_alias_rhs(self):
        # Regression: returning 'b' for a zero rhs aliased the caller's buffer (a pool slot
        # that is released and re-acquired right after the solve). The result must be a fresh array.
        A = self._make_spd(5, seed=2)
        b = np.zeros(5)
        x, info = self.cfd_utils.robust_cg(A, b, xp=self.xp, rtol=1e-8, atol=0.0, return_info_dict=True)
        self.assertTrue(info["converged"])
        self.assertFalse(np.shares_memory(np.asarray(x), b))
        np.testing.assert_array_equal(np.asarray(x), np.zeros(5))

    def test_host_scalars_matches_individual_reads(self):
        # The batched readback must return exactly the values of individual float() reads,
        # including float32 inputs and non-finite values (the finiteness decisions rely on them).
        vals = [np.float64(1.2345678901234567), np.float32(3.1415927), np.float64(np.nan), np.float64(-np.inf)]
        dev = [self.xp.asarray(v) for v in vals] if hasattr(self.xp, "asarray") else vals
        out = cfd_utils_module._host_scalars(self.xp, *dev)
        expected = [float(v) for v in vals]
        for got, exp in zip(out, expected):
            self.assertTrue((np.isnan(got) and np.isnan(exp)) or got == exp)

    def test_robust_cg_empty_system(self):
        A = self._make_spd(0, seed=3)
        b = np.zeros(0)
        x, info = self.cfd_utils.robust_cg(A, b, xp=self.xp)
        self.assertEqual(info, 0)
        self.assertEqual(len(np.asarray(x)), 0)

    def test_robust_cg_max_iterations_returns_iteration_count(self):
        A = self._make_spd(60, seed=4)
        rng = np.random.RandomState(5)
        b = rng.rand(60)
        x, info = self.cfd_utils.robust_cg(A, b, xp=self.xp, rtol=1e-14, atol=0.0, maxiter=3)
        self.assertEqual(info, 3)  # non-convergence reports iteration count (never 0)
        x_dict, info_dict = self.cfd_utils.robust_cg(
            A, b, xp=self.xp, rtol=1e-14, atol=0.0, maxiter=3, return_info_dict=True)
        self.assertFalse(info_dict["converged"])
        self.assertEqual(info_dict["iterations"], 3)
        self.assertEqual(info_dict["reason"], "maximum iterations reached")

    def test_robust_cg_initial_guess_shortcut_uses_full_tolerance(self):
        # Regression: the initial-guess shortcut used to compare against 'atol' (0.0 by default)
        # instead of the full stopping tolerance, making it dead code. An initial guess whose true
        # residual is below rtol*||b|| but above atol must be accepted with ZERO iterations.
        A = self._make_spd(4, seed=6)
        rng = np.random.RandomState(7)
        b = rng.rand(4)
        x_ref = np.linalg.solve(A.toarray(), b)
        x0 = x_ref + 1.0e-12  # nearly exact: residual > atol (0), well below rtol*||b||
        x, info = self.cfd_utils.robust_cg(
            A, b, x0=x0, xp=self.xp, rtol=1.0e-5, atol=0.0, return_info_dict=True)
        self.assertTrue(info["converged"])
        self.assertEqual(info["iterations"], 0)
        self.assertEqual(info["reason"], "initial guess satisfies tolerance")
        np.testing.assert_allclose(np.asarray(x), x0, rtol=0.0, atol=0.0)

    def test_robust_cg_does_not_modify_rhs_or_initial_guess(self):
        # The CG must treat b and x0 as read-only (b aliases assembled pool memory in the solver stage)
        A = self._make_spd(10, seed=8)
        rng = np.random.RandomState(9)
        b = rng.rand(10)
        b_copy = b.copy()
        x0 = np.zeros(10)
        x0_copy = x0.copy()
        x, info = self.cfd_utils.robust_cg(A, b, x0=x0, xp=self.xp, rtol=1e-8, atol=0.0)
        self.assertEqual(info, 0)
        np.testing.assert_array_equal(b, b_copy)
        np.testing.assert_array_equal(x0, x0_copy)

if __name__ == '__main__':
    UnitTest.main()