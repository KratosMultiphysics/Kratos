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

    def test_Compute_N_DN(self):
        n_dn = self.cfd_utils.Compute_N_DN(self.N, self.DN, self.pelem)
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

    def test_Compute_DN_N(self):
        dn_n = self.cfd_utils.Compute_DN_N(self.N, self.DN, self.pelem)
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
        a_gauss = self.cfd_utils.InterpolateValue(self.Ngauss, self.velem)
        
        conv_contrib = self.cfd_utils.ComputeConvectiveContribution(self.Ngauss, grad_u, a_gauss)
        ref_conv_contrib = np.array([[[[0.333801136363636, 0.0831875        , 0.21725          ],
         [0.333801136363636, 0.0831875        , 0.21725          ],
         [0.333801136363636, 0.0831875        , 0.21725          ],
         [0.333801136363636, 0.0831875        , 0.21725          ]]],
       [[[0.358247596153846, 0.21775          , 0.1373125        ],
         [0.358247596153846, 0.21775          , 0.1373125        ],
         [0.358247596153846, 0.21775          , 0.1373125        ],
         [0.358247596153846, 0.21775          , 0.1373125        ]]],
       [[[0.0723515625     , 0.1784140625     , 0.418876358695652],
         [0.0723515625     , 0.1784140625     , 0.418876358695652],
         [0.0723515625     , 0.1784140625     , 0.418876358695652],
         [0.0723515625     , 0.1784140625     , 0.418876358695652]]],
       [[[0.0723515625     , 0.285432942708333, 0.24328125       ],
         [0.0723515625     , 0.285432942708333, 0.24328125       ],
         [0.0723515625     , 0.285432942708333, 0.24328125       ],
         [0.0723515625     , 0.285432942708333, 0.24328125       ]]],
       [[[0.135015625      , 0.0831875        , 0.3967703125     ],
         [0.135015625      , 0.0831875        , 0.3967703125     ],
         [0.135015625      , 0.0831875        , 0.3967703125     ],
         [0.135015625      , 0.0831875        , 0.3967703125     ]]],
       [[[0.15471875       , 0.323561079545455, 0.1373125        ],
         [0.15471875       , 0.323561079545455, 0.1373125        ],
         [0.15471875       , 0.323561079545455, 0.1373125        ],
         [0.15471875       , 0.323561079545455, 0.1373125        ]]]])
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
        v_el_gauss = self.cfd_utils.InterpolateValue(self.Ngauss, self.velem)
        tmp = self.rho * self.cfd_utils.ComputeConvectiveContribution(self.Ngauss, grad_v, v_el_gauss)
        convective = np.tensordot(tmp, self.w_int_order, axes=(1, 0))
        convective *= self.elemental_volumes[:, None, None]

        pi_conv = np.zeros((len(self.model_part.Nodes), self.dim), dtype=np.float64)
        self.cfd_utils.AssembleVector(self.connectivity, convective, pi_conv)
        pi_conv = pi_conv.ravel() * self.Minv
        pi_conv = pi_conv.reshape((len(self.model_part.Nodes), self.dim))
        
        conv_proj_el = self.ElemData(pi_conv, self.connectivity)
        proj_el_gauss = self.cfd_utils.InterpolateValue(self.Ngauss, conv_proj_el)

        stab_mom = self.cfd_utils.ComputeMomentumStabilization(self.Ngauss, self.DN, self.velem, v_el_gauss, proj_el_gauss, self.rho)
        # import sys
        # np.set_printoptions(threshold=sys.maxsize, precision=15, suppress=True)
        # print(repr(stab_mom))
        ref_stab_mom = np.array([[[[ 0.050596189905262, -0.20832754990043 , -0.144868568241734],
            [ 0.054494858342217, -0.224380142889461, -0.156031355707913],
            [ 0.044856311968763, -0.184693859113093, -0.128433973065697],
            [-0.149947360216242,  0.617401551902984,  0.429333897015344]]],
          [[[-0.231654375683096,  0.17712339862571 ,  0.605183295920595],
            [ 0.080819286464665, -0.061794587954208, -0.211135585125084],
            [ 0.070166956060473, -0.053649794837692, -0.183307004705507],
            [ 0.080668133157958, -0.061679015833809, -0.210740706090004]]],
          [[[-0.262040837789505, -0.128842037752164,  0.097860279702437],
            [-0.192091896104556, -0.094449061980442,  0.071737546101359],
            [-0.16838545520599 , -0.082792916400279,  0.06288427362424 ],
            [ 0.622518189100052,  0.306084016132886, -0.232482099428036]]],
          [[[ 0.576861502480971, -0.008252176665896,  0.352386389773009],
            [-0.156483634675516,  0.002238545288791, -0.095590887664899],
            [-0.245278326174092,  0.003508779960526, -0.149832747511013],
            [-0.175099541631363,  0.002504851416579, -0.106962754597097]]],
          [[[ 0.542507971090554,  0.580843603236126, -0.162912525002362],
            [-0.246538213440828, -0.263959521078567,  0.074034235442621],
            [-0.136962806043036, -0.146641107616258,  0.041129269527646],
            [-0.15900695160669 , -0.170242974541302,  0.047749020032096]]],
          [[[-0.204071235621765,  0.067305510572914, -0.253005676772859],
            [-0.117826073328397,  0.038860665492651, -0.146079702673982],
            [-0.177827213572146,  0.058649869820043, -0.220468575011921],
            [ 0.499724522522309, -0.164816045885608,  0.619553954458762]]]])
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

    def test_ComputePressureStabilization_ProjectionTerm(self):
        pi_press_el = self.cfd_utils.Compute_N_DN(self.N, self.DN, self.pelem)
        pi_press_el *= self.elemental_volumes[:, np.newaxis, np.newaxis]
        pi_press = np.zeros((len(self.model_part.Nodes), self.dim), dtype=np.float64)
        self.cfd_utils.AssembleVector(self.connectivity, pi_press_el, pi_press)
        pi_press = pi_press.ravel() * self.Minv
        pi_press = pi_press.reshape((len(self.model_part.Nodes), self.dim))
        
        pres_proj_el = self.ElemData(pi_press, self.connectivity)

        stab_press = self.cfd_utils.ComputePressureStabilization_ProjectionTerm(self.N, self.DN, pres_proj_el)
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