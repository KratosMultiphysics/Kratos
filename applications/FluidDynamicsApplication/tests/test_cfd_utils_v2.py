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

        model_part_io = KM.ModelPartIO("small_cube")
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
        a_gauss = self.cfd_utils.InterpolateValue(self.Ngauss, self.velem)
        ref_a_gauss = np.array([[[1.088125, 0.3025  , 0.725   ]],
                               [[1.140625, 0.725   , 0.4225  ]],
                               [[0.275625, 0.633125, 1.2725  ]],
                               [[0.275625, 0.993125, 0.7825  ]],
                               [[0.525625, 0.3025  , 1.215   ]],
                               [[0.578125, 1.085   , 0.4225  ]]])
        np.testing.assert_allclose(a_gauss, ref_a_gauss, rtol=self.rtol, atol=self.atol)

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
        
        ref_stab_mom = np.array([[[[ 0.489258731865069, -0.05842702811148 ,  0.173285926323712],
         [ 0.501237452436738, -0.069339076864483,  0.169898292853235],
         [ 0.471622865920407, -0.042361587078928,  0.178273424729533],
         [-0.126914504767669,  0.502877692054891,  0.34754235609352 ]]],
       [[[-0.199482133861028,  0.138283474021625,  0.490049969991781],
         [ 0.552827641046223,  0.245474181373478,  0.014249896282616],
         [ 0.527181152207667,  0.241820014737088,  0.030470078030581],
         [ 0.552463725222524,  0.245422329867809,  0.014480055695022]]],
       [[[-0.108974137986125,  0.153589758607585,  0.68741288766921 ],
         [-0.060571221264133,  0.160216335594546,  0.615729992117872],
         [-0.044166957569104,  0.162462153071081,  0.591435895797163],
         [ 0.503118566819363,  0.237388002726788, -0.219073340801637]]],
       [[[ 0.474275865313851,  0.006552606931646,  0.298988954304771],
         [-0.036677338714999,  0.361084045085899,  0.228169574127112],
         [-0.098544432151892,  0.404011321865176,  0.219594641675281],
         [-0.04964784444696 ,  0.370083796950612,  0.226371829892835]]],
       [[[ 0.438214998612516,  0.471293858838921, -0.146650804173918],
         [-0.002770796343348, -0.093184191167652,  0.643723497319555],
         [ 0.058469216074463, -0.01479472105911 ,  0.533963653370334],
         [ 0.046149081656369, -0.030564946612159,  0.556044903484029]]],
       [[[ 0.051749226383825,  0.511550348392397, -0.012760831328323],
         [ 0.095266497441111,  0.432101792344222,  0.050663584057229],
         [ 0.064991339805582,  0.487374503611141,  0.006538938207743],
         [ 0.406867936369482, -0.136782326165941,  0.504808309063351]]]])
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