//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Mohammad R. Hashemi
//
//

// System includes


// External includes


// Project includes
#include "elements/levelset_convection_element_simplex_algebraic_stabilization.h"


namespace Kratos
{
    template< unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::LevelSetConvectionElementSimplexAlgebraicStabilization() 
    : LevelSetConvectionElementSimplex<TDim, TNumNodes>::LevelSetConvectionElementSimplex()
    {}

    template< unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::LevelSetConvectionElementSimplexAlgebraicStabilization(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : LevelSetConvectionElementSimplex<TDim, TNumNodes>::LevelSetConvectionElementSimplex(NewId, pGeometry)
    {}

    template< unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::LevelSetConvectionElementSimplexAlgebraicStabilization(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : LevelSetConvectionElementSimplex<TDim, TNumNodes>::LevelSetConvectionElementSimplex(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    template< unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::~LevelSetConvectionElementSimplexAlgebraicStabilization() {}

    template< unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new LevelSetConvectionElementSimplexAlgebraicStabilization(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    template< unsigned int TDim, unsigned int TNumNodes>
    void LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!

        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;
        const double theta = rCurrentProcessInfo.Has(TIME_INTEGRATION_THETA) ? rCurrentProcessInfo[TIME_INTEGRATION_THETA] : 0.5;

        auto p_conv_diff_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const auto& rUnknownVar = p_conv_diff_settings->GetUnknownVariable();
        const auto& rConvVar = p_conv_diff_settings->GetConvectionVariable();
        const auto& rGradVar = p_conv_diff_settings->GetGradientVariable();

        //getting data for the given geometry
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        double Volume;

        const auto& r_geom = this->GetGeometry();
        GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, Volume);

        //here we get all the variables we will need
        array_1d<double,TNumNodes> phi, phi_old;
        array_1d< array_1d<double,3 >, TNumNodes> v, vold;

        array_1d<double,3 > X_mean_tmp = ZeroVector(3);
        array_1d<double,3 > grad_phi_mean_tmp = ZeroVector(3);
        array_1d< array_1d<double,3 >, TNumNodes> X_node;
        double phi_mean_old = 0.0;
        double phi_mean = 0.0;

        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            const auto& r_node = r_geom[i];

            phi[i] = r_node.FastGetSolutionStepValue(rUnknownVar);
            phi_old[i] = r_node.FastGetSolutionStepValue(rUnknownVar,1);

            v[i] = r_node.FastGetSolutionStepValue(rConvVar);
            vold[i] = r_node.FastGetSolutionStepValue(rConvVar,1);

            X_mean_tmp += r_node.Coordinates();
            X_node[i] = r_node.Coordinates();

            grad_phi_mean_tmp += r_node.GetValue(rGradVar);

            phi_mean_old += r_node.FastGetSolutionStepValue(rUnknownVar,1);
            phi_mean += r_node.FastGetSolutionStepValue(rUnknownVar);
        }

        phi_mean /= static_cast<double>(TNumNodes);
        phi_mean_old /= static_cast<double>(TNumNodes);

        array_1d<double,TDim> X_mean, grad_phi_mean;
        for(unsigned int k = 0; k < TDim; k++)
        {
            grad_phi_mean[k] = grad_phi_mean_tmp[k]/static_cast<double>(TNumNodes);
            X_mean[k] = X_mean_tmp[k]/static_cast<double>(TNumNodes);
        }

        array_1d<double,TDim> grad_phi_diff = prod(trans(DN_DX), phi_old) - grad_phi_mean;

        BoundedMatrix<double,TNumNodes, TNumNodes> K_matrix = ZeroMatrix(TNumNodes, TNumNodes); // convection
        BoundedMatrix<double,TNumNodes, TNumNodes> S_matrix = ZeroMatrix(TNumNodes, TNumNodes); // LHS stabilization
        Vector S_vector = ZeroVector(TNumNodes); // RHS stabilization
        BoundedMatrix<double,TNumNodes, TNumNodes> Mc_matrix = ZeroMatrix(TNumNodes, TNumNodes); // consistent mass matrix
        BoundedMatrix<double,TNumNodes, TNumNodes> Ml_matrix = IdentityMatrix(TNumNodes, TNumNodes); // lumped mass matrix

        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        this->GetShapeFunctionsOnGauss(Ncontainer);

        array_1d<double, TDim > vel_gauss;
        array_1d<double, TDim > X_gauss;
        array_1d<double, TNumNodes > v_dot_grad_N;

        for(unsigned int igauss=0; igauss<TDim+1; ++igauss)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity/coordinate at the gauss point
            vel_gauss = ZeroVector(TDim);
            X_gauss = ZeroVector(TDim);
            double phi_gauss = 0.0;
            double phi_gauss_old = 0.0;

            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                for(unsigned int k=0; k<TDim; ++k)
                {
                    vel_gauss[k] += N[i]*v[i][k];
                    X_gauss[k] += N[i]*X_node[i][k];
                }
                phi_gauss += N[i]*phi[i];
                phi_gauss_old += N[i]*phi_old[i];
            }

            v_dot_grad_N = prod(DN_DX, vel_gauss);

            noalias(Mc_matrix) += outer_prod(N, N);
            noalias(K_matrix) += outer_prod(N, v_dot_grad_N);

            for (unsigned int i = 0; i < TNumNodes; ++i){
                S_vector[i] += ( (phi_gauss_old - phi_mean_old) - inner_prod( grad_phi_mean, (X_gauss - X_mean) ) )*N[i];
            }
        }

        noalias(S_matrix) = (1.0/static_cast<double>(TNumNodes))*(Ml_matrix-Mc_matrix);

        double nu_e = 0.0;
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                if (i != j){
                    nu_e = std::max( nu_e, std::max(0.0, K_matrix(i, j))/std::abs(S_matrix(i, j)) );
                }
            }
        }

        const double limiter = this->GetValue(LIMITER_COEFFICIENT);
        if (limiter > 1.0e-15){
            noalias(rLeftHandSideMatrix)  = dt_inv*((1.0-limiter)*Ml_matrix + limiter*Mc_matrix) + theta*(K_matrix + (1.0-limiter)*nu_e*S_matrix);
            noalias(rRightHandSideVector) = prod( dt_inv*((1.0-limiter)*Ml_matrix + limiter*Mc_matrix) - (1.0 - theta)*(K_matrix + (1.0-limiter)*nu_e*S_matrix) , phi_old) - limiter*nu_e*S_vector;
        } else {
            noalias(rLeftHandSideMatrix)  = dt_inv*Ml_matrix + theta*(K_matrix + nu_e*S_matrix);
            noalias(rRightHandSideVector) = prod( dt_inv*Ml_matrix - (1.0 - theta)*(K_matrix + nu_e*S_matrix) , phi_old);
        }

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Levelset Element")
    }

    // template< unsigned int TDim, unsigned int TNumNodes>
    // void LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::EquationIdVector(
    //     EquationIdVectorType& rResult,
    //     const ProcessInfo& rCurrentProcessInfo) const
    // {
    //     KRATOS_TRY

    //     ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    //     const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    //     if (rResult.size() != TNumNodes)
    //         rResult.resize(TNumNodes, false);

    //     for (unsigned int i = 0; i < TNumNodes; i++)
    //     {
    //         rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
    //     }
    //     KRATOS_CATCH("")
    // }

    // template< unsigned int TDim, unsigned int TNumNodes>
    // void LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::GetDofList(
    //     DofsVectorType& ElementalDofList,
    //     const ProcessInfo& rCurrentProcessInfo) const
    // {
    //     KRATOS_TRY

    //     ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    //     const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    //     if (ElementalDofList.size() != TNumNodes)
    //         ElementalDofList.resize(TNumNodes);

    //     for (unsigned int i = 0; i < TNumNodes; i++)
    //     {
    //         ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

    //     }
    //     KRATOS_CATCH("");
    // }

    // //gauss points for the 3D case
    // template< unsigned int TDim, unsigned int TNumNodes>
    // void LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::GetShapeFunctionsOnGauss(BoundedMatrix<double,4,4>& Ncontainer)
    // {
    //     Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
    //     Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
    //     Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
    //     Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    // }

    // //gauss points for the 2D case
    // template< unsigned int TDim, unsigned int TNumNodes>
    // void LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::GetShapeFunctionsOnGauss(BoundedMatrix<double,3,3>& Ncontainer)
    // {
    //     const double one_sixt = 1.0/6.0;
    //     const double two_third = 2.0/3.0;
    //     Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
    //     Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
    //     Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    // }

    template class LevelSetConvectionElementSimplexAlgebraicStabilization<2, 3>;
    template class LevelSetConvectionElementSimplexAlgebraicStabilization<3, 4>;

} // namespace Kratos.
