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

    template <unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new LevelSetConvectionElementSimplexAlgebraicStabilization(NewId, pGeom, pProperties));
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
        const auto& r_unknown_var = p_conv_diff_settings->GetUnknownVariable();
        const auto& r_conv_var = p_conv_diff_settings->GetConvectionVariable();
        const auto& r_grad_var = p_conv_diff_settings->GetGradientVariable();

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
        array_1d< array_1d<double,3 >, TNumNodes> X_node;
        double phi_mean_old = 0.0;
        double phi_mean = 0.0;

        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            const auto& r_node = r_geom[i];

            phi[i] = r_node.FastGetSolutionStepValue(r_unknown_var);
            phi_old[i] = r_node.FastGetSolutionStepValue(r_unknown_var,1);

            v[i] = r_node.FastGetSolutionStepValue(r_conv_var);
            vold[i] = r_node.FastGetSolutionStepValue(r_conv_var,1);

            X_mean_tmp += r_node.Coordinates();
            X_node[i] = r_node.Coordinates();

            phi_mean_old += r_node.FastGetSolutionStepValue(r_unknown_var,1);
            phi_mean += r_node.FastGetSolutionStepValue(r_unknown_var);
        }

        const double aux_weight = 1.0/static_cast<double>(TNumNodes);

        phi_mean *= aux_weight;
        phi_mean_old *= aux_weight;

        array_1d<double,TDim> X_mean;
        for(unsigned int k = 0; k < TDim; k++)
        {
            X_mean[k] = aux_weight*X_mean_tmp[k];
        }

        BoundedMatrix<double,TNumNodes, TNumNodes> K_matrix = ZeroMatrix(TNumNodes, TNumNodes); // convection
        BoundedMatrix<double,TNumNodes, TNumNodes> S_matrix = ZeroMatrix(TNumNodes, TNumNodes); // LHS stabilization
        Vector S_vector = ZeroVector(TNumNodes); // RHS stabilization
        BoundedMatrix<double,TNumNodes, TNumNodes> Mc_matrix = ZeroMatrix(TNumNodes, TNumNodes); // consistent mass matrix
        BoundedMatrix<double,TNumNodes, TNumNodes> Ml_matrix = IdentityMatrix(TNumNodes, TNumNodes); // lumped mass matrix

        BoundedMatrix<double,TNumNodes, TNumNodes> M_matrix;
        BoundedMatrix<double,TNumNodes, TNumNodes> L_matrix;

        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        this->GetShapeFunctionsOnGauss(Ncontainer);

        array_1d<double, TDim > vel_gauss;
        array_1d<double, TDim > X_gauss;
        array_1d<double, TNumNodes > v_dot_grad_N;

        const double limiter = this->GetValue(LIMITER_COEFFICIENT);

        for(unsigned int igauss=0; igauss<TDim+1; ++igauss)
        {
            noalias(N) = row(Ncontainer,igauss);

            // Obtaining the velocity/coordinate at the gauss point
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

            // Consistent and lumped mass matrices
            noalias(Mc_matrix) += outer_prod(N, N);
            noalias(K_matrix) += outer_prod(N, v_dot_grad_N);

            // If the high-order terms are necessary
            if (limiter > 1.0e-15){
                array_1d<double,3 > grad_phi_mean_tmp = ZeroVector(3);

                for (unsigned int i = 0; i < TNumNodes; ++i){
                    const auto& r_node = r_geom[i];
                    grad_phi_mean_tmp += r_node.GetValue(r_grad_var);
                }

                array_1d<double,TDim> grad_phi_mean;
                for(unsigned int k = 0; k < TDim; k++) {
                    grad_phi_mean[k] = aux_weight*grad_phi_mean_tmp[k];
                }

                for (unsigned int i = 0; i < TNumNodes; ++i){
                    S_vector[i] += ( (phi_gauss_old - phi_mean_old) - inner_prod( grad_phi_mean, (X_gauss - X_mean) ) )*N[i];
                }
            }
        }

        // Coefficient of the artificial viscosity
        noalias(S_matrix) = aux_weight*(Ml_matrix-Mc_matrix);

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

        // Determining the necessity to add higher-order terms
        if (limiter > 1.0e-15){
            noalias(M_matrix)  = dt_inv*((1.0-limiter)*Ml_matrix + limiter*Mc_matrix);
            noalias(L_matrix) = K_matrix + (1.0-limiter)*nu_e*S_matrix;
        } else {
            noalias(M_matrix)  = dt_inv*Ml_matrix;
            noalias(L_matrix) = K_matrix + nu_e*S_matrix;
        }

        noalias(rLeftHandSideMatrix)  = M_matrix + theta*L_matrix;
        noalias(rRightHandSideVector) = prod( M_matrix - (1.0 - theta)*L_matrix , phi_old) - limiter*nu_e*S_vector;

        // Taking out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

        const double gauss_pt_weigth = Volume * aux_weight;

        rRightHandSideVector *= gauss_pt_weigth;
        rLeftHandSideMatrix *= gauss_pt_weigth;

        KRATOS_CATCH("Error in Levelset Element")
    }

    template class LevelSetConvectionElementSimplexAlgebraicStabilization<2, 3>;
    template class LevelSetConvectionElementSimplexAlgebraicStabilization<3, 4>;

} // namespace Kratos.
