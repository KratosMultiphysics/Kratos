//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Uxue Chasco
//
//


// System includes

// External includes

// Project includes
#include "elements/levelset_convection_element_simplex_bdf.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{

   template< unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::LevelSetConvectionElementSimplexBDF()
    : LevelSetConvectionElementSimplex<TDim, TNumNodes>::LevelSetConvectionElementSimplex()
    {}

    template< unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::LevelSetConvectionElementSimplexBDF(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : LevelSetConvectionElementSimplex<TDim, TNumNodes>::LevelSetConvectionElementSimplex(NewId, pGeometry)
    {}

    template< unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::LevelSetConvectionElementSimplexBDF(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : LevelSetConvectionElementSimplex<TDim, TNumNodes>::LevelSetConvectionElementSimplex(NewId, pGeometry, pProperties)
    {}
    /// Destructor.
    template <unsigned int TDim, unsigned int TNumNodes>
    LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::~LevelSetConvectionElementSimplexBDF()
    {}

    template <unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::Create(
        IndexType NewId,
        NodesArrayType const &ThisNodes,
        PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new LevelSetConvectionElementSimplexBDF(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer
    LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new LevelSetConvectionElementSimplexBDF(NewId, pGeom, pProperties));
        KRATOS_CATCH("");
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!
        rLeftHandSideMatrix.clear();

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!
        rRightHandSideVector.clear();

        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;
        const Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];

        const ConvectionDiffusionSettings::Pointer& my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        // const Variable<double> &rVolumeSourceVar = my_settings->GetVolumeSourceVariable(); // To be activated for the a_n source term
       
        const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();
        //getting data for the given geometry
        double Volume;
        array_1d<double, TNumNodes > N;
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);
        double h = this->ComputeH(DN_DX, Volume);
        //here we get all the variables we will need
        array_1d<double, TNumNodes> phi, phi_old, phi_old_1, a_n, proj_oss, phi_frac, phi_frac_old;
        array_1d< array_1d<double,3 >, TNumNodes> v;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            // phi_frac[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            // phi_frac_old[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
            phi[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);



            phi_old[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar, 1);
            phi_old_1[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar, 2);
            a_n[i] = this->GetGeometry()[i].FastGetSolutionStepValue(HEAT_FLUX);
            v[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rConvVar);


            // proj_oss[i] = this->GetGeometry()[i].GetValue(VOLUMETRIC_STRAIN_PROJECTION);
        }

        // compute the gradient of the unknown variable
        array_1d<double,TDim> grad_phi = prod(trans(DN_DX), phi);

        //compute the divergence of v
        double div_v = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
            for(unsigned int k=0; k<TDim; k++)
                div_v += DN_DX(i,k)*v[i][k];



        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        this->GetShapeFunctionsOnGauss(Ncontainer);
        for(unsigned int igauss=0; igauss<TDim+1; igauss++)
        {



            noalias(N) = row(Ncontainer,igauss);


            //velocity interpolation
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                 for(unsigned int k=0; k<TDim; k++){
                    vel_gauss[k] += N[i]*v[i][k];
                 }
            }

            const double norm_vel = norm_2(vel_gauss);

            // convection operator
            array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

                // calculate stabilization operator
                double tau = 0.0;
            if (this->IsElementTauNodal())
            {

                for (unsigned int i = 0; i < TNumNodes; i++)
                {
                    tau += N[i] * this->GetGeometry()[i].GetValue(TAU);
                }
            }
            else{
                const double tau_denom = std::max(dyn_st_beta * dt_inv + 2.0 * norm_vel / h, 1e-2);
                tau = 1.0 / (tau_denom);
            }

            // double proj_oss=0.0;

            // for (unsigned int i = 0; i < TNumNodes; i++)
            // {
            //     proj_oss += N[i]*this->GetGeometry()[i].GetValue(VOLUMETRIC_STRAIN_PROJECTION);
            // }

            // assemble RHS and LHS Gaus point contribution

             for (unsigned int i = 0; i < TNumNodes; i++)
                {
                    for (unsigned int j = 0; j < TNumNodes; j++)
                    {
                    // // // Mass terms
                    rRightHandSideVector[i] -= BDFcoeffs[0] * N[i] * N[j] * phi[j];
                    rRightHandSideVector[i] -= BDFcoeffs[1] * N[i] * N[j] * phi_old[j];
                    rRightHandSideVector[i] -= BDFcoeffs[2] * N[i] * N[j] * phi_old_1[j];
                    rLeftHandSideMatrix(i,j) += BDFcoeffs[0] * N[i] * N[j];
                    // // // Convective term
                    rRightHandSideVector[i] -= N[i] * a_dot_grad[j] * phi[j];
                    rLeftHandSideMatrix(i,j) += N[i] * a_dot_grad[j];

                    // Fractional acceleration
                    rRightHandSideVector[i] += N[i] * N[j] * a_n[j];

                    // Stabilization terms

                        for (unsigned int k = 0; k < TDim; ++k)
                        {
                            // ASGS stabilization
                            rRightHandSideVector[i] += tau * (DN_DX(i, k) * vel_gauss[k]) * -(BDFcoeffs[0] * N[j] * phi[j] + BDFcoeffs[1] * N[j] * phi_old[j] + BDFcoeffs[2] * N[j]* phi_old_1[j]);
                            rLeftHandSideMatrix(i, j) -= tau * (DN_DX(i, k) * vel_gauss[k]) * -((BDFcoeffs[0] * N[j]));

                            for (unsigned int l = 0; l < TDim; ++l)
                            {
                                // OSS stabilization
                                rRightHandSideVector[i] += tau * (DN_DX(i, k) * vel_gauss[k]) * -(DN_DX(j, l) * phi[j] * vel_gauss[l]);
                                // rRightHandSideVector[i] -= tau * (DN_DX(i, k) * vel_gauss[k]) * (DN_DX(j, l) * phi[j] * vel_gauss[l] - proj_oss);
                                rLeftHandSideMatrix(i, j) -= tau * (DN_DX(i, k) * vel_gauss[k]) * +(DN_DX(j, l) * vel_gauss[l]);
                            }

                        }

                        // rRightHandSideVector[i] += N[i] * div_v * phi_s_rhs;
                        // rLeftHandSideMatrix(i, j) -= N[i] * div_v * phi_s_lhs;

                    }
                }


        }
        // multiply by the integration weight (constant in the element)
        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

            KRATOS_CATCH("Error in Levelset Element")
        }
    // template <unsigned int TDim, unsigned int TNumNodes>
    // void LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::InitializeSolutionStep(
    //     const ProcessInfo &rCurrentProcessInfo)
    // {
    //     const ConvectionDiffusionSettings::Pointer &my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    //     BoundedMatrix<double, TNumNodes, TNumNodes> Ncontainer;
    //     this->GetShapeFunctionsOnGauss(Ncontainer);
    //     double Volume;
    //     array_1d<double, TNumNodes> N;
    //     BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    //     GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);
    //     array_1d<double, TNumNodes> phi;
    //     array_1d<array_1d<double, 3>, TNumNodes>v;
    //     const Variable<double> &rUnknownVar = my_settings->GetUnknownVariable();
    //     const Variable<array_1d<double, 3>> &rConvVar = my_settings->GetConvectionVariable();


    //     for (unsigned int i = 0; i < TNumNodes; i++)
    //     {
    //         phi[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
    //         v[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
    //     }

    //     array_1d<double, TDim> grad_phi = prod(trans(DN_DX), phi);
    //     // array_1d<double, TNumNodes> proj;


    //     double a_dot_grad_phi_gauss = 0.0;
    //     for (unsigned int igauss = 0; igauss < TDim + 1; igauss++)
    //     {
    //         noalias(N) = row(Ncontainer, igauss);
    //         // velocity interpolation
    //         array_1d<double, TDim> vel_gauss = ZeroVector(TDim);
    //         for (unsigned int i = 0; i < TNumNodes; i++)
    //         {
    //             for (unsigned int k = 0; k < TDim; k++)
    //                 vel_gauss[k] += N[i] * v[i][k];
    //         }
    //         // convection operator
    //         a_dot_grad_phi_gauss = inner_prod(vel_gauss, grad_phi);
    //         KRATOS_WATCH(a_dot_grad_phi_gauss)
    //         a_dot_grad_phi_gauss *= Volume / static_cast<double>(TNumNodes);
    //         auto &r_geometry = this->GetGeometry();
    //         for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    //         {
    //             AtomicAdd(r_geometry[i_node].GetValue(VOLUMETRIC_STRAIN_PROJECTION), N[i_node] * a_dot_grad_phi_gauss);
    //         }
    //     }
    //     KRATOS_WATCH(a_dot_grad_phi_gauss)
    // }

    template <unsigned int TDim, unsigned int TNumNodes>
    void LevelSetConvectionElementSimplexBDF<TDim, TNumNodes>::FinalizeSolutionStep(
        const ProcessInfo &rCurrentProcessInfo)
    {
            const ConvectionDiffusionSettings::Pointer &my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
            BoundedMatrix<double, TNumNodes, TNumNodes> Ncontainer;
            this->GetShapeFunctionsOnGauss(Ncontainer);
            double Volume;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);
            array_1d<double, TNumNodes> phi, phi_old, phi_old_1;
            array_1d<array_1d<double, 3>, TNumNodes>v;
            const Variable<double> &rUnknownVar = my_settings->GetUnknownVariable();
            const Variable<array_1d<double, 3>> &rConvVar = my_settings->GetConvectionVariable();
            // const Variable<double> &rVolumeSourceVar = my_settings->GetVolumeSourceVariable();
            const Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                phi[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,0);
                phi_old[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
                phi_old_1[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,2);
                v[i] = this->GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
            }
            // velocity interpolation
            array_1d<double, TDim> vel_gauss = ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                for (unsigned int k = 0; k < TDim; k++)
                {
                    vel_gauss[k] += N[i] * v[i][k];
                }
            }
            array_1d<double, TDim> grad_phi = prod(trans(DN_DX), phi);
            double fractional_acceleration_gauss = 0.0;
            for(unsigned int igauss=0; igauss<TDim+1; igauss++){

                    for (unsigned int j = 0; j < TNumNodes; j++)
                    {
                        fractional_acceleration_gauss = BDFcoeffs[0] * N[j] * phi[j] + BDFcoeffs[1] * N[j] * phi_old[j] + BDFcoeffs[2] * N[j] * phi_old_1[j];
                        
                        


                        for (unsigned int d = 0; d < TDim; d++)
                        {
                            fractional_acceleration_gauss += DN_DX(j, d) * phi[j] * vel_gauss[d];

                        }
                    }
                        fractional_acceleration_gauss *= Volume / static_cast<double>(TNumNodes);
                        // KRATOS_WATCH(fractional_acceleration_gauss)
                        auto &r_geometry = this->GetGeometry();
                        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
                        {
                            AtomicAdd(r_geometry[i_node].FastGetSolutionStepValue(HEAT_FLUX), N[i_node] * fractional_acceleration_gauss);
                        }
                    }
    }
    template class LevelSetConvectionElementSimplexBDF<2, 3>;
    template class LevelSetConvectionElementSimplexBDF<3, 4>;

    } // namespace Kratos.

            // a_dot_grad_phi_gauss *= Volume / static_cast<double>(TNumNodes);

            // Assembly the projection contributions
            // auto &r_geometry = this->GetGeometry();
            // for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            // {
            //     AtomicAdd(r_geometry[i_node].GetValue(VOLUMETRIC_STRAIN_PROJECTION),  N[i_node] * a_dot_grad_phi_gauss);
            // }
            // for (IndexType i_node = 0; i_node < TNumNodes; ++i_node){
            //     proj[i_node] += N[i_node] * a_dot_grad_phi_gauss;
            // }

                                // Convective term
                                // rRightHandSideVector[i] -= N[i] * DN_DX(j, l) * phi[j] * vel_gauss[l];
                                // rLeftHandSideMatrix(i, j) += N[i] * DN_DX(j, l) * vel_gauss[l];