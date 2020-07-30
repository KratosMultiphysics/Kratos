//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

// Application includes
#include "custom_elements/transient_convection_diffusion_FIC_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientConvectionDiffusionFICElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new TransientConvectionDiffusionFICElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientConvectionDiffusionFICElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new TransientConvectionDiffusionFICElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int TransientConvectionDiffusionFICElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    int ierr = SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    return ierr;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim,TNumNodes>::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    GeometryData::IntegrationMethod ThisIntegrationMethod = SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    //Containers of variables at all integration points
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

    //Initializing variables
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    array_1d<double,TNumNodes> NodalPhi;
    BoundedMatrix<double,TNumNodes,TDim> GradNT;
    array_1d<double,TDim> GradPhi;

    //Nodal Variables
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        NodalPhi[i] = Geom[i].FastGetSolutionStepValue(TEMPERATURE,1);
    }

    //Extrapolation variables
    Matrix GradPhiContainer(NumGPoints,TDim);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        noalias(GradNT) = DN_DXContainer[GPoint];

        noalias(GradPhi) = prod(trans(GradNT), NodalPhi);

        this->SaveGPGradPhi(GradPhiContainer,GradPhi,GPoint);
    }

    this->ExtrapolateGPValues(GradPhiContainer);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim, TNumNodes>::CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                        VectorType& rRightHandSideVector,
                        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes;

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != element_size )
        rLeftHandSideMatrix.resize( element_size, element_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    GeometryData::IntegrationMethod ThisIntegrationMethod = SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( ThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,Geom,Prop,rCurrentProcessInfo);

    Variables.IterationNumber = rCurrentProcessInfo[NL_ITERATION_NUMBER];

    BoundedMatrix<double,TNumNodes,TNumNodes> MMatrixAux = ZeroMatrix( TNumNodes, TNumNodes );

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNT
        noalias(Variables.GradNT) = DN_DXContainer[GPoint];

        //Compute N and Interpolated velocity
        noalias(Variables.N) = row(NContainer,GPoint);

        Variables.QSource = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            Variables.QSource += Variables.N[i]*Variables.NodalQSource[i];
        }

        noalias(Variables.VelInter) = ZeroVector(TDim);
        ElementUtilities::InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);

        noalias(Variables.DifMatrix) = ZeroMatrix( TDim, TDim );

        this->CalculateDiffusivityVariables(Variables,Prop,rCurrentProcessInfo);
        this->CalculateHVector(Variables,Prop,rCurrentProcessInfo);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        array_1d <double, TNumNodes> AuxMVector;
        array_1d <double, TNumNodes> AuxMVector1;
        array_1d <double, TNumNodes> AuxMVector2;

        noalias(AuxMVector1) = Variables.rho_dot_c * Variables.N;
        noalias(AuxMVector2) = Variables.rho_dot_c * 0.5 * prod(Variables.GradNT,Variables.HVector);

        //// M matrix
        BoundedMatrix<double,TNumNodes,TNumNodes> MMatrixAux1 = ZeroMatrix( TNumNodes, TNumNodes );
        BoundedMatrix<double,TNumNodes,TNumNodes> MMatrixAux2 = ZeroMatrix( TNumNodes, TNumNodes );

        MMatrixAux1 = outer_prod(AuxMVector1,Variables.N) * Variables.IntegrationCoefficient;
        MMatrixAux2 = outer_prod(AuxMVector2,Variables.N) * Variables.IntegrationCoefficient;

        // We are not considering MMatrixAux2, which is the h term
        MMatrixAux += MMatrixAux1 ;

        /////////////////////////////////////////////////////////////////////////
        // To make the M matrix consistent uncomment the next 2 lines and comment the for loop next

        // noalias(AuxMVector) = AuxMVector1 + AuxMVector2;
        // noalias(rLeftHandSideMatrix) += outer_prod(AuxMVector,Variables.N) * Variables.IntegrationCoefficient;
    }

    for (unsigned int i = 0 ; i < TNumNodes ; i++ )
    {
        for (unsigned int j = 0 ; j < TNumNodes ; j ++ )
        {
            // LHS = Md lumped
            rLeftHandSideMatrix (i,i) += MMatrixAux(i,j);
        }
    }

    // We are using a lumped M matrix for the cases in which we reach a stationary. In the case of fully transient
    // problems we use a consistent matrix with the temporal stabilization term AuxMVector2

    //RightHandSideVector

    const double& Theta = rCurrentProcessInfo[THETA];
    const double& DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    array_1d<double,TNumNodes> NodalPhi0;

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        NodalPhi0[i] = Geom[i].FastGetSolutionStepValue(TEMPERATURE,1);
    }

    array_1d<double,TNumNodes> aux_vector;

    noalias(aux_vector) = Variables.NodalPhi - NodalPhi0;

    noalias(rRightHandSideVector) -= 1.0 / (Theta*DeltaTime) * prod(rLeftHandSideMatrix, aux_vector);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					      ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "TransientConvectionDiffusionFICElement::CalculateFirstDerivativesLHS not implemented" << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
					      ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int element_size = TNumNodes;

    MatrixType AuxLHS;

    //Resetting the LHS
    if ( AuxLHS.size1() != element_size )
        AuxLHS.resize( element_size, element_size, false );
    noalias( AuxLHS ) = ZeroMatrix( element_size, element_size );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    this-> CalculateFirstDerivativesContributions(AuxLHS,rRightHandSideVector, rCurrentProcessInfo);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateDiffusivityVariables(ElementVariables& rVariables, const PropertiesType& Prop,
                                                                                        const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    // GeometryType& rGeom = this->GetGeometry();
    const Geometry<Node<3> >& rGeom = this->GetGeometry();

    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    //Properties variables
    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    double conductivity = Prop[rDiffusionVar];

    const Variable<double>& rReactionVar = my_settings->GetReactionVariable();
    rVariables.absorption = Prop[rReactionVar];

    ///////////////////////////////////////////////////////////////////////////
    // Adding transient absorption
    ///////////////////////////////////////////////////////////////////////////

    array_1d<double,TNumNodes> NodalPhi0;
    array_1d<double,TNumNodes> PrevNodalPhi;

    double sum_phi = 0.0;
    double sustr_phi = 0.0;
    double theta = CurrentProcessInfo.GetValue(THETA);

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        NodalPhi0[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE,1);
        PrevNodalPhi[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE,2);

        // double aux_var = rVariables.NodalPhi[i] - NodalPhi0[i];
        // double aux_var2 = rVariables.NodalPhi[i] + NodalPhi0[i];

        double aux_var = NodalPhi0[i] - PrevNodalPhi[i];
        double aux_var2 = NodalPhi0[i] + PrevNodalPhi[i];

        // TODO: try std::abs
        if(aux_var > sustr_phi)
        {
            sustr_phi = aux_var;
        }

        if(aux_var2 > sum_phi)
        {
            sum_phi = aux_var2;
        }

        if(std::abs(sum_phi) < 1e-5)
        {
            sum_phi = 1e-5;
        }
    }

    double delta_time = CurrentProcessInfo.GetValue(DELTA_TIME);

    //TODO
    double previous_absorption = rVariables.absorption;

    // beta = 1
    // 2.626 = 2 / tanh(1)

    double fk = 2.626 * tanh(1.0 * (sustr_phi) / (sum_phi));

    double st = rVariables.rho_dot_c / (theta * delta_time) * fk;

    rVariables.TransientAbsorption = rVariables.absorption + st;

    if (std::abs(st) > 0.1 * std::abs(previous_absorption))
    {
       rVariables.TransientAbsorption = previous_absorption + 0.1 * previous_absorption;
    }

    // If absorption = 0; no transient absorption is added
    if (std::abs(previous_absorption) < rVariables.LowTolerance)
    {
        rVariables.TransientAbsorption = previous_absorption;
    }

    //rVariables.TransientAbsorption = previous_absorption;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    double NormVel = norm_2(rVariables.VelInter);

    noalias(rVariables.GradPhi) = ZeroVector(TDim);

    rVariables.GradPhi = prod(trans(rVariables.GradNT), NodalPhi0);

    for(unsigned int i = 0; i < TDim; i++)
    {
        noalias(rVariables.SecondGradPhi[i]) = prod(trans(rVariables.GradNT), rVariables.NodalPhiGradient[i]);
    }

    rVariables.NormGradPhi = norm_2(rVariables.GradPhi);

    // Unitary velocity vector
    if (NormVel < rVariables.HighTolerance)
    {
        for (unsigned int i = 0; i < TDim; i++)
        {
            rVariables.VelInterHat[i] = 0.0;
        }
    }
    else
    {
        for (unsigned int i = 0; i < TDim; i++)
        {
            rVariables.VelInterHat[i] = rVariables.VelInter[i]/NormVel;
        }
    }

    //////////////////////////////////////////////////////
    // Calculate Dv
    //////////////////////////////////////////////////////

    rVariables.AuxDiffusion = inner_prod(prod(trans(rVariables.VelInterHat), rVariables.DifMatrixK), rVariables.VelInterHat);

    if (NormVel < rVariables.HighTolerance)
    {
        rVariables.AuxDiffusion = conductivity;
    }

    double Domain = rGeom.DomainSize();

    if (TDim == 2)
    {
        rVariables.lv = std::sqrt(2.0*Domain);
        rVariables.lsc = rVariables.lv;

        if (TNumNodes == 3)
        {
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                array_1d <double, 3> AuxNodeNormal = rGeom[i].FastGetSolutionStepValue(NORMAL);
                double NormAuxNodeNormal = norm_2 (AuxNodeNormal);

                if (NormAuxNodeNormal > rVariables.LowTolerance)
                {
                    rVariables.lsc = std::sqrt(2.0) * rVariables.lv;
                }
            }
        }
    }
    else
    {
        rVariables.lv = pow( (6.0*Domain/Globals::Pi) , (1.0/3.0) );
        rVariables.lsc = rVariables.lv;
    }

    if (NormVel > rVariables.HighTolerance)
    {
        array_1d <double, 3> Velocity3;
        ElementUtilities::FillArray1dOutput (Velocity3, rVariables.VelInterHat);

        rVariables.lv = ElementSizeCalculator<TDim,TNumNodes>::ProjectedElementSize(rGeom,Velocity3);
    }

    this->CalculateBoundaryLv(rVariables);

    rVariables.AlphaVBar = 0.0;
    rVariables.OmegaV = 0.0;
    rVariables.OmegaVt = 0.0;
    rVariables.OmegaVBarra = 0.0;
    rVariables.SigmaV = 0.0;
    rVariables.SigmaVt = 0.0;
    rVariables.SigmaVBarra = 0.0;
    rVariables.Peclet = 0.0;

    if (NormVel < rVariables.HighTolerance)
    {

        // TODO: S'ha de posar OmegaV = 0 si v = 0??
        rVariables.OmegaV = rVariables.absorption * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / conductivity;
        rVariables.OmegaVt = st * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / conductivity;

        rVariables.OmegaVBarra = rVariables.OmegaV + rVariables.OmegaVt;

        rVariables.SigmaV = rVariables.OmegaV / (2.0 * rVariables.HighTolerance);
        rVariables.SigmaVt = rVariables.rho_dot_c * rVariables.lv / (theta * delta_time * NormVel) * fk;

        rVariables.SigmaVBarra = rVariables.SigmaV + rVariables.SigmaVt;

        rVariables.Peclet = NormVel * rVariables.lv * rVariables.rho_dot_c / (2.0 * rVariables.AuxDiffusion);
    }
    else
    {
        if (conductivity < rVariables.HighTolerance)
        {
            rVariables.Peclet = 1.0;
        }
        else
        {
            rVariables.Peclet = NormVel * rVariables.lv * rVariables.rho_dot_c / (2.0 * rVariables.AuxDiffusion);
        }

        rVariables.OmegaV = rVariables.absorption * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / rVariables.AuxDiffusion;
        rVariables.OmegaVt = st * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / rVariables.AuxDiffusion;

        rVariables.OmegaVBarra = rVariables.OmegaV + rVariables.OmegaVt;

        rVariables.SigmaV = rVariables.OmegaV / (2.0 * rVariables.Peclet);
        rVariables.SigmaVt = rVariables.rho_dot_c * rVariables.lv / (theta * delta_time * NormVel) * fk;

        rVariables.SigmaVBarra = rVariables.SigmaV + rVariables.SigmaVt;

        rVariables.AlphaVBar = 1.0 / tanh(rVariables.Peclet) - 1.0 / rVariables.Peclet;
    }

    rVariables.LambdaV = std::sqrt(rVariables.Peclet * rVariables.Peclet + rVariables.OmegaVBarra);

    double epsilon = rVariables.LambdaV - rVariables.Peclet;

    // rVariables.XiV = (cosh(rVariables.LambdaV) / cosh(rVariables.Peclet));

    rVariables.XiV = (exp(epsilon) + exp(- 2 * rVariables.Peclet - epsilon)) / (1.0 + exp(- 2 * rVariables.Peclet));

    if(rVariables.SigmaVBarra < 0.00024414) // 2^-12
    {
        rVariables.AlphaV = rVariables.SigmaVBarra / 3.0 + rVariables.AlphaVBar * (1.0 - rVariables.SigmaVBarra / rVariables.Peclet);
    }
    else
    {
        rVariables.AlphaV = 2.0 / rVariables.SigmaVBarra * (1.0 - (rVariables.SigmaVBarra * tanh (rVariables.Peclet)) / (rVariables.XiV - 1.0));
    }

    if (std::abs(rVariables.TransientAbsorption) < rVariables.HighTolerance)
    {
        rVariables.AlphaV = 1.0 / tanh(rVariables.Peclet) - 1.0 / rVariables.Peclet;
    }

    if (rVariables.Peclet < rVariables.LowTolerance)
    {
        rVariables.AlphaV = 0.0;
    }

    if (conductivity < rVariables.HighTolerance)
    {
        rVariables.AlphaV = 1.0;
    }

    //////////////////////////////////////////////////////
    // Calculate Ds
    //////////////////////////////////////////////////////

    BoundedMatrix<double,TDim,TDim> BaricenterMatrix;
    noalias(BaricenterMatrix) = ZeroMatrix(TDim,TDim);


    for(unsigned int i=0; i<TNumNodes; i++)
    {
        array_1d <double, TDim> BarNodeVector;
        array_1d <double, 3> BarNodeVector3D;

        noalias(BarNodeVector3D) = rGeom[i] - rGeom.Center();

        for(unsigned int j=0; j<TDim; j++)
        {
            BarNodeVector [j] = BarNodeVector3D [j];
        }

        noalias(BaricenterMatrix) += outer_prod(BarNodeVector, BarNodeVector);
    }

    noalias(rVariables.DifMatrixS) = (rVariables.absorption / (TNumNodes + 1)) * BaricenterMatrix;

    //////////////////////////////////////////////////////
    // Calculate residuals
    //////////////////////////////////////////////////////

    this->CalculateDifTerm(rVariables);

    // noalias(rVariables.DifMatrixAux) = prod(rVariables.GradNT,rVariables.DifMatrixK);
    // noalias(rVariables.MatrixAux) = prod(rVariables.DifMatrixAux,trans(rVariables.GradNT));

    // array_1d<double,TNumNodes> NormAux1 = ZeroVector(TNumNodes);

    // for (unsigned int i = 0 ; i < TNumNodes ; i++ )
    // {
    //     NormAux1 [i] = rVariables.MatrixAux (i,i) * NodalPhi0 [i];
    // }

    // double NormAux2 = norm_2(NormAux1);

    rVariables.Residual = rVariables.rho_dot_c * inner_prod (rVariables.VelInter , rVariables.GradPhi)
                            - rVariables.NormAux2 * conductivity
                            + rVariables.absorption * inner_prod(rVariables.N, NodalPhi0)
                            - rVariables.QSource;

    rVariables.TransientResidual = rVariables.Residual + rVariables.rho_dot_c * inner_prod(rVariables.N, NodalPhi0 - PrevNodalPhi) * 1.0 / (theta * delta_time);

    //////////////////////////////////////////////////////
    // Calculate Dr
    //////////////////////////////////////////////////////

    // phi = 2 in 2D and 3D

    //TODO until +40 lines
    if (NormVel < rVariables.HighTolerance)
    {

        // TODO: S'ha de posar OmegaV = 0 si v = 0??
        rVariables.OmegaV = rVariables.absorption * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / conductivity;
        rVariables.OmegaVt = st * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / conductivity;

        rVariables.OmegaVBarra = rVariables.OmegaV + rVariables.OmegaVt;

        rVariables.SigmaV = rVariables.OmegaV / (2.0 * rVariables.HighTolerance);
        rVariables.SigmaVt = rVariables.rho_dot_c * rVariables.lv / (theta * delta_time * NormVel) * fk;

        rVariables.SigmaVBarra = rVariables.SigmaV + rVariables.SigmaVt;

    }
    else
    {
        rVariables.OmegaV = rVariables.absorption * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / rVariables.AuxDiffusion;
        rVariables.OmegaVt = st * rVariables.lv * rVariables.lv * rVariables.rho_dot_c / rVariables.AuxDiffusion;

        rVariables.OmegaVBarra = rVariables.OmegaV + rVariables.OmegaVt;

        rVariables.SigmaV = rVariables.OmegaV / (2.0 * rVariables.Peclet);
        rVariables.SigmaVt = rVariables.rho_dot_c * rVariables.lv / (theta * delta_time * NormVel) * fk;

        rVariables.SigmaVBarra = rVariables.SigmaV + rVariables.SigmaVt;
    }

    rVariables.AlphaR = rVariables.Peclet * (0.5 * rVariables.SigmaV * ((rVariables.XiV + 1.0) / (rVariables.XiV - 1.0)) - rVariables.AlphaV)
                        - 1.0 - (1.0 / rVariables.AuxDiffusion) * inner_prod(rVariables.VelInterHat, prod(rVariables.DifMatrixS, rVariables.VelInterHat));

    if (rVariables.absorption < rVariables.HighTolerance)
    {
        rVariables.AlphaR = 0.0;
    }
    else
    {
        if (rVariables.Peclet < rVariables.LowTolerance)
        {
            rVariables.AlphaR = rVariables.OmegaV / (4.0 * sinh(sqrt(rVariables.OmegaV) / 2.0) * sinh(sqrt(rVariables.OmegaV) / 2.0)) +
                                rVariables.OmegaV / 6.0 - 1.0;

            if (rVariables.OmegaV < rVariables.LowTolerance)
            {
                rVariables.AlphaR = 1.0 / 6.0;
            }
        }
        else if (conductivity < rVariables.HighTolerance)
        {
            rVariables.AlphaR = rVariables.absorption * rVariables.lv * rVariables.lv / 6.0;
        }
    }

    if(rVariables.Residual < rVariables.LowTolerance)
    {
        noalias(rVariables.DifMatrixR) =  rVariables.AlphaR * rVariables.AuxDiffusion
                                    * outer_prod(rVariables.VelInterHat , rVariables.VelInterHat);
    }
    else
    {
        //TODO
        noalias(rVariables.DifMatrixR) =  std::abs(rVariables.TransientResidual / rVariables.Residual) * rVariables.AlphaR * rVariables.AuxDiffusion
                                    * outer_prod(rVariables.VelInterHat , rVariables.VelInterHat);
        // noalias(rVariables.DifMatrixR) =  rVariables.AlphaR * rVariables.AuxDiffusion
        //                             * outer_prod(rVariables.VelInterHat , rVariables.VelInterHat);
    }

    //////////////////////////////////////////////////////
    // Calculate Dsc
    //////////////////////////////////////////////////////

    // Identity matrix
    noalias(rVariables.IdentityMatrix) = ZeroMatrix(TDim,TDim);
    for(unsigned int i = 0; i < TDim; i++)
    {
        rVariables.IdentityMatrix(i,i) = 1.0;
    }

    this->CalculateFICBeta(rVariables);

    //KRATOS_ERROR_IF(rVariables.Beta < 1.0) << this->Id() << rVariables.Beta << std::endl;

    BoundedMatrix<double,TDim,TDim> AuxMatrix;
    BoundedMatrix<double,TDim,TDim> AuxMatrix2;
    BoundedMatrix<double,TDim,TDim> AuxMatrix3;
    BoundedMatrix<double,TDim,TDim> AuxMatrix4;

    noalias(AuxMatrix2) = outer_prod(rVariables.VelInterHat , rVariables.VelInterHat);

    noalias (AuxMatrix) = rVariables.IdentityMatrix - AuxMatrix2;

    // Double dot product
    AuxMatrix4 = prod((rVariables.DifMatrixK + rVariables.DifMatrixS), trans(AuxMatrix));
    double DoubleDotScalar = 0.0;

    for (unsigned int i = 0 ; i < TDim ; i++ )
    {
        DoubleDotScalar += AuxMatrix4 (i,i);
    }

    rVariables.DifSC = (0.5 * rVariables.lsc * std::abs (rVariables.TransientResidual) / rVariables.NormGradPhi
                                        - DoubleDotScalar) * (1.0 - rVariables.Beta * rVariables.Beta);

    // TODO Check LowTolerance
    if (rVariables.NormGradPhi < rVariables.LowTolerance)
    {
        rVariables.DifSC = (0.5 * rVariables.lsc * std::abs (rVariables.TransientResidual) / rVariables.LowTolerance
                                        - DoubleDotScalar) * (1.0 - rVariables.Beta * rVariables.Beta);
    }

    if (rVariables.DifSC < 0.0)
    {
        rVariables.DifSC = 0.0;
    }

    // On the first iteration, SC can't be used
    if (rVariables.IterationNumber == 1)
    {
        rVariables.DifSC = 0.0;
    }

    noalias(rVariables.DifMatrixSC) = rVariables.DifSC * rVariables.IdentityMatrix;

    // double producte = inner_prod (rVariables.VelInter , rVariables.GradPhi);

    //Calculate Dv and new lsc term for Dv
    this->CalculateNormalsAngle(rVariables);

    if (std::abs(std::abs(rVariables.Beta) - 1.0) > rVariables.LowTolerance)
    {
        rVariables.CosinusNormals = 1.0;
    }

    // TODO: Calculate properly
    if (std::abs(rVariables.GradPhi[0]) < rVariables.LowTolerance || std::abs(rVariables.GradPhi[1]) < rVariables.LowTolerance)
    {
        rVariables.CosinusGradPhi = 1.0;
    }
    else
    {
        rVariables.CosinusGradPhi = 0.0;
    }

    // On the first iteration, SC can't be used
    if (rVariables.IterationNumber == 1)
    {
        rVariables.CosinusGradPhi = 1.0;
        rVariables.CosinusNormals = 1.0;
    }

    noalias(rVariables.DifMatrixV) = 0.5 * (rVariables.lv + rVariables.lsc * (1.0 - rVariables.CosinusGradPhi) * (1.0 - rVariables.CosinusNormals * rVariables.CosinusNormals))
                                                 * rVariables.AlphaV * outer_prod(rVariables.VelInterHat , rVariables.VelInter);

KRATOS_CATCH("")
}


//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateHVector(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo)
{
KRATOS_TRY

    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    //Properties variables
    double NormVel = rVariables.VelInter[0]*rVariables.VelInter[0];
    for (unsigned int d = 1; d < TDim; d++)
        NormVel += rVariables.VelInter[d]*rVariables.VelInter[d];
    NormVel = std::sqrt(NormVel);

    // Compute HvVector
    for (unsigned int i = 0; i < TDim; i++)
    {
        rVariables.HvVector[i] = rVariables.AlphaV * rVariables.lv * rVariables.VelInterHat[i];
    }

    // Compute HrVector
    BoundedMatrix<double,TDim,TDim> AuxMatrix;
    BoundedMatrix<double,TDim,TDim> AuxMatrix2;
    BoundedMatrix<double,TDim,TDim> AuxMatrix3;
    BoundedMatrix<double,TDim,TDim> AuxMatrix4;

    if (std::abs(rVariables.Residual) < rVariables.LowTolerance)
    {
        for (unsigned int i = 0 ; i < TDim ; i++ )
        {
            rVariables.HrVector [i] = 0.0;
        }
    }
    else
    {
        if(std::abs(rVariables.TransientResidual) < rVariables.LowTolerance)
        {
            AuxMatrix3 = (2.0 / rVariables.Residual ) * (rVariables.DifMatrixS
                            + rVariables.AlphaR * rVariables.AuxDiffusion * outer_prod(rVariables.VelInterHat, rVariables.VelInterHat));
        }
        else
        {
            // TODO
            AuxMatrix3 = (2.0 * (rVariables.TransientResidual / std::abs(rVariables.TransientResidual)) / rVariables.Residual ) * (rVariables.DifMatrixS
                            + rVariables.AlphaR * rVariables.AuxDiffusion * outer_prod(rVariables.VelInterHat, rVariables.VelInterHat));

            // AuxMatrix3 = (2.0 / rVariables.Residual ) * (rVariables.DifMatrixS
            //                 + rVariables.AlphaR * rVariables.AuxDiffusion * outer_prod(rVariables.VelInterHat, rVariables.VelInterHat));
        }

        rVariables.HrVector = prod(AuxMatrix3, rVariables.GradPhi);

    }

    // Compute HscVector
    if (std::abs(rVariables.NormGradPhi) < rVariables.LowTolerance)
    {
        for (unsigned int i = 0 ; i < TDim ; i++ )
        {
            rVariables.HscVector [i] = 0.0;
        }
    }
    else
    {
        noalias(AuxMatrix2) = outer_prod(rVariables.VelInterHat , rVariables.VelInterHat);

        for (unsigned int i = 0 ; i < TDim ; i++ )
        {
            for (unsigned int j = 0 ; j < TDim ; j++ )
            {
                AuxMatrix(i,j) = rVariables.IdentityMatrix(i,j) - AuxMatrix2(i,j);
            }
        }

        // Double dot product
        AuxMatrix4 = prod((rVariables.DifMatrixK + rVariables.DifMatrixS), AuxMatrix);
        double DoubleDotScalar = 0.0;

        for (unsigned int i = 0 ; i < TDim ; i++ )
        {
            DoubleDotScalar += AuxMatrix4 (i,i);
        }


        double AuxScalar = (rVariables.lsc * (rVariables.Residual / std::abs(rVariables.Residual)) - 2.0 * rVariables.NormGradPhi / rVariables.Residual
                            * DoubleDotScalar) * (1.0 - rVariables.Beta * rVariables.Beta);

        if (std::abs(rVariables.Residual) < rVariables.LowTolerance)
        {
            AuxScalar = 0.0;
        }

        // On the first iteration, SC can't be used
        if (rVariables.IterationNumber == 1)
        {
            AuxScalar = 0.0;
        }

        rVariables.HscVector = AuxScalar * rVariables.GradPhi / rVariables.NormGradPhi;
    }

    // Compute HVector
    rVariables.HVector = rVariables.HvVector + rVariables.HrVector + rVariables.HscVector;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim,TNumNodes>::GetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& Geom = this->GetGeometry();
    GeometryData::IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    if ( rVariable == FIC_BETA || rVariable == PECLET || rVariable == CFL_NUMBER)
    {
        if ( rValues.size() != NumGPoints )
            rValues.resize(NumGPoints);

        this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
        if ( rValues.size() != NumGPoints )
            rValues.resize(NumGPoints);

        for ( unsigned int i = 0;  i < NumGPoints; i++ )
        {
            rValues[i] = 0.0;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rOutput,
                                                                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if(rVariable == PECLET)
    {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        GeometryData::IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
        const unsigned int NumGPoints = integration_points.size();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( ThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Vector detJContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
        {
            //Compute GradNT
            noalias(Variables.GradNT) = DN_DXContainer[GPoint];
            noalias(Variables.VelInter) = ZeroVector(TDim);

            ElementUtilities::InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);

            double NormVel = norm_2(Variables.VelInter);

            noalias(Variables.GradPhi) = ZeroVector(TDim);

            Variables.GradPhi = prod(trans(Variables.GradNT), Variables.NodalPhi);

            Variables.NormGradPhi = norm_2(Variables.GradPhi);

            // Unitary velocity vector
            if (NormVel < Variables.LowTolerance)
            {
                for (unsigned int i = 0; i < TDim; i++)
                {
                    Variables.VelInterHat[i] = 0.0;
                }
            }
            else
            {
                for (unsigned int i = 0; i < TDim; i++)
                {
                    Variables.VelInterHat[i] = Variables.VelInter[i]/NormVel;
                }
            }

            this->CalculatePeclet(Variables, Geom, NormVel, rCurrentProcessInfo, Prop);
            rOutput[GPoint] = Variables.Peclet;
        }
    }
    else if(rVariable == FIC_BETA)
    {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        GeometryData::IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
        const unsigned int NumGPoints = integration_points.size();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( ThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Vector detJContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
        {
            //Compute GradNT
            noalias(Variables.GradNT) = DN_DXContainer[GPoint];

            noalias(Variables.VelInter) = ZeroVector(TDim);
            ElementUtilities::InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);

            double NormVel = norm_2(Variables.VelInter);

            noalias(Variables.GradPhi) = ZeroVector(TDim);

            Variables.GradPhi = prod(trans(Variables.GradNT), Variables.NodalPhi);

            Variables.NormGradPhi = norm_2(Variables.GradPhi);

            // Unitary velocity vector
            if (NormVel < Variables.LowTolerance)
            {
                for (unsigned int i = 0; i < TDim; i++)
                {
                    Variables.VelInterHat[i] = 0.0;
                }
            }
            else
            {
                for (unsigned int i = 0; i < TDim; i++)
                {
                    Variables.VelInterHat[i] = Variables.VelInter[i]/NormVel;
                }
            }

            this->CalculateFICBeta(Variables);
            rOutput[GPoint] = Variables.Beta;
        }
    }
    else if(rVariable == CFL_NUMBER)
    {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        GeometryData::IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
        const unsigned int NumGPoints = integration_points.size();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( ThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Vector detJContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,Geom,Prop,rCurrentProcessInfo);
        double delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);


        //Loop over integration points
        for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
        {
            //Compute GradNT
            noalias(Variables.GradNT) = DN_DXContainer[GPoint];

            noalias(Variables.VelInter) = ZeroVector(TDim);
            ElementUtilities::InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);

            double NormVel = norm_2(Variables.VelInter);

            if (NormVel < Variables.HighTolerance)
            {
                for (unsigned int i = 0; i < TDim; i++)
                {
                    Variables.VelInterHat[i] = 0.0;
                }
            }
            else
            {
                for (unsigned int i = 0; i < TDim; i++)
                {
                    Variables.VelInterHat[i] = Variables.VelInter[i]/NormVel;
                }
            }

            double Domain = Geom.DomainSize();

            if (TDim == 2)
            {
                Variables.lv = std::sqrt(2.0*Domain);
            }
            else
            {
                Variables.lv = pow( (6.0*Domain/Globals::Pi) , (1.0/3.0) );
            }

            if (NormVel > Variables.HighTolerance)
            {
                array_1d <double, 3> Velocity3;
                ElementUtilities::FillArray1dOutput (Velocity3, Variables.VelInterHat);

                Variables.lv = ElementSizeCalculator<TDim,TNumNodes>::ProjectedElementSize(Geom,Velocity3);
            }

            this->CalculateBoundaryLv(Variables);

            Variables.Courant = NormVel * delta_time / Variables.lv;

            rOutput[GPoint] = Variables.Courant;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template class TransientConvectionDiffusionFICElement<2,3>;
template class TransientConvectionDiffusionFICElement<2,4>;
template class TransientConvectionDiffusionFICElement<3,4>;
template class TransientConvectionDiffusionFICElement<3,8>;

} // Namespace Kratos