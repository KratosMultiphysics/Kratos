// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "custom_elements/eulerian_conv_diff_epoxy.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        m_degree_of_cure_vector = ZeroVector(GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2));
        m_glass_transition_temperature = ZeroVector(GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2));
        m_heat_of_reaction = ZeroVector(GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2));
        m_pre_strain_vector = ZeroVector(GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2));

        auto& r_geometry = GetGeometry();
        auto integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
        IndexType number_of_integration_points = integration_points.size();
        Matrix Ncontainer = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        double specific_heat = 0;
        double conductivity = 0;
        for (IndexType i = 0; i < number_of_integration_points; ++i)
        {
            m_glass_transition_temperature[0] = ComputeGlassTransitionTemperature(
                m_degree_of_cure_vector[0]);

            double temperature = 0;
            for (IndexType j = 0; j < TNumNodes; j++)
            {
                const double& temp = r_geometry[j].GetSolutionStepValue(TEMPERATURE, 0);
                temperature += Ncontainer(i, j) * temp;
            }

            specific_heat += this->ComputeSpecificHeatCapacity(
                temperature, m_glass_transition_temperature[i],
                m_degree_of_cure_vector[i]);

            conductivity += this->ComputeThermalConductivity(
                temperature,m_degree_of_cure_vector[i]);
        }
        specific_heat /= number_of_integration_points;
        m_specific_heat_capacity = specific_heat;

        conductivity /= number_of_integration_points;
        m_thermal_conductivity = conductivity;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement<TDim,TNumNodes>::CalculateLocalSystem(
        Matrix& rLeftHandSideMatrix,
        Vector& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Resize of the Left and Right Hand side
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!
        rLeftHandSideMatrix = ZeroMatrix(TNumNodes, TNumNodes);

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!
        rRightHandSideVector = ZeroVector(TNumNodes);


        /// SHANELLE
        //Vector heat_flux_vector = ZeroVector(TNumNodes);

        //Element variables
        ElementVariables Variables;
        this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

        // Compute the geometry
        BoundedMatrix<double,TNumNodes, TDim> DN_DX;
        array_1d<double,TNumNodes > N;
        double Volume;
        this-> CalculateGeometry(DN_DX,Volume);

        // Getting the values of shape functions on Integration Points
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        const GeometryType& Geom = this->GetGeometry();
        Ncontainer = Geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

        // Getting the values of Current Process Info and computing the value of h
        this-> GetNodalValues(Variables,rCurrentProcessInfo);

        Variables.specific_heat = m_specific_heat_capacity;
        Variables.conductivity = m_thermal_conductivity;

        double h = this->ComputeH(DN_DX);

        //Computing the divergence
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for(unsigned int k=0; k<TDim; k++)
            {
                Variables.div_v += DN_DX(i,k)*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
        }

        //Some auxilary definitions
        BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
        BoundedMatrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
        bounded_matrix<double,TNumNodes, TDim> tmp;

        /// SHANELLE
        auto& integration_points = Geom.IntegrationPoints(GeometryData::GI_GAUSS_2);

        // Gauss points and Number of nodes coincides in this case.
        for(unsigned int igauss=0; igauss<TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity in the middle of the tiem step
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                 for(unsigned int k=0; k<TDim; k++)
                    vel_gauss[k] += N[i]*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
            const double norm_vel = norm_2(vel_gauss);
            array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

            const double tau = this->CalculateTau(Variables,norm_vel,h);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);

            //terms which multiply the gradient of phi
            noalias(aux2) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);


            /*
            /// SHANELLE: computation of current degree of cure
            double temperature = 0;
            for (IndexType j = 0; j < TNumNodes; j++)
            {
                const double& temp = Geom[j].GetSolutionStepValue(TEMPERATURE, 0);

                temperature += N[j] * temp;
            }

            double degree_of_cure = this->ComputeDegreeOfCure(m_degree_of_cure_vector[igauss], temperature);

            double heat_flux = this->ComputeHeatOfReaction(
                m_degree_of_cure_vector[igauss],
                degree_of_cure,
                rCurrentProcessInfo[DELTA_TIME]);

            KRATOS_WATCH(heat_flux)

            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                heat_flux_vector[i] += N[i] * heat_flux * integration_points[igauss].Weight()* Volume;
            }*/
        }

        //KRATOS_WATCH(Variables.specific_heat)

        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (Variables.dt_inv*Variables.density*Variables.specific_heat + Variables.theta*Variables.beta*Variables.div_v)*aux1;
        noalias(rRightHandSideVector) = (Variables.dt_inv*Variables.density*Variables.specific_heat - (1.0-Variables.theta)*Variables.beta*Variables.div_v)*prod(aux1,Variables.phi_old);

        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(DN_DX, trans(DN_DX))),Variables.phi_old)*static_cast<double>(TNumNodes) ;


        //terms in aux2
        noalias(rLeftHandSideMatrix) += Variables.density*Variables.specific_heat*Variables.theta*aux2;
        noalias(rRightHandSideVector) -= Variables.density*Variables.specific_heat*(1.0-Variables.theta)*prod(aux2,Variables.phi_old);

        // volume source terms (affecting the RHS only)
        noalias(rRightHandSideVector) += prod(aux1, Variables.volumetric_source);

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Variables.phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Eulerian ConvDiff Element")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement<TDim,TNumNodes>::InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rVariables.theta = rCurrentProcessInfo[THETA]; //Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        rVariables.dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        rVariables.dt_inv = 1.0 / delta_t;
		rVariables.lumping_factor = 1.00 / double(TNumNodes);

        rVariables.conductivity = 0.0;
        rVariables.specific_heat = 0.0;
        rVariables.density = 0.0;
        rVariables.beta = 0.0;
        rVariables.div_v = 0.0;


        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement<TDim,TNumNodes>::CalculateGeometry(BoundedMatrix<double,TNumNodes,TDim>& rDN_DX, double& rVolume)
    {

        const GeometryType& Geom = this->GetGeometry();

        // We select GI_GAUSS_1 due to we are computing at the barycenter.
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( GeometryData::GI_GAUSS_1 );
        const unsigned int NumGPoints = integration_points.size();
        rVolume = Geom.Area();
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,GeometryData::GI_GAUSS_1);

        noalias( rDN_DX ) = DN_DXContainer[0];

    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::ComputeH(BoundedMatrix<double,TNumNodes,TDim >& DN_DX)
    {
        double h=0.0;

        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement<TDim,TNumNodes>::GetNodalValues(ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
    {
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

        //////storing locally the flags to avoid repeated check in the nodal loops
        const bool IsDefinedVelocityVariable = my_settings->IsDefinedVelocityVariable();
        const bool IsDefinedMeshVelocityVariable = my_settings->IsDefinedMeshVelocityVariable();
        const bool IsDefinedDensityVariable = my_settings->IsDefinedDensityVariable();
        const bool IsDefinedSpecificHeatVariableVariable = my_settings->IsDefinedSpecificHeatVariable();
        const bool IsDefinedDiffusionVariable = my_settings->IsDefinedDiffusionVariable();
        const bool IsDefinedVolumeSourceVariable = my_settings->IsDefinedVolumeSourceVariable();

        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            rVariables.phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
            //dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];

			rVariables.v[i]=ZeroVector(3);
			rVariables.vold[i]=ZeroVector(3);
            rVariables.volumetric_source[i] = 0.0;
            if (IsDefinedVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
				  rVariables.v[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
				  rVariables.vold[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar,1);
				  //active_convection=true;
			}

			if (IsDefinedMeshVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
				  rVariables.v[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
				  rVariables.vold[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar,1);
				  //active_convection=true;
			}

			if (IsDefinedDensityVariable)
			{
				const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
				rVariables.density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
			}
			else
				rVariables.density += 1.0;

			if (IsDefinedSpecificHeatVariableVariable)
			{
				const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
				rVariables.specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);
			}
			else
				rVariables.specific_heat += 1.0;

			if (IsDefinedDiffusionVariable)
			{
				const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
				rVariables.conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
			}
			//if not, then the conductivity = 0

            if (IsDefinedVolumeSourceVariable)
            {
                const Variable<double>& rVolumeSourceVar = my_settings->GetVolumeSourceVariable();
                rVariables.volumetric_source[i] += GetGeometry()[i].FastGetSolutionStepValue(rVolumeSourceVar);
            }
        }

        //array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
        //const double norm_grad = norm_2(grad_phi_halfstep);

        rVariables.conductivity *= rVariables.lumping_factor;
        rVariables.density *= rVariables.lumping_factor;
        rVariables.specific_heat *= rVariables.lumping_factor;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement<TDim,TNumNodes>::CalculateTau(const ElementVariables& rVariables, double norm_vel, double h)
    {
        // Dynamic part
        double inv_tau = rVariables.dyn_st_beta * rVariables.dt_inv;

        // Convection
        inv_tau += 2.0 * norm_vel / h + rVariables.beta*rVariables.div_v;

        // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
        inv_tau *= rVariables.density * rVariables.specific_heat;

        // Diffusion
        inv_tau += 4.0 * rVariables.conductivity / (h*h);

        // Limiting
        inv_tau = std::max(inv_tau, 1e-2);

        return (rVariables.density*rVariables.specific_heat) / inv_tau;
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::FinalizeSolutionStep(
        const ProcessInfo& rCurrentProcessInfo)
    {
        auto& r_geometry = GetGeometry();
        auto integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);

        IndexType number_of_integration_points = integration_points.size();

        Matrix Ncontainer = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        for (IndexType i = 0; i < number_of_integration_points; ++i)
        {
            double temperature = 0;
            for (IndexType j = 0; j < TNumNodes; j++)
            {
                const double& temp = r_geometry[j].GetSolutionStepValue(TEMPERATURE, 0);

                temperature += Ncontainer(i, j) * temp;
            }

            //KRATOS_WATCH(temperature)

            double degree_of_cure = this->ComputeDegreeOfCure(m_degree_of_cure_vector[i], temperature);
            m_heat_of_reaction[i] = ComputeHeatOfReaction(m_degree_of_cure_vector[i], degree_of_cure, rCurrentProcessInfo[DELTA_TIME]);
            m_degree_of_cure_vector[i] = degree_of_cure;
        }

        for (IndexType j = 0; j < TNumNodes; j++)
        {
            //KRATOS_WATCH(j)
            //KRATOS_WATCH(r_geometry[j].GetSolutionStepValue(TEMPERATURE, 0));
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                //KRATOS_WATCH(integration_points[i].Weight());
            }
        }
       // KRATOS_WATCH(m_heat_of_reaction)
            //KRATOS_WATCH(GetGeometry().DomainSize())

        for (IndexType i = 0; i < number_of_integration_points; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; j++)
            {
                #pragma omp critical
                r_geometry[j].GetSolutionStepValue(TEMPERATURE , 0) += (m_heat_of_reaction[i] * (integration_points[i].Weight()/ number_of_integration_points))/m_specific_heat_capacity;
            }
        }
        for (IndexType j = 0; j < TNumNodes; j++)
        {
            //KRATOS_WATCH(j)
                //KRATOS_WATCH(r_geometry[j].GetSolutionStepValue(TEMPERATURE, 0));
        }
        for (IndexType i = 0; i < number_of_integration_points; ++i)
        {
            double glass_transition_temperature = this->ComputeGlassTransitionTemperature(
                m_degree_of_cure_vector[i]);
            m_glass_transition_temperature[i] = glass_transition_temperature;
        }

        double specific_heat = 0;
        double conductivity = 0;
        for (IndexType i = 0; i < number_of_integration_points; ++i)
        {
            double temperature = 0;
            for (IndexType j = 0; j < TNumNodes; j++)
            {
                const double& temp = r_geometry[j].GetSolutionStepValue(TEMPERATURE, 0);

                temperature += Ncontainer(i, j) * temp;
            }

            double previous_temperature = 0;
            for (IndexType j = 0; j < TNumNodes; j++)
            {
                const double& temp = r_geometry[j].GetSolutionStepValue(TEMPERATURE, 1);

                previous_temperature += Ncontainer(i, j) * temp;
            }

            specific_heat += this->ComputeSpecificHeatCapacity(
                temperature, m_glass_transition_temperature[i],
                m_degree_of_cure_vector[i]);
            m_pre_strain_vector[i] = this->ComputePreStrainFactor(
                m_degree_of_cure_vector[i], temperature, previous_temperature);

            conductivity += this->ComputeThermalConductivity(
                temperature, m_degree_of_cure_vector[i]);
            //KRATOS_WATCH(temperature)
        }
        specific_heat /= number_of_integration_points;
        conductivity /= number_of_integration_points;
        //KRATOS_WATCH(m_glass_transition_temperature)
        //KRATOS_WATCH(m_degree_of_cure_vector)
        //KRATOS_WATCH(specific_heat)
        m_specific_heat_capacity = specific_heat;
        m_thermal_conductivity = conductivity;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        Matrix LeftHandSide;
        this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::GetValueOnIntegrationPoints(
            const Variable<double>& rVariable,
            std::vector<double>& rOutput,
            const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_geometry = GetGeometry();
        auto integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);

        IndexType number_of_integration_points = integration_points.size();

        if (rOutput.size() != number_of_integration_points)
            rOutput.resize(number_of_integration_points);

        Matrix Ncontainer = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        if (rVariable == DEGREE_OF_CURE)
        {
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                rOutput[i] = m_degree_of_cure_vector[i];
            }
        }
        else if (rVariable == GLASS_TRANSITION_TEMPERATURE)
        {
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                rOutput[i] = m_glass_transition_temperature[i];
            }
        }
        else if (rVariable == HEAT_OF_REACTION)
        {
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                rOutput[i] = m_heat_of_reaction[i];
            }
        }
        else if (rVariable == PRE_STRAIN_FACTOR)
        {
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                rOutput[i] = m_pre_strain_vector[i];
            }
        }
        else if (rVariable == SPECIFIC_HEAT_CAPACITY)
        {
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                rOutput[i] = m_specific_heat_capacity;
            }
        }
        else if (rVariable == THERMAL_CONDUCTIVITY)
        {
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                rOutput[i] = m_thermal_conductivity;
            }
        }
        else if (rVariable == CONDUCTIVITY)
        {
            for (IndexType i = 0; i < number_of_integration_points; ++i)
            {
                rOutput[i] = GetProperties()[CONDUCTIVITY];
            }
        }
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::ComputeDegreeOfCure(double DegreeOfCure, double Temperature)
    {

        double A_1 = 1184;
        double A_2 = 68380;
        double E_1 = 60700;
        double E_2 = 48250;
        double m_factor = 0.314;
        double n_factor = 1.664;

        double rate_of_conversion = ((A_1 * exp(-E_1 / (8.3145 * Temperature)))
            + ((A_2 * exp(-E_2 / (8.3145 * Temperature))) * pow(DegreeOfCure, m_factor))) * pow((1 - DegreeOfCure), n_factor);

        return DegreeOfCure + rate_of_conversion;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::ComputeGlassTransitionTemperature(
        double DegreeOfCure)
    {
        double Tg0 = -44.143;
        double Tginf = 67.688;
        double lamda = 0.288;

        double glass_transition_temperature = ((Tg0 + 273.15) + (Tginf - Tg0) * ((lamda * DegreeOfCure) / (1 - (1 - lamda) * DegreeOfCure)))-273.15;

        return glass_transition_temperature;

    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::ComputeSpecificHeatCapacity(
        double Temperature, double GlassTransitionTemperature, double DegreeOfCure)
    {
        double Crub = 2.10;
        double Crub_alpha = 0.16;
        double Crub_T = 0.000419;
        double Cglass = 1.55;
        double Cglass_T = 0.00244;
        double Cw = 0.474;
        double sigma = 107.6;
        double sigma_T = -0.454;

        double specific_heat_capacity = Crub + Crub_alpha * DegreeOfCure + Crub_T * Temperature + (Cglass + Cglass_T * Temperature - Crub - Crub_alpha * DegreeOfCure - Crub_T * Temperature) /
            (1 + exp(Cw * (Temperature - GlassTransitionTemperature - sigma - sigma_T * Temperature)));

        return specific_heat_capacity;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::ComputeThermalConductivity(
        double Temperature, double DegreeOfCure)
    {
        double theta_conductivity = 0.44;
        double beta_conductivity = -12.1;
        double gamma_conductivity = 0.061;
        

        double thermal_conductivity = (1  + theta_conductivity * DegreeOfCure) / (beta_conductivity + gamma_conductivity * Temperature);

        return thermal_conductivity;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::ComputeHeatOfReaction(
        double degree_of_cure_previous,
        double degree_of_cure_current,
        double delta_time)
    {
        double total_heat_of_reaction = 117/1000;

        const double density = GetProperties()[DENSITY];
        const double volume = GetGeometry().DomainSize();
        const double specific_heat = GetProperties()[SPECIFIC_HEAT];

        double heat_flux = ((degree_of_cure_current - degree_of_cure_previous)/delta_time)
            * total_heat_of_reaction * density * volume;

        return heat_flux;


   //Calculating shrinkage and CTE effects
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionEpoxyElement< TDim, TNumNodes >::ComputePreStrainFactor(
        double DegreeOfCureCurrent,
        double temperature_current,
        double temperature_previous)
    {
        const double CTE = 5;

        if (DegreeOfCureCurrent < 0.77) {
            double pre_strain_factor = (CTE * (temperature_current - temperature_previous)) - (8.233 * DegreeOfCureCurrent - 0.4199)/100;
            return pre_strain_factor;
        }
            double pre_strain_factor = (CTE * (temperature_current - temperature_previous)) - (18.75 * DegreeOfCureCurrent - 8.7915)/100;
                return pre_strain_factor;
    }
//----------------------------------------------------------------------------------------


template class EulerianConvectionDiffusionEpoxyElement<2,3>;
template class EulerianConvectionDiffusionEpoxyElement<2,4>;
template class EulerianConvectionDiffusionEpoxyElement<3,4>;
template class EulerianConvectionDiffusionEpoxyElement<3,8>;

}