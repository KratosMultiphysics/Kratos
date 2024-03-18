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
#include "includes/define.h"
#include "custom_elements/conv_diff_IGA_element.h"
// #include "convection_diffusion_application.h"  // Is this essential?
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------

    
    void ConvDiffIGAElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        unsigned int NumNodes = this->GetGeometry().size();
        if (rResult.size() != NumNodes)
            rResult.resize(NumNodes, false);

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    
    void ConvDiffIGAElement::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        unsigned int NumNodes = this->GetGeometry().size();
        if (ElementalDofList.size() != NumNodes)
            ElementalDofList.resize(NumNodes);

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    
    void ConvDiffIGAElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                        VectorType& rRightHandSideVector,
                        const ProcessInfo& rCurrentProcessInfo)
    {
        const GeometryType& Geom = this->GetGeometry();

        // Read the refinements.iga.json
        const Parameters refinements_parameters = ReadParamatersFile("refinements.iga.json");
        int insertions = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
        int basisFunctionsOrder = refinements_parameters["refinements"][0]["parameters"]["increase_degree_u"].GetInt()+1;

        const unsigned int NumNodes = this->GetGeometry().size();
        const unsigned int TDim = 2;
        // Resize of the Left and Right Hand side
        if (rLeftHandSideMatrix.size1() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes, NumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes, false); //false says not to preserve existing storage!!

        //Element variables
        ElementVariables Variables;
        this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

        // Compute the geometry
        Matrix DN_DX(NumNodes,TDim);
        // array_1d<double,NumNodes > N;
        const Matrix& N_gausspoint = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());
        auto N = row(N_gausspoint, 0); // these are the N which correspond to the gauss point "i_point"
        double Volume;
        this-> CalculateGeometry(DN_DX,Volume);

        // Getting the values of Current Process Info and computing the value of h
        this-> GetNodalValues(Variables,rCurrentProcessInfo);

        // KRATOS_WATCH(Variables.theta) // == 1

        // Variables.vold --> vold
        std::vector<std::vector<int>> vold(NumNodes, std::vector<int>(TDim)); 
        std::vector<std::vector<int>> v(NumNodes, std::vector<int>(TDim)); 
        Vector phi = ZeroVector(NumNodes) ;
        Vector phi_old = ZeroVector(NumNodes) ;
        Vector volumetric_source = ZeroVector(NumNodes) ;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            phi[i] = Geom[i].FastGetSolutionStepValue(TEMPERATURE);
            phi_old[i] = Geom[i].FastGetSolutionStepValue(TEMPERATURE,1);
            volumetric_source[i] = Geom[i].FastGetSolutionStepValue(HEAT_FLUX) ;
            for(unsigned int k=0; k<TDim; k++)
            {
                vold[i][k]=0.0;
                v[i][k]=0.0;
            }
        }
        // double h = this->ComputeH(DN_DX);
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( GeometryData::IntegrationMethod::GI_GAUSS_1 );
        double h = sqrt(integration_points[0].Weight());

        //Computing the divergence
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            for(unsigned int k=0; k<TDim; k++)
            {
                Variables.div_v += DN_DX(i,k)*(vold[i][k]*Variables.theta + vold[i][k]*(1.0-Variables.theta));
            }
        }


        // Get all the derivatives we need
        const Matrix DDN_DDe = Geom.ShapeFunctionDerivatives(2, 0, Geom.GetDefaultIntegrationMethod());
        const Matrix DDDN_DDDe = Geom.ShapeFunctionDerivatives(3, 0, Geom.GetDefaultIntegrationMethod());
        const Matrix DDDDN_DDDDe = Geom.ShapeFunctionDerivatives(4, 0, Geom.GetDefaultIntegrationMethod());
        // KRATOS_WATCH(DDN_DDe)
        // KRATOS_WATCH(DDDN_DDDe)


        // laplacian N
        Vector laplacian_of_N = ZeroVector(NumNodes);
        // Vector laplacian_of_laplacian_of_N = ZeroVector(NumNodes);
        for(unsigned int i=0; i<NumNodes; i++) {
            if (basisFunctionsOrder > 1) {
                laplacian_of_N[i] =  DDN_DDe(i,0) + DDN_DDe(i,2) ;
            }
            // if (basisFunctionsOrder > 3) {
            //     laplacian_of_laplacian_of_N[i] =  DDDDN_DDDDe(i,0) + 2 * DDDDN_DDDDe(i,2) + DDDDN_DDDDe(i,4) ;
            // }
            
        }
        //Some auxilary definitions
        Matrix aux1 = ZeroMatrix(NumNodes, NumNodes); //terms multiplying dphi/dt
        Matrix aux2 = ZeroMatrix(NumNodes, NumNodes); //terms multiplying phi
        Matrix aux3 = ZeroMatrix(NumNodes, NumNodes); //terms multiplying laplacian of phi
        Matrix aux4 = ZeroMatrix(NumNodes, NumNodes); //terms multiplying laplacian of f
        Matrix aux5 = ZeroMatrix(NumNodes, NumNodes); //terms multiplying laplacian of laplacian of N

        Matrix tmp = ZeroMatrix(NumNodes, TDim);
        Matrix grad_of_v_dot_grad = ZeroMatrix(NumNodes, TDim);
        Matrix grad_of_laplacian = ZeroMatrix(NumNodes, TDim);
        Matrix DN_DX_tau = ZeroMatrix(NumNodes, TDim);


            //obtain the velocity in the middle of the time step
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            for (unsigned int igauss = 0; igauss < NumNodes; igauss++)
            {
                 for(unsigned int k=0; k<TDim; k++)
                    vel_gauss[k] += N[igauss]*(vold[igauss][k]*Variables.theta + vold[igauss][k]*(1.0-Variables.theta));
            }
            const double norm_vel = norm_2(vel_gauss);

            Vector a_dot_grad = ZeroVector(NumNodes);
            a_dot_grad = prod(DN_DX, vel_gauss);

            const double tau = this->CalculateTau(Variables,norm_vel,h);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);

            //terms which multiply the gradient of phi
            noalias(aux2) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);

            //terms which multiply the gradient of phi
            // noalias(aux3) += - tau*outer_prod(a_dot_grad, laplacian_of_N);

            //terms which multiply the laplacian of f
            noalias(aux4) += - tau*outer_prod(N, laplacian_of_N);

            //terms which multiply the laplacian_of_laplacian_of_phi
            noalias(aux5) += - tau*outer_prod(laplacian_of_N, laplacian_of_N);

            // // grad(\laplacian N)        
            // if (basisFunctionsOrder>2) {
            //     grad_of_laplacian(i,0) = -tau * (DDDN_DDDe(i,0)+DDDN_DDDe(i,2)) ;
            //     grad_of_laplacian(i,1) = -tau * (DDDN_DDDe(i,1)+DDDN_DDDe(i,3)) ;
            // }
            // grad(c \cot grad N)
            // if (basisFunctionsOrder>1) {
            //     grad_of_v_dot_grad(i,0) = +tau * (vold[i][0]*DDN_DDe(i,0) + vold[i][1]*DDN_DDe(i,1)) ;
            //     grad_of_v_dot_grad(i,1) = +tau * (vold[i][0]*DDN_DDe(i,1) + vold[i][1]*DDN_DDe(i,2)) ;
            // }
            // DN_DX_tau(i,0) = +tau * DN_DX(i,0);
            // DN_DX_tau(i,1) = +tau * DN_DX(i,1); 


        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  += (Variables.dt_inv*Variables.density*Variables.specific_heat + Variables.theta*Variables.beta*Variables.div_v)*aux1;
        noalias(rRightHandSideVector) += (Variables.dt_inv*Variables.density*Variables.specific_heat - (1.0-Variables.theta)*Variables.beta*Variables.div_v)*prod(aux1,phi_old);

        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(DN_DX, trans(DN_DX))); // *static_cast<double>(NumNodes);
        noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(DN_DX, trans(DN_DX))),phi_old); //*static_cast<double>(NumNodes) ;

        //terms in aux2
        noalias(rLeftHandSideMatrix) += Variables.density*Variables.specific_heat*Variables.theta*aux2;
        noalias(rRightHandSideVector) -= Variables.density*Variables.specific_heat*(1.0-Variables.theta)*prod(aux2,phi_old);

        //adding in aux3   --->     STABILIZATION laplacian term: tau*(v \cdot grad w, laplacian phi)
        // noalias(rLeftHandSideMatrix)  += Variables.conductivity *Variables.theta* aux3;
        // noalias(rRightHandSideVector) -= Variables.conductivity *(1.0-Variables.theta) * prod(aux3,phi_old);

        //adding the diffusion from STABILIZATION grad_of_v_dot_grad: \grad w . \grad(v . \grad phi)
        // noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(grad_of_v_dot_grad, trans(DN_DX))); 
        // noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(grad_of_v_dot_grad, trans(DN_DX))),phi_old); 

        // //adding the diffusion from STABILIZATION grad_of_laplacian: \grad w . \grad(\laplacian phi)
        // noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(DN_DX, trans(grad_of_laplacian))); 
        // noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(DN_DX, trans(grad_of_laplacian))),phi_old); 
        // Alternative WITHOUT integration by part the stabilization term
        noalias(rLeftHandSideMatrix) += Variables.conductivity*Variables.theta*aux5;
        noalias(rRightHandSideVector) -= Variables.conductivity*(1.0-Variables.theta)*prod(aux5,phi_old);

        
        //adding the diffusion from STABILIZATION d/dt grad phi
        // noalias(rLeftHandSideMatrix)  += (Variables.dt_inv*Variables.density*Variables.specific_heat )* prod(DN_DX_tau, trans(DN_DX));
        // noalias(rRightHandSideVector) -= prod(Variables.dt_inv*Variables.density*Variables.specific_heat * prod(DN_DX_tau, trans(DN_DX)),phi_old);

        // volume source terms (affecting the RHS only)
        noalias(rRightHandSideVector) += prod(aux1, volumetric_source);
        noalias(rRightHandSideVector) += prod(aux4, volumetric_source);

        // DA TOGLIERE:
            // const double qgauss = inner_prod(N, volumetric_source);
            // noalias(rRightHandSideVector) +=  qgauss * N;

        rRightHandSideVector *= Volume;
        rLeftHandSideMatrix *= Volume;

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);
        
        
        // KRATOS_WATCH(rRightHandSideVector)
        // KRATOS_WATCH(rLeftHandSideMatrix)
        // KRATOS_WATCH(phi)
        // exit(0);
    }

    // void ConvDiffIGAElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
    //                                             VectorType& rRightHandSideVector,
    //                                             const ProcessInfo& rCurrentProcessInfo)
    // {
    //     KRATOS_TRY

    //     const ProcessInfo& r_process_info = rCurrentProcessInfo;
    //     ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    //     auto& r_settings = *p_settings;

    //     const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable(); // Temperature
    //     const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // Conductivity
    //     const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable(); // HeatFlux

    //     const auto& r_geometry = this->GetGeometry();
    //     const unsigned int number_of_points = r_geometry.size();
    //     const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 3
        
    //     //resizing as needed the LHS
    //     if(rLeftHandSideMatrix.size1() != number_of_points)
    //         rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    //     noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS
        
    //     // resizing as needed the RHS
    //     if(rRightHandSideVector.size() != number_of_points)
    //         rRightHandSideVector.resize(number_of_points,false);
    //     noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    //     // reading integration points and local gradients
    //     const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    //     const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    //     const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    //     // Initialize DN_DX
    //     Matrix DN_DX(number_of_points,2);
    //     Matrix InvJ0(2,2);
    //     Vector temp(number_of_points);

    //     // Initialize Jacobian
    //     GeometryType::JacobiansType J0;
    //     r_geometry.Jacobian(J0,this->GetIntegrationMethod());


    //     Vector heat_flux_local(number_of_points);
    //     Vector nodal_conductivity(number_of_points);
    //     for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    //     {
    //         heat_flux_local[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_volume_source_var);
    //         nodal_conductivity[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_var);
    //     }

    //     double DetJ0;
    //     Matrix Jacobian = ZeroMatrix(2,2);
    //     Jacobian(0,0) = J0[0](0,0);
    //     Jacobian(0,1) = J0[0](0,1);
    //     Jacobian(1,0) = J0[0](1,0);
    //     Jacobian(1,1) = J0[0](1,1);
    //     // Calculating inverse jacobian and jacobian determinant
    //     MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    //     // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    //     noalias(DN_DX) = prod(DN_De[0],InvJ0);

    //     auto N = row(N_gausspoint,0); // these are the N which correspond to the gauss point "i_point"
    //     const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0);

    //     // Also the conductivity is multiplied by the value of the shape function at the GaussPoint
    //     const double conductivity_gauss = inner_prod(N, nodal_conductivity);
    //     noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); //

    //     // KRATOS_WATCH(IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)))
    //     // KRATOS_WATCH(IntToReferenceWeight)
    //     // KRATOS_WATCH(conductivity_gauss)
    //     // KRATOS_WATCH(prod(DN_DX, trans(DN_DX)))

    //     // Calculating the local RHS
    //     const double qgauss = inner_prod(N, heat_flux_local);
    //     noalias(rRightHandSideVector) += IntToReferenceWeight * qgauss * N;

    //     // RHS = ExtForces - K*temp;
    //     for (unsigned int i = 0; i < number_of_points; i++) {
    //         temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);

    //         // std::ofstream outputFile("txt_files/Id_active_control_points.txt", std::ios::app);
    //         // outputFile << r_geometry[i].GetId() << "  " <<r_geometry[i].GetDof(r_unknown_var).EquationId() <<"\n";
    //         // outputFile.close();
    //     }
    //     // KRATOS_WATCH(rRightHandSideVector)
    //     // KRATOS_WATCH(heat_flux_local)
    //     // RHS -= K*temp
    //     noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    //     // KRATOS_WATCH(rRightHandSideVector)
    //     // KRATOS_WATCH(rLeftHandSideMatrix)
    //     // KRATOS_WATCH(temp)
    //     // exit(0);

    //     KRATOS_CATCH("")
    // }

//----------------------------------------------------------------------------------------

    
    void ConvDiffIGAElement::InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        const unsigned int NumNodes = this->GetGeometry().size();

        rVariables.theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA]; //Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        rVariables.dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        rVariables.dt_inv = 1.0 / delta_t;
		rVariables.lumping_factor = 1.00 / double(NumNodes);

        rVariables.conductivity = 0.0;
        rVariables.specific_heat = 0.0;
        rVariables.density = 0.0;
        rVariables.beta = 0.0;
        rVariables.div_v = 0.0;


        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------

    void ConvDiffIGAElement::CalculateGeometry(Matrix& rDN_DX, double& rVolume)
    {

        const GeometryType& Geom = this->GetGeometry();

        const GeometryType::ShapeFunctionsGradientsType& DN_De = Geom.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        GeometryType::JacobiansType J0;
        Geom.Jacobian(J0,this->GetIntegrationMethod());
        Matrix InvJ0(2,2);
        double DetJ0;
        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);
        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // We select GI_GAUSS_1 due to we are computing at the barycenter.
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( GeometryData::IntegrationMethod::GI_GAUSS_1 );
        const unsigned int NumGPoints = integration_points.size();
        // rVolume = Geom.Area();
        rVolume = integration_points[0].Weight() ;
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);

        // noalias(DN_DXContainer) = prod(DN_De[0],InvJ0);
        // Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,GeometryData::IntegrationMethod::GI_GAUSS_1);

        noalias( rDN_DX ) = prod(DN_De[0],InvJ0);
    }

//----------------------------------------------------------------------------------------

    
    // double ConvDiffIGAElement::ComputeH(Matrix& DN_DX)
    // {
    //     double h=0.0;
    //     const unsigned int NumNodes = this->GetGeometry().size();
    //     const unsigned int TDim = 2;
    //     for(unsigned int i=0; i<NumNodes; i++)
    //     {
    //         double h_inv = 0.0;
    //         for(unsigned int k=0; k<TDim; k++)
    //         {
    //             h_inv += DN_DX(i,k)*DN_DX(i,k);
    //         }
    //         h += 1.0/h_inv;
    //     }
    //     h = sqrt(h)/static_cast<double>(NumNodes);
    //     return h;
    // }

//----------------------------------------------------------------------------------------

    
    void ConvDiffIGAElement::GetNodalValues(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo) const
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

        const unsigned int NumNodes = this->GetGeometry().size();

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            // rVariables.phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            // rVariables.phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
            
            ////dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];
			// rVariables.v[i]=ZeroVector(3);
			// rVariables.vold[i]=ZeroVector(3);
            // rVariables.volumetric_source[i] = 0.0;
            // if (IsDefinedVelocityVariable)
            // {
            //     //   const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
            //     //   KRATOS_WATCH(rVelocityVar)
			// 	   //   rVariables.v[i] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
			// 	   //   rVariables.vold[i] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
			// 	   //    //active_convection=true;
			// }

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

            // if (IsDefinedVolumeSourceVariable)
            // {
            //     const Variable<double>& rVolumeSourceVar = my_settings->GetVolumeSourceVariable();
            //     rVariables.volumetric_source[i] += GetGeometry()[i].FastGetSolutionStepValue(rVolumeSourceVar);
            //     // KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(rVolumeSourceVar))
            //     // KRATOS_WATCH(rVariables.volumetric_source)
                
            //     // rVariables.volumetric_source[i] += (1-rVariables.theta) * GetGeometry()[i].FastGetSolutionStepValue(HEAT_FLUX,1) + rVariables.theta * GetGeometry()[i].FastGetSolutionStepValue(HEAT_FLUX);

            // }
        }


        //array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
        //const double norm_grad = norm_2(grad_phi_halfstep);

        rVariables.conductivity *= rVariables.lumping_factor;
        rVariables.density *= rVariables.lumping_factor;
        rVariables.specific_heat *= rVariables.lumping_factor;
    }

    
    double ConvDiffIGAElement::CalculateTau(const ElementVariables& rVariables, double norm_vel, double h)
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

    
    void ConvDiffIGAElement::CalculateRightHandSide(
        VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType LeftHandSide;
        KRATOS_WATCH('LefthandSIDE \n \n')
        this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);
    }

//----------------------------------------------------------------------------------------


    void ConvDiffIGAElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_geometry = GetGeometry();
        const SizeType nb_nodes = r_geometry.size();

        // Integration Points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        // Shape function values
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        const ProcessInfo& r_process_info = rCurrentProcessInfo;
        ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
        auto& r_settings = *p_settings;
        const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

        GeometryType::JacobiansType J0;
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());

        double rOutput = 0;
        for (IndexType i = 0; i < nb_nodes; ++i)
        {
            double output_solution_step_value = r_geometry[i].GetSolutionStepValue(r_unknown_var);
            rOutput += r_N(0, i) * output_solution_step_value;
        }        

        std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
        if (output_file.is_open()) {
            output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file << rOutput << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << integration_points[0].Weight() << std::endl;
            output_file.close();
        }
    }

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters ConvDiffIGAElement::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        std::ifstream infile(rDataFileName);

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    };


}
