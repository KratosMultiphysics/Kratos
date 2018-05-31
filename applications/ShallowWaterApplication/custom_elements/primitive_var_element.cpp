//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_elements/primitive_var_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    int PrimitiveVarElement<TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // Base class checks for positive Jacobian and Id > 0
        int ierr = Element::Check(rCurrentProcessInfo);
        if(ierr != 0) return ierr;

        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(HEIGHT)
        KRATOS_CHECK_VARIABLE_KEY(PROJECTED_SCALAR1)
        KRATOS_CHECK_VARIABLE_KEY(PROJECTED_VECTOR1)
        KRATOS_CHECK_VARIABLE_KEY(BATHYMETRY)
        KRATOS_CHECK_VARIABLE_KEY(RAIN)
        KRATOS_CHECK_VARIABLE_KEY(MANNING)
        KRATOS_CHECK_VARIABLE_KEY(GRAVITY)
        KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME)
        KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_TAU)
        KRATOS_CHECK_VARIABLE_KEY(WATER_HEIGHT_UNIT_CONVERTER)

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            Node<3> &rnode = this->GetGeometry()[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BATHYMETRY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN,rnode)

            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rnode)
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rnode)
            KRATOS_CHECK_DOF_IN_NODE(HEIGHT,rnode)
        }

        return ierr;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        unsigned int element_size = TNumNodes*3;
        if(rResult.size() != element_size)
            rResult.resize(element_size,false);                         // False says not to preserve existing storage!!
        
        GeometryType& rGeom = GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[counter++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
            rResult[counter++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
            rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
        }
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        const unsigned int element_size = TNumNodes*3;
        if(rElementalDofList.size() != element_size)
            rElementalDofList.resize(element_size);
        
        GeometryType& rGeom = GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_X);
            rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_Y);
            rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
        }
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Resize of the Left and Right Hand side
        unsigned int element_size = TNumNodes*3;
        if(rLeftHandSideMatrix.size1() != element_size)
            rLeftHandSideMatrix.resize(element_size,element_size,false); // False says not to preserve existing storage!!

        if(rRightHandSideVector.size() != element_size)
            rRightHandSideVector.resize(element_size,false);             // False says not to preserve existing storage!!

        // Initialize element variables
        ElementVariables variables;
        this-> InitializeElement(variables, rCurrentProcessInfo);

        // Compute the geometry
        BoundedMatrix<double,TNumNodes, 2> DN_DX;
        array_1d<double,TNumNodes> N;
        double Area;
        this-> CalculateGeometry(DN_DX,Area);
        double elem_length = this->ComputeElemSize(DN_DX);

        // Getting the values of shape functions on Integration Points
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;  // In this case, number of Gauss points and number of nodes coincides
        const GeometryType& rGeom = this->GetGeometry();
        Ncontainer = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

        // Get nodal values for current step and projected variables (this function inlcudes the units conversion)
        this-> GetNodalValues(variables);

        // Get element values (this function inlcudes the units conversion)
        this-> GetElementValues(DN_DX, variables );
        double abs_vel = norm_2(variables.velocity );
        double height43 = std::pow(variables.height, 1.33333 );

        // Compute stabilization and discontinuity capturing parameters
        double tau_u;
        double tau_h;
        double k_dc;
        this-> ComputeStabilizationParameters(variables, elem_length, tau_u, tau_h, k_dc);

        // Some auxilary definitions
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_q= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_w= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_u  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_u_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);

        this-> ComputeAuxMatrices(Ncontainer, DN_DX, variables, mass_matrix_q, mass_matrix_w, aux_w_grad_h, aux_q_div_u, aux_h_diffus, aux_u_diffus);

        noalias(mass_matrix) = mass_matrix_w + mass_matrix_q;

        // Build LHS
        // Cross terms
        noalias(rLeftHandSideMatrix)  = variables.height * aux_q_div_u;      // Add <q*h*div(u)> to Mass Eq.
        noalias(rLeftHandSideMatrix) += variables.gravity * aux_w_grad_h;    // Add <w*g*grad(h)> to Momentum Eq.

        // Inertia terms
        noalias(rLeftHandSideMatrix) += variables.dt_inv * mass_matrix;      // Add <N,N> to both Eq's

        // Stabilization terms
        noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;    // Add art. diff. to Mass Eq.
        noalias(rLeftHandSideMatrix) +=         tau_u  * aux_u_diffus;    // Add art. diff. to Momentum Eq.

        // Friction term
        noalias(rLeftHandSideMatrix) += variables.gravity * variables.manning2 * abs_vel / height43 * mass_matrix_w;

        // Build RHS
        // Source term (bathymetry contribution)
        noalias(rRightHandSideVector)  = -variables.gravity * prod(aux_w_grad_h, variables.depth);
        
        // Source term (rain contribution)
        noalias(rRightHandSideVector) += prod(mass_matrix, variables.rain);

        // Inertia terms
        noalias(rRightHandSideVector) += variables.dt_inv * prod(mass_matrix, variables.proj_unk);

        // Substracting the botton diffusion due to stabilization (eta = h + H)
        noalias(rRightHandSideVector) -= (k_dc + tau_h) * prod(aux_h_diffus, variables.depth);

        // Substracting the Dirichlet term (since we use a residualbased approach)
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, variables.unknown);

        rRightHandSideVector *= Area * variables.lumping_factor;
        rLeftHandSideMatrix  *= Area * variables.lumping_factor;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM || rVariable == MIU)
        {
            for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
                rValues[PointNumber] = double(this->GetValue(rVariable));
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::InitializeElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        rVariables.dt_inv = 1.0 / delta_t;
        rVariables.lumping_factor = 1.0 / static_cast<double>(TNumNodes);
        rVariables.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        rVariables.gravity = rCurrentProcessInfo[GRAVITY_Z];
        rVariables.manning2 = std::pow( GetProperties()[MANNING], 2);
        rVariables.height_units = rCurrentProcessInfo[WATER_HEIGHT_UNIT_CONVERTER];
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::CalculateGeometry(BoundedMatrix<double, TNumNodes, 2>& rDN_DX, double& rArea)
    {
        const GeometryType& rGeom = this->GetGeometry();

        // We select GI_GAUSS_1 due to we are computing at the barycenter.
        const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);
        const unsigned int NumGPoints = integration_points.size();
        rArea = rGeom.Area();
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer( NumGPoints );
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, GeometryData::GI_GAUSS_1);

        noalias( rDN_DX ) = DN_DXContainer[0];

    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    double PrimitiveVarElement<TNumNodes>::ComputeElemSize(const BoundedMatrix<double, TNumNodes, 2>& rDN_DX)
    {
        double l = 0.0;

        for(unsigned int i = 0; i < TNumNodes; i++)
        {
            double l_inv = 0.0;
            for(unsigned int k = 0; k < 2; k++)
            {
                l_inv += rDN_DX(i,k) * rDN_DX(i,k);
            }
            l += 1.0 / l_inv;
        }
        l = std::sqrt(l) / static_cast<double>(TNumNodes);
        return l;
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::GetNodalValues(ElementVariables& rVariables)
    {
        GeometryType& rGeom = GetGeometry();
        unsigned int counter = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.depth[counter] = 0;
            rVariables.rain[counter]  = 0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);
            rVariables.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_X);
            counter++;

            rVariables.depth[counter] = 0;
            rVariables.rain[counter]  = 0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
            rVariables.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Y);
            counter++;

            rVariables.depth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY) / rVariables.height_units;
            rVariables.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
            rVariables.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_SCALAR1);
            counter++;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables)
    {
        // Initialize outputs
        rVariables.height = 0;
        rVariables.velocity = ZeroVector(2);
        rVariables.height_grad = ZeroVector(2);

        // integrate over the element
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.velocity[0] += rVariables.unknown[  + 3*i];
            rVariables.velocity[1] += rVariables.unknown[1 + 3*i];
            rVariables.height += rVariables.unknown[2 + 3*i];
            rVariables.height_grad[0] += rDN_DX(i,0) * rVariables.unknown[2 + 3*i];
            rVariables.height_grad[1] += rDN_DX(i,1) * rVariables.unknown[2 + 3*i];
        }

        rVariables.velocity *= rVariables.lumping_factor;
        rVariables.height *= rVariables.lumping_factor * rVariables.height_units;
        rVariables.height_grad *= rVariables.height_units;
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::ComputeStabilizationParameters(const ElementVariables& rVariables,
                                                                        const double& rElemSize,
                                                                        double& rTauU,
                                                                        double& rTauH,
                                                                        double& rKdc)
    {
        // Initialize outputs
        rTauU = 0;
        rTauH = 0;
        rKdc  = 0;

        // Get element values
        double height_grad_norm = norm_2(rVariables.height_grad);

        // Compute stabilization parameters
        bool stabilization = true;
        double Ctau = rVariables.dyn_tau;  // 0.005 ~ 0.002  // 0.002; //
        double abs_height = std::abs(rVariables.height);
        if (stabilization && abs_height > 1e-6)
        {
            rTauU = Ctau/rElemSize*std::sqrt(rVariables.gravity/abs_height);
            rTauH = Ctau/rElemSize*std::sqrt(abs_height/rVariables.gravity);
        }

        // Compute discontinuity capturing parameters
        bool discontinuity_capturing = true;
        double gradient_threshold = 1e-6;
        //~ double residual;
        if (discontinuity_capturing && height_grad_norm > gradient_threshold)
        {
            rKdc = 0.5*0.1*rElemSize*height_grad_norm;  // Residual formulation
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::ComputeAuxMatrices(
            const BoundedMatrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff )
    {
        // Initialize solution
        noalias(rMassMatrixVector) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rMassMatrixScalar) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rScalarGrad) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rVectorDiv) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rScalarDiff) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rVectorDiff) = ZeroMatrix(TNumNodes*3, TNumNodes*3);

        // Some auxilary definitions
        array_1d<double,TNumNodes> N;
        BoundedMatrix<double,2,TNumNodes*3> N_vel        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
        BoundedMatrix<double,1,TNumNodes*3> N_height     = ZeroMatrix(1,TNumNodes*3);  // Shape functions vector (for height unknown)
        BoundedMatrix<double,1,TNumNodes*3> DN_DX_vel    = ZeroMatrix(1,TNumNodes*3);  // Shape functions gradients vector (for velocity unknown)
        BoundedMatrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)

        // Loop on Gauss points. In this case, number of Gauss points and number of nodes coincides
        for(unsigned int igauss = 0; igauss < TNumNodes; igauss++)
        {
            noalias(N) = row(rNcontainer, igauss);
            
            // Build shape and derivatives functions at Gauss points
            for(unsigned int nnode = 0; nnode < TNumNodes; nnode++)
            {
                // Height gradient
                DN_DX_height(0, 2+nnode*3) = rDN_DX(nnode,0);
                DN_DX_height(1, 2+nnode*3) = rDN_DX(nnode,1);
                // Velocity divergence
                DN_DX_vel(0,   nnode*3) = rDN_DX(nnode,0);
                DN_DX_vel(0, 1+nnode*3) = rDN_DX(nnode,1);
                // Height shape funtions
                N_height(0, 2+nnode*3) = N[nnode];
                // Velocity shape functions
                N_vel(0,   nnode*3) = N[nnode];
                N_vel(1, 1+nnode*3) = N[nnode];
            }
            N_height     *= rVariables.height_units;
            DN_DX_height *= rVariables.height_units;

            noalias(rMassMatrixVector) += prod(trans(N_vel),N_vel);       // q * h
            noalias(rMassMatrixScalar) += prod(trans(N_height),N_height); // w * u

            noalias(rVectorDiv)  += prod(trans(N_height),DN_DX_vel);       // q * div_u
            noalias(rScalarGrad) += prod(trans(N_vel),DN_DX_height);       // w * grad_h

            noalias(rVectorDiff) += prod(trans(DN_DX_vel),DN_DX_vel);       // div_w * div_u
            noalias(rScalarDiff) += prod(trans(DN_DX_height),DN_DX_height); // grad_q * grad_h
        }

    }

//----------------------------------------------------------------------

template class PrimitiveVarElement<3>;
template class PrimitiveVarElement<4>;

} // namespace Kratos
