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

// Project includes
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_elements/euler_cons_var_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    int EulerConsVarElement<TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // Base class checks for positive Jacobian and Id > 0
        int ierr = Element::Check(rCurrentProcessInfo);
        if(ierr != 0) return ierr;

        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(MOMENTUM)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(HEIGHT)
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
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BATHYMETRY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN,rnode)

            KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X,rnode)
            KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y,rnode)
            KRATOS_CHECK_DOF_IN_NODE(HEIGHT,rnode)
        }

        return ierr;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        unsigned int element_size = TNumNodes*3;
        if(rResult.size() != element_size)
            rResult.resize(element_size,false);                         // False says not to preserve existing storage!!

        GeometryType& rGeom = this->GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[counter++] = rGeom[i].GetDof(MOMENTUM_X).EquationId();
            rResult[counter++] = rGeom[i].GetDof(MOMENTUM_Y).EquationId();
            rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
        }
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        unsigned int element_size = TNumNodes*3;
        if(rElementalDofList.size() != element_size)
            rElementalDofList.resize(element_size);

        GeometryType& rGeom = this->GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_X);
            rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_Y);
            rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
        }
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        // Resize of the Left and Right Hand side
        unsigned int element_size = TNumNodes*3;
        if(rLeftHandSideMatrix.size1() != element_size)
            rLeftHandSideMatrix.resize(element_size,element_size,false); // Resizing the system in case it does not have the right size
        
        if(rRightHandSideVector.size() != element_size)
            rRightHandSideVector.resize(element_size,false);

        // Initialize element variables
        ElementVariables variables;
        this-> InitializeElement(variables, rCurrentProcessInfo);

        // Compute the geometry
        BoundedMatrix<double,TNumNodes, 2> DN_DX;
        array_1d<double,TNumNodes> N;
        double Area;
        this-> CalculateGeometry(DN_DX, Area);
        double elem_length = this->ComputeElemSize(DN_DX);

        // Getting the values of shape functions on Integration Points
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;  // In this case, number of Gauss points and number of nodes coincides
        const GeometryType& rGeom = this->GetGeometry();
        Ncontainer = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

        // Get nodal values for current step and projected variables (this function inlcudes the units conversion)
        this-> GetNodalValues(variables);

        // Get element values (this function inlcudes the units conversion)
        this-> GetElementValues(DN_DX, variables );
        double abs_mom = norm_2(variables.momentum );
        double height73 = std::pow(variables.height, 2.33333 );

        // Compute stabilization and discontinuity capturing parameters
        double tau_m;
        double tau_h;
        double k_dc;
        this-> ComputeStabilizationParameters(variables, elem_length, tau_m, tau_h, k_dc);

        // Some auxilary definitions
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_q= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_w= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_convect  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_non_lin  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_m  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_m_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);

        // ComputeAuxMatrices: loops the Gauss points and build the matrices
        this-> ComputeAuxMatrices(Ncontainer, DN_DX, variables, mass_matrix_q, mass_matrix_w, aux_w_grad_h, aux_q_div_m, aux_h_diffus, aux_m_diffus, aux_convect, aux_non_lin);

        noalias(mass_matrix) = mass_matrix_w + mass_matrix_q;

        // Build LHS
        // Inertia terms
        noalias(rLeftHandSideMatrix)  = variables.dt_inv * mass_matrix;     // Add <N,N> to both Eq's

        // Convective and non linear terms
        noalias(rLeftHandSideMatrix) += aux_convect;                    // Add <w,u*grad(hu)> to Momentum Eq.
        noalias(rLeftHandSideMatrix) += aux_non_lin;                    // Add <w,grad(u)*hu> to Momentum Eq.

        // Cross terms
        noalias(rLeftHandSideMatrix) += aux_q_div_m;                     // Add <q,div(hu)> to Mass Eq.
        noalias(rLeftHandSideMatrix) += variables.gravity * variables.height * aux_w_grad_h; // Add <w,g*h*grad(h)> to Momentum Eq.

        // Stabilization terms
        noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;  // Add art. diff. to Mass Eq.
        noalias(rLeftHandSideMatrix) +=         tau_m  * aux_m_diffus;  // Add art. diff. to Momentum Eq.

        // Friction term
        noalias(rLeftHandSideMatrix) += variables.gravity * variables.manning2 * abs_mom / height73 * mass_matrix_w;

        // Build RHS
        // Source terms (bathymetry contribution)
        noalias(rRightHandSideVector)  = -variables.gravity * variables.height * prod(aux_w_grad_h, variables.depth); // Add <w,-g*h*grad(H)> to RHS (Momentum Eq.)

        // Source terms (rain contribution)
        noalias(rRightHandSideVector) += prod(mass_matrix, variables.rain);

        // Inertia terms
        noalias(rRightHandSideVector) += variables.dt_inv * prod(mass_matrix, variables.proj_unk);

        // Subtracting the dirichlet term (since we use a residualbased approach)
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, variables.unknown);

        rRightHandSideVector *= Area * variables.lumping_factor;
        rLeftHandSideMatrix *= Area * variables.lumping_factor;

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::GetNodalValues(ElementVariables& rVariables)
    {
        GeometryType& rGeom = this->GetGeometry();
        unsigned int counter = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.depth[counter] = 0;
            rVariables.rain[counter] =  0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X,1);
            counter++;

            rVariables.depth[counter] = 0;
            rVariables.rain[counter] =  0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y,1);
            counter++;

            rVariables.depth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY);
            rVariables.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(HEIGHT,1);
            counter++;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables)
    {
        // Initialize outputs
        rVariables.height = 0;
        rVariables.momentum = ZeroVector(2);
        rVariables.velocity = ZeroVector(2);
        rVariables.velocity_grad = ZeroMatrix(2,2);

        // check if the element is close to dry
        bool near_dry = false;
        double threshold = 1e-1;       
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if (rVariables.unknown[2 + 3*i] < threshold)
                near_dry = true;
        }

        // integrate over the element
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.momentum[0] += rVariables.unknown[  + 3*i];
            rVariables.momentum[1] += rVariables.unknown[1 + 3*i];
            rVariables.height += rVariables.unknown[2 + 3*i];
            if (near_dry)
            {
                rVariables.velocity[0] += rVariables.unknown[  + 3*i];
                rVariables.velocity[1] += rVariables.unknown[1 + 3*i];
                rVariables.velocity_grad(0,0) += rDN_DX(i,0) * rVariables.unknown[  + 3*i];
                rVariables.velocity_grad(0,1) += rDN_DX(i,0) * rVariables.unknown[1 + 3*i];
                rVariables.velocity_grad(1,0) += rDN_DX(i,1) * rVariables.unknown[  + 3*i];
                rVariables.velocity_grad(1,1) += rDN_DX(i,1) * rVariables.unknown[1 + 3*i];
            }
            else
            {
                rVariables.velocity[0] += rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity[1] += rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_grad(0,0) += rDN_DX(i,0) * rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_grad(0,1) += rDN_DX(i,0) * rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_grad(1,0) += rDN_DX(i,1) * rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_grad(1,1) += rDN_DX(i,1) * rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
            }
        }

        rVariables.momentum *= rVariables.lumping_factor;
        rVariables.height *= rVariables.lumping_factor * rVariables.height_units;
        if (near_dry)
        {
            rVariables.velocity = rVariables.momentum / rVariables.height;
            rVariables.velocity_grad /= rVariables.height;
        }
        else
        {
            rVariables.velocity *= rVariables.lumping_factor / rVariables.height_units;
            rVariables.velocity_grad /= rVariables.height_units;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::ComputeAuxMatrices(
            const BoundedMatrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rConvection,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rNonLinear )
    {
        // Initialize solution
        noalias(rMassMatrixVector) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rMassMatrixScalar) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rScalarGrad) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rVectorDiv) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rScalarDiff) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rVectorDiff) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rConvection) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rNonLinear) = ZeroMatrix(TNumNodes*3, TNumNodes*3);

        // Some auxilary definitions
        array_1d<double,TNumNodes> N;
        BoundedMatrix<double,2,TNumNodes*3> N_mom        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
        array_1d<double,TNumNodes*3> N_height             = ZeroVector(TNumNodes*3);    // Shape functions vector (for height unknown)
        array_1d<double,TNumNodes*3> DN_DX_mom            = ZeroVector(TNumNodes*3);    // Shape functions gradients vector (for velocity unknown)
        BoundedMatrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)
        BoundedMatrix<double,2,TNumNodes*3> Grad_mom_1   = ZeroMatrix(2,TNumNodes*3);
        BoundedMatrix<double,2,TNumNodes*3> Grad_mom_2   = ZeroMatrix(2,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> temp_non_lin;

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
                // Momentum divergence
                DN_DX_mom[  nnode*3] = rDN_DX(nnode,0);
                DN_DX_mom[1+nnode*3] = rDN_DX(nnode,1);
                // Height shape funtions
                N_height[2+nnode*3] = N[nnode];
                // Momentum shape functions
                N_mom(0,   nnode*3) = N[nnode];
                N_mom(1, 1+nnode*3) = N[nnode];
                // Momentum gradient
                Grad_mom_1(0,   nnode*3) = rDN_DX(nnode,0);
                Grad_mom_1(1, 1+nnode*3) = rDN_DX(nnode,0);
                Grad_mom_2(0,   nnode*3) = rDN_DX(nnode,1);
                Grad_mom_2(1, 1+nnode*3) = rDN_DX(nnode,1);
            }
            N_height     *= rVariables.height_units;
            DN_DX_height *= rVariables.height_units;

            noalias(rMassMatrixVector) += prod(trans(N_mom),N_mom);         // w * hu
            noalias(rMassMatrixScalar) += outer_prod(N_height,N_height);    // q * h

            noalias(rVectorDiv)  += outer_prod(N_height,DN_DX_mom);         // q * div(hu)
            noalias(rScalarGrad) += prod(trans(N_mom),DN_DX_height);        // w * grad(h)

            noalias(rVectorDiff) += outer_prod(DN_DX_mom,DN_DX_mom);        // div_w * div(hu)
            noalias(rScalarDiff) += prod(trans(DN_DX_height),DN_DX_height); // grad_q * grad(hu)

            noalias(rConvection) += rVariables.velocity[0] * prod(trans(N_mom), Grad_mom_1);     // w * u * grad(hu)
            noalias(rConvection) += rVariables.velocity[1] * prod(trans(N_mom), Grad_mom_2);     // w * u * grad(hu)

            temp_non_lin = prod(rVariables.velocity_grad, N_mom);
            rNonLinear += prod(trans(N_mom), temp_non_lin);                 // w * grad_u * hu
        }
    }

//----------------------------------------------------------------------

template class EulerConsVarElement<3>;
template class EulerConsVarElement<4>;

} // namespace Kratos
