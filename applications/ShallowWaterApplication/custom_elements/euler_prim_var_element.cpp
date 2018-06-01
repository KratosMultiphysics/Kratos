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
#include "custom_elements/euler_prim_var_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    int EulerPrimVarElement<TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // Base class checks for positive Jacobian and Id > 0
        int ierr = Element::Check(rCurrentProcessInfo);
        if(ierr != 0) return ierr;

        // Check that all required variables have been registered
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
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT,rnode)
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
    void EulerPrimVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        unsigned int element_size = TNumNodes*3;
        if(rResult.size() != element_size)
            rResult.resize(element_size,false);                         // False says not to preserve existing storage!!

        GeometryType& rGeom = this->GetGeometry();
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
    void EulerPrimVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int element_size = TNumNodes*3;
        if(rElementalDofList.size() != element_size)
            rElementalDofList.resize(element_size);

        GeometryType& rGeom = this->GetGeometry();
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
    void EulerPrimVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
        this-> CalculateGeometry(DN_DX, Area);
        double elem_length = this->ComputeElemSize(DN_DX);

        // Getting the values of shape functions on Integration Points
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;  // In this case, number of Gauss points and number of nodes coincides
        const GeometryType& rGeom = this->GetGeometry();
        Ncontainer = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

        // Get nodal values for current step and previous step (this function inlcudes the units conversion)
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
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_convect  = ZeroMatrix(TNumNodes*3,TNumNodes*3);

        // ComputeAuxMatrices loops the Gauss points and build the matrices
        this-> ComputeAuxMatrices(Ncontainer, DN_DX, variables, mass_matrix_q, mass_matrix_w, aux_w_grad_h, aux_q_div_u, aux_h_diffus, aux_u_diffus, aux_convect);

        noalias(mass_matrix) = mass_matrix_w + mass_matrix_q;


        // Build LHS
        // Cross terms
        noalias(rLeftHandSideMatrix)  = variables.height * aux_q_div_u;     // Add <q*h*div(u)> to Mass Eq.
        noalias(rLeftHandSideMatrix) += variables.gravity * aux_w_grad_h;   // Add <w*g*grad(h)> to Momentum Eq.

        // Convective term
        noalias(rLeftHandSideMatrix) += aux_convect;                    // Add <q*u*grad(h)> and <w*u*grad(u)> to both Eq.'s

        // Inertia terms
        noalias(rLeftHandSideMatrix) += variables.dt_inv * mass_matrix; // Add <N,N> to both Eq.'s

        // Stabilization term
        noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;  // Add art. diff. to Mass Eq.
        noalias(rLeftHandSideMatrix) +=         tau_u  * aux_u_diffus;  // Add art. diff. to Momentum Eq.

        // Friction term
        noalias(rLeftHandSideMatrix) += variables.gravity * variables.manning2 * abs_vel / height43 * mass_matrix_w;

        // Build RHS
        // Source term (bathymetry contribution)
        noalias(rRightHandSideVector)  = -variables.gravity * prod(aux_w_grad_h, variables.depth);  // Add <w,-g*grad(H)> to RHS (Momentum Eq.)

        // Source terms (rain contribution)
        noalias(rRightHandSideVector) += prod(mass_matrix, variables.rain);

        // Inertia terms
        noalias(rRightHandSideVector) += variables.dt_inv * prod(mass_matrix, variables.proj_unk);

        // Substracting the Dirichlet term (since we use a residualbased approach)
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, variables.unknown);

        rRightHandSideVector *= Area * variables.lumping_factor;
        rLeftHandSideMatrix *= Area * variables.lumping_factor;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::GetNodalValues(ElementVariables& rVariables)
    {
        GeometryType& rGeom = this->GetGeometry();
        unsigned int counter = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.depth[counter] = 0;
            rVariables.rain[counter]  = 0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,1);
            counter++;

            rVariables.depth[counter] = 0;
            rVariables.rain[counter]  = 0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,1);
            counter++;

            rVariables.depth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY) / rVariables.height_units;
            rVariables.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
            rVariables.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT,1);
            counter++;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables)
    {
        // Initialize outputs
        rVariables.height = 0;
        rVariables.velocity = ZeroVector(2);

        // integrate over the element
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.velocity[0] += rVariables.unknown[  + 3*i];
            rVariables.velocity[1] += rVariables.unknown[1 + 3*i];
            rVariables.height += rVariables.unknown[2 + 3*i];
        }

        rVariables.velocity *= rVariables.lumping_factor;
        rVariables.height *= rVariables.lumping_factor * rVariables.height_units;
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::ComputeAuxMatrices(
            const BoundedMatrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rConvection )
    {
        // Initialize solution
        noalias(rMassMatrixVector) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rMassMatrixScalar) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rScalarGrad) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rVectorDiv) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rScalarDiff) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rVectorDiff) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rConvection) = ZeroMatrix(TNumNodes*3, TNumNodes*3);

        // Some auxilary definitions
        array_1d<double,TNumNodes> N;
        BoundedMatrix<double,2,TNumNodes*3> N_vel        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
        array_1d<double,TNumNodes*3> N_height             = ZeroVector(TNumNodes*3);    // Shape functions vector (for height unknown)
        array_1d<double,TNumNodes*3> DN_DX_vel            = ZeroVector(TNumNodes*3);    // Shape functions gradients vector (for velocity unknown)
        BoundedMatrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)
        BoundedMatrix<double,2,TNumNodes*3> Grad_vel_1   = ZeroMatrix(2,TNumNodes*3);
        BoundedMatrix<double,2,TNumNodes*3> Grad_vel_2   = ZeroMatrix(2,TNumNodes*3);
        array_1d<double,TNumNodes*3> temp_convect;

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
                DN_DX_vel[  nnode*3] = rDN_DX(nnode,0);
                DN_DX_vel[1+nnode*3] = rDN_DX(nnode,1);
                // Height shape funtions
                N_height[2+nnode*3] = N[nnode];
                // Velocity shape functions
                N_vel(0,   nnode*3) = N[nnode];
                N_vel(1, 1+nnode*3) = N[nnode];
                // Velocity gradient
                Grad_vel_1(0,   nnode*3) = rDN_DX(nnode,0);
                Grad_vel_1(1, 1+nnode*3) = rDN_DX(nnode,0);
                Grad_vel_2(0,   nnode*3) = rDN_DX(nnode,1);
                Grad_vel_2(1, 1+nnode*3) = rDN_DX(nnode,1);
            }
            N_height     *= rVariables.height_units;
            DN_DX_height *= rVariables.height_units;

            noalias(rMassMatrixVector) += prod(trans(N_vel),N_vel);         // q * h
            noalias(rMassMatrixScalar) += outer_prod(N_height,N_height);    // w * u

            noalias(rVectorDiv)  += outer_prod(N_height,DN_DX_vel);         // q * div_u
            noalias(rScalarGrad) += prod(trans(N_vel),DN_DX_height);        // w * grad_h

            noalias(rVectorDiff) += outer_prod(DN_DX_vel,DN_DX_vel);        // div_w * div_u
            noalias(rScalarDiff) += prod(trans(DN_DX_height),DN_DX_height); // grad_q * grad_h

            temp_convect = prod(rVariables.velocity, DN_DX_height);
            noalias(rConvection) += outer_prod(N_height, temp_convect);     // q * u * grad_h
            noalias(rConvection) += rVariables.velocity[0] * prod(trans(N_vel), Grad_vel_1);     // w * u * grad_u
            noalias(rConvection) += rVariables.velocity[1] * prod(trans(N_vel), Grad_vel_2);     // w * u * grad_u
        }
    }

//----------------------------------------------------------------------

template class EulerPrimVarElement<3>;
template class EulerPrimVarElement<4>;

} // namespace Kratos
