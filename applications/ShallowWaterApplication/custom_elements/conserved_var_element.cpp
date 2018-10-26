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
#include "custom_elements/conserved_var_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    int ConservedVarElement<TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // Base class checks for positive Jacobian and Id > 0
        int ierr = Element::Check(rCurrentProcessInfo);
        if(ierr != 0) return ierr;

        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(MOMENTUM)
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
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1,rnode)
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
    void ConservedVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
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
    void ConservedVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
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
    void ConservedVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) 
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

        // Some auxiliar definitions
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_q= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_w= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_m  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_m_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        BoundedMatrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);

        // ComputeAuxMatrices loops the Gauss points and build the matrices
        this-> ComputeAuxMatrices(Ncontainer, DN_DX, variables, mass_matrix_q, mass_matrix_w, aux_w_grad_h, aux_q_div_m, aux_h_diffus, aux_m_diffus);

        noalias(mass_matrix) = mass_matrix_w + mass_matrix_q;

        // Build LHS 
        // Cross terms 
        //~ noalias(rLeftHandSideMatrix)  = aux_q_div_m;                                         // Add <q*div(hu)> to Mass Eq.
        noalias(rLeftHandSideMatrix)  = variables.gravity * variables.height * aux_w_grad_h; // Add <w,g*h*grad(h)> to Momentum Eq.

        // Inertia terms 
        noalias(rLeftHandSideMatrix) += variables.dt_inv * mass_matrix;        // Add <N,N> to both Eq's

        // Non linear terms 
        noalias(rLeftHandSideMatrix) += variables.velocity_div * mass_matrix;    // Add <q,div(u)*h> to Mass Eq. and <w,div(u)*hu> to Momentum Eq.
        //~ noalias(rLeftHandSideMatrix) += variables.velocity_div * mass_matrix_w;  // Add <w,div(u)*hu> to Momentum Eq.
        //~ noalias(rLeftHandSideMatrix) += aux_non_linear;                          // Add  and <w,hu*grad(u)> to Momentum Eq.

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
    void ConservedVarElement<TNumNodes>::GetNodalValues(ElementVariables& rVariables)
    {
        GeometryType& rGeom = this->GetGeometry();
        unsigned int counter = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.depth[counter] = 0;
            rVariables.rain[counter]  = 0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_X);
            counter++;

            rVariables.depth[counter] = 0;
            rVariables.rain[counter]  = 0;
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Y);
            counter++;

            rVariables.depth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY) / rVariables.height_units;
            rVariables.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
            rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
            rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(PROJECTED_SCALAR1);
            //~ rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(HEIGHT,1);
            counter++;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes > 
    void ConservedVarElement<TNumNodes>::GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables)
    {
        // Initialize outputs
        rVariables.height = 0;
        rVariables.momentum = ZeroVector(2);
        rVariables.height_grad = ZeroVector(2);
        rVariables.velocity_grad = ZeroMatrix(2,2);
        rVariables.velocity_div = 0;

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
            rVariables.height_grad[0] += rDN_DX(i,0) * rVariables.unknown[2 + 3*i];
            rVariables.height_grad[1] += rDN_DX(i,1) * rVariables.unknown[2 + 3*i];
            if (near_dry)
            {
                rVariables.velocity_grad(0,0) += rDN_DX(i,0) * rVariables.unknown[  + 3*i];
                rVariables.velocity_grad(0,1) += rDN_DX(i,0) * rVariables.unknown[1 + 3*i];
                rVariables.velocity_grad(1,0) += rDN_DX(i,1) * rVariables.unknown[  + 3*i];
                rVariables.velocity_grad(1,1) += rDN_DX(i,1) * rVariables.unknown[1 + 3*i];
                rVariables.velocity_div += rDN_DX(i,0) * rVariables.unknown[  + 3*i];
                rVariables.velocity_div += rDN_DX(i,1) * rVariables.unknown[1 + 3*i];
            }
            else
            {
                rVariables.velocity_grad(0,0) += rDN_DX(i,0) * rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_grad(0,1) += rDN_DX(i,0) * rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_grad(1,0) += rDN_DX(i,1) * rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_grad(1,1) += rDN_DX(i,1) * rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_div += rDN_DX(i,0) * rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.velocity_div += rDN_DX(i,1) * rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
            }
        }

        rVariables.momentum *= rVariables.lumping_factor;
        rVariables.height *= rVariables.lumping_factor * rVariables.height_units;
        rVariables.height_grad *= rVariables.height_units;
        if (near_dry)
        {
            rVariables.velocity_grad /= rVariables.height;
            rVariables.velocity_div  /= rVariables.height;
        }
        else
        {
            rVariables.velocity_grad /= rVariables.height_units;
            rVariables.velocity_div  /= rVariables.height_units;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void ConservedVarElement<TNumNodes>::CalculateLumpedMassMatrix(BoundedMatrix<double, TNumNodes*3, TNumNodes*3>& rMassMatrix) 
    {
        const unsigned int element_size = 3*TNumNodes;
        rMassMatrix  = IdentityMatrix(element_size, element_size);
        rMassMatrix /= static_cast<double>(TNumNodes);
    }

//----------------------------------------------------------------------

template class ConservedVarElement<3>;
template class ConservedVarElement<4>;

} // namespace Kratos
