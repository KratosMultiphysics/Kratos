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
        bounded_matrix<double,TNumNodes, 2> DN_DX;
        array_1d<double,TNumNodes> N;
        double Area;
        this-> CalculateGeometry(DN_DX, Area);
        double elem_length = this->ComputeElemSize(DN_DX);

        // Getting the values of shape functions on Integration Points
        bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;  // In this case, number of Gauss points and number of nodes coincides
        const GeometryType& rGeom = this->GetGeometry();
        Ncontainer = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

        // Get nodal values for current step and projected variables (this function inlcudes the units conversion)
        this-> GetNodalValues(variables);

        // Get element values (this function inlcudes the units conversion)
        this-> GetElementValues(DN_DX, variables );
        double abs_mom = norm_2(variables.vector );
        double height73 = std::pow(variables.scalar, 2.33333 );

        // Compute stabilization and discontinuity capturing parameters
        double tau_m;
        double tau_h;
        double k_dc;
        this-> ComputeStabilizationParameters(variables, elem_length, tau_m, tau_h, k_dc);




        //~ // Getting gravity
        //~ double gravity = rCurrentProcessInfo[GRAVITY_Z];
        //~ 
        //~ // Getting properties
        //~ double manning2 = pow( GetProperties()[MANNING], 2 );
        //~ 
        //~ // Getting the time step (not fixed to allow variable time step)
        //~ const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        //~ double BDFcoeffs[2] = {1.0/delta_t, 1.0/delta_t};
        //~ 
        //~ // Compute the geometry
        //~ boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2> DN_DX;
        //~ array_1d<double,TNumNodes> N;
        //~ double Area;
        //~ this-> CalculateGeometry(DN_DX,Area);
        //~ double elem_length = this->ComputeElemSize(DN_DX);
        //~ 
        //~ // Getting the values of shape functions on Integration Points
        //~ boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;  // In this case, number of Gauss points and number of nodes coincides
        //~ const GeometryType& rGeom = this->GetGeometry();
        //~ Ncontainer = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );
        //~ 
        //~ // Get nodal values for current and previous step
        //~ array_1d<double, TNumNodes*3> v_depth;
        //~ array_1d<double, TNumNodes*3> v_rain;
        //~ array_1d<double, TNumNodes*3> v_unknown;
        //~ array_1d<double, TNumNodes*3> v_unknown_n;
        //~ GetNodalValues(v_depth, v_unknown, v_unknown_n );
        //~ 
        //~ // Get element values
        //~ double height;
        //~ array_1d<double,2> momentum;
        //~ array_1d<double,2> velocity;
        //~ boost::numeric::ublas::bounded_matrix<double,2,2> grad_u;
        //~ GetElementValues(DN_DX, v_unknown, momentum, height, velocity, grad_u);
        //~ double height73 = pow(height, 2.33333);
        //~ double abs_mom = norm_2(momentum);
        //~ 
        //~ // Some auxilary definitions
        //~ boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> N_mom        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for momentum unknown)
        //~ boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> N_height     = ZeroMatrix(1,TNumNodes*3);  // Shape functions vector (for height unknown)
        //~ boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> DN_DX_mom    = ZeroMatrix(1,TNumNodes*3);  // Shape functions for vector divergence (for momentum unknown)
        //~ boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions for scalar gradient (for height unknown)
        //~ boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> grad_mom_1   = ZeroMatrix(2,TNumNodes*3);  // Shape functions for vector gradient (for momentum unknown)
        //~ boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> grad_mom_2   = ZeroMatrix(2,TNumNodes*3);  // Shape functions for vector gradient (for momentum unknown)
        //
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_q= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_w= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_convect  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_non_lin  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_m  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_m_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);

        //~ // Loop on Gauss points. In this case, number of Gauss points are equal to number of nodes
        //~ for (unsigned int igauss = 0; igauss < TNumNodes; igauss++)
        //~ {
            //~ noalias(N) = row(Ncontainer, igauss);
            //~ 
            //~ // Build shape and derivatives functions at Gauss points
            //~ for (unsigned int nnode = 0; nnode < TNumNodes; nnode++)
            //~ {
                //~ // Height gradient
                //~ DN_DX_height(0, 2+nnode*3) = DN_DX(nnode,0);
                //~ DN_DX_height(1, 2+nnode*3) = DN_DX(nnode,1);
                //~ // Momentum divergence
                //~ DN_DX_mom(0,   nnode*3) = DN_DX(nnode,0);
                //~ DN_DX_mom(0, 1+nnode*3) = DN_DX(nnode,1);
                //~ // Momentum gradient
                //~ grad_mom_1(0,   nnode*3) = DN_DX(nnode,0);
                //~ grad_mom_1(0, 1+nnode*3) = DN_DX(nnode,0);
                //~ grad_mom_2(1,   nnode*3) = DN_DX(nnode,1);
                //~ grad_mom_2(1, 1+nnode*3) = DN_DX(nnode,1);
                //~ // Height shape funtions
                //~ N_height(0, 2+nnode*3) = N[nnode];
                //~ // Momentum shape functions
                //~ N_mom(0,   nnode*3) = N[nnode];
                //~ N_mom(1, 1+nnode*3) = N[nnode];
            //~ }
            //~ noalias(mass_matrix_w)+= prod(trans(N_mom),N_mom);
            //~ noalias(mass_matrix_q)+= prod(trans(N_height),N_height);
            //~ 
            //~ noalias(aux_convect)  += velocity[0] * prod(trans(N_mom),grad_mom_1);
            //~ noalias(aux_convect)  += velocity[1] * prod(trans(N_mom),grad_mom_2);
            //~ 
            //~ noalias(aux_non_lin)  += prod(trans(N_mom),Matrix(prod(grad_u,N_mom)));
            //~ 
            //~ noalias(aux_q_div_m)  += prod(trans(N_height),DN_DX_mom);
            //~ noalias(aux_w_grad_h) += prod(trans(N_mom),DN_DX_height);
            //~ 
            //~ noalias(aux_m_diffus) += prod(trans(DN_DX_mom),DN_DX_mom);
            //~ noalias(aux_h_diffus) += prod(trans(DN_DX_height),DN_DX_height);
        //~ }

        // ComputeAuxMatrices: loops the Gauss points and build the matrices
        this-> ComputeAuxMatrices(Ncontainer, DN_DX, variables, mass_matrix_q, mass_matrix_w, aux_w_grad_h, aux_q_div_m, aux_h_diffus, aux_m_diffus, aux_convect, aux_non_lin);

        noalias(mass_matrix) = mass_matrix_w + mass_matrix_q;

        // Build LHS
        // Inertia terms
        noalias(rLeftHandSideMatrix)  = variables.dt_inv * mass_matrix;     // Add <N,N> to both Eq's

        // Convective and non linear terms
        //~ noalias(rLeftHandSideMatrix) += aux_convect;                    // Add <w,u*grad(hu)> to Momentum Eq.
        noalias(rLeftHandSideMatrix) += aux_non_lin;                    // Add <w,grad(u)*hu> to Momentum Eq.

        // Cross terms
        noalias(rLeftHandSideMatrix) += aux_q_div_m;                     // Add <q,div(hu)> to Mass Eq.
        noalias(rLeftHandSideMatrix) += variables.gravity * variables.scalar * aux_w_grad_h; // Add <w,g*h*grad(h)> to Momentum Eq.

        // Stabilization terms
        noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;  // Add art. diff. to Mass Eq.
        noalias(rLeftHandSideMatrix) +=         tau_m  * aux_m_diffus;  // Add art. diff. to Momentum Eq.

        // Friction term
        noalias(rLeftHandSideMatrix) += variables.gravity * variables.manning2 * abs_mom / height73 * mass_matrix_w;

        // Build RHS
        // Source terms (bathymetry contribution)
        noalias(rRightHandSideVector)  = -variables.gravity * variables.scalar * prod(aux_w_grad_h, variables.depth); // Add <w,-g*h*grad(H)> to RHS (Momentum Eq.)

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
    void EulerConsVarElement<TNumNodes>::GetElementValues(const bounded_matrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables)
    {
        // Initialize outputs
        rVariables.scalar = 0;
        rVariables.vector = ZeroVector(2);
        rVariables.scalar_grad = ZeroVector(2);
        rVariables.vector_grad = ZeroMatrix(2,2);
        rVariables.vector_div = 0;

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
            rVariables.vector[0] += rVariables.unknown[  + 3*i];
            rVariables.vector[1] += rVariables.unknown[1 + 3*i];
            rVariables.scalar += rVariables.unknown[2 + 3*i];
            rVariables.scalar_grad[0] += rDN_DX(i,0) * rVariables.unknown[2 + 3*i];
            rVariables.scalar_grad[1] += rDN_DX(i,1) * rVariables.unknown[2 + 3*i];
            if (near_dry)
            {
                //~ rVel[0] += rVariables.unknown[  + 3*i];
                //~ rVel[1] += rVariables.unknown[1 + 3*i];
                rVariables.vector_grad(0,0) += rDN_DX(i,0) * rVariables.unknown[  + 3*i];
                rVariables.vector_grad(1,1) += rDN_DX(i,0) * rVariables.unknown[1 + 3*i];
                rVariables.vector_grad(0,0) += rDN_DX(i,1) * rVariables.unknown[  + 3*i];
                rVariables.vector_grad(1,1) += rDN_DX(i,1) * rVariables.unknown[1 + 3*i];
            }
            else
            {
                //~ rVel[0] += rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                //~ rVel[1] += rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.vector_grad(0,0) += rDN_DX(i,0) * rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.vector_grad(0,1) += rDN_DX(i,0) * rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.vector_grad(1,0) += rDN_DX(i,1) * rVariables.unknown[  + 3*i] / rVariables.unknown[2 + 3*i];
                rVariables.vector_grad(1,1) += rDN_DX(i,1) * rVariables.unknown[1 + 3*i] / rVariables.unknown[2 + 3*i];
            }
        }

        rVariables.vector *= rVariables.lumping_factor;
        rVariables.scalar *= rVariables.lumping_factor * rVariables.height_units;
        //~ rVel    *= rVariables.lumping_factor;
        if (near_dry)
        {
            //~ rVel /= rVariables.scalar;
            rVariables.vector_grad /= rVariables.scalar;
        }
        else
        {
            rVariables.vector_grad /= rVariables.height_units;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::ComputeAuxMatrices(
            const bounded_matrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const bounded_matrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rConvection,
            bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rNonLinear )
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
        bounded_matrix<double,2,TNumNodes*3> N_mom        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
        array_1d<double,TNumNodes*3> N_height             = ZeroVector(TNumNodes*3);    // Shape functions vector (for height unknown)
        array_1d<double,TNumNodes*3> DN_DX_mom            = ZeroVector(TNumNodes*3);    // Shape functions gradients vector (for velocity unknown)
        bounded_matrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)
        bounded_matrix<double,2,TNumNodes*3> Grad_mom_1   = ZeroMatrix(2,TNumNodes*3);
        bounded_matrix<double,2,TNumNodes*3> Grad_mom_2   = ZeroMatrix(2,TNumNodes*3);
        array_1d<double,TNumNodes*3> temp_convect;
        bounded_matrix<double,TNumNodes*3,TNumNodes*3> temp_non_lin;

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

            noalias(rVectorDiv)  += outer_prod(N_height,DN_DX_mom);         // q * div_u
            noalias(rScalarGrad) += prod(trans(N_mom),DN_DX_height);        // w * grad_h

            noalias(rVectorDiff) += outer_prod(DN_DX_mom,DN_DX_mom);        // div_w * div_u
            noalias(rScalarDiff) += prod(trans(DN_DX_height),DN_DX_height); // grad_q * grad_h

            //~ temp_convect = prod(rVariables.velocity, DN_DX_height);
            //~ noalias(rConvection) += outer_prod(N_height, temp_convect);     // q * u * grad_h
            //~ noalias(rConvection) += rVariables.velocity[0] * prod(trans(N_mom), Grad_mom_1);     // w * u * grad_u
            //~ noalias(rConvection) += rVariables.velocity[1] * prod(trans(N_mom), Grad_mom_2);     // w * u * grad_u

            temp_non_lin = prod(rVariables.vector_grad, N_mom);
            rNonLinear += prod(trans(N_mom), temp_non_lin);                  // w * grad_u * hu
        }
    }

//----------------------------------------------------------------------

template class EulerConsVarElement<3>;
template class EulerConsVarElement<4>;

} // namespace Kratos
