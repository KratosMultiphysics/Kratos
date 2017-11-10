/*
==============================================================================
KratosShallowWaterApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//   
//   Project Name:        Kratos       
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                June 28th 2017
//   Revision:            1.4
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/cfd_variables.h"
#include "custom_elements/primitive_var_element.hpp"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    int PrimitiveVarElement<TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        
        const GeometryType& rGeom = this->GetGeometry();
        const PropertiesType& rProp = this->GetProperties();
        
        // verify nodal variables and dofs
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            // Verify basic variables
            if (rGeom[i].SolutionStepsDataHas( HEIGHT ) == false)
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable HEIGHT on node ", rGeom[i].Id() )
            
            if ( rGeom[i].SolutionStepsDataHas( VELOCITY ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable VELOCITY on node ", rGeom[i].Id() )
            
            // Verify auxiliar variables
            if (rGeom[i].SolutionStepsDataHas( BATHYMETRY ) == false)
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable BATHYMETRY on node ", rGeom[i].Id() )
            
            if (rGeom[i].SolutionStepsDataHas( RAIN ) == false)
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable RAIN on node ", rGeom[i].Id() )
            
            // Verify degrees of freedom
            if (rGeom[i].HasDofFor( HEIGHT ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable HEIGHT on node ", rGeom[i].Id() )
            
            if (rGeom[i].HasDofFor( VELOCITY_X ) == false ||
                rGeom[i].HasDofFor( VELOCITY_Y ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable VELOCITY on node ", rGeom[i].Id() )
        }
        
        // Verify properties
        if (MANNING.Key() == 0 ||
            rProp.Has( MANNING ) == false ||
            rProp[MANNING] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,"MANNING has Key zero, is not defined or has an invalid value at element", this->Id() )
        
        return 0;
        
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
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2> DN_DX;
        array_1d<double,TNumNodes> N;
        double Area;
        this-> CalculateGeometry(DN_DX,Area);
        double elem_length = this->ComputeElemSize(DN_DX);
        
        // Getting the values of shape functions on Integration Points
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;  // In this case, number of Gauss points and number of nodes coincides
        const GeometryType& rGeom = this->GetGeometry();
        Ncontainer = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );
        
        // Get nodal values for current step and projected variables (this function inlcudes the units conversion)
        this-> GetNodalValues(variables);
        
        // Get element values (this function inlcudes the units conversion)
        this-> GetElementValues(DN_DX, variables );
        double abs_vel = norm_2(variables.vector );
        double height43 = pow(variables.scalar, 1.33333 );
        
        // Compute stabilization and discontinuity capturing parameters
        double tau_u;
        double tau_h;
        double k_dc;
        this-> ComputeStabilizationParameters(variables, elem_length, tau_u, tau_h, k_dc);
        
        // Some auxilary definitions
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_q= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_w= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_u  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_u_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        
        this-> ComputeAuxMatrices(Ncontainer, DN_DX, variables, mass_matrix_q, mass_matrix_w, aux_w_grad_h, aux_q_div_u, aux_h_diffus, aux_u_diffus);
        
        noalias(mass_matrix) = mass_matrix_w + mass_matrix_q;
        
        // Build LHS
        // Cross terms
        noalias(rLeftHandSideMatrix)  = variables.scalar * aux_q_div_u;      // Add <q*h*div(u)> to Mass Eq.
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
        rVariables.manning2 = pow( GetProperties()[MANNING], 2 );
        rVariables.height_units = rCurrentProcessInfo[WATER_HEIGHT_UNIT_CONVERTER];
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea)
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
    double PrimitiveVarElement<TNumNodes>::ComputeElemSize(const boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX)
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
        l = sqrt(l) / static_cast<double>(TNumNodes);
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
    void PrimitiveVarElement<TNumNodes>::GetElementValues(const boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2>& rDN_DX,
                                                          ElementVariables& rVariables)
    {
        // Initialize outputs
        rVariables.scalar = 0;
        rVariables.vector = ZeroVector(2);
        rVariables.scalar_grad = ZeroVector(2);

        // integrate over the element
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.vector[0] += rVariables.unknown[  + 3*i];
            rVariables.vector[1] += rVariables.unknown[1 + 3*i];
            rVariables.scalar += rVariables.unknown[2 + 3*i];
            rVariables.scalar_grad[0] += rDN_DX(i,0) * rVariables.unknown[2 + 3*i];
            rVariables.scalar_grad[1] += rDN_DX(i,1) * rVariables.unknown[2 + 3*i];
        }

        rVariables.vector *= rVariables.lumping_factor;
        rVariables.scalar *= rVariables.lumping_factor * rVariables.height_units;
        rVariables.scalar_grad *= rVariables.height_units;
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
        double height_grad_norm = norm_2(rVariables.scalar_grad);
        
        // Compute stabilization parameters
        bool stabilization = true;
        double Ctau = rVariables.dyn_tau;  // 0.005 ~ 0.002  // 0.002; //
        double fheight = fabs(rVariables.scalar);
        if (stabilization && fheight > 1e-6)
        {
            rTauU = Ctau/rElemSize*pow(rVariables.gravity/fheight,0.5);
            rTauH = Ctau/rElemSize*pow(fheight/rVariables.gravity,0.5);
        }
        
        // Compute discontinuity capturing parameters
        bool discontinuity_capturing = true;
        double gradient_threshold = 1e-6;
        //~ double residual;
        if (discontinuity_capturing && height_grad_norm > gradient_threshold)
        {
            rKdc = 0.5*0.4*rElemSize*height_grad_norm;  // Residual formulation
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::ComputeAuxMatrices(
            const boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const boost::numeric::ublas::bounded_matrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff )
    {
        // Initialize solution
        noalias(rMassMatrixVector) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        noalias(rMassMatrixScalar) = ZeroMatrix(TNumNodes*3, TNumNodes*3);
        
        // Some auxilary definitions
        array_1d<double,TNumNodes> N;
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> N_vel        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> N_height     = ZeroMatrix(1,TNumNodes*3);  // Shape functions vector (for height unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> DN_DX_vel    = ZeroMatrix(1,TNumNodes*3);  // Shape functions gradients vector (for velocity unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)
        
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
