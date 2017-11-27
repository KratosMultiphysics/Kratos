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
#include "includes/define.h"
#include "custom_elements/euler_prim_var_element.hpp"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
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
    void EulerPrimVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
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
    void EulerPrimVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        // Resize of the Left and Right Hand side
        unsigned int element_size = TNumNodes*3;
        if(rLeftHandSideMatrix.size1() != element_size)
            rLeftHandSideMatrix.resize(element_size,element_size,false); // False says not to preserve existing storage!!
        
        if(rRightHandSideVector.size() != element_size)
            rRightHandSideVector.resize(element_size,false);             // False says not to preserve existing storage!!
        
        // Getting gravity
        //~ array_1d<double,3> v_gravity = rCurrentProcessInfo[GRAVITY];
        double gravity = 9.8; //-v_gravity[2];
        
        // Getting the time step (not fixed to allow variable time step)
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        double BDFcoeffs[2] = {1.0/delta_t, 1.0/delta_t};
        
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
        
        // Get nodal values for current and previous step
        array_1d<double, TNumNodes*3> v_depth;
        array_1d<double, TNumNodes*3> v_unknown;
        array_1d<double, TNumNodes*3> v_unknown_n;
        boost::numeric::ublas::bounded_matrix<double,1,2> velocity;
        double height;
        GetNodalValues(v_depth,v_unknown,v_unknown_n);
        GetElementValues(DN_DX,v_unknown,height,velocity);
        
        // Some auxilary definitions
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> N_vel        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> N_height     = ZeroMatrix(1,TNumNodes*3);  // Shape functions vector (for height unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> DN_DX_vel    = ZeroMatrix(1,TNumNodes*3);  // Shape functions for vector divergence (for velocity unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions for scalar gradients (for height unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> grad_vel_1   = ZeroMatrix(2,TNumNodes*3);  // Shape functions for vector gradient (for velocity unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> grad_vel_2   = ZeroMatrix(2,TNumNodes*3);  // Shape functions for vector gradient (for velocity unknown)
        //
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_convect  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_u  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_u_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        
        // Loop on Gauss points. In this case, number of Gauss points and number of nodes coincides
        for(unsigned int igauss = 0; igauss < TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer, igauss);
            
            // Build shape and derivatives functions at Gauss points
            for(unsigned int nnode = 0; nnode < TNumNodes; nnode++)
            {
                // Height gradient
                DN_DX_height(0, 2+nnode*3) = DN_DX(nnode,0);
                DN_DX_height(1, 2+nnode*3) = DN_DX(nnode,1);
                // Velocity divergence
                DN_DX_vel(0,  nnode*3) = DN_DX(nnode,0);
                DN_DX_vel(0,1+nnode*3) = DN_DX(nnode,1);
                // Height shape funtions
                N_height(0,2+nnode*3) = N[nnode];
                // Velocity shape functions
                N_vel(0,  nnode*3) = N[nnode];
                N_vel(1,1+nnode*3) = N[nnode];
                // Velocity gradient
                grad_vel_1(0,  nnode*3) = DN_DX(nnode,0);
                grad_vel_1(1,1+nnode*3) = DN_DX(nnode,0);
                grad_vel_2(0,  nnode*3) = DN_DX(nnode,1);
                grad_vel_2(0,1+nnode*3) = DN_DX(nnode,1);
            }
            
            noalias(mass_matrix)  += prod(trans(N_vel),N_vel);
            noalias(mass_matrix)  += prod(trans(N_height),N_height);
            
            noalias(aux_convect)  += prod(trans(N_height),Matrix(prod(velocity,DN_DX_height)));
            noalias(aux_convect)  += velocity(0,0) * prod(trans(N_vel),grad_vel_1);
            noalias(aux_convect)  += velocity(0,1) * prod(trans(N_vel),grad_vel_2);
            
            noalias(aux_q_div_u)  += prod(trans(N_height),DN_DX_vel);
            noalias(aux_w_grad_h) += prod(trans(N_vel),DN_DX_height);
            
            noalias(aux_u_diffus) += prod(trans(DN_DX_vel),DN_DX_vel);
            noalias(aux_h_diffus) += prod(trans(DN_DX_height),DN_DX_height);
        }
        
        // Copmute stabilization parameters
        bool stabilization = true;
        double Ctau = 0.005;
        double tau_u = 0, tau_h = 0;
        double fheight = fabs(height);
        if (stabilization && fheight > 1e-6)
        {
            tau_u = Ctau/elem_length*pow(gravity/fheight,0.5);
            tau_h = Ctau/elem_length*pow(fheight/gravity,0.5);
        }
        // Compute discontinuity capturing parameters
        bool discontinuity_capturing = true;
        double gradient_threshold = 1e-6;
        //~ double residual;
        double height_grad_norm = norm_2(prod(DN_DX_height,v_unknown));
        double k_dc = 0;
        if( discontinuity_capturing && height_grad_norm > gradient_threshold )
        {
            k_dc = 0.5*0.4*elem_length*height_grad_norm;  // Residual formulation
        }
        
        // Build LHS
        // Cross terms
        noalias(rLeftHandSideMatrix)  = height * aux_q_div_u;           // Add <q*h*div(u)> to Mass Eq.
        noalias(rLeftHandSideMatrix) += gravity * aux_w_grad_h;         // Add <w*g*grad(h)> to Momentum Eq.
        
        // Convective term
        noalias(rLeftHandSideMatrix) += aux_convect;                    // Add <q*u*grad(h)> and <w*u*grad(u)> to both Eq.'s
        
        // Inertia terms
        noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * mass_matrix;     // Add <N,N> to both Eq.'s
        
        // Stabilization term
        noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;  // Add art. diff. to Mass Eq.
        noalias(rLeftHandSideMatrix) +=         tau_u  * aux_u_diffus;  // Add art. diff. to Momentum Eq.
        
        // Build RHS
        // Source term (bathymetry contribution)
        noalias(rRightHandSideVector)  = -gravity * prod(aux_w_grad_h, v_depth);
        
        // Inertia terms
        noalias(rRightHandSideVector) += BDFcoeffs[1] * prod(mass_matrix, v_unknown_n);
        
        // Substracting the Dirichlet term (since we use a residualbased approach)
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, v_unknown);
        
        rRightHandSideVector *= Area / static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Area  / static_cast<double>(TNumNodes);
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM || rVariable == MIU)
        {
            for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
                rValues[PointNumber] = double(this->GetValue(rVariable));
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea)
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
    double EulerPrimVarElement<TNumNodes>::ComputeElemSize(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX)
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
    void EulerPrimVarElement<TNumNodes>::GetNodalValues(array_1d<double, TNumNodes*3>& rdepth,
                                                        array_1d<double, TNumNodes*3>& runkn, 
                                                        array_1d<double, TNumNodes*3>& rprev
                                                       )
    {
        GeometryType& rGeom = GetGeometry();
        unsigned int counter = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rdepth[counter] = 0;
            runkn[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);
            rprev[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,1);
            counter++;

            rdepth[counter] = 0;
            runkn[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
            rprev[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,1);
            counter++;

            rdepth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY);
            runkn[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
            rprev[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT,1);
            counter++;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerPrimVarElement<TNumNodes>::GetElementValues(boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2>& rDN_DX,
                                                          array_1d<double, TNumNodes*3>& r_nodal_var,
                                                          double& rheight,
                                                          boost::numeric::ublas::bounded_matrix<double,1, 2>& rvel
                                                         )
    {
        double lumping_factor = 1 / double(TNumNodes);
        
        rheight = 0;
        rvel = ZeroMatrix(1,2);
        
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rvel(0,0) += r_nodal_var[  + 3*i];
            rvel(0,1) += r_nodal_var[1 + 3*i];
            rheight   += r_nodal_var[2 + 3*i];
        }

        rheight *= lumping_factor;
        rvel    *= lumping_factor;
    }

//----------------------------------------------------------------------

template class EulerPrimVarElement<3>;
template class EulerPrimVarElement<4>;

} // namespace Kratos
