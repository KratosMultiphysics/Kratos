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
#include "includes/define.h"
#include "custom_elements/euler_cons_var_element.hpp"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        unsigned int element_size = TNumNodes*3;
        if(rResult.size() != element_size)
            rResult.resize(element_size,false);                         // False says not to preserve existing storage!!

        GeometryType& rGeom = GetGeometry();
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

        GeometryType& rGeom = GetGeometry();
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

        // Getting gravity
        double gravity = rCurrentProcessInfo[GRAVITY_Z];
        
        // Getting properties
        double manning2 = pow( GetProperties()[MANNING], 2 );
        
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
        array_1d<double, TNumNodes*3> v_rain;
        array_1d<double, TNumNodes*3> v_unknown;
        array_1d<double, TNumNodes*3> v_unknown_n;
        GetNodalValues(v_depth, v_unknown, v_unknown_n );
        
        // Get element values
        double height;
        array_1d<double,2> momentum;
        array_1d<double,2> velocity;
        boost::numeric::ublas::bounded_matrix<double,2,2> grad_u;
        GetElementValues(DN_DX, v_unknown, momentum, height, velocity, grad_u);
        double height73 = pow(height, 2.33333);
        double abs_mom = norm_2(momentum);
        
        // Some auxilary definitions
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> N_mom        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for momentum unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> N_height     = ZeroMatrix(1,TNumNodes*3);  // Shape functions vector (for height unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> DN_DX_mom    = ZeroMatrix(1,TNumNodes*3);  // Shape functions for vector divergence (for momentum unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions for scalar gradient (for height unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> grad_mom_1   = ZeroMatrix(2,TNumNodes*3);  // Shape functions for vector gradient (for momentum unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> grad_mom_2   = ZeroMatrix(2,TNumNodes*3);  // Shape functions for vector gradient (for momentum unknown)
        //
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_q= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix_w= ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_convect  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_non_lin  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_m  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_m_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);

        // Loop on Gauss points. In this case, number of Gauss points are equal to number of nodes
        for (unsigned int igauss = 0; igauss < TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer, igauss);
            
            // Build shape and derivatives functions at Gauss points
            for (unsigned int nnode = 0; nnode < TNumNodes; nnode++)
            {
                // Height gradient
                DN_DX_height(0, 2+nnode*3) = DN_DX(nnode,0);
                DN_DX_height(1, 2+nnode*3) = DN_DX(nnode,1);
                // Momentum divergence
                DN_DX_mom(0,   nnode*3) = DN_DX(nnode,0);
                DN_DX_mom(0, 1+nnode*3) = DN_DX(nnode,1);
                // Momentum gradient
                grad_mom_1(0,   nnode*3) = DN_DX(nnode,0);
                grad_mom_1(0, 1+nnode*3) = DN_DX(nnode,0);
                grad_mom_2(1,   nnode*3) = DN_DX(nnode,1);
                grad_mom_2(1, 1+nnode*3) = DN_DX(nnode,1);
                // Height shape funtions
                N_height(0, 2+nnode*3) = N[nnode];
                // Momentum shape functions
                N_mom(0,   nnode*3) = N[nnode];
                N_mom(1, 1+nnode*3) = N[nnode];
            }
            noalias(mass_matrix_w)+= prod(trans(N_mom),N_mom);
            noalias(mass_matrix_q)+= prod(trans(N_height),N_height);
            
            noalias(aux_convect)  += velocity[0] * prod(trans(N_mom),grad_mom_1);
            noalias(aux_convect)  += velocity[1] * prod(trans(N_mom),grad_mom_2);
            
            noalias(aux_non_lin)  += prod(trans(N_mom),Matrix(prod(grad_u,N_mom)));
            
            noalias(aux_q_div_m)  += prod(trans(N_height),DN_DX_mom);
            noalias(aux_w_grad_h) += prod(trans(N_mom),DN_DX_height);
            
            noalias(aux_m_diffus) += prod(trans(DN_DX_mom),DN_DX_mom);
            noalias(aux_h_diffus) += prod(trans(DN_DX_height),DN_DX_height);
        }
        noalias(mass_matrix) = mass_matrix_w + mass_matrix_q;

        // Copmute stabilization parameters
        bool stabilization = true;
        double height_threshold = 1e-6;
        double Ctau = 0.01;       // Stabilization parameter >0.005 (R.Codina, CMAME 197, 2008, 1305-1322)
        double tau_h = 0, tau_m = 0;
        if (stabilization && height > height_threshold)
        {
            tau_m = Ctau/elem_length*pow(gravity/height,0.5);
            tau_h = Ctau/elem_length*pow(height/gravity,0.5);
        }
        // Compute discontinuity capturing parameters
        bool discontinuity_capturing = true;
        double gradient_threshold = 1e-6;    // Shock capturing parameters
        //~ double residual = norm_1(prod(msN_height,ms_unknown)) * norm_1(prod(msDN_DX_mom,ms_unknown)) + dt_inv*norm_1(prod(msN_height, (ms_unknown - ms_proj_unknown)));
        double height_grad_norm = norm_2(prod(DN_DX_height,v_unknown));
        double k_dc = 0;
        if (discontinuity_capturing && height_grad_norm > gradient_threshold)
        {
            // Residual formulation
            k_dc = 0.5*0.4*elem_length*height_grad_norm;
        }

        // Build LHS
        // Inertia terms
        noalias(rLeftHandSideMatrix)  = BDFcoeffs[0] * mass_matrix;     // Add <N,N> to both Eq's

        // Convective and non linear terms
        noalias(rLeftHandSideMatrix) += aux_convect;                    // Add <w,u*grad(hu)> to Momentum Eq.
        noalias(rLeftHandSideMatrix) += aux_non_lin;                    // Add <w,drad(u)*hu> to Momentum Eq.

        // Cross terms
        noalias(rLeftHandSideMatrix) += aux_q_div_m;                     // Add <q,div(hu)> to Mass Eq.
        noalias(rLeftHandSideMatrix) += gravity * height * aux_w_grad_h; // Add <w,g*h*grad(h)> to Momentum Eq.

        // Stabilization terms
        noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;  // Add art. diff. to Mass Eq.
        noalias(rLeftHandSideMatrix) +=         tau_m  * aux_m_diffus;  // Add art. diff. to Momentum Eq.

        // Friction term
        noalias(rLeftHandSideMatrix) += gravity * manning2 * abs_mom / height73 * mass_matrix_w;

        // Build RHS
        // Source terms (bathymetry contribution)
        noalias(rRightHandSideVector)  = -gravity * height * prod(aux_w_grad_h, v_depth); // Add <w,-g*h*grad(H)> to RHS (Momentum Eq.)

        // Source terms (rain contribution)
        noalias(rRightHandSideVector) += prod(mass_matrix, v_rain);

        // Inertia terms
        noalias(rRightHandSideVector) += BDFcoeffs[1] * prod(mass_matrix, v_unknown_n);

        // Subtracting the dirichlet term (since we use a residualbased approach)
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, v_unknown);

        rRightHandSideVector *= Area / static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Area  / static_cast<double>(TNumNodes);

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM || rVariable == MIU)
        {
            for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
                rValues[PointNumber] = double(this->GetValue(rVariable));
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea)
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
    double EulerConsVarElement<TNumNodes>::ComputeElemSize(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX)
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
    void EulerConsVarElement<TNumNodes>::GetNodalValues(array_1d<double, TNumNodes*3>& rdepth,
                                                        array_1d<double, TNumNodes*3>& runkn, 
                                                        array_1d<double, TNumNodes*3>& rprev
                                                       )
    {
        GeometryType& rGeom = GetGeometry();
        unsigned int counter = 0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rdepth[counter] = 0;
            runkn[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X);
            rprev[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X,1);
            counter++;

            rdepth[counter] = 0;
            runkn[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y);
            rprev[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y,1);
            counter++;

            rdepth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY);
            runkn[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
            rprev[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT,1);
            counter++;
        }
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::GetElementValues(boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2>& rDN_DX,
                                                          array_1d<double,TNumNodes*3>& rNodalVar,
                                                          array_1d<double,2>& rMom,
                                                          double& rHeight,
                                                          array_1d<double,2>& rVel,
                                                          boost::numeric::ublas::bounded_matrix<double,2,2>& rGradU
                                                         )
    {
        double lumping_factor = 1 / double(TNumNodes);
        bool near_dry = false;
        
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if (rNodalVar[2 + 3*i] < 1e-1)
                near_dry = true;
        }
        
        rMom = ZeroVector(2);
        rHeight = 0;
        rVel = ZeroVector(2);
        rGradU = ZeroMatrix(2,2);
        
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rMom[0] += rNodalVar[  + 3*i];
            rMom[1] += rNodalVar[1 + 3*i];
            rHeight += rNodalVar[2 + 3*i];
            if (near_dry)
            {
                rVel[0] += rNodalVar[  + 3*i];
                rVel[1] += rNodalVar[1 + 3*i];
                rGradU(0,0) += rDN_DX(i,0) * rNodalVar[  + 3*i];
                rGradU(1,1) += rDN_DX(i,0) * rNodalVar[1 + 3*i];
                rGradU(0,0) += rDN_DX(i,1) * rNodalVar[  + 3*i];
                rGradU(1,1) += rDN_DX(i,1) * rNodalVar[1 + 3*i];
            }
            else
            {
                rVel[0] += rNodalVar[  + 3*i] / rNodalVar[2 + 3*i];
                rVel[1] += rNodalVar[1 + 3*i] / rNodalVar[2 + 3*i];
                rGradU(0,0) += rDN_DX(i,0) * rNodalVar[  + 3*i] / rNodalVar[2 + 3*i];
                rGradU(0,1) += rDN_DX(i,0) * rNodalVar[1 + 3*i] / rNodalVar[2 + 3*i];
                rGradU(1,0) += rDN_DX(i,1) * rNodalVar[  + 3*i] / rNodalVar[2 + 3*i];
                rGradU(1,1) += rDN_DX(i,1) * rNodalVar[1 + 3*i] / rNodalVar[2 + 3*i];
            }
        }
        if (near_dry)
        {
            rVel /= rHeight;
            rGradU /= rHeight;
        }
        
        rMom    *= lumping_factor;
        rHeight *= lumping_factor;
        rVel    *= lumping_factor;
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void EulerConsVarElement<TNumNodes>::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double, TNumNodes*3, TNumNodes*3>& rMassMatrix) 
    {
        const unsigned int element_size = 3*TNumNodes;
        rMassMatrix  = IdentityMatrix(element_size, element_size);
        rMassMatrix /= static_cast<double>(TNumNodes);
    }

//----------------------------------------------------------------------

template class EulerConsVarElement<3>;
template class EulerConsVarElement<4>;

} // namespace Kratos
