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
//   Date:                July 3rd 2017 
//   Revision:            1.2 
// 
// 
 
// Project includes 
#include "includes/define.h" 
#include "custom_elements/conserved_var_element.hpp" 
#include "shallow_water_application.h" 
#include "utilities/math_utils.h" 
#include "utilities/geometry_utilities.h" 
 
namespace Kratos 
{ 
 
//---------------------------------------------------------------------- 
 
    template< unsigned int TNumNodes > 
    void ConservedVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) 
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
    void ConservedVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo) 
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
    void ConservedVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) 
    { 
        KRATOS_TRY 
         
        // Resize of the Left and Right Hand side 
        unsigned int element_size = TNumNodes*3; 
        if(rLeftHandSideMatrix.size1() != element_size) 
            rLeftHandSideMatrix.resize(element_size,element_size,false); // Resizing the system in case it does not have the right size 
         
        if(rRightHandSideVector.size() != element_size) 
            rRightHandSideVector.resize(element_size,false); 
 
        // Getting gravity 
        mGravity = rCurrentProcessInfo[GRAVITY_Z];
        
        // Getting water height unit converter
        mHeightUnitConvert = rCurrentProcessInfo[WATER_HEIGHT_UNIT_CONVERTER];
     
        // Getting the time step (not fixed to allow variable time step) 
        const double delta_t = rCurrentProcessInfo[DELTA_TIME]; 
        const double dt_inv = 1.0 / delta_t; 
        
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
        array_1d<double, TNumNodes*3> v_depth;
        array_1d<double, TNumNodes*3> v_rain;
        array_1d<double, TNumNodes*3> v_unknown;
        array_1d<double, TNumNodes*3> v_proj_unknown;
        GetNodalValues(v_depth, v_rain, v_unknown, v_proj_unknown);
        
        // Get element values (this function inlcudes the units conversion)
        double height;
        array_1d<double,2> height_grad;
        double div_u;
        GetElementValues(DN_DX, v_unknown, height, height_grad, div_u);
        
        // Compute stabilization and discontinuity capturing parameters
        double tau_m;
        double tau_h;
        double k_dc;
        ComputeStabilizationParameters(height, height_grad, elem_length, tau_m, tau_h, k_dc);
        
        // Some auxilary definitions
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> N_mom        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for momentum unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> N_height     = ZeroMatrix(1,TNumNodes*3);  // Shape functions vector (for height unknown)
        boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> DN_DX_mom    = ZeroMatrix(1,TNumNodes*3);  // Shape functions divergence vector (for momentum unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradient matrix (for height unknown)
        boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> grad_mom     = ZeroMatrix(2,TNumNodes*3);  // Shaoe functions gradient vector (for momentum unknown)
        //
        boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
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
                grad_mom(0,   nnode*3) = DN_DX(nnode,0);
                grad_mom(0, 1+nnode*3) = DN_DX(nnode,0);
                grad_mom(1,   nnode*3) = DN_DX(nnode,1);
                grad_mom(1, 1+nnode*3) = DN_DX(nnode,1);
                // Height shape funtions
                N_height(0, 2+nnode*3) = N[nnode];
                // Momentum shape functions
                N_mom(0,   nnode*3) = N[nnode];
                N_mom(1, 1+nnode*3) = N[nnode];
            }
            N_height     *= mHeightUnitConvert;
            DN_DX_height *= mHeightUnitConvert;
            
            noalias(mass_matrix)  += prod(trans(N_mom),N_mom);
            noalias(mass_matrix)  += prod(trans(N_height),N_height);
         
            noalias(aux_q_div_m)  += prod(trans(N_height),DN_DX_mom);
            noalias(aux_w_grad_h) += prod(trans(N_mom),DN_DX_height);
        
            noalias(aux_m_diffus) += prod(trans(DN_DX_mom),DN_DX_mom);
            noalias(aux_h_diffus) += prod(trans(DN_DX_height),DN_DX_height);
        }
        
        //~ CalculateLumpedMassMatrix(mass_matrix);
        
        
        // Build LHS 
        // Cross terms 
        noalias(rLeftHandSideMatrix)  = mGravity * height * aux_w_grad_h; // Add <w,g*h*grad(h)> to Momentum Eq. 
        
        // Inertia terms 
        noalias(rLeftHandSideMatrix) += dt_inv * mass_matrix;           // Add <N,N> to both Eq's 
        
        // Non linear terms 
        noalias(rLeftHandSideMatrix) += div_u * mass_matrix;             // Add <q,div(u)*h> to Mass Eq. and <w,div(u)*hu> to Momentum Eq. 
        
        // Stabilization terms 
        noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;  // Add art. diff. to Mass Eq. 
        noalias(rLeftHandSideMatrix) +=         tau_m  * aux_m_diffus;  // Add art. diff. to Momentum Eq. 
        
        // Build RHS 
        // Source terms (bathymetry contribution) 
        noalias(rRightHandSideVector)  = -mGravity * prod(aux_w_grad_h, v_depth); // Add <w,-g*h*grad(H)> to RHS (Momentum Eq.) 
        
        // Source terms (rain contribution)
        noalias(rRightHandSideVector) += prod(mass_matrix, v_rain);
        
        // Inertia terms 
        noalias(rRightHandSideVector) += dt_inv * prod(mass_matrix, v_proj_unknown); 
     
        // Subtracting the dirichlet term (since we use a residualbased approach) 
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, v_unknown); 
     
        rRightHandSideVector *= Area / static_cast<double>(TNumNodes); 
        rLeftHandSideMatrix *= Area  / static_cast<double>(TNumNodes); 
     
        KRATOS_CATCH(""); 
    } 
 
//---------------------------------------------------------------------- 
 
    template< unsigned int TNumNodes > 
    void ConservedVarElement<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) 
    { 
        KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , ""); 
    } 
 
//---------------------------------------------------------------------- 
 
    template< unsigned int TNumNodes > 
    void ConservedVarElement<TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) 
    { 
        if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM || rVariable == MIU) 
        { 
            for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
                rValues[PointNumber] = double(this->GetValue(rVariable)); 
        } 
    }
 
//---------------------------------------------------------------------- 
 
    template< unsigned int TNumNodes > 
    void ConservedVarElement<TNumNodes>::CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea) 
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
    double ConservedVarElement<TNumNodes>::ComputeElemSize(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX) 
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
    void ConservedVarElement<TNumNodes>::GetNodalValues(array_1d<double, TNumNodes*3>& rDepth,
                                                        array_1d<double, TNumNodes*3>& rRain,
                                                        array_1d<double, TNumNodes*3>& rUnkn,
                                                        array_1d<double, TNumNodes*3>& rProj
                                                       ) 
    { 
        GeometryType& rGeom = GetGeometry(); 
        unsigned int counter = 0; 
        for (unsigned int i = 0; i < TNumNodes; i++) 
        { 
            rDepth[counter] = 0; 
            rRain[counter]  = 0; 
            rUnkn[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X); 
            rProj[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_X); 
            counter++; 
 
            rDepth[counter] = 0; 
            rRain[counter]  = 0; 
            rUnkn[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y); 
            rProj[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Y); 
            counter++; 
 
            rDepth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY) / mHeightUnitConvert;
            rRain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN); 
            rUnkn[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT); 
            rProj[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_SCALAR1); 
            counter++; 
        } 
    } 
 
//---------------------------------------------------------------------- 
 
    template< unsigned int TNumNodes > 
    void ConservedVarElement<TNumNodes>::GetElementValues(const boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2>& rDN_DX, 
                                                          const array_1d<double, TNumNodes*3>& rNodalVar, 
                                                          double& rHeight,
                                                          array_1d<double,2>& rHeightGrad,
                                                          double& rDivU 
                                                          )
    {
        double lumping_factor = 1 / double(TNumNodes);
        bool near_dry = false;
        
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            if (rNodalVar[2 + 3*i] < 1e-1)
                near_dry = true;
        }
        
        // Initialize outputs
        rHeight = 0;
        rHeightGrad = ZeroVector(2);
        rDivU  = 0;
        
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rHeight += rNodalVar[2 + 3*i]; 
            rHeightGrad[0] += rDN_DX(i,0) * rNodalVar[2 + 3*i];
            rHeightGrad[1] += rDN_DX(i,1) * rNodalVar[2 + 3*i];
            if (near_dry)
            {
                rDivU  += rDN_DX(i,0) * rNodalVar[  + 3*i];
                rDivU  += rDN_DX(i,1) * rNodalVar[1 + 3*i];
            }
            else
            {
                rDivU  += rDN_DX(i,0) * rNodalVar[  + 3*i] / rNodalVar[2 + 3*i];
                rDivU  += rDN_DX(i,1) * rNodalVar[1 + 3*i] / rNodalVar[2 + 3*i];
            }
        }
        if (near_dry)
            rDivU /= rHeight;
        
        rHeight *= lumping_factor * mHeightUnitConvert;
        rHeightGrad *= mHeightUnitConvert;
        rDivU /= (mHeightUnitConvert);
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void ConservedVarElement<TNumNodes>::ComputeStabilizationParameters(const double& rHeight,
                                                                        const array_1d<double,2>& rHeightGrad,
                                                                        const double& rElemSize,
                                                                        double& rTauU,
                                                                        double& rTauH,
                                                                        double& rKdc
                                                                        )
    {
        // Initialize outputs
        rTauU = 0;
        rTauH = 0;
        rKdc  = 0;
        
        // Get element values
        double height_grad_norm = norm_2(rHeightGrad);
        
        // Compute stabilization parameters
        bool stabilization = true;
        double Ctau = 0.01;       // Stabilization parameter >0.005 (R.Codina, CMAME 197, 2008, 1305-1322) 
        double fheight = fabs(rHeight);
        if (stabilization && fheight > 1e-6)
        {
            rTauU = Ctau/rElemSize*pow(mGravity/fheight,0.5);
            rTauH = Ctau/rElemSize*pow(fheight/mGravity,0.5);
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
    void ConservedVarElement<TNumNodes>::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double, TNumNodes*3, TNumNodes*3>& rMassMatrix) 
    {
        const unsigned int element_size = 3*TNumNodes;
        rMassMatrix  = IdentityMatrix(element_size, element_size);
        rMassMatrix /= static_cast<double>(TNumNodes);
    }
 
//---------------------------------------------------------------------- 
 
template class ConservedVarElement<3>; 
template class ConservedVarElement<4>; 
 
} // namespace Kratos 
