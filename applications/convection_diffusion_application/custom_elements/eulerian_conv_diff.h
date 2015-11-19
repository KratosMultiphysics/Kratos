/*
==============================================================================
KratosIncompressibleFluidApplication
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




#if !defined(KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED )
#define  KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED


// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/serializer.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 



namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

///formulation described in https://docs.google.com/document/d/13a_zGLj6xORDuLgoOG5LwHI6BwShvfO166opZ815zLY/edit?usp=sharing
template< unsigned int TDim, unsigned int TNumNodes>
class EulerianConvectionDiffusionElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of 
    KRATOS_CLASS_POINTER_DEFINITION(EulerianConvectionDiffusionElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    EulerianConvectionDiffusionElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    EulerianConvectionDiffusionElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~EulerianConvectionDiffusionElement() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new EulerianConvectionDiffusionElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
       
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!


//         noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
//         noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

//         //Crank-Nicholson factor
//         const double cr_nk = 0.5;

        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;
		const double lumping_factor = 1.00 / double(TNumNodes);

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        //getting data for the given geometry
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
        double h = ComputeH(DN_DX, Volume);
        

        //here we get all the variables we will need
        array_1d<double,TNumNodes> phi, phi_old;
        array_1d< array_1d<double,3 >, TNumNodes> v, vold;
        //bool active_convection=false;    // to kill some terms in case active_convection=false. For the moment it is inactive.
        //using only one Gauss Point for the material properties and volumetric heat flux
        double conductivity = 0.0;
        double specific_heat = 0.0;
        double density = 0.0;

		//storing locally the flags to avoid repeated check in the nodal loops
        const bool IsDefinedVelocityVariable = my_settings->IsDefinedVelocityVariable();
        const bool IsDefinedMeshVelocityVariable = my_settings->IsDefinedMeshVelocityVariable();
        const bool IsDefinedDensityVariable = my_settings->IsDefinedDensityVariable();
        const bool IsDefinedSpecificHeatVariableVariable = my_settings->IsDefinedSpecificHeatVariable();
        const bool IsDefinedDiffusionVariable = my_settings->IsDefinedDiffusionVariable();


        for (unsigned int i = 0; i < TNumNodes; i++)
        {	
            phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
//             dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];

			v[i]=ZeroVector(3);
			vold[i]=ZeroVector(3);
            if (IsDefinedVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
				  v[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
				  vold[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar,1);
				  //active_convection=true;
			}
			
			if (IsDefinedMeshVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
				  v[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
				  vold[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar,1);
				  //active_convection=true;
			}
			
			if (IsDefinedDensityVariable)
			{
				const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
				density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
			}
			else
				density += 1.0;
				
			if (IsDefinedSpecificHeatVariableVariable)
			{
				const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
				specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);
			}
			else
				specific_heat += 1.0;			
				
			if (IsDefinedDiffusionVariable)
			{
				const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
				conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
			}
			//if not, then the conductivity = 0   
        }
        array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
        //const double norm_grad = norm_2(grad_phi_halfstep);

        conductivity *= lumping_factor;
        density *= lumping_factor;
        specific_heat *= lumping_factor;
        //heat_flux *= lumping_factor;



        //here we use a term beta which takes into account a reaction term of the type "beta*div_v"
        //compute the divergence of v
        double div_v = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
            for(unsigned int k=0; k<TDim; k++)
                div_v += 0.5*DN_DX(i,k)*(v[i][k] + vold[i][k]);
            
        double beta = 0.0; //1.0;
        
//         unsigned int nneg=0;
//         for(unsigned int i=0; i<TNumNodes; i++) if(phi[i] < 0.0) nneg++;
//         if(nneg > 0) beta = 1.0; //beta = 0.1;

        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
        bounded_matrix<double,TNumNodes, TDim> tmp;

            
        bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);
        for(unsigned int igauss=0; igauss<TDim+1; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity in the middle of the tiem step
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                 for(unsigned int k=0; k<TDim; k++)
                    vel_gauss[k] += 0.5*N[i]*(v[i][k]+vold[i][k]);
            }
            const double norm_vel = norm_2(vel_gauss);
            array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

            const double tau_denom = std::max(dyn_st_beta *dt_inv + 4.0 * conductivity/(h*h) + 2.0 *density*specific_heat* norm_vel / h + beta*div_v,  1e-2);
            const double tau = 1.0 / (tau_denom);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*beta*div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);
            
            //terms which multiply the gradient of phi
            noalias(aux2) += (1.0+tau*beta*div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);
            
            //cross-wind term  
//             if(norm_grad > 1e-3 && norm_vel > 1e-9)
//             {
//                 const double C = 0.7;
//                 const double time_derivative = dt_inv*(inner_prod(N,phi)-inner_prod(N,phi_old));
//                 const double res = -time_derivative -inner_prod(vel_gauss, grad_phi_halfstep);
//                 
//                 const double disc_capturing_coeff = 0.5*C*h*fabs(res/norm_grad);
//                 bounded_matrix<double,TDim,TDim> D = disc_capturing_coeff*( IdentityMatrix(TDim,TDim));
//                 const double norm_vel_squared = norm_vel*norm_vel;
//                 D += (std::max( disc_capturing_coeff - tau*norm_vel_squared , 0.0) - disc_capturing_coeff)/(norm_vel_squared) * outer_prod(vel_gauss,vel_gauss);
// 
//                 noalias(tmp) = prod(DN_DX,D);
//                 noalias(aux2) += prod(tmp,trans(DN_DX));
//             }
        }
        
        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (dt_inv*density*specific_heat + 0.5*beta*div_v)*aux1; //the 0.5 comes from the use of Crank Nichlson
        noalias(rRightHandSideVector) = (dt_inv*density*specific_heat - 0.5*beta*div_v)*prod(aux1,phi_old); //the 0.5 comes from the use of Crank Nichlson
        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (conductivity *0.5 * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) -= prod((conductivity *0.5 * prod(DN_DX, trans(DN_DX))),phi_old)*static_cast<double>(TNumNodes) ;
        
        //terms in aux2
        noalias(rLeftHandSideMatrix) += 0.5*aux2; //the 0.5 comes from the use of Crank Nichlson
        noalias(rRightHandSideVector) -= 0.5*prod(aux2,phi_old); //the 0.5 comes from the use of Crank Nichlson
        
        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Eulerian ConvDiff Element")
    }


    
    
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "CalculateRightHandSide not implemented","");
    }
    
    
    
    

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
        }
        KRATOS_CATCH("")

    }
    
    
    
    

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

        }
        KRATOS_CATCH("");
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    virtual std::string Info() const
    {
        return "LevelSetConvectionElementSimplex #";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    EulerianConvectionDiffusionElement() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{
    double ComputeH(boost::numeric::ublas::bounded_matrix<double,TNumNodes, TDim>& DN_DX, const double Volume)
    {
        double h=0.0;
                 for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }
    
    //gauss points for the 3D case
    void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,4, 4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
        Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;	
        Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
        Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    }

    //gauss points for the 2D case
    void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,3,3>& Ncontainer)
    {
        const double one_sixt = 1.0/6.0;
        const double two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third; 
        Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt; 
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt; 
    }



    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}




private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;
    //         ASGS2D() : Element()
    //         {
    //         }

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}

    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{





    ///@}

}; 

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                                    Fluid2DASGS& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                                    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

} // namespace Kratos.

#endif // KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED  defined 


