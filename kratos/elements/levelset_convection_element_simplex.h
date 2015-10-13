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




#if !defined(KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_INCLUDED )
#define  KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_INCLUDED


// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"



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

template< unsigned int TDim, unsigned int TNumNodes>
class LevelSetConvectionElementSimplex
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Fluid2DASGS
    KRATOS_CLASS_POINTER_DEFINITION(LevelSetConvectionElementSimplex);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    LevelSetConvectionElementSimplex(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    LevelSetConvectionElementSimplex(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~LevelSetConvectionElementSimplex() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new LevelSetConvectionElementSimplex(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const double lumping_factor = 1.00 / static_cast<double>(TNumNodes);

        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!


        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

        //Crank-Nicholson factor
        const double cr_nk = 0.5;

        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;


        //getting data for the given geometry
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
        double h = ComputeH(DN_DX, Volume);

        //INERTIA CONTRIBUTION - does not depend on gauss point
        boost::numeric::ublas::bounded_matrix<double, TNumNodes,TNumNodes > msMassFactors = lumping_factor* IdentityMatrix(TNumNodes,TNumNodes);
        noalias(rLeftHandSideMatrix) = dt_inv * msMassFactors;


        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();

        array_1d<double,TNumNodes> phi, phi_old;
        array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);

            const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
            const array_1d<double,3>& vold = GetGeometry()[i].FastGetSolutionStepValue(rConvVar,1);
            for(unsigned int k=0; k<TDim; k++)
                vel_gauss[k] += 0.5*N[i]*(v[k]+vold[k]);
        }
        const double norm_vel = norm_2(vel_gauss);

        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        const double tau = 1.0 / (dyn_st_beta *dt_inv + 2.0 * norm_vel / h);

        //Advective term
        array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);
        boost::numeric::ublas::bounded_matrix<double, TNumNodes,TNumNodes > Advective_Matrix = outer_prod(N, a_dot_grad);
        noalias(rLeftHandSideMatrix) += (1.0 - cr_nk) * Advective_Matrix;

        //stabilization terms
        array_1d<double, TNumNodes > a_dot_grad_and_mass;
        a_dot_grad_and_mass = dt_inv * N + (1.0 - cr_nk) * a_dot_grad;
        noalias(rLeftHandSideMatrix) += tau * outer_prod(a_dot_grad, a_dot_grad_and_mass);

        //compute shock capturing term
        array_1d<double,TDim> grad_g = prod(trans(DN_DX),phi);

        double res = inner_prod(vel_gauss,grad_g);

        double dphi_dt = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
            dphi_dt += N[i]*(phi[i] - phi_old[i]);
        res += dt_inv * dphi_dt;

        //Add all n_step terms
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes > old_step_matrix = dt_inv*msMassFactors;
        old_step_matrix -= (cr_nk * Advective_Matrix);
        noalias(rRightHandSideVector) = prod(old_step_matrix, phi_old);

        //Add n_Stabilization terms
        a_dot_grad_and_mass = dt_inv * N - cr_nk * a_dot_grad;
        double old_res = inner_prod(a_dot_grad_and_mass, phi_old);
        /*	old_res += heat_source;*/
        noalias(rRightHandSideVector) += tau * a_dot_grad * old_res;

        //subtracting the dirichlet term
        // RHS -= LHS*temperatures
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);


        rRightHandSideVector *= Volume;
        rLeftHandSideMatrix *= Volume;

        KRATOS_CATCH("Error in Levelset Element")
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

    LevelSetConvectionElementSimplex() : Element()
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

#endif // KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_INCLUDED  defined 


