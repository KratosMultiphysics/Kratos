//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//


#if !defined(KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_INCLUDED )
#define  KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
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
class LevelSetConvectionElementSimplex
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
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
    ~LevelSetConvectionElementSimplex() override {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new LevelSetConvectionElementSimplex(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
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

        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();
        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];


        //getting data for the given geometry
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
        double h = ComputeH(DN_DX, Volume);


        //here we get all the variables we will need
        array_1d<double,TNumNodes> phi, phi_old;
        array_1d< array_1d<double,3 >, TNumNodes> v, vold;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
//             dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];

            v[i] = GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
            vold[i] = GetGeometry()[i].FastGetSolutionStepValue(rConvVar,1);
        }
        array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
        const double norm_grad = norm_2(grad_phi_halfstep);

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

        BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
        BoundedMatrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
        BoundedMatrix<double,TNumNodes, TDim> tmp;


        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
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

            const double tau_denom = std::max(dyn_st_beta *dt_inv + 2.0 * norm_vel / h + std::abs(/*beta**/div_v),  1e-2); //the term std::abs(div_v) is added following Pablo Becker's suggestion
            const double tau = 1.0 / (tau_denom);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*beta*div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);

            //terms which multiply the gradient of phi
            noalias(aux2) += (1.0+tau*beta*div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);

            //cross-wind term
            if(norm_grad > 1e-3 && norm_vel > 1e-9)
            {
                const double C = rCurrentProcessInfo.GetValue(CROSS_WIND_STABILIZATION_FACTOR);
                const double time_derivative = dt_inv*(inner_prod(N,phi)-inner_prod(N,phi_old));
                const double res = -time_derivative -inner_prod(vel_gauss, grad_phi_halfstep);

                const double disc_capturing_coeff = 0.5*C*h*fabs(res/norm_grad);
                BoundedMatrix<double,TDim,TDim> D = disc_capturing_coeff*( IdentityMatrix(TDim));
                const double norm_vel_squared = norm_vel*norm_vel;
                D += (std::max( disc_capturing_coeff - tau*norm_vel_squared , 0.0) - disc_capturing_coeff)/(norm_vel_squared) * outer_prod(vel_gauss,vel_gauss);

                noalias(tmp) = prod(DN_DX,D);
                noalias(aux2) += prod(tmp,trans(DN_DX));
            }
        }

        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (dt_inv + 0.5*beta*div_v)*aux1; //the 0.5 comes from the use of Crank Nichlson
        noalias(rRightHandSideVector) = (dt_inv - 0.5*beta*div_v)*prod(aux1,phi_old); //the 0.5 comes from the use of Crank Nichlson

        //terms in aux2
        noalias(rLeftHandSideMatrix) += 0.5*aux2; //the 0.5 comes from the use of Crank Nichlson
        noalias(rRightHandSideVector) -= 0.5*prod(aux2,phi_old); //the 0.5 comes from the use of Crank Nichlson

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Levelset Element")
    }





    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_THROW_ERROR(std::runtime_error, "CalculateRightHandSide not implemented","");
    }





    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
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





    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override
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

    std::string Info() const override
    {
        return "LevelSetConvectionElementSimplex #";
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
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
    double ComputeH(BoundedMatrix<double,TNumNodes, TDim>& DN_DX, const double Volume)
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
    void GetShapeFunctionsOnGauss(BoundedMatrix<double,4, 4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
        Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
        Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
        Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    }

    //gauss points for the 2D case
    void GetShapeFunctionsOnGauss(BoundedMatrix<double,3,3>& Ncontainer)
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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
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


