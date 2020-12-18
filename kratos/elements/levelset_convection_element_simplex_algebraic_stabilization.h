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
//                   Mohammad R. Hashemi
//


#if !defined(KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_ALGEBRAIC_STABILIZATION )
#define  KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_ALGEBRAIC_STABILIZATION

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "includes/convection_diffusion_settings.h"
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

/// Formulation is based on the works by Kuzmin et al. (especifically see Comput. Methods Appl. Mech. Engrg. 322 (2017) 23–41)
/// Dirichlet boundary condition for rUnknownVar at velocity inlets/outles is essential to be set for this solver since it is based on the flux
template< unsigned int TDim, unsigned int TNumNodes>
class LevelSetConvectionElementSimplexAlgebraicStabilization
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LevelSetConvectionElementSimplexAlgebraicStabilization);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LevelSetConvectionElementSimplexAlgebraicStabilization() : Element()
    {}

    LevelSetConvectionElementSimplexAlgebraicStabilization(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    LevelSetConvectionElementSimplexAlgebraicStabilization(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~LevelSetConvectionElementSimplexAlgebraicStabilization() override {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new LevelSetConvectionElementSimplexAlgebraicStabilization(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!

        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;
        //const double theta = rCurrentProcessInfo.Has(TIME_INTEGRATION_THETA) ? rCurrentProcessInfo[TIME_INTEGRATION_THETA] : 0.5;

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();

        //getting data for the given geometry
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
        //double h = ComputeH(DN_DX, Volume);

        //here we get all the variables we will need
        array_1d<double,TNumNodes> phi, phi_old;
        array_1d< array_1d<double,3 >, TNumNodes> v, vold;

        array_1d<double,3 > X_mean_tmp = ZeroVector(3);
        array_1d<double,3 > grad_phi_mean_tmp = ZeroVector(3);
        array_1d< array_1d<double,3 >, TNumNodes> X_node;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);

            v[i] = GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
            vold[i] = GetGeometry()[i].FastGetSolutionStepValue(rConvVar,1);

            X_mean_tmp += GetGeometry()[i].Coordinates();
            X_node[i] = GetGeometry()[i].Coordinates();

            grad_phi_mean_tmp += GetGeometry()[i].FastGetSolutionStepValue(DISTANCE_GRADIENT);
         }

        array_1d<double,TDim> X_mean, grad_phi_mean;
        for(unsigned int k = 0; k < TDim; k++)
        {
            grad_phi_mean[k] = grad_phi_mean_tmp[k]/static_cast<double>(TNumNodes);
            X_mean[k] = X_mean_tmp[k]/static_cast<double>(TNumNodes);
        }

        array_1d<double,TDim> grad_phi_diff = prod(trans(DN_DX), phi_old) - grad_phi_mean;

        BoundedMatrix<double,TNumNodes, TNumNodes> K_matrix = ZeroMatrix(TNumNodes, TNumNodes); // convection
        BoundedMatrix<double,TNumNodes, TNumNodes> S_matrix = ZeroMatrix(TNumNodes, TNumNodes); // LHS stabilization
        Vector S_vector = ZeroVector(TNumNodes); // RHS stabilization
        BoundedMatrix<double,TNumNodes, TNumNodes> Mc_matrix = ZeroMatrix(TNumNodes, TNumNodes); // consistent mass matrix
        BoundedMatrix<double,TNumNodes, TNumNodes> Ml_matrix = ZeroMatrix(TNumNodes, TNumNodes); // lumped mass matrix

        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);
        for(unsigned int igauss=0; igauss<TDim+1; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity/coordinate at the gauss point
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            array_1d<double, TDim > X_gauss=ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                for(unsigned int k=0; k<TDim; k++)
                {
                    vel_gauss[k] += N[i]*v[i][k]; //0.5*N[i]*(v[i][k]+vold[i][k]);
                    X_gauss[k] += N[i]*X_node[i][k];
                }
            }

            array_1d<double, TNumNodes > v_dot_grad_N = prod(DN_DX, vel_gauss);

            noalias(Mc_matrix) += outer_prod(N, N);
            noalias(K_matrix) += outer_prod(N, v_dot_grad_N);

            for (unsigned int i = 0; i < TNumNodes; i++){
                S_vector[i] += inner_prod(grad_phi_diff,(X_gauss - X_mean))*N[i];
            }
        }

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Ml_matrix(i,i) = 1.0;
        }

        noalias(S_matrix) = (1.0/(static_cast<double>(TDim)+1.0))*(Ml_matrix-Mc_matrix);

        double nu_e = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                if (i == j)
                    continue;

                nu_e = std::max( nu_e, std::max(0.0, K_matrix(i, j))/std::abs(S_matrix(i, j)) );
            }
        }

        const double limiter = GetValue(LIMITER_COEFFICIENT);
        noalias(rLeftHandSideMatrix)  = dt_inv*((1.0-limiter)*Ml_matrix + limiter*Mc_matrix) + K_matrix + (1.0-limiter)*nu_e*S_matrix;
        noalias(rRightHandSideVector) = dt_inv*prod((1.0-limiter)*Ml_matrix + limiter*Mc_matrix, phi_old) - limiter*nu_e*S_vector;

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Levelset Element")
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_THROW_ERROR(std::runtime_error, "CalculateRightHandSide not implemented","");
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override
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

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override
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
        return "LevelSetConvectionElementSimplexAlgebraicStabilization #";
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
            h = std::max(h, 1.0 / h_inv);
        }
        h = std::sqrt(h);
        return h;
    }

    //gauss points for the 3D case
    void GetShapeFunctionsOnGauss(BoundedMatrix<double,4,4>& Ncontainer)
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

#endif // KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_ALGEBRAIC_STABILIZATION  defined


