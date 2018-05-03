//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: KratosAdjointFluidApplication/license.txt
//
//  Main author:    Suneth Warnakulasriya, https://github.com/sunethwarna
//


#if !defined(KRATOS_NUMERICAL_DIFFUSION)
#define KRATOS_NUMERICAL_DIFFUSION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Eigen library includes
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>


// Application includes
#include "custom_elements/vms_adjoint_element.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A Numerical Diffusion Calculation class.
class NumericalDiffusion
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(NumericalDiffusion);

    typedef Element::IndexType IndexType;

    typedef Element::SizeType SizeType;

    typedef Element::MatrixType MatrixType;
    
    typedef Element::GeometryType GeometryType;

    typedef Element::VectorType VectorType;

    enum class NumericalDiffusionMethod
    {
        singularValuePressureCoupled,
        singularValuePressureDecoupled,
        eigenValueFullMatrix
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    ~NumericalDiffusion()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    void SetNumericalDiffusionParameters(Parameters DiffusionParameters)
    {
        KRATOS_TRY;

        Parameters default_params(R"(
        {
            "method": "singular_value_pressure_coupled",
            "beta": 0.0
        })");

        DiffusionParameters.ValidateAndAssignDefaults(default_params);

        mBeta = DiffusionParameters["beta"].GetDouble();

        std::string method_name = DiffusionParameters["method"].GetString();
        
        if (method_name.compare("singular_value_pressure_coupled")==0)
            mNumericalDiffusionMethod = NumericalDiffusionMethod::singularValuePressureCoupled;
        else if (method_name.compare("singular_value_pressure_decoupled")==0)
            mNumericalDiffusionMethod = NumericalDiffusionMethod::singularValuePressureDecoupled;
        else if (method_name.compare("eigen_value_full_element_matrix")==0)
            mNumericalDiffusionMethod = NumericalDiffusionMethod::eigenValueFullMatrix;
        else
            KRATOS_ERROR<<"numerical diffusion method only supports singular_value_pressure_coupled or singular_value_pressure_decoupled or eigen_value_full_element_matrix"<<DiffusionParameters.PrettyPrintJsonString();
            
        KRATOS_CATCH("");
        
    }

    void CalculateNumericalDiffusion(
        Element::Pointer pCurrentElement,
        MatrixType& rAdjointMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const unsigned int domain_size = static_cast<unsigned int>(rCurrentProcessInfo[DOMAIN_SIZE]);

        double numerical_diffusion = 0.0;

        if (domain_size == 2)
        {
            InitializeMatrix<2>(rAdjointMatrix);
            switch (mNumericalDiffusionMethod)
            {
                case NumericalDiffusionMethod::singularValuePressureCoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureCoupled<2>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case NumericalDiffusionMethod::singularValuePressureDecoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureDecoupled<2>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case NumericalDiffusionMethod::eigenValueFullMatrix:
                    numerical_diffusion = CalculateNumericalDiffusionEigenFullMatrix<2>(pCurrentElement, rCurrentProcessInfo);
                    break;                    
            }

            BoundedMatrix<double, 3, 2> dn_dx;
            array_1d< double, 3 > n;
            double volume;
    
            GeometryUtils::CalculateGeometryData(pCurrentElement->GetGeometry(),dn_dx,n,volume);

            numerical_diffusion = mBeta * numerical_diffusion * volume;

            AddNumericalDiffusionTerm<2>(rAdjointMatrix, dn_dx, numerical_diffusion);            
        }
        else if (domain_size == 3)
        {
            InitializeMatrix<3>(rAdjointMatrix);

            switch (mNumericalDiffusionMethod)
            {
                case NumericalDiffusionMethod::singularValuePressureCoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureCoupled<3>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case NumericalDiffusionMethod::singularValuePressureDecoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureDecoupled<3>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case NumericalDiffusionMethod::eigenValueFullMatrix:
                    numerical_diffusion = CalculateNumericalDiffusionEigenFullMatrix<3>(pCurrentElement, rCurrentProcessInfo);
                    break;                                        
            }

            BoundedMatrix<double, 4, 3> dn_dx;
            array_1d< double, 4 > n;
            double volume;
    
            GeometryUtils::CalculateGeometryData(pCurrentElement->GetGeometry(),dn_dx,n,volume);

            numerical_diffusion = mBeta * numerical_diffusion * volume;

            AddNumericalDiffusionTerm<3>(rAdjointMatrix, dn_dx, numerical_diffusion);
        }
        else
        {
            KRATOS_ERROR<<"numerical diffusion method only supports 2D or 3D elements";
        }        

        pCurrentElement->SetValue(NUMERICAL_DIFFUSION, numerical_diffusion);

        KRATOS_CATCH("");
    }
    ///@}

private:
    ///@name Member Variables
    ///@{
    NumericalDiffusionMethod mNumericalDiffusionMethod = NumericalDiffusionMethod::singularValuePressureCoupled;
    double mBeta = 0.0;
    ///@}
    ///@name Private Operators
    ///@{
    template<unsigned int TDim>
    void InitializeMatrix(MatrixType& rAdjointMatrix)
    {
        KRATOS_TRY;

        constexpr unsigned int TNumNodes = TDim + 1;        
        constexpr unsigned int TBlockSize = TDim + 1;    
        constexpr unsigned int TFluidLocalSize = TBlockSize * TNumNodes;    
        
        if (rAdjointMatrix.size1() != TFluidLocalSize || rAdjointMatrix.size2() != TFluidLocalSize)
            rAdjointMatrix.resize(TFluidLocalSize,TFluidLocalSize,false);

        rAdjointMatrix.clear();

        KRATOS_CATCH("")
    }

    template<unsigned int TDim>
    void CalculateVelocityGradientTensor(
        Element::Pointer const pElement,
        MatrixType& rMatrix
    )
    {
        constexpr unsigned int TNumNodes = TDim + 1;        

        if (rMatrix.size1() != TDim || rMatrix.size2() != TDim)
            rMatrix.resize(TDim,TDim,false);
        
        rMatrix.clear();

        BoundedMatrix<double, TNumNodes, TDim> dn_dx;
        array_1d< double, TNumNodes > n;
        double volume;
        GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,n,volume);

        BoundedVector<array_1d<double, TDim>, TNumNodes> r_nodal_velocity_vectors;

        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
            r_nodal_velocity_vectors[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                for (unsigned int k = 0; k < TNumNodes; k++)
                    rMatrix(i,j) += dn_dx(k,i)*r_nodal_velocity_vectors[k][j];

    }

    
    template<unsigned int TSize, unsigned int TDim>
    Eigen::Matrix<double, TSize, TSize>  CalculateSVMethodCharacteristicMatrix(Element::Pointer const pElement)
    {
        MatrixType velocity_gradient;

        CalculateVelocityGradientTensor<TDim>(pElement, velocity_gradient);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<TDim; i++)
            velocity_divergence += velocity_gradient(i,i);
        
        Eigen::Matrix<double, TSize, TSize>  characteristic_matrix;
        for (IndexType i=0; i < TDim; i++)
            characteristic_matrix(i,i) = 0.5 * velocity_divergence - velocity_gradient(i,i);
        for (IndexType i=0; i < TDim; i++)
            for (IndexType j=i+1; j < TDim; j++)
            {
                characteristic_matrix(i,j) =  velocity_gradient(i,j);
                characteristic_matrix(j,i) =  velocity_gradient(j,i);
            }
        
        return characteristic_matrix;
    }

    template<unsigned int TDim>
    double CalculateVelocityDivergence(Element::Pointer const pElement)
    {
        MatrixType velocity_gradient;

        CalculateVelocityGradientTensor<TDim>(pElement, velocity_gradient);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<TDim; i++)
            velocity_divergence += velocity_gradient(i,i);
        
        return velocity_divergence;
    }

    template<unsigned int TDim>
    double CalculateNumericalDiffusionSVMethodPressureCoupled(Element::Pointer const pElement, const ProcessInfo& rCurrentProcessInfo)
    {

        Eigen::Matrix<double, TDim+1, TDim+1>  characteristic_matrix = 
                            CalculateSVMethodCharacteristicMatrix<TDim+1, TDim>(pElement);

        for (IndexType i=0; i < TDim + 1; i++)
        {
            characteristic_matrix(TDim, i) = 0.0;
            characteristic_matrix(i, TDim) = 0.0;
        }
        
        double velocity_divergence = CalculateVelocityDivergence<TDim>(pElement);
        characteristic_matrix(TDim, TDim) = 0.5 * velocity_divergence;

        Eigen::JacobiSVD<Eigen::Matrix<double, TDim+1, TDim+1>> svd(
                characteristic_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV
                );

        const auto& S = svd.singularValues();

        double numerical_viscosity;
        numerical_viscosity = S[0];

        return numerical_viscosity;
    }

    template<unsigned int TDim>
    double CalculateNumericalDiffusionSVMethodPressureDecoupled(Element::Pointer const pElement, const ProcessInfo& rCurrentProcessInfo)
    {
        Eigen::Matrix<double, TDim, TDim>  characteristic_matrix = 
                CalculateSVMethodCharacteristicMatrix<TDim, TDim>(pElement);

        Eigen::JacobiSVD<Eigen::Matrix<double, TDim, TDim>> svd(
            characteristic_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV
            );

        const auto& S = svd.singularValues();

        double numerical_viscosity;
        numerical_viscosity = S[0];
        
        return numerical_viscosity;    
    }

    template<unsigned int TDim>
    double CalculateNumericalDiffusionEigenFullMatrix(Element::Pointer const pElement, const ProcessInfo& rCurrentProcessInfo)
    {
        constexpr unsigned int TNumNodes = TDim + 1;

        BoundedMatrix<double, TNumNodes, TDim> dn_dx;
        array_1d< double, TNumNodes > n;
        double volume;

        GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,n,volume);

        MatrixType vms_steady_term_primal_gradient;

        pElement->Calculate(
                            VMS_STEADY_TERM_PRIMAL_GRADIENT_MATRIX,
                            vms_steady_term_primal_gradient,
                            rCurrentProcessInfo);

        Vector adjoint_values_vector;
        Vector temp;

        pElement->GetValuesVector(adjoint_values_vector, 1);

        noalias(temp) = prod(vms_steady_term_primal_gradient, adjoint_values_vector);
        double adjoint_energy = inner_prod(temp, adjoint_values_vector);
        
        MatrixType numerical_diffusion_matrix;
        InitializeMatrix<TDim>(numerical_diffusion_matrix);

        AddNumericalDiffusionTerm<TDim>(numerical_diffusion_matrix, dn_dx, volume);

        noalias(temp) = prod(numerical_diffusion_matrix, adjoint_values_vector);
        double diffusion_energy = inner_prod(temp,adjoint_values_vector);

        double numerical_diffusion = 0.0;

        if (adjoint_energy > 0.0 and diffusion_energy > 0.0)
            numerical_diffusion = adjoint_energy/diffusion_energy;

        return numerical_diffusion;
    }    

    template<IndexType TDim>
    void AddNumericalDiffusionTerm(
        MatrixType& rResult,
        const BoundedMatrix<double, TDim + 1, TDim>& rDN_DX,
        const double Weight)
    {
        constexpr SizeType NumNodes = TDim + 1;
        constexpr SizeType TBlockSize = TDim + 1;

        for (IndexType a = 0; a < NumNodes; ++a)
        {
            for (IndexType b = 0; b < NumNodes; ++b)
            {
                // (dN_a/dx_k dN_b/dx_k)
                double value = 0.0;
                for (IndexType k = 0; k < TDim; k++)
                    value += rDN_DX(a,k) * rDN_DX(b,k);
                
                value *= Weight;

                for (IndexType i=0; i < TBlockSize; i++)
                    rResult(a*TBlockSize+i,b*TBlockSize+i) += value;
            }
        }
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_NUMERICAL_DIFFUSION defined */
