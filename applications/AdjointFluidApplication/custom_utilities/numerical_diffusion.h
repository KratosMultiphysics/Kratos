//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Suneth Warnakulasriya, https://github.com/sunethwarna
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
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Eigen library includes
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>


// Application includes
#include "../custom_elements/vms_adjoint_element.h"

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

    enum class numericalDiffusionMethod
    {
        singularValuePressureCoupled,
        singularValuePressureDecoupled,
        eigenValueFullMatrix
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NumericalDiffusion()
    {
    }

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
    void SetNumericalDiffusionParamters(Parameters& rParameters)
    {
        KRATOS_TRY;

        Parameters default_params(R"(
        {
            "method": "singular_value_pressure_coupled",
            "beta": 0.0
        })");

        rParameters.ValidateAndAssignDefaults(default_params);

        mBeta = rParameters["beta"].GetDouble();

        std::string method_name = rParameters["method"].GetString();
        
        if (method_name.compare("singular_value_pressure_coupled")==0)
            mNumericalDiffusionMethod = numericalDiffusionMethod::singularValuePressureCoupled;
        else if (method_name.compare("singular_value_pressure_decoupled")==0)
            mNumericalDiffusionMethod = numericalDiffusionMethod::singularValuePressureDecoupled;
        else if (method_name.compare("eigen_value_full_matrix")==0)
            mNumericalDiffusionMethod = numericalDiffusionMethod::eigenValueFullMatrix;
        else
            KRATOS_THROW_ERROR(std::runtime_error,
                "numerical diffusion method only supports singular_value_pressure_coupled or singular_value_pressure_decoupled or eigen_values",
                rParameters.PrettyPrintJsonString())
            
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
                case numericalDiffusionMethod::singularValuePressureCoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureCoupled<2>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case numericalDiffusionMethod::singularValuePressureDecoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureDecoupled<2>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case numericalDiffusionMethod::eigenValueFullMatrix:
                    numerical_diffusion = CalculateNumericalDiffusionEigenFullMatrix<2>(pCurrentElement, rCurrentProcessInfo);
                    break;
            }

            boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
            array_1d< double, 3 > N;
            double Volume;
    
            GeometryUtils::CalculateGeometryData(pCurrentElement->GetGeometry(),DN_DX,N,Volume);

            numerical_diffusion = mBeta * numerical_diffusion;

            AddNumericalDiffusionTerm<2>(rAdjointMatrix, DN_DX, numerical_diffusion);            
        }
        else if (domain_size == 3)
        {
            InitializeMatrix<3>(rAdjointMatrix);

            switch (mNumericalDiffusionMethod)
            {
                case numericalDiffusionMethod::singularValuePressureCoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureCoupled<3>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case numericalDiffusionMethod::singularValuePressureDecoupled:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureDecoupled<3>(pCurrentElement, rCurrentProcessInfo);
                    break;
                case numericalDiffusionMethod::eigenValueFullMatrix:
                    numerical_diffusion = CalculateNumericalDiffusionEigenFullMatrix<3>(pCurrentElement, rCurrentProcessInfo);
                    break;
            }

            boost::numeric::ublas::bounded_matrix<double, 4, 3> DN_DX;
            array_1d< double, 4 > N;
            double Volume;
    
            GeometryUtils::CalculateGeometryData(pCurrentElement->GetGeometry(),DN_DX,N,Volume);

            numerical_diffusion = mBeta * numerical_diffusion;

            AddNumericalDiffusionTerm<3>(rAdjointMatrix, DN_DX, numerical_diffusion);
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                "numerical diffusion method only supports 2D or 3D elements",
                "")            
        }        

        pCurrentElement->SetValue(NUMERICAL_DIFFUSION, numerical_diffusion);

        KRATOS_CATCH("");
    }
    ///@}

private:
    ///@name Member Variables
    ///@{
    numericalDiffusionMethod mNumericalDiffusionMethod = numericalDiffusionMethod::singularValuePressureCoupled;
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

        for (IndexType i=0; i < TFluidLocalSize; i++)
            for (IndexType j=0; j < TFluidLocalSize; j++)
                rAdjointMatrix(i,j) = 0.0;

        KRATOS_CATCH("")
    }

    template<unsigned int TDim>
    void CalculateVelocityGradientTensor(
        Element::Pointer pElement,
        const ProcessInfo& rCurrentProcessInfo,
        MatrixType& rMatrix
    )
    {
        constexpr unsigned int TNumNodes = TDim + 1;        

        if (rMatrix.size1() != TDim || rMatrix.size2() != TDim)
            rMatrix.resize(TDim,TDim,false);
        
        rMatrix.clear();

        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),DN_DX,N,Volume);

        std::vector<array_1d<double, TDim>> rNodalVelocityVectors;
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
            rNodalVelocityVectors.push_back(pElement->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY));

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                for (unsigned int k = 0; k < TNumNodes; k++)
                    rMatrix(i,j) += DN_DX(k,i)*rNodalVelocityVectors[k][j];

    }

    template<unsigned int TDim>
    double CalculateNumericalDiffusionSVMethodPressureCoupled(Element::Pointer pCurrentElement, const ProcessInfo& rCurrentProcessInfo)
    {
        constexpr unsigned int TNumNodes = TDim + 1;        

        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;

        GeometryType& rGeom = pCurrentElement->GetGeometry();

        GeometryUtils::CalculateGeometryData(rGeom,DN_DX,N,Volume);

        MatrixType GradVel;

        CalculateVelocityGradientTensor<TDim>(pCurrentElement, rCurrentProcessInfo, GradVel);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<TDim; i++)
            velocity_divergence += GradVel(i,i);
        
        Eigen::Matrix<double, TDim+1, TDim+1>  M;
        for (IndexType i=0; i < TDim; i++)
            M(i,i) = 0.5 * velocity_divergence - GradVel(i,i);
        for (IndexType i=0; i < TDim; i++)
            for (IndexType j=i+1; j < TDim; j++)
            {
                M(i,j) = -GradVel(i,j);
                M(j,i) = -GradVel(j,i);
            }
        for (IndexType i=0; i < TDim + 1; i++)
        {
            M(TDim, i) = 0.0;
            M(i, TDim) = 0.0;
        }
        M(TDim, TDim) = 0.5 * velocity_divergence;

        Eigen::JacobiSVD<Eigen::Matrix<double, TDim+1, TDim+1>> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);

        const auto& S = svd.singularValues();

        double numerical_viscosity;
        numerical_viscosity = S[0];

        return numerical_viscosity;
    }

    template<unsigned int TDim>
    double CalculateNumericalDiffusionSVMethodPressureDecoupled(Element::Pointer pCurrentElement, const ProcessInfo& rCurrentProcessInfo)
    {
        constexpr unsigned int TNumNodes = TDim + 1;

        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;

        GeometryType& rGeom = pCurrentElement->GetGeometry();

        GeometryUtils::CalculateGeometryData(rGeom,DN_DX,N,Volume);

        MatrixType GradVel;
        CalculateVelocityGradientTensor<TDim>(pCurrentElement, rCurrentProcessInfo, GradVel);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<TDim; i++)
            velocity_divergence += GradVel(i,i);
        
        Eigen::Matrix<double, TDim+1, TDim+1>  M;
        for (IndexType i=0; i < TDim; i++)
            M(i,i) = 0.5 * velocity_divergence - GradVel(i,i);
        for (IndexType i=0; i < TDim; i++)
            for (IndexType j=i+1; j < TDim; j++)
            {
                M(i,j) = -GradVel(i,j);
                M(j,i) = -GradVel(j,i);
            }

        Eigen::JacobiSVD<Eigen::Matrix<double, TDim+1, TDim+1>> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);

        const auto& S = svd.singularValues();

        double numerical_viscosity;
        numerical_viscosity = S[0];
        
        return numerical_viscosity;    
    }

    template<IndexType TDim>
    void AddNumericalDiffusionTerm(
        MatrixType& rResult,
        const boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim>& rDN_DX,
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

    template<unsigned int TDim>
    double CalculateNumericalDiffusionEigenFullMatrix(Element::Pointer pCurrentElement, const ProcessInfo& rCurrentProcessInfo)
    {
        // constexpr unsigned int TNumNodes = TDim + 1;

        // boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        // array_1d< double, TNumNodes > N;
        // double Volume;

        // GeometryType& rGeom = pCurrentElement->GetGeometry();

        // GeometryUtils::CalculateGeometryData(rGeom,DN_DX,N,Volume);

        // MatrixType vms_steady_term_primal_gradient;

        // pCurrentElement->Calculate(
        //                     VMS_STEADY_TERM_PRIMAL_GRADIENT_MATRIX,
        //                     vms_steady_term_primal_gradient,
        //                     rCurrentProcessInfo);

        // boost::numeric::ublas::vector<double> adjoint_values_vector;
        // boost::numeric::ublas::vector<double> temp;

        // adjoint_values_vector.resize(TNumNodes*(TDim+1));
        // temp.resize(TNumNodes*(TDim+1));

        // pCurrentElement->GetValuesVector(adjoint_values_vector, 1);

        // noalias(temp) = prod(vms_steady_term_primal_gradient, adjoint_values_vector);
        
        // double adjoint_energy = 0.0;
        // for (IndexType i = 0; i < temp.size(); i++)
        //     adjoint_energy += temp[i]*adjoint_values_vector[i];
        

        // MatrixType numerical_diffusion_matrix;
        // InitializeMatrix<TDim>(numerical_diffusion_matrix);

        // AddNumericalDiffusionTerm<TDim>(numerical_diffusion_matrix, DN_DX, Volume);

        // double diffusion_energy = 0.0;
        // noalias(temp) = prod(numerical_diffusion_matrix, adjoint_values_vector);
        // for (IndexType i = 0; i < temp.size(); i++)
        //     diffusion_energy += temp[i]*adjoint_values_vector[i];

        double numerical_diffusion = 0.0;

        // std::cout<<adjoint_energy<<", "<<diffusion_energy<<std::endl;

        // if (adjoint_energy < 0.0 and diffusion_energy != 0.0)
        //     numerical_diffusion = -adjoint_energy/diffusion_energy;

        return numerical_diffusion;
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
