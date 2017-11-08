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

// // Boost includes
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/lapack/lapack_names.h>
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

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

    enum class NUMERICAL_DIFFUSION_METHOD
    {
        SINGULAR_VALUE_PRESSURE_COUPLED,
        SINGULAR_VALUE_PRESSURE_DECOUPLED,
        EIGEN_VALUE
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
    static void SetNumericalDiffusionParamters(Parameters& rParameters)
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
            mNumericalDiffusionMethod = NUMERICAL_DIFFUSION_METHOD::SINGULAR_VALUE_PRESSURE_COUPLED;
        else if (method_name.compare("singular_value_pressure_decoupled")==0)
            mNumericalDiffusionMethod = NUMERICAL_DIFFUSION_METHOD::SINGULAR_VALUE_PRESSURE_DECOUPLED; //type: camel case, values: camel case or lower case 
        else if (method_name.compare("eigen_values")==0)
            mNumericalDiffusionMethod = NUMERICAL_DIFFUSION_METHOD::EIGEN_VALUE;
        else
            KRATOS_THROW_ERROR(std::runtime_error,
                "numerical diffusion method only supports singular_value_pressure_coupled or singular_value_pressure_decoupled or eigen_values",
                rParameters.PrettyPrintJsonString())
            
        KRATOS_CATCH("");
        
    }

    static void CalculateNumericalDiffusion(
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
            VMSAdjointElement<2> *pVMSAdjointElement = dynamic_cast<VMSAdjointElement<2>*>(pCurrentElement.get());

            InitializeMatrix<2>(rAdjointMatrix);

            switch (mNumericalDiffusionMethod)
            {
                case NUMERICAL_DIFFUSION_METHOD::SINGULAR_VALUE_PRESSURE_COUPLED:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureCoupled<2>(pVMSAdjointElement);
                    break;
                case NUMERICAL_DIFFUSION_METHOD::SINGULAR_VALUE_PRESSURE_DECOUPLED:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureDecoupled<2>(pVMSAdjointElement);
                    break;
                case NUMERICAL_DIFFUSION_METHOD::EIGEN_VALUE:
                    break;
            }

            boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX; DN_DX;
            array_1d< double, 3 > N;
            double Volume;
    
            GeometryUtils::CalculateGeometryData(pCurrentElement->GetGeometry(),DN_DX,N,Volume);

            AddNumericalDiffusionTerm2D<2>(rAdjointMatrix, DN_DX, mBeta * numerical_diffusion * Volume);            
        }
        else if (domain_size == 3)
        {
            VMSAdjointElement<3> *pVMSAdjointElement = static_cast<VMSAdjointElement<3>*>(pCurrentElement.get());

            InitializeMatrix<3>(rAdjointMatrix);

            switch (mNumericalDiffusionMethod)
            {
                case NUMERICAL_DIFFUSION_METHOD::SINGULAR_VALUE_PRESSURE_COUPLED:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureCoupled<3>(pVMSAdjointElement);
                    break;
                case NUMERICAL_DIFFUSION_METHOD::SINGULAR_VALUE_PRESSURE_DECOUPLED:
                    numerical_diffusion = CalculateNumericalDiffusionSVMethodPressureDecoupled<3>(pVMSAdjointElement);
                    break;
                case NUMERICAL_DIFFUSION_METHOD::EIGEN_VALUE:
                    break;
            }

            boost::numeric::ublas::bounded_matrix<double, 4, 3> DN_DX; DN_DX;
            array_1d< double, 4 > N;
            double Volume;
    
            GeometryUtils::CalculateGeometryData(pCurrentElement->GetGeometry(),DN_DX,N,Volume);
            
            AddNumericalDiffusionTerm3D<3>(rAdjointMatrix, DN_DX, mBeta * numerical_diffusion * Volume);
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                "numerical diffusion method only supports 2D or 3D elements",
                "")            
        }

        KRATOS_CATCH("");
    }
    ///@}

private:
    ///@name Member Variables
    ///@{
    static NUMERICAL_DIFFUSION_METHOD mNumericalDiffusionMethod;
    static double mBeta;
    ///@}
    ///@name Private Operators
    ///@{
    template<unsigned int TDim>
    static void InitializeMatrix(MatrixType& rAdjointMatrix)
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
    double static CalculateNumericalDiffusionSVMethodPressureCoupled(VMSAdjointElement<TDim> *pCurrentElement)
    {
        constexpr unsigned int TNumNodes = TDim + 1;        

        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;

        GeometryType& rGeom = pCurrentElement->GetGeometry();

        GeometryUtils::CalculateGeometryData(rGeom,DN_DX,N,Volume);

        boost::numeric::ublas::bounded_matrix< double, TDim, TDim > GradVel;
        pCurrentElement->CalculateVelocityGradient(GradVel,DN_DX);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<TDim; i++)
            velocity_divergence += GradVel(i,i);
        
        boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> M( TDim + 1, TDim + 1);
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

        boost::numeric::ublas::vector<double> S(TDim);
        int ierr = boost::numeric::bindings::lapack::gesvd(M, S);

        double numerical_viscosity;
        if (ierr == 0) {
            numerical_viscosity = S(0);
        }
        else
            numerical_viscosity = 0.0;

        return numerical_viscosity;
    }

    template<unsigned int TDim>
    double static CalculateNumericalDiffusionSVMethodPressureDecoupled(VMSAdjointElement<TDim> *pCurrentElement)
    {
        constexpr unsigned int TNumNodes = TDim + 1;

        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;

        GeometryType& rGeom = pCurrentElement->GetGeometry();

        GeometryUtils::CalculateGeometryData(rGeom,DN_DX,N,Volume);

        boost::numeric::ublas::bounded_matrix< double, TDim, TDim > GradVel;
        pCurrentElement->CalculateVelocityGradient(GradVel,DN_DX);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<TDim; i++)
            velocity_divergence += GradVel(i,i);
        
        boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> M( TDim , TDim);
        for (IndexType i=0; i < TDim; i++)
            M(i,i) = 0.5 * velocity_divergence - GradVel(i,i);
        for (IndexType i=0; i < TDim; i++)
            for (IndexType j=i+1; j < TDim; j++)
            {
                M(i,j) = -GradVel(i,j);
                M(j,i) = -GradVel(j,i);
            }

        boost::numeric::ublas::vector<double> S(TDim);
        int ierr = boost::numeric::bindings::lapack::gesvd(M, S);

        double numerical_viscosity;
        if (ierr == 0) {
            numerical_viscosity = S(0);
        }
        else
            numerical_viscosity = 0.0;

        return numerical_viscosity;    
    }

    template<unsigned int TDim>
    void static AddNumericalDiffusionTerm2D(
        MatrixType& rResult,
        const boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim>& rDN_DX,
        const double Weight)
    {
      const SizeType NumNodes = 3;
    
      IndexType FirstRow(0), FirstCol(0);
    
      for (IndexType j = 0; j < NumNodes; ++j)
      {
        for (IndexType i = 0; i < NumNodes; ++i)
        {
          // (dN_i/dx_k dN_j/dx_k)
          const double Diag = rDN_DX(i,0) * rDN_DX(j,0)
          + rDN_DX(i,1) * rDN_DX(j,1);
    
          // First Row
          rResult(FirstRow,FirstCol) += Weight
              * (Diag);
    
          // Second Row
          rResult(FirstRow+1,FirstCol+1) += Weight
              * (Diag);
    
          // Third Row
          rResult(FirstRow+2,FirstCol+2) += Weight
          * (Diag);
    
          // Update Counter
          FirstRow += 3;
        }
        FirstRow = 0;
        FirstCol += 3;
      }
    }
    
    template<unsigned int TDim>
    void static AddNumericalDiffusionTerm3D(
        MatrixType& rResult,
        const boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim>& rDN_DX,
        const double Weight)
    {
      const unsigned int NumNodes = 4;
    
      unsigned int FirstRow(0), FirstCol(0);
    
      for (unsigned int j = 0; j < NumNodes; ++j)
      {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
          // (dN_i/dx_k dN_j/dx_k)
          const double Diag = rDN_DX(i,0) * rDN_DX(j,0)
              + rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,2) * rDN_DX(j,2);
    
          // First Row
          rResult(FirstRow,FirstCol) += Weight
              * (Diag);
          // Second Row
          rResult(FirstRow + 1, FirstCol + 1) += Weight
              * (Diag);
    
          // Third Row
          rResult(FirstRow+2,FirstCol+2) += Weight
              * (Diag);
    
          // Fourth Row
          rResult(FirstRow+3,FirstCol+3) += Weight
          * (Diag);
    
          // Update Counter
          FirstRow += 4;
        }
        FirstRow = 0;
        FirstCol += 4;
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
