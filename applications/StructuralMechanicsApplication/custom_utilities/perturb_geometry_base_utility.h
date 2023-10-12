// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Manuel Messmer
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "spaces/ublas_space.h"


namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class PerturbGeometryBaseUtility
 * @ingroup StructuralMechanicsApplication
 * @brief Base class for geometry perturbation utilities.
 * @author Manuel Messmer
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PerturbGeometryBaseUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef TDenseSpaceType::MatrixPointerType DenseMatrixPointerType;

    typedef TDenseSpaceType::VectorType DenseVectorType;

    typedef TDenseSpaceType::MatrixType DenseMatrixType;

    /// Pointer definition of PerturbGeometryBaseUtility
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometryBaseUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PerturbGeometryBaseUtility( ModelPart& rInitialModelPart, Parameters Settings);

    /// Destructor.
    virtual ~PerturbGeometryBaseUtility() {}

    ///@}
    ///@name Operations
    ///@{

    virtual int CreateRandomFieldVectors() = 0;

    /**
     * @brief Assemble random field and apply to initial geometry
     * @param rPerturbationMatrix Perturbation matrix. Stores eigenvectors of correlation matrix (colum-wise).
     * @param random_field Random field vector. Stores nodal deviations.
     */
    void ApplyRandomFieldVectorsToGeometry(ModelPart& rThisModelPart, const std::vector<double>& variables );

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PerturbGeometryBaseUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PerturbGeometryBaseUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}

protected:

    ///@name Protected member Variables
    ///@{

    DenseMatrixPointerType mpPerturbationMatrix;

    ModelPart& mrInitialModelPart;

    double mCorrelationLength;

    double mTruncationError;

    int mEchoLevel;

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Correlation function
     * @return Correlation value of two nodes
     */
    double CorrelationFunction( ModelPart::NodeType& itNode1, ModelPart::NodeType& itNode2, double CorrelationLenth);

    ///@}

private:

    ///@name Member Variables
    ///@{

    double mMaximalDisplacement;

    /// Assignment operator.
    PerturbGeometryBaseUtility& operator=(PerturbGeometryBaseUtility const& rOther) = delete;

    /// Copy constructor.
    PerturbGeometryBaseUtility(PerturbGeometryBaseUtility const& rOther) = delete;

    ///@}

    }; // Class PerturbGeometryBaseUtility

///@}

///@} addtogroup block
}