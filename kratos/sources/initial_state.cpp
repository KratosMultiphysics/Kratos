//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//


// System includes

// External includes

// Project includes
#include "includes/initial_state.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{


    /// Only defining Dimension constructor.
    InitialState::InitialState(const SizeType Dimension)
        : mReferenceCounter(0)
    {
        const SizeType voigt_size = (Dimension == 3) ? 6 : 3;
        mInitialStressVector.resize(voigt_size, false);
        mInitialStrainVector.resize(voigt_size, false);
        mInitialDeformationGradientMatrix.resize(Dimension, Dimension, false);

        noalias(mInitialStressVector) = ZeroVector(voigt_size);
        noalias(mInitialStrainVector) = ZeroVector(voigt_size);
        noalias(mInitialDeformationGradientMatrix) = ZeroMatrix(Dimension, Dimension);
    }

    // Full constructor
    InitialState::InitialState(const Vector& rInitialStrainVector,
                    const Vector& rInitialStressVector,
                    const Matrix& rInitialDeformationGradientMatrix)
        : mReferenceCounter(0)
    {
        const SizeType voigt_size = rInitialStrainVector.size();
        const SizeType dimension  = rInitialDeformationGradientMatrix.size1();

        KRATOS_ERROR_IF(voigt_size <= 0 || dimension <= 0) << "The imposed vector or matrix is null..." << std::endl;

        mInitialStressVector.resize(voigt_size, false);
        mInitialStrainVector.resize(voigt_size, false);
        mInitialDeformationGradientMatrix.resize(dimension, dimension, false);

        noalias(mInitialStressVector) = rInitialStressVector;
        noalias(mInitialStrainVector) = rInitialStrainVector;
        noalias(mInitialDeformationGradientMatrix) = rInitialDeformationGradientMatrix;
    }

    // Selective constructor for vectors
    InitialState::InitialState(const Vector& rImposingEntity,
                    const InitialImposingType InitialImposition)
        : mReferenceCounter(0)
    {
        const SizeType voigt_size = rImposingEntity.size();
        const SizeType dimension = (voigt_size == 6) ? 3 : 2;

        mInitialStrainVector.resize(voigt_size, false);
        mInitialStressVector.resize(voigt_size, false);
        mInitialDeformationGradientMatrix.resize(dimension, dimension, false);

        noalias(mInitialDeformationGradientMatrix) = ZeroMatrix(dimension, dimension);
        noalias(mInitialStrainVector) = ZeroVector(voigt_size);
        noalias(mInitialStressVector) = ZeroVector(voigt_size);

        if (InitialImposition == InitialImposingType::STRAIN_ONLY) {
            noalias(mInitialStrainVector) = rImposingEntity;
        } else if (InitialImposition == InitialImposingType::STRESS_ONLY) {
            noalias(mInitialStressVector) = rImposingEntity;
        }
    }

    // Selective constructor for vectors (E, S)
    InitialState::InitialState(const Vector& rInitialStrainVector,
                    const Vector& rInitialStressVector)
        : mReferenceCounter(0)
    {
        const SizeType voigt_size_1 = rInitialStrainVector.size();
        const SizeType voigt_size_2 = rInitialStressVector.size();
        const SizeType dimension = (voigt_size_1 == 6) ? 3 : 2;
        KRATOS_ERROR_IF(voigt_size_1 <= 0 || voigt_size_2 <= 0) << "The imposed vector is null..." << std::endl;

        mInitialStressVector.resize(voigt_size_1, false);
        mInitialStrainVector.resize(voigt_size_1, false);
        mInitialDeformationGradientMatrix.resize(dimension, dimension, false);

        noalias(mInitialDeformationGradientMatrix) = ZeroMatrix(dimension, dimension);
        noalias(mInitialStressVector) = rInitialStressVector;
        noalias(mInitialStrainVector) = rInitialStrainVector;
    }

    // Selective constructor for Deformation Gradient only
    InitialState::InitialState(const Matrix& rInitialDeformationGradientMatrix)
        : mReferenceCounter(0)
    {
        const SizeType dimension = rInitialDeformationGradientMatrix.size1();
        const SizeType voigt_size = (dimension == 3) ? 6 : 3;
        KRATOS_ERROR_IF(dimension <= 0) << "The imposed Matrix is null..." << std::endl;

        mInitialDeformationGradientMatrix.resize(dimension, dimension, false);
        noalias(mInitialDeformationGradientMatrix) = rInitialDeformationGradientMatrix;

        mInitialStressVector.resize(voigt_size, false);
        mInitialStrainVector.resize(voigt_size, false);
        noalias(mInitialStrainVector) = ZeroVector(voigt_size);
        noalias(mInitialStressVector) = ZeroVector(voigt_size);
    }


    /**
     * @brief This method sets the initial strain vector
     * @param rInitialStrainVector The vector to be set
     */
    void InitialState::SetInitialStrainVector(const Vector& rInitialStrainVector) {
        const SizeType voigt_size = rInitialStrainVector.size();
        KRATOS_ERROR_IF(voigt_size <= 0) << "The imposed vector is null..." << std::endl;

        mInitialStrainVector.resize(voigt_size, false);
        noalias(mInitialStrainVector) = rInitialStrainVector;
    }

    /**
     * @brief This method sets the initial stress vector
     * @param rInitialStressVector The vector to be set
     */
    void InitialState::SetInitialStressVector(const Vector& rInitialStressVector) {
        const SizeType voigt_size = rInitialStressVector.size();
        KRATOS_ERROR_IF(voigt_size <= 0) << "The imposed vector is null..." << std::endl;

        mInitialStressVector.resize(voigt_size, false);
        noalias(mInitialStressVector) = rInitialStressVector;
    }

    /**
     * @brief This method sets the initial deformation gradient matrix
     * @param rInitialDeformationGradientMatrix The vector to be set
     */
    void InitialState::SetInitialDeformationGradientMatrix(const Matrix& rInitialDeformationGradientMatrix) {
        const SizeType dimension = rInitialDeformationGradientMatrix.size1();
        KRATOS_ERROR_IF(dimension <= 0) << "The imposed Matrix is null..." << std::endl;

        mInitialDeformationGradientMatrix.resize(dimension, dimension, false);
        noalias(mInitialDeformationGradientMatrix) = rInitialDeformationGradientMatrix;
    }

    /**
     * @brief This method returns the initial strain vector if was set before
     */
    const Vector& InitialState::GetInitialStrainVector() const
    {
        return mInitialStrainVector;
    }

    /**
     * @brief This method returns the initial stress vector if was set before
     */
    const Vector& InitialState::GetInitialStressVector() const
    {
        return mInitialStressVector;
    }

    /**
     * @brief This method returns the initial stress vector if was set before
     */
    const Matrix& InitialState::GetInitialDeformationGradientMatrix() const
    {
        return mInitialDeformationGradientMatrix;
    }

} // namespace Kratos
