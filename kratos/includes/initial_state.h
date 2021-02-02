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

#if !defined(KRATOS_INITIAL_STATE_H_INCLUDED)
#define  KRATOS_INITIAL_STATE_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes

namespace Kratos 
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class InitialState
 * @ingroup StructuralMechanicsApplication
 * @brief Define the initial state of the material in terms of initial stress/strain/F
 * @details Storages the information regarding initial stresses/strains/F
 * @author Alejandro Cornejo
 */
class KRATOS_API(KRATOS_CORE) InitialState
{
    public:
        ///@name Type Definitions
        ///@{
        typedef std::size_t SizeType;

        enum class InitialImposingType {StrainOnly = 0, StressOnly = 1, DeformationGradientOnly = 2, StrainAndStress = 3, DeformationGradientAndStress = 4};

        /// Pointer definition of NodeSearchUtility
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(InitialState);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        InitialState()
        {}

        /// Only defining Dimension constructor.
        InitialState(const SizeType Dimension)
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
        InitialState(const Vector& rInitialStrainVector,  
                     const Vector& rInitialStressVector, 
                     const Matrix& rInitialDeformationGradientMatrix)
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
        InitialState(const Vector& rImposingEntity,
                     const int InitialImposingType = 0)
        {
            const SizeType voigt_size = rImposingEntity.size();
            const SizeType dimension = (voigt_size == 6) ? 3 : 2;

            mInitialStrainVector.resize(voigt_size, false);
            mInitialStressVector.resize(voigt_size, false);
            mInitialDeformationGradientMatrix.resize(dimension, dimension, false);

            noalias(mInitialDeformationGradientMatrix) = ZeroMatrix(dimension, dimension);
            noalias(mInitialStrainVector) = ZeroVector(voigt_size);
            noalias(mInitialStressVector) = ZeroVector(voigt_size);

            if (InitialImposingType == static_cast<int>(InitialImposingType::StrainOnly)) {
                noalias(mInitialStrainVector) = rImposingEntity;
            } else if (InitialImposingType == static_cast<int>(InitialImposingType::StressOnly)) {
                noalias(mInitialStressVector) = rImposingEntity;       
            }
        }

        // Selective constructor for vectors (E, S)
        InitialState(const Vector& rInitialStrainVector,
                     const Vector& rInitialStressVector)
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
        InitialState(const Matrix& rInitialDeformationGradientMatrix)
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

        /// Destructor.
        ~InitialState() override {}


        ///@}
        ///@name Operations
        ///@{
        
        //*********************************************
        //public API of intrusive_ptr
        unsigned int use_count() const noexcept
        {
            return mReferenceCounter;
        }
        //*********************************************

        /**
         * @brief This method sets the initial strain vector
         * @param rInitialStrainVector The vector to be set
         */
        void SetInitialStrainVector(const Vector& rInitialStrainVector) {
            const SizeType voigt_size = rInitialStrainVector.size();
            KRATOS_ERROR_IF(voigt_size <= 0) << "The imposed vector is null..." << std::endl;
            
            mInitialStrainVector.resize(voigt_size, false);
            noalias(mInitialStrainVector) = rInitialStrainVector;
        }

        /**
         * @brief This method sets the initial stress vector
         * @param rInitialStressVector The vector to be set
         */
        void SetInitialStressVector(const Vector& rInitialStressVector) {
            const SizeType voigt_size = rInitialStressVector.size();
            KRATOS_ERROR_IF(voigt_size <= 0) << "The imposed vector is null..." << std::endl;
            
            mInitialStressVector.resize(voigt_size, false);
            noalias(mInitialStressVector) = rInitialStressVector;
        }

        /**
         * @brief This method sets the initial deformation gradient matrix
         * @param rInitialDeformationGradientMatrix The vector to be set
         */
        void SetInitialDeformationGradientMatrix(const Matrix& rInitialDeformationGradientMatrix) {
            const SizeType dimension = rInitialDeformationGradientMatrix.size1();
            KRATOS_ERROR_IF(dimension <= 0) << "The imposed Matrix is null..." << std::endl;
            
            mInitialDeformationGradientMatrix.resize(dimension, dimension, false);
            noalias(mInitialDeformationGradientMatrix) = rInitialDeformationGradientMatrix;
        }

        /**
         * @brief This method returns the initial strain vector if was set before
         */
        const Vector& GetInitialStrainVector() const
        {
            return mInitialStrainVector;
        }

        /**
         * @brief This method returns the initial stress vector if was set before
         */
        const Vector& GetInitialStressVector() const
        {
            return mInitialStressVector;
        }

        /**
         * @brief This method returns the initial stress vector if was set before
         */
        const Matrix& GetInitialDeformationGradientMatrix() const
        {
            return mInitialDeformationGradientMatrix;
        }


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "InitialState" ;

            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "InitialState";}

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const {}

        ///@}

    private:

        ///@name Member Variables
        ///@{
        Vector mInitialStrainVector;
        Vector mInitialStressVector;
        Matrix mInitialDeformationGradientMatrix;


    //*********************************************
        //this block is needed for refcounting
        mutable std::atomic<int> mReferenceCounter;
        friend void intrusive_ptr_add_ref(const InitialState* x)
        {
            x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
        }

        friend void intrusive_ptr_release(const InitialState* x)
        {
            if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
                std::atomic_thread_fence(std::memory_order_acquire);
                delete x;
            }
        }
    //*********************************************

}; // class
///@}

///@} addtogroup block

} // namespace Kratos

#endif // KRATOS_INITIAL_STATE_H_INCLUDED  defined