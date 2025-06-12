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

# pragma once

// System includes
#include <atomic>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

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

        // Size type definition
        using SizeType = std::size_t;

        /// Pointer definition of NodeSearchUtility
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(InitialState);

        ///@}
        ///@name  Enum's
        ///@{

        enum class InitialImposingType
        {
            STRAIN_ONLY = 0,
            STRESS_ONLY = 1,
            DEFORMATION_GRADIENT_ONLY = 2,
            STRAIN_AND_STRESS = 3,
            DEFORMATION_GRADIENT_AND_STRESS = 4
        };

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        InitialState()
        {}

        /// Only defining Dimension constructor.
        InitialState(const SizeType Dimension);

        // Full constructor
        InitialState(const Vector& rInitialStrainVector,
                     const Vector& rInitialStressVector,
                     const Matrix& rInitialDeformationGradientMatrix);

        // Selective constructor for vectors
        InitialState(const Vector& rImposingEntity,
                     const InitialImposingType InitialImposition = InitialImposingType::STRAIN_ONLY);

        // Selective constructor for vectors (E, S)
        InitialState(const Vector& rInitialStrainVector,
                     const Vector& rInitialStressVector);

        // Selective constructor for Deformation Gradient only
        InitialState(const Matrix& rInitialDeformationGradientMatrix);

        /// Destructor.
        virtual ~InitialState() {}


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
        void SetInitialStrainVector(const Vector& rInitialStrainVector);

        /**
         * @brief This method sets the initial stress vector
         * @param rInitialStressVector The vector to be set
         */
        void SetInitialStressVector(const Vector& rInitialStressVector);

        /**
         * @brief This method sets the initial deformation gradient matrix
         * @param rInitialDeformationGradientMatrix The vector to be set
         */
        void SetInitialDeformationGradientMatrix(const Matrix& rInitialDeformationGradientMatrix);

        /**
         * @brief This method returns the initial strain vector if was set before
         */
        const Vector& GetInitialStrainVector() const;

        /**
         * @brief This method returns the initial stress vector if was set before
         */
        const Vector& GetInitialStressVector() const;

        /**
         * @brief This method returns the initial stress vector if was set before
         */
        const Matrix& GetInitialDeformationGradientMatrix() const;


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
        mutable std::atomic<int> mReferenceCounter{0};
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


        friend class Serializer;

        void save(Serializer& rSerializer) const
        {
            rSerializer.save("InitialStrainVector",mInitialStrainVector);
            rSerializer.save("InitialStressVector",mInitialStressVector);
            rSerializer.save("InitialDeformationGradientMatrix",mInitialDeformationGradientMatrix);
        }

        void load(Serializer& rSerializer)
        {
            rSerializer.load("InitialStrainVector",mInitialStrainVector);
            rSerializer.load("InitialStressVector",mInitialStressVector);
            rSerializer.load("InitialDeformationGradientMatrix",mInitialDeformationGradientMatrix);
        }


}; // class
///@}

///@} addtogroup block

} // namespace Kratos
