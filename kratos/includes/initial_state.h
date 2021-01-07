//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_INITIAL_STATE_H_INCLUDED)
#define  KRATOS_INITIAL_STATE_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes

namespace Kratos {
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

class InitialState
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of NodeSearchUtility
        KRATOS_CLASS_POINTER_DEFINITION(InitialState);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        InitialState(const Vector& rInitialStrain,  const Vector& rInitialStress, const Matrix& rInitialDeformationGradient) {
            // todo...
        }

        /// Destructor.
        ~InitialState() {}


        ///@}
        ///@name Operations
        ///@{


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
        virtual void PrintData(std::ostream& rOStream) const  {}

        ///@}

    private:

        ///@name Member Variables
        ///@{
        Vector mInitialStrainVector;
        Vector mInitialStressVector;
        Matrix mInitialDeformationGradientMatrix;

}; // class
///@}

///@} addtogroup block

} // namespace Kratos

#endif // KRATOS_INITIAL_STATE_H_INCLUDED  defined