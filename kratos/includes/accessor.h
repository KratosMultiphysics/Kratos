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
//                   Riccardo Rossi
//                   Carlos Roig
//
//

# pragma once

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
 * @class Accessor
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines the way a certain property is accessed to
 * @author Alejandro Cornejo, Riccardo Rossi and Carlos Roig
 */
class KRATOS_API(KRATOS_CORE) Accessor
{
    public:
        ///@name Type Definitions
        ///@{


        /// Pointer definition of NodeSearchUtility
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Accessor);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        Accessor()
        {}

        /// Destructor.
        virtual ~Accessor() {}

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
            buffer << "Accessor" ;

            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "Accessor";}

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const {}

        ///@}

    private:

        ///@name Member Variables
        ///@{


        friend class Serializer;

        void save(Serializer& rSerializer) const
        {
        }

        void load(Serializer& rSerializer)
        {
        }


}; // class
///@}

///@} addtogroup block

} // namespace Kratos
