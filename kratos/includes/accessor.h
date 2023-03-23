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
// #include "includes/define.h"
// #include "includes/variables.h"
#include "geometries/geometry.h"
#include "includes/process_info.h"
#include "includes/node.h"
// #include "includes/properties.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class Accessor
 * @ingroup Kratos Core
 * @brief This class defines the way a certain property is accessed to
 * @author Alejandro Cornejo, Riccardo Rossi and Carlos Roig
 */
class KRATOS_API(KRATOS_CORE) Accessor
{
    public:
        ///@name Type Definitions
        ///@{
            typedef Node<3> NodeType;
            typedef Geometry<NodeType> GeometryType;
            

        /// Pointer definition of NodeSearchUtility
        KRATOS_CLASS_POINTER_DEFINITION(Accessor);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        Accessor(){};

        Accessor(const double Value)
        {
            mValue = Value;
        };

        /// Destructor.
        virtual ~Accessor() {}

        ///@}
        ///@name Operations
        ///@{

        virtual double GetProperty(
            const Variable<double> &rVariable,
            const Properties::Pointer pProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo)
        {
            return 0.0;
        }
        virtual double GetDoubleTypeProperty() { return 0.0; }

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
        double mValue = 0.0;

        friend class Serializer;

        void save(Serializer& rSerializer) const
        {
            // to be added
        }

        void load(Serializer& rSerializer)
        {
        }


}; // class
///@}

///@} addtogroup block

} // namespace Kratos
