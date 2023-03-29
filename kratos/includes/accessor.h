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
#include "geometries/geometry.h"
#include "includes/process_info.h"
#include "includes/node.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

class Properties;

/**
 * @class Accessor
 * @ingroup Kratos Core
 * @brief This class defines the way a certain property is accessed.
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

        /// Destructor.
        virtual ~Accessor() {}

        ///@}
        ///@name Operations
        ///@{
        
        /**
         * Custom method to retrieve double type properties
         */
        virtual double const& GetProperty(
            const Variable<double> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        /**
         * Custom method to retrieve Vector type properties
         */
        virtual Vector const& GetProperty(
            const Variable<Vector> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        /**
         * Custom method to retrieve bool type properties
         */
        virtual bool const& GetProperty(
            const Variable<bool> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        /**
         * Custom method to retrieve int type properties
         */
        virtual int const& GetProperty(
            const Variable<int> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        /**
         * Custom method to retrieve Matrix type properties
         */
        virtual Matrix const& GetProperty(
            const Variable<Matrix> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        /**
         * Custom method to retrieve array_1d<double, 3 > type properties
         */
        virtual array_1d<double, 3 > const& GetProperty(
            const Variable<array_1d<double, 3 >> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        /**
         * Custom method to retrieve array_1d<double, 6 > type properties
         */
        virtual array_1d<double, 6 > const& GetProperty(
            const Variable<array_1d<double, 6 >> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        /**
         * Custom method to retrieve string type properties
         */
        virtual std::string const& GetProperty(
            const Variable<std::string> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const;

        // Getting a pointer tot he class
        virtual Accessor::Pointer Clone() const;

        /**
         * Copy constructor.
         */
        Accessor(const Accessor &rOther);

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
        {}

        void load(Serializer& rSerializer)
        {}


}; // class
///@}

///@} addtogroup block

} // namespace Kratos
