// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#if !defined(KRATOS_GEO_THERMAL_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_THERMAL_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/stress_strain_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

    template< unsigned int TDim, unsigned int TNumNodes >
    class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientThermalElement :
        public Element
    {

    public:

        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientThermalElement);
        //
        typedef Element BaseType;


        /// Default Constructor
        TransientThermalElement(IndexType NewId = 0) : BaseType(NewId) {}

        /// Constructor using an array of nodes
        TransientThermalElement(IndexType NewId, const NodesArrayType& ThisNodes) : BaseType(NewId, ThisNodes) {}

        /// Constructor using Geometry
        TransientThermalElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry) {}

        /// Constructor using Properties
        TransientThermalElement(IndexType NewId, GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties) {}

        /// Destructor
        ~TransientThermalElement() override {}

        ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        int Check(const ProcessInfo& rCurrentProcessInfo) const override;



    }; // Class TransientThermalElement

} // namespace Kratos

#endif // KRATOS_GEO_THERMAL_ELEMENT_H_INCLUDED  defined