//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_SMALL_STRAIN_U_PW_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_U_PW_INTERFACE_ELEMENT_H_INCLUDED

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/poromechanics_math_utilities.hpp"

namespace Kratos
{

class SmallStrainUPwInterfaceElement : public Element
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SmallStrainUPwInterfaceElement );
    
    typedef GeometryData::IntegrationMethod IntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SmallStrainUPwInterfaceElement();
    
    // Constructor 1
    SmallStrainUPwInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry);
    
    // Constructor 2
    SmallStrainUPwInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    
    // Destructor
    virtual ~SmallStrainUPwInterfaceElement();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    int Check(const ProcessInfo& rCurrentProcessInfo);


}; // Class SmallStrainUPwInterfaceElement

} // namespace Kratos

#endif // KRATOS_SMALL_STRAIN_U_PW_INTERFACE_ELEMENT_H_INCLUDED  defined 
