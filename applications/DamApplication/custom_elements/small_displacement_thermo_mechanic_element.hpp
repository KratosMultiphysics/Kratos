//
//   Project Name:        
//   Last modified by:    $Author:      
//   Date:                $Date:          
//   Revision:            $Revision:        
//

#if !defined(KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED

/* Project includes */
#include "includes/serializer.h"
#include "custom_elements/solid_elements/small_displacement_element.hpp"

#include "custom_utilities/element_utilities.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

class SmallDisplacementThermoMechanicElement : public SmallDisplacementElement
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SmallDisplacementThermoMechanicElement );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SmallDisplacementThermoMechanicElement();
    
    // Constructor 1
    SmallDisplacementThermoMechanicElement(IndexType NewId, GeometryType::Pointer pGeometry);
    
    // Constructor 2
    SmallDisplacementThermoMechanicElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    
    // Destructor
    virtual ~SmallDisplacementThermoMechanicElement();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);
    
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    void SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint);

    void ExtrapolateGPStress(const Matrix& StressContainer, const unsigned int& Dim, const unsigned int& VoigtSize);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementElement )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementElement )
    }
    
    
}; // Class SmallDisplacementThermoMechanicElement

} // namespace Kratos

#endif // KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED  defined 
