//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//


#if !defined(KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED

/* Project includes */
#include "includes/serializer.h"
#include "custom_elements/small_displacement_element.hpp"

#include "custom_utilities/poro_element_utilities.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

class SmallDisplacementThermoMechanicElement : public SmallDisplacementElement
{

  
  
public:
  
    ///Type for element variables
    typedef SmallDisplacementElement::ElementDataType ElementDataType;


  
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

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint);

    void ExtrapolateGPStress(const Matrix& StressContainer, const unsigned int& Dim, const unsigned int& VoigtSize);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementElement )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementElement )
    }


}; // Class SmallDisplacementThermoMechanicElement

} // namespace Kratos

#endif // KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED  defined
