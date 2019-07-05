//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined (KRATOS_HISTORY_LINEAR_ELASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_HISTORY_LINEAR_ELASTIC_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) HistoryLinearElastic3DLaw : public LinearElastic3DLaw
{

public:

    /// Definition of the base class
    typedef LinearElastic3DLaw BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(HistoryLinearElastic3DLaw);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    HistoryLinearElastic3DLaw() {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<HistoryLinearElastic3DLaw>(HistoryLinearElastic3DLaw(*this));
    }

    /// Copy Constructor
    HistoryLinearElastic3DLaw(const HistoryLinearElastic3DLaw& rOther)
        : BaseType(rOther),
            mInitialStressVector(rOther.mInitialStressVector) {}

    /// Destructor
    ~HistoryLinearElastic3DLaw() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) override;

    void SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    Vector mInitialStressVector;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateStress( const Vector &rStrainVector,
                          const Matrix &rConstitutiveMatrix,
                          Vector& rStressVector) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
        rSerializer.save("InitialStressVector", mInitialStressVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
        rSerializer.load("InitialStressVector", mInitialStressVector);
    }

}; // Class HistoryLinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_HISTORY_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined
