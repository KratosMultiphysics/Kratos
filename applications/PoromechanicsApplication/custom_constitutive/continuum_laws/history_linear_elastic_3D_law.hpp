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
#include "includes/checks.h"

// Application includes
#include "utilities/math_utils.h"
#include "custom_constitutive/continuum_laws/linear_elastic_3D_law.hpp"
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
        : BaseType(rOther) {}

    /// Destructor
    ~HistoryLinearElastic3DLaw() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddInitialStresses( Parameters& rValues, Vector& rStressVector);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

}; // Class HistoryLinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_HISTORY_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined
