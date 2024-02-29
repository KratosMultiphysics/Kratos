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

#if !defined (KRATOS_HISTORY_LINEAR_ELASTIC_PLANE_STRESS_2D_LAW_H_INCLUDED)
#define  KRATOS_HISTORY_LINEAR_ELASTIC_PLANE_STRESS_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/continuum_laws/history_linear_elastic_plane_strain_2D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) HistoryLinearElasticPlaneStress2DLaw : public HistoryLinearElasticPlaneStrain2DLaw
{

public:

    typedef HistoryLinearElasticPlaneStrain2DLaw BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(HistoryLinearElasticPlaneStress2DLaw);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    HistoryLinearElasticPlaneStress2DLaw() {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<HistoryLinearElasticPlaneStress2DLaw>(HistoryLinearElasticPlaneStress2DLaw(*this));
    }

    /// Copy Constructor
    HistoryLinearElasticPlaneStress2DLaw (const HistoryLinearElasticPlaneStress2DLaw& rOther)
        : BaseType(rOther) {}

    /// Destructor
    ~HistoryLinearElasticPlaneStress2DLaw() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient ) override;

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

}; // Class HistoryLinearElasticPlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_HISTORY_LINEAR_ELASTIC_PLANE_STRESS_2D_LAW_H_INCLUDED  defined
