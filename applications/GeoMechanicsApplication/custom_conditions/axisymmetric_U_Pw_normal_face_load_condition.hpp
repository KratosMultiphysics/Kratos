// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Vahid Galavi
//


#if !defined(KRATOS_GEO_AXISYMMETRIC_U_PW_NORMAL_FACE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_AXISYMMETRIC_U_PW_NORMAL_FACE_LOAD_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_normal_face_load_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION)
    AxisymmetricUPwNormalFaceLoadCondition : public UPwNormalFaceLoadCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AxisymmetricUPwNormalFaceLoadCondition );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwNormalFaceLoadCondition<TDim,TNumNodes>::mThisIntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    AxisymmetricUPwNormalFaceLoadCondition() : UPwNormalFaceLoadCondition<TDim,TNumNodes>() {}

    // Constructor 1
    AxisymmetricUPwNormalFaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwNormalFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry) {}

    // Constructor 2
    AxisymmetricUPwNormalFaceLoadCondition( IndexType NewId,
                                            GeometryType::Pointer pGeometry,
                                            PropertiesType::Pointer pProperties )
                                            : UPwNormalFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    // Destructor
    ~AxisymmetricUPwNormalFaceLoadCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateIntegrationCoefficient( const IndexType PointNumber,
                                            const GeometryType::IntegrationPointsArrayType& IntegrationPoints ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }

}; // class AxisymmetricUPwNormalFaceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_AXISYMMETRIC_U_PW_NORMAL_FACE_LOAD_CONDITION_H_INCLUDED defined
