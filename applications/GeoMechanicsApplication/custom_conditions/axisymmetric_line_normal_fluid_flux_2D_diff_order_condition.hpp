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


#if !defined(KRATOS_GEO_AXISYMMETRIC_LINE_NORMAL_FLUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_AXISYMMETRIC_LINE_NORMAL_FLUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION)
    AxisymmetricLineNormalFluidFlux2DDiffOrderCondition : public LineNormalFluidFlux2DDiffOrderCondition
{

public:

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AxisymmetricLineNormalFluidFlux2DDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    AxisymmetricLineNormalFluidFlux2DDiffOrderCondition();
    
    // Constructor 1
    AxisymmetricLineNormalFluidFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    AxisymmetricLineNormalFluidFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    ~AxisymmetricLineNormalFluidFlux2DDiffOrderCondition() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateIntegrationCoefficient(const IndexType PointNumber,
                                           const GeometryType::JacobiansType& JContainer,
                                           const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineNormalFluidFlux2DDiffOrderCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineNormalFluidFlux2DDiffOrderCondition )
    }
    
}; // class AxisymmetricLineNormalFluidFlux2DDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_AXISYMMETRIC_LINE_NORMAL_FLUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED defined 
