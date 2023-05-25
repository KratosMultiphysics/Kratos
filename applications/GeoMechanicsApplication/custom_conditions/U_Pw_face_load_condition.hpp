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
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//


#if !defined(KRATOS_GEO_U_PW_FACE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_U_PW_FACE_LOAD_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwFaceLoadCondition
    : public UPwCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwFaceLoadCondition );
    
    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwFaceLoadCondition() : UPwCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwFaceLoadCondition( IndexType NewId,
                          GeometryType::Pointer pGeometry ) :
                        UPwCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwFaceLoadCondition( IndexType NewId,
                          GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties ) :
                        UPwCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPwFaceLoadCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateRHS(VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo) override;

    virtual double CalculateIntegrationCoefficient(const Matrix& Jacobian, const double& Weight);

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
    
}; // class UPwFaceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_U_PW_FACE_LOAD_CONDITION_H_INCLUDED defined 
