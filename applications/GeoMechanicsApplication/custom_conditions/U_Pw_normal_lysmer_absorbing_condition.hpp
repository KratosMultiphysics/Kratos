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
//  Main authors:    Aron Noordam
//


#if !defined(KRATOS_GEO_U_PW_NORMAL_LYSMER_ABSORBING_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_U_PW_NORMAL_LYSMER_ABSORBING_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_conditions/U_Pw_face_load_condition.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwLysmerAbsorbingCondition : public UPwFaceLoadCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwLysmerAbsorbingCondition);
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwLysmerAbsorbingCondition() : UPwFaceLoadCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwLysmerAbsorbingCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwLysmerAbsorbingCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPwLysmerAbsorbingCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct NormalLysmerAbsorbingVariables
    {
        //array_1d<double, TNumNodes> NormalStressVector;
        //array_1d<double, TNumNodes> TangentialStressVector;
        double NormalAbsorb;
        array_1d<double, TDim> localRelVelVector;
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> Np;
        BoundedMatrix<double, TDim, TNumNodes* TDim> Nu;
        array_1d<double, TNumNodes* TDim> UVector;
    };
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                

    //void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo);

    void CalculateRHS(VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo) override;
    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalLysmerAbsorbingVariables& rVariables);

    void CalculateLocalVelocityVector(array_1d<double, TDim>& rLocalVelocityVector,
        const Matrix& Jacobian,
        const Matrix& NContainer,
        const NormalLysmerAbsorbingVariables& Variables,
 
        const unsigned int& GPoint);

    void CalculateRotationMatrix(BoundedMatrix<double, TDim, TDim>& rRotationMatrix, const Element::GeometryType& Geom);

    void CalculateAbsorbingMatrix(MatrixType& rAbsorbingMatrix,
        const ProcessInfo& rCurrentProcessInfo);
    
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
    
}; // class UPwLysmerAbsorbingCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_U_PW_NORMAL_LYSMER_ABSORBING_CONDITION_H_INCLUDED defined 
