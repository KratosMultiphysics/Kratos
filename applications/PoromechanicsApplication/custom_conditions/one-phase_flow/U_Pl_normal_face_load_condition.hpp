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


#if !defined(KRATOS_U_PL_NORMAL_FACE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_NORMAL_FACE_LOAD_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/one-phase_flow/U_Pl_condition.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlNormalFaceLoadCondition : public UPlCondition<TDim,TNumNodes>
{

public:

<<<<<<< HEAD:applications/PoromechanicsApplication/custom_conditions/U_Pw_normal_face_load_condition.hpp
    KRATOS_CLASS_POINTER_DEFINITION( UPwNormalFaceLoadCondition );

=======
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlNormalFaceLoadCondition );
    
>>>>>>> master:applications/PoromechanicsApplication/custom_conditions/one-phase_flow/U_Pl_normal_face_load_condition.hpp
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
<<<<<<< HEAD:applications/PoromechanicsApplication/custom_conditions/U_Pw_normal_face_load_condition.hpp
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwNormalFaceLoadCondition() : UPwCondition<TDim,TNumNodes>() {}

    // Constructor 1
    UPwNormalFaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwCondition<TDim,TNumNodes>(NewId, pGeometry) {}

=======
    using UPlCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPlNormalFaceLoadCondition() : UPlCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPlNormalFaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPlCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
>>>>>>> master:applications/PoromechanicsApplication/custom_conditions/one-phase_flow/U_Pl_normal_face_load_condition.hpp
    // Constructor 2
    UPlNormalFaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPlCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPlNormalFaceLoadCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct NormalFaceLoadVariables
    {
        array_1d<double,TNumNodes> NormalStressVector;
        array_1d<double,TNumNodes> TangentialStressVector;
    };

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
<<<<<<< HEAD:applications/PoromechanicsApplication/custom_conditions/U_Pw_normal_face_load_condition.hpp

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

=======
    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;
    
>>>>>>> master:applications/PoromechanicsApplication/custom_conditions/one-phase_flow/U_Pl_normal_face_load_condition.hpp
    void InitializeConditionVariables(NormalFaceLoadVariables& rVariables, const GeometryType& Geom);

    void CalculateTractionVector(array_1d<double,TDim>& rTractionVector,const Matrix& Jacobian,const Matrix& NContainer,
                                    const NormalFaceLoadVariables& Variables,const unsigned int& GPoint);

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& Weight, const ProcessInfo& CurrentProcessInfo);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        typedef UPlCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseCondition )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseCondition )
    }
<<<<<<< HEAD:applications/PoromechanicsApplication/custom_conditions/U_Pw_normal_face_load_condition.hpp

}; // class UPwNormalFaceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_U_PW_NORMAL_FACE_LOAD_CONDITION_H_INCLUDED defined
=======
    
}; // class UPlNormalFaceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_NORMAL_FACE_LOAD_CONDITION_H_INCLUDED defined 
>>>>>>> master:applications/PoromechanicsApplication/custom_conditions/one-phase_flow/U_Pl_normal_face_load_condition.hpp
