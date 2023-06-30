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
//                   Danilo Cavalcanti
//


#if !defined(KRATOS_U_PW_PG_NORMAL_GAS_FLUX_CONDITION_H_INCLUDED )
#define  KRATOS_U_PW_PG_NORMAL_GAS_FLUX_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/multiphase_flow/U_Pw_Pg_condition.hpp"
#include "custom_conditions/multiphase_flow/U_Pw_Pg_face_load_condition.hpp"
#include "custom_utilities/poro_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwPgNormalGasFluxCondition : public UPwPgFaceLoadCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPwPgNormalGasFluxCondition );
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwPgNormalGasFluxCondition() : UPwPgFaceLoadCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwPgNormalGasFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwPgFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwPgNormalGasFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwPgFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPwPgNormalGasFluxCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    struct NormalFluxVariables
    {
        double NormalFlux;
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> Np;
        array_1d<double,TNumNodes> PVector;
    };
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;
    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables);
    
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
    
}; // class UPwPgNormalGasFluxCondition.

} // namespace Kratos.

#endif // KRATOS_U_PW_PG_NORMAL_GAS_FLUX_CONDITION_H_INCLUDED defined 
