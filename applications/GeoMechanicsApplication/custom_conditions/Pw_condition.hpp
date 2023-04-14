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
//                   Vahid Galavi,
//                   Aron Noordam
//


#if !defined(KRATOS_GEO_PW_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_PW_CONDITION_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) PwCondition : public Condition
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( PwCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    PwCondition() : PwCondition(0, nullptr, nullptr) {}

    PwCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : PwCondition(NewId, pGeometry, nullptr)
    {}

    PwCondition( IndexType               NewId,
                 GeometryType::Pointer   pGeometry,
                 PropertiesType::Pointer pProperties )
        : Condition(NewId, pGeometry, pProperties)
    {}

    ~PwCondition() override = default;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties ) const override;
 
    void GetDofList(DofsVectorType& rConditionDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateRHS(VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
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
    
}; // class PwCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_PW_CONDITION_H_INCLUDED defined 
