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


#if !defined(KRATOS_GEO_U_PW_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_U_PW_CONDITION_H_INCLUDED

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
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwCondition : public Condition
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    UPwCondition() : UPwCondition(0, nullptr, nullptr){}

    UPwCondition( IndexType               NewId,
                  GeometryType::Pointer   pGeometry )
        : UPwCondition(NewId, pGeometry, nullptr)
    {}

    UPwCondition( IndexType               NewId,
                  GeometryType::Pointer   pGeometry,
                  PropertiesType::Pointer pProperties )
        : Condition(NewId, pGeometry, pProperties)
    {}

    ~UPwCondition() override = default;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties ) const override;

    void GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const override;

    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;
    }

    void SetIntegrationMethod(IntegrationMethod method)
    {
        mThisIntegrationMethod = method;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateRHS(VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

    [[nodiscard]] DofsVectorType GetDofs() const;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    GeometryData::IntegrationMethod mThisIntegrationMethod{ Condition::GetIntegrationMethod() };

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
    
}; // class UPwCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_U_PW_CONDITION_H_INCLUDED defined 
