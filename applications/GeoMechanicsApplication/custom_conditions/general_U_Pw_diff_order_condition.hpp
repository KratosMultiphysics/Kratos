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

#pragma once

// System includes
#include <cmath>

// Project includes
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/serializer.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeneralUPwDiffOrderCondition : public Condition
{
public:
    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeneralUPwDiffOrderCondition);

    GeneralUPwDiffOrderCondition() : GeneralUPwDiffOrderCondition(0, nullptr, nullptr) {};

    GeneralUPwDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : GeneralUPwDiffOrderCondition(NewId, pGeometry, nullptr)
    {
    }

    GeneralUPwDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const override;

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    std::string Info() const override;

protected:
    struct ConditionVariables {
        // Variables at all integration points
        Matrix                      NuContainer;
        Matrix                      NpContainer;
        GeometryType::JacobiansType JContainer;

        // Variables at each integration point
        Vector Nu; // Contains the displacement shape functions at every node
        Vector Np; // Contains the pressure shape functions at every node
        double IntegrationCoefficient;

        // Imposed condition at all nodes
        Vector ConditionVector;
    };

    // Member Variables
    Geometry<Node>::Pointer mpPressureGeometry;

    void CalculateAll(const Matrix&      ,
                      Vector&            rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      bool               ,
                      bool               CalculateResidualVectorFlag);

    void InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateKinematics(ConditionVariables& rVariables, unsigned int PointNumber);

    virtual void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber);

    virtual double CalculateIntegrationCoefficient(IndexType                          PointNumber,
                                                   const GeometryType::JacobiansType& JContainer,
                                                   const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const;

    void CalculateAndAddRHS(Vector& rRightHandSideVector, ConditionVariables& rVariables);

    virtual void CalculateAndAddConditionForce(Vector& rRightHandSideVector, ConditionVariables& rVariables);

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    }

}; // class GeneralUPwDiffOrderCondition.

} // namespace Kratos.
