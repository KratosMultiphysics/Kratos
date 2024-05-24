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


#if !defined(KRATOS_GEO_GENERAL_U_PW_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_GENERAL_U_PW_DIFF_ORDER_CONDITION_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeneralUPwDiffOrderCondition : public Condition
{

public:

    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( GeneralUPwDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    GeneralUPwDiffOrderCondition() : GeneralUPwDiffOrderCondition(0, nullptr, nullptr) {};

    GeneralUPwDiffOrderCondition( IndexType               NewId,
                                  GeometryType::Pointer   pGeometry )
        : GeneralUPwDiffOrderCondition(NewId, pGeometry, nullptr)
    {}

    GeneralUPwDiffOrderCondition( IndexType               NewId,
                                  GeometryType::Pointer   pGeometry,
                                  PropertiesType::Pointer pProperties )
        : Condition(NewId, pGeometry, pProperties)
    {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ConditionVariables
    {
        //Variables at all integration points
        Matrix NuContainer;
        Matrix NpContainer;
        GeometryType::JacobiansType JContainer;

        //Variables at each integration point
        Vector Nu; //Contains the displacement shape functions at every node
        Vector Np; //Contains the pressure shape functions at every node
        double IntegrationCoefficient;

        //Imposed condition at all nodes
        Vector ConditionVector;
    };

    // Member Variables
    Geometry< Node >::Pointer mpPressureGeometry;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      bool CalculateLHSMatrixFlag,
                      bool CalculateResidualVectorFlag);

    void InitializeConditionVariables(ConditionVariables& rVariables,
                                      const ProcessInfo& rCurrentProcessInfo);

    void CalculateKinematics(ConditionVariables& rVariables, unsigned int PointNumber);

    virtual void CalculateConditionVector(ConditionVariables& rVariables,
                                          unsigned int PointNumber);

    virtual double CalculateIntegrationCoefficient(const IndexType PointNumber,
                                                   const GeometryType::JacobiansType& JContainer,
                                                   const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const;


    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables);

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables);

    virtual void CalculateAndAddConditionForce(VectorType& rRightHandSideVector,
                                                ConditionVariables& rVariables);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

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

}; // class GeneralUPwDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_GENERAL_U_PW_DIFF_ORDER_CONDITION_H_INCLUDED defined
