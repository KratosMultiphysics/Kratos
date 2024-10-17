// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_STEADY_STATE_PW_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_STEADY_STATE_PW_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/transient_Pw_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) SteadyStatePwElement : public TransientPwElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SteadyStatePwElement);

    using BaseType = TransientPwElement<TDim, TNumNodes>;

    using IndexType            = std::size_t;
    using PropertiesType       = Properties;
    using NodeType             = Node;
    using GeometryType         = Geometry<NodeType>;
    using NodesArrayType       = GeometryType::PointsArrayType;
    using VectorType           = Vector;
    using MatrixType           = Matrix;
    using DofsVectorType       = Element::DofsVectorType;
    using EquationIdVectorType = Element::EquationIdVectorType;

    /// The definition of the sizetype
    using SizeType = std::size_t;
    using BaseType::mConstitutiveLawVector;
    using BaseType::mIsInitialised;
    using BaseType::mRetentionLawVector;

    using ElementVariables = typename BaseType::ElementVariables;

    explicit SteadyStatePwElement(IndexType NewId = 0) : BaseType(NewId) {}

    /// Constructor using an array of nodes
    SteadyStatePwElement(IndexType NewId, const NodesArrayType& ThisNodes, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : BaseType(NewId, ThisNodes, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Geometry
    SteadyStatePwElement(IndexType NewId, GeometryType::Pointer pGeometry, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : BaseType(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Properties
    SteadyStatePwElement(IndexType                          NewId,
                         GeometryType::Pointer              pGeometry,
                         PropertiesType::Pointer            pProperties,
                         std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : BaseType(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
    }

    ~SteadyStatePwElement() = default;

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    // Turn back information as a string.
    std::string Info() const override
    {
        const std::string retention_info =
            !mRetentionLawVector.empty() ? mRetentionLawVector[0]->Info() : "not defined";
        return "steady-state Pw flow Element #" + std::to_string(this->Id()) + "\nRetention law: " + retention_info;
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

protected:
    /// Member Variables

    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      bool               CalculateStiffnessMatrixFlag,
                      bool               CalculateResidualVectorFlag) override;

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint) override;

private:
    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override{KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)}

    // Assignment operator.
    SteadyStatePwElement&
    operator=(SteadyStatePwElement const& rOther);

    // Copy constructor.
    SteadyStatePwElement(SteadyStatePwElement const& rOther);

}; // Class SteadyStatePwElement

} // namespace Kratos

#endif // KRATOS_GEO_STEADY_STATE_PW_ELEMENT_H_INCLUDED  defined
