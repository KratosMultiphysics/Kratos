// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#if !defined(KRATOS_GEO_STEADY_STATE_PW_PIPING_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_STEADY_STATE_PW_PIPING_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/steady_state_Pw_interface_element.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) SteadyStatePwPipingElement
    : public SteadyStatePwInterfaceElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SteadyStatePwPipingElement);

    using BaseType = SteadyStatePwInterfaceElement<TDim, TNumNodes>;

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;

    using DofsVectorType       = Element::DofsVectorType;
    using EquationIdVectorType = Element::EquationIdVectorType;

    /// The definition of the sizetype
    using SizeType = std::size_t;

    using BaseType::CalculateRetentionResponse;
    using BaseType::mRetentionLawVector;
    using BaseType::mThisIntegrationMethod;

    using InterfaceElementVariables = typename BaseType::InterfaceElementVariables;
    using SFGradAuxVariables        = typename BaseType::SFGradAuxVariables;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    SteadyStatePwPipingElement(IndexType NewId = 0)
        : SteadyStatePwInterfaceElement<TDim, TNumNodes>(NewId)
    {
    }

    /// Constructor using an array of nodes
    SteadyStatePwPipingElement(IndexType                          NewId,
                               const NodesArrayType&              ThisNodes,
                               std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : SteadyStatePwInterfaceElement<TDim, TNumNodes>(NewId, ThisNodes, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Geometry
    SteadyStatePwPipingElement(IndexType                          NewId,
                               GeometryType::Pointer              pGeometry,
                               std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : SteadyStatePwInterfaceElement<TDim, TNumNodes>(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Properties
    SteadyStatePwPipingElement(IndexType                          NewId,
                               GeometryType::Pointer              pGeometry,
                               PropertiesType::Pointer            pProperties,
                               std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : SteadyStatePwInterfaceElement<TDim, TNumNodes>(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
    }

    /// Destructor
    ~SteadyStatePwPipingElement() override {}

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    bool InEquilibrium(const PropertiesType& Prop, const GeometryType& Geom);

    double CalculateHeadGradient(const PropertiesType& Prop, const GeometryType& Geom, double pipe_length);

    double CalculateEquilibriumPipeHeight(const PropertiesType& Prop, const GeometryType& Geom, double dx);

    void CalculateLength(const GeometryType& Geom);

protected:
    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      const bool         CalculateStiffnessMatrixFlag,
                      const bool         CalculateResidualVectorFlag) override;

    using BaseType::CalculateOnIntegrationPoints;
    void CalculateOnIntegrationPoints(const Variable<bool>& rVariable,
                                      std::vector<bool>&    rValues,
                                      const ProcessInfo&    rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    double CalculateParticleDiameter(const PropertiesType& Prop);

    double pipe_initialised = false;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    /// Member Variables

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Assignment operator.
    SteadyStatePwPipingElement& operator=(SteadyStatePwPipingElement const& rOther);

    /// Copy constructor.
    SteadyStatePwPipingElement(SteadyStatePwPipingElement const& rOther);

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

}; // Class SteadyStatePwInterfaceElement

} // namespace Kratos

#endif // KRATOS_GEO_STEADY_STATE_PW_PIPING_ELEMENT_H_INCLUDED  defined
