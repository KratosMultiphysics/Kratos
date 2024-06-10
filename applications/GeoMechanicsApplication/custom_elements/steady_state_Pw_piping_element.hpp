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
#include "custom_elements/transient_Pw_line_element.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) SteadyStatePwPipingElement
    : public TransientPwLineElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SteadyStatePwPipingElement);

    using BaseType = TransientPwLineElement<TDim, TNumNodes>;

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
    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    explicit SteadyStatePwPipingElement(IndexType NewId = 0)
        : TransientPwLineElement<TDim, TNumNodes>(NewId)
    {
    }

    SteadyStatePwPipingElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : TransientPwLineElement(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    SteadyStatePwPipingElement(IndexType                          NewId,
                               GeometryType::Pointer              pGeometry,
                               PropertiesType::Pointer            pProperties)
        : TransientPwLineElement<TDim, TNumNodes>(NewId, pGeometry, pProperties)
    {
    }

    ~SteadyStatePwPipingElement() = default;
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

}; // Class SteadyStatePwPipingElement

} // namespace Kratos

#endif // KRATOS_GEO_STEADY_STATE_PW_PIPING_ELEMENT_H_INCLUDED  defined
