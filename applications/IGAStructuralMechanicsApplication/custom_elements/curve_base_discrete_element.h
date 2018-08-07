/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//
//  Authors: Tobias Teschemacher
*/

#if !defined(KRATOS_CURVE_BASE_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_CURVE_BASE_DISCRETE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"

// External includes

// Project includes
#include "custom_elements/base_discrete_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** CurveBaseDiscreteElement deals as base class for curve structured element formulations.
*/
class  CurveBaseDiscreteElement
    : public BaseDiscreteElement
{
protected:

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "CurveBaseDiscreteElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CurveBaseDiscreteElement #" << Id();
    }

public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of CurveBaseDiscreteElement
    KRATOS_CLASS_POINTER_DEFINITION(CurveBaseDiscreteElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
     // Constructor using an array of nodes
    CurveBaseDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteElement(NewId, pGeometry)
    {};
     // Constructor using an array of nodes with properties
    CurveBaseDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : BaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    CurveBaseDiscreteElement() : BaseDiscreteElement()
    {};

    /// Destructor.
    virtual ~CurveBaseDiscreteElement() override
    {};

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_ERROR << "Trying to create a \"CurveBaseDiscreteElement\"" << std::endl;
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * Called to initialize the element.
    * Must be called before any calculation is done
    */
    void Initialize() override;


    ///@}
protected:
    ///@name Static Member Variables
    ///@{
    Vector mBaseVector0;
    ///@}
    ///@name Operations
    ///@{

    void GetBaseVector(
        Vector& rBaseVector, 
        const Matrix& rDN_De);

    /**
    * GetBoundaryEdgeBaseVector computes t3 of the boundary edge
    * @param DN_De derivatives of shape functions.
    * @param Tangents in Parameter space
    * @see rBaseVector t3 of the edge
    */
    void GetBoundaryEdgeBaseVector(const Matrix& DN_De,
        const array_1d<double, 2>& Tangents,
        Vector& rBaseVector);

    void Get1stVariationsAxialStrain(
        Vector& rEpsilon1stVariationDoF,
        const Vector& rBaseVector,
        const int& rNumberOfDoFs, 
        const Matrix& rDN_De);

    void Get2ndVariationsAxialStrain(
        Matrix& rEpsilon2ndVariationDoF,
        const int& rNumberOfDoFs, 
        const Matrix& rDN_De);

private:
    ///@name Operations
    ///@{

    ///@}

    ///@name Static Member Variables

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseDiscreteElement)
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseDiscreteElement)
    }
    ///@}
};     // Class CurveBaseDiscreteElement
///@}
}  // namespace Kratos.

#endif // KRATOS_CURVE_BASE_DISCRETE_ELEMENT_H_INCLUDED  defined