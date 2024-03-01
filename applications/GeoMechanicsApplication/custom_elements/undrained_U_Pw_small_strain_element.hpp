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

#if !defined(KRATOS_GEO_UNDRAINED_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_UNDRAINED_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UndrainedUPwSmallStrainElement
    : public UPwSmallStrainElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UndrainedUPwSmallStrainElement);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;
    /// The definition of the sizetype
    using SizeType = std::size_t;
    using UPwSmallStrainElement<TDim, TNumNodes>::mConstitutiveLawVector;
    using UPwSmallStrainElement<TDim, TNumNodes>::mRetentionLawVector;
    using UPwSmallStrainElement<TDim, TNumNodes>::mStressVector;
    using UPwSmallStrainElement<TDim, TNumNodes>::mStateVariablesFinalized;
    using ElementVariables     = typename UPwSmallStrainElement<TDim, TNumNodes>::ElementVariables;
    using EquationIdVectorType = Element::EquationIdVectorType;
    using DofsVectorType       = Element::DofsVectorType;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UndrainedUPwSmallStrainElement(IndexType NewId = 0)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId)
    {
    }

    /// Constructor using an array of nodes
    UndrainedUPwSmallStrainElement(IndexType                          NewId,
                                   const NodesArrayType&              ThisNodes,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId, ThisNodes, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Geometry
    UndrainedUPwSmallStrainElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Properties
    UndrainedUPwSmallStrainElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   PropertiesType::Pointer            pProperties,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
    }

    /// Destructor
    ~UndrainedUPwSmallStrainElement() override {}

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "undrained small strain Element #" << this->Id()
               << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "undrained small strain Element #" << this->Id()
                 << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    /// Member Variables
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint) override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    /// Assignment operator.
    UndrainedUPwSmallStrainElement& operator=(UndrainedUPwSmallStrainElement const& rOther);

    /// Copy constructor.
    UndrainedUPwSmallStrainElement(UndrainedUPwSmallStrainElement const& rOther);

}; // Class UndrainedUPwSmallStrainElement

} // namespace Kratos

#endif // KRATOS_GEO_UNDRAINED_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED  defined
