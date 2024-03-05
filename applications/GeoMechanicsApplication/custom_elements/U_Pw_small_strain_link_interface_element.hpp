// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_U_PW_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_U_PW_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_elements/U_Pw_small_strain_interface_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwSmallStrainLinkInterfaceElement
    : public UPwSmallStrainInterfaceElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwSmallStrainLinkInterfaceElement);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;
    using UPwBaseElement<TDim, TNumNodes>::mConstitutiveLawVector;
    using UPwBaseElement<TDim, TNumNodes>::mRetentionLawVector;
    using UPwBaseElement<TDim, TNumNodes>::mStressVector;
    using UPwBaseElement<TDim, TNumNodes>::mThisIntegrationMethod;

    using SFGradAuxVariables = typename UPwSmallStrainInterfaceElement<TDim, TNumNodes>::SFGradAuxVariables;
    using InterfaceElementVariables =
        typename UPwSmallStrainInterfaceElement<TDim, TNumNodes>::InterfaceElementVariables;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwSmallStrainLinkInterfaceElement() : UPwSmallStrainInterfaceElement<TDim, TNumNodes>() {}

    // Constructor 1
    UPwSmallStrainLinkInterfaceElement(IndexType                          NewId,
                                       GeometryType::Pointer              pGeometry,
                                       std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainInterfaceElement<TDim, TNumNodes>(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    // Constructor 2
    UPwSmallStrainLinkInterfaceElement(IndexType                          NewId,
                                       GeometryType::Pointer              pGeometry,
                                       PropertiesType::Pointer            pProperties,
                                       std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainInterfaceElement<TDim, TNumNodes>(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
    }

    // Destructor
    ~UPwSmallStrainLinkInterfaceElement() override {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "U-Pw small strain link interface Element #" << this->Id()
               << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "U-Pw small strain link interface Element #" << this->Id()
                 << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    // Member Variables

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      const bool         CalculateStiffnessMatrixFlag,
                      const bool         CalculateResidualVectorFlag) override;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    // Member Variables

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

}; // Class UPwSmallStrainLinkInterfaceElement

} // namespace Kratos

#endif // KRATOS_GEO_U_PW_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED defined
