// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi,
//                   Aron Noordam
//                   Mohamed Nabi
//

// Application includes
#include "geo_mechanics_application.h"

namespace Kratos
{
// We define the node type
using NodeType = Node;

KratosGeoMechanicsApplication::KratosGeoMechanicsApplication()
    : KratosApplication("GeoMechanicsApplication")
{
}

void KratosGeoMechanicsApplication::Register()
{
    KRATOS_INFO("") << " KRATOS___                             \n"
                    << "     //   ) )                          \n"
                    << "    //         ___      ___            \n"
                    << "   //  ____  //___) ) //   ) )         \n"
                    << "  //    / / //       //   / /          \n"
                    << " ((____/ / ((____   ((___/ /  MECHANICS\n"
                    << " Initializing KratosGeoMechanicsApplication..." << std::endl;

    // Register Elements
    //  transient one-phase flow elements:
    KRATOS_REGISTER_ELEMENT("TransientPwElement2D3N", mTransientPwElement2D3N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement2D4N", mTransientPwElement2D4N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement3D4N", mTransientPwElement3D4N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement3D8N", mTransientPwElement3D8N)

    KRATOS_REGISTER_ELEMENT("TransientPwElement2D6N", mTransientPwElement2D6N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement2D8N", mTransientPwElement2D8N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement2D9N", mTransientPwElement2D9N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement2D10N", mTransientPwElement2D10N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement2D15N", mTransientPwElement2D15N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement3D10N", mTransientPwElement3D10N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement3D20N", mTransientPwElement3D20N)
    KRATOS_REGISTER_ELEMENT("TransientPwElement3D27N", mTransientPwElement3D27N)

    KRATOS_REGISTER_ELEMENT("TransientPwLineElement2D2N", mTransientPwLineElement2D2N)
    KRATOS_REGISTER_ELEMENT("TransientPwLineElement2D3N", mTransientPwLineElement2D3N)
    KRATOS_REGISTER_ELEMENT("TransientPwLineElement2D4N", mTransientPwLineElement2D4N)
    KRATOS_REGISTER_ELEMENT("TransientPwLineElement2D5N", mTransientPwLineElement2D5N)
    KRATOS_REGISTER_ELEMENT("TransientPwLineElement3D2N", mTransientPwLineElement3D2N)
    KRATOS_REGISTER_ELEMENT("TransientPwLineElement3D3N", mTransientPwLineElement3D3N)

    KRATOS_REGISTER_ELEMENT("TransientPwInterfaceElement2D4N", mTransientPwInterfaceElement2D4N)
    KRATOS_REGISTER_ELEMENT("TransientPwInterfaceElement3D6N", mTransientPwInterfaceElement3D6N)
    KRATOS_REGISTER_ELEMENT("TransientPwInterfaceElement3D8N", mTransientPwInterfaceElement3D8N)

    // Steady-State one-phase flow elements:
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement2D3N", mSteadyStatePwElement2D3N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement2D4N", mSteadyStatePwElement2D4N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement3D4N", mSteadyStatePwElement3D4N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement3D8N", mSteadyStatePwElement3D8N)

    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement2D6N", mSteadyStatePwElement2D6N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement2D8N", mSteadyStatePwElement2D8N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement2D9N", mSteadyStatePwElement2D9N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement2D10N", mSteadyStatePwElement2D10N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement2D15N", mSteadyStatePwElement2D15N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement3D10N", mSteadyStatePwElement3D10N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement3D20N", mSteadyStatePwElement3D20N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwElement3D27N", mSteadyStatePwElement3D27N)

    KRATOS_REGISTER_ELEMENT("SteadyStatePwInterfaceElement2D4N", mSteadyStatePwInterfaceElement2D4N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwInterfaceElement3D6N", mSteadyStatePwInterfaceElement3D6N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwInterfaceElement3D8N", mSteadyStatePwInterfaceElement3D8N)

    KRATOS_REGISTER_ELEMENT("SteadyStatePwPipingElement2D4N", mSteadyStatePwPipingElement2D4N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwPipingElement3D6N", mSteadyStatePwPipingElement3D6N)
    KRATOS_REGISTER_ELEMENT("SteadyStatePwPipingElement3D8N", mSteadyStatePwPipingElement3D8N)

    // Small strain elements
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement2D3N", mUPwSmallStrainElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement2D4N", mUPwSmallStrainElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement3D4N", mUPwSmallStrainElement3D4N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement3D8N", mUPwSmallStrainElement3D8N)

    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement2D6N", mUPwSmallStrainElement2D6N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement2D8N", mUPwSmallStrainElement2D8N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement2D9N", mUPwSmallStrainElement2D9N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement2D10N", mUPwSmallStrainElement2D10N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement2D15N", mUPwSmallStrainElement2D15N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement3D10N", mUPwSmallStrainElement3D10N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement3D20N", mUPwSmallStrainElement3D20N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainElement3D27N", mUPwSmallStrainElement3D27N)

    // Drained small strain elements
    KRATOS_REGISTER_ELEMENT("DrainedUPwSmallStrainElement2D3N", mDrainedUPwSmallStrainElement2D3N)
    KRATOS_REGISTER_ELEMENT("DrainedUPwSmallStrainElement2D4N", mDrainedUPwSmallStrainElement2D4N)
    KRATOS_REGISTER_ELEMENT("DrainedUPwSmallStrainElement3D4N", mDrainedUPwSmallStrainElement3D4N)
    KRATOS_REGISTER_ELEMENT("DrainedUPwSmallStrainElement3D8N", mDrainedUPwSmallStrainElement3D8N)

    // Undrained small strain elements
    KRATOS_REGISTER_ELEMENT("UndrainedUPwSmallStrainElement2D3N", mUndrainedUPwSmallStrainElement2D3N)
    KRATOS_REGISTER_ELEMENT("UndrainedUPwSmallStrainElement2D4N", mUndrainedUPwSmallStrainElement2D4N)
    KRATOS_REGISTER_ELEMENT("UndrainedUPwSmallStrainElement3D4N", mUndrainedUPwSmallStrainElement3D4N)
    KRATOS_REGISTER_ELEMENT("UndrainedUPwSmallStrainElement3D8N", mUndrainedUPwSmallStrainElement3D8N)

    // Small strain FIC elements
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainFICElement2D3N", mUPwSmallStrainFICElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainFICElement2D4N", mUPwSmallStrainFICElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainFICElement3D4N", mUPwSmallStrainFICElement3D4N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainFICElement3D8N", mUPwSmallStrainFICElement3D8N)

    // Small strain different order elements
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement2D6N", mSmallStrainUPwDiffOrderElement2D6N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement2D8N", mSmallStrainUPwDiffOrderElement2D8N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement2D9N", mSmallStrainUPwDiffOrderElement2D9N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement2D10N", mSmallStrainUPwDiffOrderElement2D10N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement2D15N", mSmallStrainUPwDiffOrderElement2D15N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement3D10N", mSmallStrainUPwDiffOrderElement3D10N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement3D20N", mSmallStrainUPwDiffOrderElement3D20N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderElement3D27N", mSmallStrainUPwDiffOrderElement3D27N)

    // small strain axisymmtric elements:
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricElement2D3N", mUPwSmallStrainAxisymmetricElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricElement2D4N", mUPwSmallStrainAxisymmetricElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricElement2D6N", mUPwSmallStrainAxisymmetricElement2D6N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricElement2D8N", mUPwSmallStrainAxisymmetricElement2D8N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricElement2D9N", mUPwSmallStrainAxisymmetricElement2D9N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricElement2D10N", mUPwSmallStrainAxisymmetricElement2D10N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricElement2D15N", mUPwSmallStrainAxisymmetricElement2D15N)

    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricFICElement2D3N", mUPwSmallStrainAxisymmetricFICElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainAxisymmetricFICElement2D4N", mUPwSmallStrainAxisymmetricFICElement2D4N)

    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderAxisymmetricElement2D6N", mSmallStrainUPwDiffOrderAxisymmetricElement2D6N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderAxisymmetricElement2D8N", mSmallStrainUPwDiffOrderAxisymmetricElement2D8N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderAxisymmetricElement2D9N", mSmallStrainUPwDiffOrderAxisymmetricElement2D9N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderAxisymmetricElement2D10N", mSmallStrainUPwDiffOrderAxisymmetricElement2D10N)
    KRATOS_REGISTER_ELEMENT("SmallStrainUPwDiffOrderAxisymmetricElement2D15N", mSmallStrainUPwDiffOrderAxisymmetricElement2D15N)

    // Small strain interface elements
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainInterfaceElement2D4N", mUPwSmallStrainInterfaceElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainInterfaceElement3D6N", mUPwSmallStrainInterfaceElement3D6N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainInterfaceElement3D8N", mUPwSmallStrainInterfaceElement3D8N)

    KRATOS_REGISTER_ELEMENT("UPwSmallStrainLinkInterfaceElement2D4N", mUPwSmallStrainLinkInterfaceElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainLinkInterfaceElement3D6N", mUPwSmallStrainLinkInterfaceElement3D6N)
    KRATOS_REGISTER_ELEMENT("UPwSmallStrainLinkInterfaceElement3D8N", mUPwSmallStrainLinkInterfaceElement3D8N)

    KRATOS_REGISTER_ELEMENT("Geo_ULineInterfacePlaneStrainElement2Plus2N", mULineInterfacePlaneStrainElement2Plus2N)
    KRATOS_REGISTER_ELEMENT("Geo_ULineInterfacePlaneStrainElement3Plus3N", mULineInterfacePlaneStrainElement3Plus3N)

    // Updated-Lagrangian elements
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement2D3N", mUPwUpdatedLagrangianElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement2D4N", mUPwUpdatedLagrangianElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement3D4N", mUPwUpdatedLagrangianElement3D4N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement3D8N", mUPwUpdatedLagrangianElement3D8N)

    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement2D6N", mUPwUpdatedLagrangianElement2D6N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement2D8N", mUPwUpdatedLagrangianElement2D8N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement2D9N", mUPwUpdatedLagrangianElement2D9N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement2D10N", mUPwUpdatedLagrangianElement2D10N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement2D15N", mUPwUpdatedLagrangianElement2D15N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement3D10N", mUPwUpdatedLagrangianElement3D10N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement3D20N", mUPwUpdatedLagrangianElement3D20N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianElement3D27N", mUPwUpdatedLagrangianElement3D27N)

    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianFICElement2D3N", mUPwUpdatedLagrangianFICElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianFICElement2D4N", mUPwUpdatedLagrangianFICElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianFICElement3D4N", mUPwUpdatedLagrangianFICElement3D4N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianFICElement3D8N", mUPwUpdatedLagrangianFICElement3D8N)

    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement2D6N", mUpdatedLagrangianUPwDiffOrderElement2D6N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement2D8N", mUpdatedLagrangianUPwDiffOrderElement2D8N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement2D9N", mUpdatedLagrangianUPwDiffOrderElement2D9N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement2D10N", mUpdatedLagrangianUPwDiffOrderElement2D10N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement2D15N", mUpdatedLagrangianUPwDiffOrderElement2D15N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement3D10N", mUpdatedLagrangianUPwDiffOrderElement3D10N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement3D20N", mUpdatedLagrangianUPwDiffOrderElement3D20N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderElement3D27N", mUpdatedLagrangianUPwDiffOrderElement3D27N)

    // Updated-Lagrangian axisymmetric elements
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricElement2D3N", mUPwUpdatedLagrangianAxisymmetricElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricElement2D4N", mUPwUpdatedLagrangianAxisymmetricElement2D4N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricElement2D6N", mUPwUpdatedLagrangianAxisymmetricElement2D6N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricElement2D8N", mUPwUpdatedLagrangianAxisymmetricElement2D8N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricElement2D9N", mUPwUpdatedLagrangianAxisymmetricElement2D9N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricElement2D10N", mUPwUpdatedLagrangianAxisymmetricElement2D10N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricElement2D15N", mUPwUpdatedLagrangianAxisymmetricElement2D15N)

    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderAxisymmetricElement2D6N",
                            mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D6N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderAxisymmetricElement2D8N",
                            mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D8N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderAxisymmetricElement2D9N",
                            mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D9N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderAxisymmetricElement2D10N",
                            mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D10N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianUPwDiffOrderAxisymmetricElement2D15N",
                            mUpdatedLagrangianUPwDiffOrderAxisymmetricElement2D15N)

    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricFICElement2D3N", mUPwUpdatedLagrangianAxisymmetricFICElement2D3N)
    KRATOS_REGISTER_ELEMENT("UPwUpdatedLagrangianAxisymmetricFICElement2D4N", mUPwUpdatedLagrangianAxisymmetricFICElement2D4N)

    // Register geo structural elements
    KRATOS_REGISTER_ELEMENT("GeoTrussElement2D2N", mGeoTrussElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoTrussElement3D2N", mGeoTrussElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoLinearTrussElement2D2N", mGeoLinearTrussElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoLinearTrussElement3D2N", mGeoLinearTrussElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoCableElement2D2N", mGeoCableElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoCableElement3D2N", mGeoCableElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElement3D2N", mGeoCrBeamElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElement2D2N", mGeoCrBeamElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElementLinear2D2N", mGeoCrBeamElementLinear2D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElementLinear3D2N", mGeoCrBeamElementLinear3D2N)
    KRATOS_REGISTER_ELEMENT("GeoCurvedBeamElement2D3N", mGeoCurvedBeamElement2D3N)

    // Register thermal elements
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement2D3N", mTransientThermalElement2D3N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement2D4N", mTransientThermalElement2D4N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement3D4N", mTransientThermalElement3D4N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement3D8N", mTransientThermalElement3D8N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement2D6N", mTransientThermalElement2D6N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement2D8N", mTransientThermalElement2D8N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement2D9N", mTransientThermalElement2D9N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement2D10N", mTransientThermalElement2D10N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement2D15N", mTransientThermalElement2D15N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement3D10N", mTransientThermalElement3D10N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement3D20N", mTransientThermalElement3D20N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalElement3D27N", mTransientThermalElement3D27N)

    KRATOS_REGISTER_ELEMENT("GeoTransientThermalLineElement2D2N", mTransientThermalLineElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalLineElement2D3N", mTransientThermalLineElement2D3N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalLineElement2D4N", mTransientThermalLineElement2D4N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalLineElement2D5N", mTransientThermalLineElement2D5N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalLineElement3D2N", mTransientThermalLineElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoTransientThermalLineElement3D3N", mTransientThermalLineElement3D3N)

    // Register Conditions
    KRATOS_REGISTER_CONDITION("UPwForceCondition2D1N", mUPwForceCondition2D1N)
    KRATOS_REGISTER_CONDITION("UPwForceCondition3D1N", mUPwForceCondition3D1N)
    KRATOS_REGISTER_CONDITION("UPwFaceLoadCondition2D2N", mUPwFaceLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("UPwFaceLoadCondition2D3N", mUPwFaceLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("UPwFaceLoadCondition2D4N", mUPwFaceLoadCondition2D4N)
    KRATOS_REGISTER_CONDITION("UPwFaceLoadCondition2D5N", mUPwFaceLoadCondition2D5N)
    KRATOS_REGISTER_CONDITION("UPwFaceLoadCondition3D3N", mUPwFaceLoadCondition3D3N)
    KRATOS_REGISTER_CONDITION("UPwFaceLoadCondition3D4N", mUPwFaceLoadCondition3D4N)
    KRATOS_REGISTER_CONDITION("UPwNormalFaceLoadCondition2D2N", mUPwNormalFaceLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("UPwNormalFaceLoadCondition2D3N", mUPwNormalFaceLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("UPwNormalFaceLoadCondition2D4N", mUPwNormalFaceLoadCondition2D4N)
    KRATOS_REGISTER_CONDITION("UPwNormalFaceLoadCondition2D5N", mUPwNormalFaceLoadCondition2D5N)
    KRATOS_REGISTER_CONDITION("UPwNormalFaceLoadCondition3D3N", mUPwNormalFaceLoadCondition3D3N)
    KRATOS_REGISTER_CONDITION("UpwNormalFaceLoadCondition3D4N", mUPwNormalFaceLoadCondition3D4N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxCondition2D2N", mUPwNormalFluxCondition2D2N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxCondition2D3N", mUPwNormalFluxCondition2D3N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxCondition2D4N", mUPwNormalFluxCondition2D4N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxCondition2D5N", mUPwNormalFluxCondition2D5N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxCondition3D3N", mUPwNormalFluxCondition3D3N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxCondition3D4N", mUPwNormalFluxCondition3D4N)
    KRATOS_REGISTER_CONDITION("PwNormalFluxCondition2D2N", mPwNormalFluxCondition2D2N)
    KRATOS_REGISTER_CONDITION("PwNormalFluxCondition2D3N", mPwNormalFluxCondition2D3N)
    KRATOS_REGISTER_CONDITION("PwNormalFluxCondition2D4N", mPwNormalFluxCondition2D4N)
    KRATOS_REGISTER_CONDITION("PwNormalFluxCondition2D5N", mPwNormalFluxCondition2D5N)
    KRATOS_REGISTER_CONDITION("PwNormalFluxCondition3D3N", mPwNormalFluxCondition3D3N)
    KRATOS_REGISTER_CONDITION("PwNormalFluxCondition3D4N", mPwNormalFluxCondition3D4N)

    KRATOS_REGISTER_CONDITION("PwPointFluxCondition2D1N", mPwPointFluxCondition2D1N)
    KRATOS_REGISTER_CONDITION("PwPointFluxCondition3D1N", mPwPointFluxCondition3D1N)

    KRATOS_REGISTER_CONDITION("UPwFaceLoadInterfaceCondition2D2N", mUPwFaceLoadInterfaceCondition2D2N)
    KRATOS_REGISTER_CONDITION("UPwFaceLoadInterfaceCondition3D4N", mUPwFaceLoadInterfaceCondition3D4N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxInterfaceCondition2D2N", mUPwNormalFluxInterfaceCondition2D2N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxInterfaceCondition3D4N", mUPwNormalFluxInterfaceCondition3D4N)

    KRATOS_REGISTER_CONDITION("UPwNormalFluxFICCondition2D2N", mUPwNormalFluxFICCondition2D2N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxFICCondition3D3N", mUPwNormalFluxFICCondition3D3N)
    KRATOS_REGISTER_CONDITION("UPwNormalFluxFICCondition3D4N", mUPwNormalFluxFICCondition3D4N)

    KRATOS_REGISTER_CONDITION("LineLoadDiffOrderCondition2D3N", mLineLoadDiffOrderCondition2D3N)
    KRATOS_REGISTER_CONDITION("LineLoadDiffOrderCondition2D4N", mLineLoadDiffOrderCondition2D4N)
    KRATOS_REGISTER_CONDITION("LineLoadDiffOrderCondition2D5N", mLineLoadDiffOrderCondition2D5N)
    KRATOS_REGISTER_CONDITION("LineNormalLoadDiffOrderCondition2D3N", mLineNormalLoadDiffOrderCondition2D3N)
    KRATOS_REGISTER_CONDITION("LineNormalLoadDiffOrderCondition2D4N", mLineNormalLoadDiffOrderCondition2D4N)
    KRATOS_REGISTER_CONDITION("LineNormalLoadDiffOrderCondition2D5N", mLineNormalLoadDiffOrderCondition2D5N)
    KRATOS_REGISTER_CONDITION("LineNormalFluidFluxDiffOrderCondition2D3N", mLineNormalFluidFluxDiffOrderCondition2D3N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadDiffOrderCondition3D6N", mSurfaceLoadDiffOrderCondition3D6N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadDiffOrderCondition3D8N", mSurfaceLoadDiffOrderCondition3D8N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadDiffOrderCondition3D9N", mSurfaceLoadDiffOrderCondition3D9N)
    KRATOS_REGISTER_CONDITION("SurfaceNormalLoadDiffOrderCondition3D6N", mSurfaceNormalLoadDiffOrderCondition3D6N)
    KRATOS_REGISTER_CONDITION("SurfaceNormalLoadDiffOrderCondition3D8N", mSurfaceNormalLoadDiffOrderCondition3D8N)
    KRATOS_REGISTER_CONDITION("SurfaceNormalLoadDiffOrderCondition3D9N", mSurfaceNormalLoadDiffOrderCondition3D9N)
    KRATOS_REGISTER_CONDITION("SurfaceNormalFluidFluxDiffOrderCondition3D6N", mSurfaceNormalFluidFluxDiffOrderCondition3D6N)
    KRATOS_REGISTER_CONDITION("SurfaceNormalFluidFluxDiffOrderCondition3D8N", mSurfaceNormalFluidFluxDiffOrderCondition3D8N)
    KRATOS_REGISTER_CONDITION("SurfaceNormalFluidFluxDiffOrderCondition3D9N", mSurfaceNormalFluidFluxDiffOrderCondition3D9N)

    KRATOS_REGISTER_CONDITION("AxisymmetricUPwNormalFaceLoadCondition2D2N", mAxisymmetricUPwNormalFaceLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("AxisymmetricUPwNormalFaceLoadCondition2D3N", mAxisymmetricUPwNormalFaceLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("AxisymmetricUPwNormalFaceLoadCondition2D4N", mAxisymmetricUPwNormalFaceLoadCondition2D4N)
    KRATOS_REGISTER_CONDITION("AxisymmetricUPwNormalFaceLoadCondition2D5N", mAxisymmetricUPwNormalFaceLoadCondition2D5N)

    KRATOS_REGISTER_CONDITION("AxisymmetricLineNormalLoadDiffOrderCondition2D3N",
                              mAxisymmetricLineNormalLoadDiffOrderCondition2D3N)
    KRATOS_REGISTER_CONDITION("AxisymmetricLineNormalLoadDiffOrderCondition2D4N",
                              mAxisymmetricLineNormalLoadDiffOrderCondition2D4N)
    KRATOS_REGISTER_CONDITION("AxisymmetricLineNormalLoadDiffOrderCondition2D5N",
                              mAxisymmetricLineNormalLoadDiffOrderCondition2D5N)

    KRATOS_REGISTER_CONDITION("AxisymmetricLineNormalFluidFluxDiffOrderCondition2D3N",
                              mAxisymmetricLineNormalFluidFluxDiffOrderCondition2D3N)

    KRATOS_REGISTER_CONDITION("UPwLysmerAbsorbingCondition2D2N", mUPwLysmerAbsorbingCondition2D2N)
    KRATOS_REGISTER_CONDITION("UPwLysmerAbsorbingCondition2D3N", mUPwLysmerAbsorbingCondition2D3N)
    KRATOS_REGISTER_CONDITION("UPwLysmerAbsorbingCondition3D3N", mUPwLysmerAbsorbingCondition3D3N)
    KRATOS_REGISTER_CONDITION("UPwLysmerAbsorbingCondition3D4N", mUPwLysmerAbsorbingCondition3D4N)

    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition2D2N", mGeoTNormalFluxCondition2D2N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition2D3N", mGeoTNormalFluxCondition2D3N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition2D4N", mGeoTNormalFluxCondition2D4N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition2D5N", mGeoTNormalFluxCondition2D5N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition3D3N", mGeoTNormalFluxCondition3D3N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition3D4N", mGeoTNormalFluxCondition3D4N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition3D6N", mGeoTNormalFluxCondition3D6N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition3D8N", mGeoTNormalFluxCondition3D8N)
    KRATOS_REGISTER_CONDITION("GeoTNormalFluxCondition3D9N", mGeoTNormalFluxCondition3D9N)

    KRATOS_REGISTER_CONDITION("GeoThermalPointFluxCondition2D1N", mGeoThermalPointFluxCondition2D1N)
    KRATOS_REGISTER_CONDITION("GeoThermalPointFluxCondition3D1N", mGeoThermalPointFluxCondition3D1N)

    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition2D2N", mGeoTMicroClimateFluxCondition2D2N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition2D3N", mGeoTMicroClimateFluxCondition2D3N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition2D4N", mGeoTMicroClimateFluxCondition2D4N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition2D5N", mGeoTMicroClimateFluxCondition2D5N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition3D3N", mGeoTMicroClimateFluxCondition3D3N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition3D4N", mGeoTMicroClimateFluxCondition3D4N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition3D6N", mGeoTMicroClimateFluxCondition3D6N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition3D8N", mGeoTMicroClimateFluxCondition3D8N)
    KRATOS_REGISTER_CONDITION("GeoTMicroClimateFluxCondition3D9N", mGeoTMicroClimateFluxCondition3D9N)

    // Register Constitutive Laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive3DLaw", mBilinearCohesive3DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive2DLaw", mBilinearCohesive2DLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStrainK02DLaw", mLinearPlaneStrainK0Law)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticK03DLaw", mElasticIsotropicK03DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("GeoLinearElasticPlaneStrain2DLaw", mLinearElasticPlaneStrain2DLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("GeoLinearElasticPlaneStress2DLaw", mLinearElasticPlaneStress2DLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM3DLaw", mSmallStrainUDSM3DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM2DPlaneStrainLaw", mSmallStrainUDSM2DPlaneStrainLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM2DInterfaceLaw", mSmallStrainUDSM2DInterfaceLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM3DInterfaceLaw", mSmallStrainUDSM3DInterfaceLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT3DLaw", mSmallStrainUMAT3DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT2DPlaneStrainLaw", mSmallStrainUMAT2DPlaneStrainLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT2DInterfaceLaw", mSmallStrainUMAT2DInterfaceLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT3DInterfaceLaw", mSmallStrainUMAT3DInterfaceLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic2DInterfaceLaw", mLinearElastic2DInterfaceLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic3DInterfaceLaw", mLinearElastic3DInterfaceLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic2DBeamLaw", mLinearElastic2DBeamLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("TrussBackboneConstitutiveLaw", mTrussBackboneConstitutiveLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("Geo_IncrementalLinearElasticInterfaceLaw", mIncrementalLinearElasticInterfaceLaw)

    // Register Variables
    KRATOS_REGISTER_VARIABLE(VELOCITY_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(DT_PRESSURE_COEFFICIENT)

    KRATOS_REGISTER_VARIABLE(DT_WATER_PRESSURE)
    KRATOS_REGISTER_VARIABLE(NORMAL_FLUID_FLUX)

    KRATOS_REGISTER_VARIABLE(HYDRAULIC_HEAD)

    KRATOS_REGISTER_VARIABLE(HYDRAULIC_DISCHARGE)

    KRATOS_REGISTER_VARIABLE(DENSITY_SOLID)
    KRATOS_REGISTER_VARIABLE(BULK_MODULUS_SOLID)
    KRATOS_REGISTER_VARIABLE(BULK_MODULUS_FLUID)

    KRATOS_REGISTER_VARIABLE(SPECIFIC_HEAT_CAPACITY_WATER)
    KRATOS_REGISTER_VARIABLE(SPECIFIC_HEAT_CAPACITY_SOLID)
    KRATOS_REGISTER_VARIABLE(THERMAL_CONDUCTIVITY_WATER)
    KRATOS_REGISTER_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(THERMAL_CONDUCTIVITY_SOLID)
    KRATOS_REGISTER_VARIABLE(SOLID_COMPRESSIBILITY)
    KRATOS_REGISTER_VARIABLE(DT_TEMPERATURE_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(DT_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(NORMAL_HEAT_FLUX)
    KRATOS_REGISTER_VARIABLE(THERMAL_LAW_NAME)

    // Variables for Micro-Climate boundary
    KRATOS_REGISTER_VARIABLE(AIR_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(SOLAR_RADIATION)
    KRATOS_REGISTER_VARIABLE(AIR_HUMIDITY)
    KRATOS_REGISTER_VARIABLE(PRECIPITATION)
    KRATOS_REGISTER_VARIABLE(WIND_SPEED)
    KRATOS_REGISTER_VARIABLE(A1_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(A2_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(A3_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(ALPHA_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(QF_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(SMIN_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(SMAX_COEFFICIENT)

    KRATOS_REGISTER_VARIABLE(K0_MAIN_DIRECTION)
    KRATOS_REGISTER_VARIABLE(K0_VALUE_XX)
    KRATOS_REGISTER_VARIABLE(K0_VALUE_YY)
    KRATOS_REGISTER_VARIABLE(K0_VALUE_ZZ)
    KRATOS_REGISTER_VARIABLE(K0_NC)
    KRATOS_REGISTER_VARIABLE(OCR)
    KRATOS_REGISTER_VARIABLE(POISSON_UNLOADING_RELOADING)
    KRATOS_REGISTER_VARIABLE(POP)

    KRATOS_REGISTER_VARIABLE(ACCUMULATED_STRAIN)

    KRATOS_REGISTER_VARIABLE(PERMEABILITY_XX)
    KRATOS_REGISTER_VARIABLE(PERMEABILITY_YY)
    KRATOS_REGISTER_VARIABLE(PERMEABILITY_ZZ)
    KRATOS_REGISTER_VARIABLE(PERMEABILITY_XY)
    KRATOS_REGISTER_VARIABLE(PERMEABILITY_YZ)
    KRATOS_REGISTER_VARIABLE(PERMEABILITY_ZX)
    KRATOS_REGISTER_VARIABLE(PERMEABILITY_CHANGE_INVERSE_FACTOR)

    KRATOS_REGISTER_VARIABLE(MINIMUM_JOINT_WIDTH)
    KRATOS_REGISTER_VARIABLE(TRANSVERSAL_PERMEABILITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FLUID_FLUX_VECTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LOCAL_FLUID_FLUX_VECTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LOCAL_STRESS_VECTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LOCAL_RELATIVE_DISPLACEMENT_VECTOR)
    KRATOS_REGISTER_VARIABLE(PERMEABILITY_MATRIX)
    KRATOS_REGISTER_VARIABLE(LOCAL_PERMEABILITY_MATRIX)

    KRATOS_REGISTER_VARIABLE(CRITICAL_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TOTAL_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(INCREMENTAL_DISPLACEMENT)

    KRATOS_REGISTER_VARIABLE(IS_CONVERGED)

    KRATOS_REGISTER_VARIABLE(TOTAL_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(TOTAL_STRESS_VECTOR)

    KRATOS_REGISTER_VARIABLE(CAUCHY_STRAIN_TENSOR)
    KRATOS_REGISTER_VARIABLE(CAUCHY_STRAIN_VECTOR)

    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE)
    KRATOS_REGISTER_VARIABLE(ARC_LENGTH_LAMBDA)
    KRATOS_REGISTER_VARIABLE(ARC_LENGTH_RADIUS_FACTOR)

    KRATOS_REGISTER_VARIABLE(TIME_UNIT_CONVERTER)

    KRATOS_REGISTER_VARIABLE(LOCAL_EQUIVALENT_STRAIN)
    KRATOS_REGISTER_VARIABLE(NONLOCAL_EQUIVALENT_STRAIN)

    KRATOS_REGISTER_VARIABLE(JOINT_WIDTH)

    KRATOS_REGISTER_VARIABLE(NODAL_SMOOTHING)
    KRATOS_REGISTER_VARIABLE(NODAL_CAUCHY_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(ENGINEERING_STRAIN_TENSOR)
    KRATOS_REGISTER_VARIABLE(ENGINEERING_STRAIN_VECTOR)
    KRATOS_REGISTER_VARIABLE(NODAL_DAMAGE_VARIABLE)
    KRATOS_REGISTER_VARIABLE(NODAL_JOINT_AREA)
    KRATOS_REGISTER_VARIABLE(NODAL_JOINT_WIDTH)
    KRATOS_REGISTER_VARIABLE(NODAL_JOINT_DAMAGE)

    KRATOS_REGISTER_VARIABLE(BIOT_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(PLATE_SHAPE_CORRECTION_FACTOR)

    KRATOS_REGISTER_VARIABLE(RESET_DISPLACEMENTS)
    KRATOS_REGISTER_VARIABLE(CONSIDER_GEOMETRIC_STIFFNESS)

    KRATOS_REGISTER_VARIABLE(CONSIDER_GAP_CLOSURE)

    KRATOS_REGISTER_VARIABLE(USE_CONSISTENT_MASS_MATRIX)

    KRATOS_REGISTER_VARIABLE(IGNORE_UNDRAINED)
    KRATOS_REGISTER_VARIABLE(USE_HENCKY_STRAIN)

    KRATOS_REGISTER_VARIABLE(MEAN_EFFECTIVE_STRESS)
    KRATOS_REGISTER_VARIABLE(MEAN_STRESS)
    KRATOS_REGISTER_VARIABLE(ENGINEERING_VOLUMETRIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(ENGINEERING_VON_MISES_STRAIN)
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_VOLUMETRIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_VON_MISES_STRAIN)

    KRATOS_REGISTER_VARIABLE(SATURATED_SATURATION)
    KRATOS_REGISTER_VARIABLE(RESIDUAL_SATURATION)
    KRATOS_REGISTER_VARIABLE(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE)
    KRATOS_REGISTER_VARIABLE(VAN_GENUCHTEN_GN)
    KRATOS_REGISTER_VARIABLE(VAN_GENUCHTEN_GL)
    KRATOS_REGISTER_VARIABLE(MINIMUM_RELATIVE_PERMEABILITY)

    KRATOS_REGISTER_VARIABLE(RETENTION_LAW)
    KRATOS_REGISTER_VARIABLE(DEGREE_OF_SATURATION)
    KRATOS_REGISTER_VARIABLE(EFFECTIVE_SATURATION)
    KRATOS_REGISTER_VARIABLE(BISHOP_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(DERIVATIVE_OF_SATURATION)
    KRATOS_REGISTER_VARIABLE(RELATIVE_PERMEABILITY)

    KRATOS_REGISTER_VARIABLE(ABSORBING_FACTORS)
    KRATOS_REGISTER_VARIABLE(VIRTUAL_THICKNESS)

    KRATOS_REGISTER_VARIABLE(CONFINED_STIFFNESS)
    KRATOS_REGISTER_VARIABLE(SHEAR_STIFFNESS)

    KRATOS_REGISTER_VARIABLE(IS_PIPING_CONVERGED)
    KRATOS_REGISTER_VARIABLE(PIPE_ETA)
    KRATOS_REGISTER_VARIABLE(PIPE_THETA)
    KRATOS_REGISTER_VARIABLE(PIPE_D_70)
    KRATOS_REGISTER_VARIABLE(PIPE_START_ELEMENT)
    KRATOS_REGISTER_VARIABLE(PIPE_ELEMENT_LENGTH)
    KRATOS_REGISTER_VARIABLE(PIPE_IN_EQUILIBRIUM)
    KRATOS_REGISTER_VARIABLE(PIPE_MODIFIED_D)
    KRATOS_REGISTER_VARIABLE(PIPE_MODEL_FACTOR)
    KRATOS_REGISTER_VARIABLE(PIPE_HEIGHT)
    KRATOS_REGISTER_VARIABLE(PREV_PIPE_HEIGHT)
    KRATOS_REGISTER_VARIABLE(DIFF_PIPE_HEIGHT)
    KRATOS_REGISTER_VARIABLE(PIPE_EROSION)
    KRATOS_REGISTER_VARIABLE(PIPE_ACTIVE)

    // UDSM
    KRATOS_REGISTER_VARIABLE(UDSM_NAME) // Also for UMAT
    KRATOS_REGISTER_VARIABLE(UDSM_NUMBER)
    KRATOS_REGISTER_VARIABLE(IS_FORTRAN_UDSM) // Also for UMAT

    KRATOS_REGISTER_VARIABLE(UMAT_PARAMETERS)

    KRATOS_REGISTER_VARIABLE(NUMBER_OF_UMAT_STATE_VARIABLES)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_UMAT_PARAMETERS)
    KRATOS_REGISTER_VARIABLE(INDEX_OF_UMAT_C_PARAMETER)
    KRATOS_REGISTER_VARIABLE(INDEX_OF_UMAT_PHI_PARAMETER)

    KRATOS_REGISTER_VARIABLE(STATE_VARIABLES)

    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_1)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_2)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_3)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_4)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_5)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_6)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_7)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_8)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_9)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_10)

    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_11)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_12)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_13)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_14)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_15)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_16)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_17)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_18)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_19)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_20)

    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_21)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_22)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_23)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_24)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_25)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_26)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_27)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_28)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_29)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_30)

    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_31)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_32)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_33)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_34)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_35)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_36)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_37)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_38)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_39)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_40)

    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_41)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_42)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_43)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_44)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_45)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_46)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_47)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_48)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_49)
    KRATOS_REGISTER_VARIABLE(STATE_VARIABLE_50)

    KRATOS_REGISTER_VARIABLE(STRAINS_OF_PIECEWISE_LINEAR_LAW)
    KRATOS_REGISTER_VARIABLE(STRESSES_OF_PIECEWISE_LINEAR_LAW)

    KRATOS_REGISTER_VARIABLE(INTERFACE_NORMAL_STIFFNESS)
    KRATOS_REGISTER_VARIABLE(INTERFACE_SHEAR_STIFFNESS)
}
} // namespace Kratos.
