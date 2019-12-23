// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

// System includes

// External includes

// Project includes
#include "includes/define.h"

/* Variables */
#include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos {
KratosContactStructuralMechanicsApplication::KratosContactStructuralMechanicsApplication():
    KratosApplication("ContactStructuralMechanicsApplication"),
    /* CONDITIONS */
    // Mesh tying mortar conditions
    // 2D
    mMeshTyingMortarCondition2D2NTriangle( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mMeshTyingMortarCondition2D2NQuadrilateral( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    // 3D
    mMeshTyingMortarCondition3D3NTetrahedron( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mMeshTyingMortarCondition3D4NHexahedron( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mMeshTyingMortarCondition3D3NTetrahedron4NHexahedron( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mMeshTyingMortarCondition3D4NHexahedron3NTetrahedron( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),

    // ALM Contact mortar conditions
    // Frictionless
    // 2D
    mALMFrictionlessMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mALMNVFrictionlessMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mALMFrictionlessAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mALMNVFrictionlessAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    // 3D
    mALMFrictionlessMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMNVFrictionlessMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMFrictionlessMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMNVFrictionlessMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMFrictionlessMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMNVFrictionlessMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMFrictionlessMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMNVFrictionlessMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    // Frictionless components
    // 2D
    mALMFrictionlessComponentsMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mALMNVFrictionlessComponentsMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    // 3D
    mALMFrictionlessComponentsMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMNVFrictionlessComponentsMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMFrictionlessComponentsMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMNVFrictionlessComponentsMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMFrictionlessComponentsMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMNVFrictionlessComponentsMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMFrictionlessComponentsMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMNVFrictionlessComponentsMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    // Frictional
    // 2D
    mALMFrictionalMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mALMNVFrictionalMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mALMFrictionalAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mALMNVFrictionalAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    // 3D
    mALMFrictionalMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMNVFrictionalMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMFrictionalMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMNVFrictionalMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMFrictionalMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMNVFrictionalMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mALMFrictionalMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mALMNVFrictionalMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    // Penalty method Contact mortar conditions
    // Frictionless
    // 2D
    mPenaltyFrictionlessMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mPenaltyNVFrictionlessMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mPenaltyFrictionlessAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mPenaltyNVFrictionlessAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    // 3D
    mPenaltyFrictionlessMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mPenaltyNVFrictionlessMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mPenaltyFrictionlessMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyNVFrictionlessMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyFrictionlessMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyNVFrictionlessMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyFrictionlessMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mPenaltyNVFrictionlessMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    // Frictional
    // 2D
    mPenaltyFrictionalMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mPenaltyNVFrictionalMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mPenaltyFrictionalAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    mPenaltyNVFrictionalAxisymMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    // 3D
    mPenaltyFrictionalMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mPenaltyNVFrictionalMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mPenaltyFrictionalMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyNVFrictionalMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyFrictionalMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyNVFrictionalMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mPenaltyFrictionalMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mPenaltyNVFrictionalMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    // 2D MPC
    mMPCMortarContactCondition2D2N( 0, GeometryPointerType(new LineType(PointsArrayType(2))), nullptr, GeometryPointerType(new LineType(PointsArrayType(2)))),
    // 3D MPC
    mMPCMortarContactCondition3D3N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    mMPCMortarContactCondition3D4N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mMPCMortarContactCondition3D3N4N( 0, GeometryPointerType(new TriangleType(PointsArrayType(3))), nullptr, GeometryPointerType(new QuadrilateralType(PointsArrayType(4)))),
    mMPCMortarContactCondition3D4N3N( 0, GeometryPointerType(new QuadrilateralType(PointsArrayType(4))), nullptr, GeometryPointerType(new TriangleType(PointsArrayType(3)))),
    // Master-Slave Constraint
    mContactMasterSlaveConstraint()
    {}

void KratosContactStructuralMechanicsApplication::Register()
{
    // Calling base class register to register Kratos components
    KratosApplication::Register();

    // VARIABLES
    // MPC Contact related variables
    KRATOS_REGISTER_VARIABLE( CONSTRAINT_POINTER )                                    // Pointer to the constraint of the condition
    KRATOS_REGISTER_VARIABLE( REACTION_CHECK_STIFFNESS_FACTOR )                       // The reaction factor to be considered on the tension check

    /* Mortar method general variables */
    KRATOS_REGISTER_VARIABLE( INNER_LOOP_ITERATION )                                  // The number of loops in the simplified semi-smooth inner iteration
    KRATOS_REGISTER_VARIABLE( INTEGRATION_ORDER_CONTACT )                             // The integration order computed in the contact
    KRATOS_REGISTER_VARIABLE( DISTANCE_THRESHOLD )                                    // The distance threshold considered
    KRATOS_REGISTER_VARIABLE( ZERO_TOLERANCE_FACTOR )                                 // The epsilon factor considered
    KRATOS_REGISTER_VARIABLE( ACTIVE_CHECK_FACTOR )                                   // The factor employed to serach an active/inactive node
    KRATOS_REGISTER_VARIABLE( SLIP_THRESHOLD )                                        // The threshold employed to search an slip/stick node
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUXILIAR_COORDINATES )               // Auxiliar coordinates used to map
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DELTA_COORDINATES )                  // Delta coordinates used to map


    /* Weighted values */
    KRATOS_REGISTER_VARIABLE( WEIGHTED_GAP )                                          // The integrated gap employed in mortar formulation
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_SLIP )                      // The integrated slip employed in mortar formulation
    KRATOS_REGISTER_VARIABLE( WEIGHTED_SCALAR_RESIDUAL )                              // The integrated scalar residual
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL )           // The integrated vector residual
    KRATOS_REGISTER_VARIABLE( NORMAL_GAP )                                            // The normal gap employed in contact formulation
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( TANGENT_SLIP )                       // The tangent slip employed in contact formulation

    /* For ALM mortar condition */
    KRATOS_REGISTER_VARIABLE( ACTIVE_SET_COMPUTED )                                   // To know if the active set has been computed
    KRATOS_REGISTER_VARIABLE( ACTIVE_SET_CONVERGED )                                  // To know if the active set has converged
    KRATOS_REGISTER_VARIABLE( SLIP_SET_CONVERGED )                                    // To know if the slip set has converged
    KRATOS_REGISTER_VARIABLE( OPERATOR_THRESHOLD )                                    // Consider objetive/non-objetive formulation threshold
    KRATOS_REGISTER_VARIABLE( SLIP_AUGMENTATION_COEFFICIENT )                         // Coefficient to improve the slip computation convergence (augmented part related)
    KRATOS_REGISTER_VARIABLE( DYNAMIC_FACTOR )                                        // The factor considered for dynamic problems (in order to take intro account the gap evolution)
    KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_CONTACT_PRESSURE )                  // The lagrange multiplier for normal contact pressure
    KRATOS_REGISTER_VARIABLE( AUGMENTED_NORMAL_CONTACT_PRESSURE )                     // The resultant augmented pressure in the normal direction
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUGMENTED_TANGENT_CONTACT_PRESSURE ) // The resultant augmented pressure in the tangent direction
    KRATOS_REGISTER_VARIABLE( CONSIDER_NORMAL_VARIATION )                             // A value used to check if consider normal variation or not
    KRATOS_REGISTER_VARIABLE( ADAPT_PENALTY )                                         // To set if the penalty is recalculated or not
    KRATOS_REGISTER_VARIABLE( MAX_GAP_FACTOR )                                        // The factor between the nodal H and the max gap considered to recalculate the penalty

    /* For mesh tying mortar condition */
    KRATOS_REGISTER_VARIABLE( TYING_VARIABLE )                                        // The variable name for the mesh tying

    /* For mesh tying mortar condition */
    KRATOS_REGISTER_VARIABLE( MAX_GAP_THRESHOLD )                                     // The gap considered as threshold to rescale penalty

    /* Frictional laws */
    KRATOS_REGISTER_VARIABLE( FRICTIONAL_LAW )                                        // The frictional law considered
    KRATOS_REGISTER_VARIABLE( TRESCA_FRICTION_THRESHOLD )                             // The threshold value for Tresca frictional contact

    // CONDITIONS
    // Mesh tying mortar condition
    KRATOS_REGISTER_CONDITION( "MeshTyingMortarCondition2D2NTriangle", mMeshTyingMortarCondition2D2NTriangle );
    KRATOS_REGISTER_CONDITION( "MeshTyingMortarCondition2D2NQuadrilateral", mMeshTyingMortarCondition2D2NQuadrilateral );
    KRATOS_REGISTER_CONDITION( "MeshTyingMortarCondition3D3NTetrahedron", mMeshTyingMortarCondition3D3NTetrahedron );
    KRATOS_REGISTER_CONDITION( "MeshTyingMortarCondition3D4NHexahedron", mMeshTyingMortarCondition3D4NHexahedron );
    KRATOS_REGISTER_CONDITION( "MeshTyingMortarCondition3D3NTetrahedron4NHexahedron", mMeshTyingMortarCondition3D3NTetrahedron4NHexahedron );
    KRATOS_REGISTER_CONDITION( "MeshTyingMortarCondition3D4NHexahedron3NTetrahedron", mMeshTyingMortarCondition3D4NHexahedron3NTetrahedron );

    // Mortar contact condition
    // Frictionless cases
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessMortarContactCondition2D2N", mALMFrictionlessMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessMortarContactCondition2D2N", mALMNVFrictionlessMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessAxisymMortarContactCondition2D2N", mALMFrictionlessAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessAxisymMortarContactCondition2D2N", mALMNVFrictionlessAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessMortarContactCondition3D3N", mALMFrictionlessMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessMortarContactCondition3D3N", mALMNVFrictionlessMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessMortarContactCondition3D4N", mALMFrictionlessMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessMortarContactCondition3D4N", mALMNVFrictionlessMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessMortarContactCondition3D3N4N", mALMFrictionlessMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessMortarContactCondition3D3N4N", mALMNVFrictionlessMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessMortarContactCondition3D4N3N", mALMFrictionlessMortarContactCondition3D4N3N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessMortarContactCondition3D4N3N", mALMNVFrictionlessMortarContactCondition3D4N3N );
    // Frictionless components cases
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessComponentsMortarContactCondition2D2N", mALMFrictionlessComponentsMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessComponentsMortarContactCondition2D2N", mALMNVFrictionlessComponentsMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessComponentsMortarContactCondition3D3N", mALMFrictionlessComponentsMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessComponentsMortarContactCondition3D3N", mALMNVFrictionlessComponentsMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessComponentsMortarContactCondition3D4N", mALMFrictionlessComponentsMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessComponentsMortarContactCondition3D4N", mALMNVFrictionlessComponentsMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessComponentsMortarContactCondition3D3N4N", mALMFrictionlessComponentsMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessComponentsMortarContactCondition3D3N4N", mALMNVFrictionlessComponentsMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionlessComponentsMortarContactCondition3D4N3N", mALMFrictionlessComponentsMortarContactCondition3D4N3N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionlessComponentsMortarContactCondition3D4N3N", mALMNVFrictionlessComponentsMortarContactCondition3D4N3N );
    // Frictional cases
    KRATOS_REGISTER_CONDITION( "ALMFrictionalMortarContactCondition2D2N", mALMFrictionalMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionalMortarContactCondition2D2N", mALMNVFrictionalMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionalAxisymMortarContactCondition2D2N", mALMFrictionalAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionalAxisymMortarContactCondition2D2N", mALMNVFrictionalAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionalMortarContactCondition3D3N", mALMFrictionalMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionalMortarContactCondition3D3N", mALMNVFrictionalMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionalMortarContactCondition3D4N", mALMFrictionalMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionalMortarContactCondition3D4N", mALMNVFrictionalMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionalMortarContactCondition3D3N4N", mALMFrictionalMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionalMortarContactCondition3D3N4N", mALMNVFrictionalMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "ALMFrictionalMortarContactCondition3D4N3N", mALMFrictionalMortarContactCondition3D4N3N );
    KRATOS_REGISTER_CONDITION( "ALMNVFrictionalMortarContactCondition3D4N3N", mALMNVFrictionalMortarContactCondition3D4N3N );
    // Frictionless penalty cases
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionlessMortarContactCondition2D2N", mPenaltyFrictionlessMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionlessMortarContactCondition2D2N", mPenaltyNVFrictionlessMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionlessAxisymMortarContactCondition2D2N", mPenaltyFrictionlessAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionlessAxisymMortarContactCondition2D2N", mPenaltyNVFrictionlessAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionlessMortarContactCondition3D3N", mPenaltyFrictionlessMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionlessMortarContactCondition3D3N", mPenaltyNVFrictionlessMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionlessMortarContactCondition3D4N", mPenaltyFrictionlessMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionlessMortarContactCondition3D4N", mPenaltyNVFrictionlessMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionlessMortarContactCondition3D3N4N", mPenaltyFrictionlessMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionlessMortarContactCondition3D3N4N", mPenaltyNVFrictionlessMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionlessMortarContactCondition3D4N3N", mPenaltyFrictionlessMortarContactCondition3D4N3N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionlessMortarContactCondition3D4N3N", mPenaltyNVFrictionlessMortarContactCondition3D4N3N );
    // Frictional penalty cases
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionalMortarContactCondition2D2N", mPenaltyFrictionalMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionalMortarContactCondition2D2N", mPenaltyNVFrictionalMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionalAxisymMortarContactCondition2D2N", mPenaltyFrictionalAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionalAxisymMortarContactCondition2D2N", mPenaltyNVFrictionalAxisymMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionalMortarContactCondition3D3N", mPenaltyFrictionalMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionalMortarContactCondition3D3N", mPenaltyNVFrictionalMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionalMortarContactCondition3D4N", mPenaltyFrictionalMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionalMortarContactCondition3D4N", mPenaltyNVFrictionalMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionalMortarContactCondition3D3N4N", mPenaltyFrictionalMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionalMortarContactCondition3D3N4N", mPenaltyNVFrictionalMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "PenaltyFrictionalMortarContactCondition3D4N3N", mPenaltyFrictionalMortarContactCondition3D4N3N );
    KRATOS_REGISTER_CONDITION( "PenaltyNVFrictionalMortarContactCondition3D4N3N", mPenaltyNVFrictionalMortarContactCondition3D4N3N );
    // MPC conditions
    KRATOS_REGISTER_CONDITION( "MPCMortarContactCondition2D2N", mMPCMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "MPCMortarContactCondition3D3N", mMPCMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "MPCMortarContactCondition3D4N", mMPCMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "MPCMortarContactCondition3D3N4N", mMPCMortarContactCondition3D3N4N );
    KRATOS_REGISTER_CONDITION( "MPCMortarContactCondition3D4N3N", mMPCMortarContactCondition3D4N3N );

    // Register constraints
    KRATOS_REGISTER_CONSTRAINT("ContactMasterSlaveConstraint",mContactMasterSlaveConstraint);
}

}  // namespace Kratos.
