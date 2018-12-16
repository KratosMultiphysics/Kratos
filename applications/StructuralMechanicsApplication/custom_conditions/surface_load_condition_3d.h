// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_SURFACE_LOAD_CONDITION_3D_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_CONDITION_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/base_load_condition.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SurfaceLoadCondition3D
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the responsible to add the contributions of the RHS and LHS of the surface loads of the structure
 * @details It allows to consider different types of pressure and surface loads
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  SurfaceLoadCondition3D
    : public BaseLoadCondition
{
public:

    ///@name Type Definitions
    ///@{

    // Counted pointer of SurfaceLoadCondition3D
    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoadCondition3D );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    SurfaceLoadCondition3D();

    // Constructor using an array of nodes
    SurfaceLoadCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    // Constructor using an array of nodes with properties
    SurfaceLoadCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    // Destructor
    ~SurfaceLoadCondition3D() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // Name Operations
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The local LHS contribution
     * @param rRightHandSideVector The local RHS contribution
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    /**
     * @brief This method adds the local contribution of the pressure to the LHS matrix
     * @param rK The local LHS contribution
     * @param rTangentXi The first tangent direction
     * @param rTangentEta The second tangent direction
     * @param rDN_De The local gradient of the geometry
     * @param rN The shape function of the current integration point
     * @param Pressure The pressure to be applied
     * @param Weight The integration contribution
     */
    void CalculateAndSubKp(
        Matrix& rK,
        const array_1d<double, 3>& rTangentXi,
        const array_1d<double, 3>& rTangentEta,
        const Matrix& rDN_De,
        const Vector& rN,
        const double Pressure,
        const double Weight
        );

    /**
     * @brief This method computes the cross product matrix
     * @param rM The matrix to be build
     * @param rU The vector that defines the
     */
    void MakeCrossMatrix(
        BoundedMatrix<double, 3, 3>& rM,
        const array_1d<double, 3>& rU
        );

    /**
     * @brief This method adds the pressure contribution to the RHS
     * @param rResidualVector The local contribution to the RHS
     * @param rN The corresponding shape function
     * @param rNormal The normal to the geometry surface
     * @param Pressure The pressure to be applied
     * @param Weight The integration contribution
     * @param rCurrentProcessInfo The current instance of process info
     */
    void CalculateAndAddPressureForce(
        VectorType& rResidualVector,
        const Vector& rN,
        const array_1d<double, 3 >& rNormal,
        const double Pressure,
        const double Weight,
        const ProcessInfo& rCurrentProcessInfo
        );

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseLoadCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseLoadCondition );
    }


}; // class SurfaceLoadCondition3D.

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_CONDITION_3D_H_INCLUDED  defined
