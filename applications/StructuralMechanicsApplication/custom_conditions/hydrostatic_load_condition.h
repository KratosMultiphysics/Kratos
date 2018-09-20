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

#if !defined(KRATOS_HYDROSTATIC_LOAD_CONDITION_H_INCLUDED)
#define KRATOS_HYDROSTATIC_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/surface_load_condition_3d.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/divide_geometry.h"
#include "utilities/divide_triangle_2d_3.h"
#include "includes/element.h"
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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) HydrostaticLoadCondition
    : public SurfaceLoadCondition3D
{
  public:
    ///@name Type Definitions
    ///@{

    // Counted pointer of HydrostaticLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION(HydrostaticLoadCondition);

    typedef DivideGeometry::IndexedPointGeometryType IndexedPointGeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    HydrostaticLoadCondition();

    // Constructor using an array of nodes
    HydrostaticLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    // Constructor using an array of nodes with properties
    HydrostaticLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties

    );

    // Destructor
    ~HydrostaticLoadCondition() override;

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
        NodesArrayType const &ThisNodes,
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
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag) override;

    void CalculateAndSubKpSym(
        Matrix &K,
        const array_1d<double, 3> &ge,
        const array_1d<double, 3> &gn,
        const Matrix &DN_De,
        const Vector &N,
        const array_1d<double, 3> &Normal,
        const double Pressure,
        const double Weight);

    void CalculateAndSubKpHydrostatic(
        Matrix &K,
        const Vector &N,
        const array_1d<double, 3> &Normal,
        const double &rSpecificWeight,
        const array_1d<double, 3> &rW,
        const double &Weight);

    void CalculateAndSubKpHydrostaticSym(
        Matrix &K,
        const Vector &N,
        const array_1d<double, 3> &Normal,
        const double &rSpecificWeight,
        const array_1d<double, 3> &rW,
        const double &Weight);

    void CalculateAndSubKpVolume(
        Matrix &K,
        const double &rSpecificWeight,
        const double &rIntersectedArea);

    void DyadicProduct(Matrix &M,
                       const array_1d<double, 3> &U,
                       const array_1d<double, 3> &V);

    unsigned int NumberOfCommonElements(NodeType &rNodeM,
                                        NodeType &rNodeN);

    /*      void IsNegativeOrSplit(
        GeometryType &r_geom,
        bool &r_is_split,
        bool &r_is_negative);  */

    void CalculateAllInSplitAndNegativeDistanceConditions(
        const GeometryType &rSubGeom,
        IntegrationMethod &rIntegrationMethod,
        const array_1d<double, 3> &rGe,
        const array_1d<double, 3> &rGn,
        const array_1d<double, 3> &rNormal,
        Vector &rPressureOnNodes,
        const array_1d<double, 3> &rW,
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag);
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

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseLoadCondition);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseLoadCondition);
    }

}; // class HydrostaticLoadCondition.

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_CONDITION_3D_H_INCLUDED  defined
