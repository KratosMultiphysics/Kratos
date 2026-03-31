//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti, Ricky Aristio
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos{
///@name Kratos Classes
///@{
/// Short class definition.
/** Base class for IGA 3D Beams.
*/
class KRATOS_API(IGA_APPLICATION) Base3DBeamElement
    : public Element
{
protected:
    /**
    * Internal variables used in the constitutive equations
    */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        
        ConstitutiveVariables(SizeType StrainSize = 3)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    static inline void ComputeSkewSymmetricMatrix(
        BoundedMatrix<double, 3, 3>& rSkewMatrix,
        const array_1d<double, 3>& rVec)
    {
        noalias(rSkewMatrix) = ZeroMatrix(3, 3);

        rSkewMatrix(0, 1) = -rVec[2];
        rSkewMatrix(0, 2) =  rVec[1];
        rSkewMatrix(1, 0) =  rVec[2];
        rSkewMatrix(1, 2) = -rVec[0];
        rSkewMatrix(2, 0) = -rVec[1];
        rSkewMatrix(2, 1) =  rVec[0];
    }
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Shell3pElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Base3DBeamElement);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    // GometryType
    typedef Geometry<Node> GeometryType;

    /// Constructor using an array of nodes
    Base3DBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    /// Constructor using an array of nodes with properties
    Base3DBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    /// Default constructor necessary for serialization
    Base3DBeamElement()
        : Element()
    {};

    /// Destructor.
    virtual ~Base3DBeamElement() = default;

    ///@}
    ///@name Life Cycle
    ///@{
    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<Base3DBeamElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< Base3DBeamElement >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * this determines the elemental equation ID vector for all elemental
    * DOFs
    * @param rResult: the elemental equation ID vector
    * @param rCurrentProcessInfo: the current process info instance
    */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override;

    /**
    * determines the elemental list of DOFs
    * @param ElementalDofList: the list of DOFs
    * @param rCurrentProcessInfo: the current process info instance
    */
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    static inline void ComputeRotationMatrix(
        MatrixType& rMatRod,
        const array_1d<double, 3>& rVec,
        const double Phi)
    {
        if (rMatRod.size1() != 3 || rMatRod.size2() != 3) {
            rMatRod.resize(3, 3, false);
        }

        const double c = std::cos(Phi);
        const double s = std::sin(Phi);

        BoundedMatrix<double, 3, 3> skew_matrix;
        ComputeSkewSymmetricMatrix(skew_matrix, rVec);

        noalias(rMatRod) = c * IdentityMatrix(3) + s * skew_matrix;
    }

    static inline void ComputeRotationMatrixDeriv(
        MatrixType& rMatRodDer,
        const array_1d<double, 3>& rVec,
        const array_1d<double, 3>& rVecDeriv,
        const double Phi,
        const double PhiDeriv)
    {
        if (rMatRodDer.size1() != 3 || rMatRodDer.size2() != 3) {
            rMatRodDer.resize(3, 3, false);
        }

        const double c = std::cos(Phi);
        const double s = std::sin(Phi);

        BoundedMatrix<double, 3, 3> skew_vec;
        BoundedMatrix<double, 3, 3> skew_vec_deriv;

        ComputeSkewSymmetricMatrix(skew_vec, rVec);
        ComputeSkewSymmetricMatrix(skew_vec_deriv, rVecDeriv);

        noalias(rMatRodDer) =
            (-PhiDeriv * s) * IdentityMatrix(3)
            + (PhiDeriv * c) * skew_vec
            + s * skew_vec_deriv;
    }

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{



    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;


}; // Class Base3DBeamElement

} // Namespace Kratos