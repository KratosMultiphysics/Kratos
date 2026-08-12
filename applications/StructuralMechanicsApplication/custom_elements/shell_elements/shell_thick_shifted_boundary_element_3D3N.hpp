// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Ricky Aristio   
//

#pragma once

// System includes

// External includes

// Project includes

// Application includes
#include "shell_thick_element_3D3N.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"

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

template <ShellKinematics TKinematics>
class ShellThickShiftedBoundaryElement3D3N : public ShellThickElement3D3N<TKinematics>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ShellThickShiftedBoundaryElement3D3N
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ShellThickShiftedBoundaryElement3D3N);

    typedef ShellThickElement3D3N<TKinematics> BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;


    static constexpr std::size_t NumNodes = 3;

    static constexpr std::size_t StrainSize = 6;

    static constexpr std::size_t BlockSize = 6; //TO DO

    static constexpr std::size_t LocalSize = NumNodes * 6; //TO DO

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with geometry
    ShellThickShiftedBoundaryElement3D3N(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry);

    /// Constructor with geometry and properties
    ShellThickShiftedBoundaryElement3D3N(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ShellThickShiftedBoundaryElement3D3N() = default;

    /// Copy constructor.
    ShellThickShiftedBoundaryElement3D3N(const ShellThickShiftedBoundaryElement3D3N& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ShellThickShiftedBoundaryElement3D3N& operator=(const ShellThickShiftedBoundaryElement3D3N& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ShellThickShiftedBoundaryElement3D3N<TKinematics>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeom,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ShellThickShiftedBoundaryElement3D3N<TKinematics>>(NewId, pGeom, pProperties);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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


    ///@}

protected:
    /**
     * Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        Vector  N;
        Matrix  B;
        double  detF;
        Matrix  F;
        double  detJ0;
        Matrix  J0;
        Matrix  InvJ0;
        Matrix  DN_DX;
        Vector Displacements;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         * @param Dimension The problem dimension: 2D or 3D
         * @param NumberOfNodes The number of nodes in the element
         */
        KinematicVariables(
            const SizeType StrainSize,
            const SizeType Dimension,
            const SizeType NumberOfNodes
            )
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
            Displacements = ZeroVector(Dimension * NumberOfNodes);
        }
    };

    /**
     * Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const SizeType StrainSize)
        {
            if (StrainVector.size() != StrainSize)
                StrainVector.resize(StrainSize);

            if (StressVector.size() != StrainSize)
                StressVector.resize(StrainSize);

            if (D.size1() != StrainSize || D.size2() != StrainSize)
                D.resize(StrainSize, StrainSize);

            noalias(StrainVector) = ZeroVector(StrainSize);
            noalias(StressVector) = ZeroVector(StrainSize);
            noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
        }
    };

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


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // Protected default constructor necessary for serialization
    ShellThickShiftedBoundaryElement3D3N() : ShellThickElement3D3N<TKinematics>()
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ShellThickElement3D3N<TKinematics>);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ShellThickElement3D3N<TKinematics>);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Returns the surrogate faces local ids
     * This method returns a list with the local ids of the surrogate faces that
     * are the element faces for which all their nodes lie in the surrogate boundary
     * @return std::vector<std::size_t> List with the surrogate faces local ids
     */
    std::vector<std::size_t> GetSurrogateFacesIds();

    /**
     * @brief Calculates the Cauchy traction vector
     * This function computes the Cauchy traction vector from the provided stress
     * vector (in Voigt notation) and the unit normal vector
     * @param rVoigtStress Reference to the stress vector in Voigt notation
     * @param rUnitNormal Reference to the unit normal vector
     * @param rCauchyTraction Output Cauchy traction vector
     */
    void CalculateCauchyTractionVector(
        const Vector& rVoigtStress,
        const array_1d<double,2>& rUnitNormal,
        array_1d<double,6>& rCauchyTraction);

    /**
     * @brief Computes the normal projection of the C times B product
     * This function calculates the normal projection of the standard
     * C (constitutive tensor) times B (strain matrix) product
     * @param rC Reference to the constituive tensor
     * @param rB Reference to the strain matrix
     * @param rUnitNormal Reference to the unit normal vector
     * @param rAuxMat Output result
     */
    void CalculateCBProjectionLinearisation(
        const Matrix& rC,
        const BoundedMatrix<double,8,LocalSize>& rB,
        const array_1d<double,2>& rUnitNormal,
        BoundedMatrix<double,6,LocalSize>& rAuxMat);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class ShellThickShiftedBoundaryElement3D3N

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    ShellThickShiftedBoundaryElement3D3N& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const ShellThickShiftedBoundaryElement3D3N& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.
