// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes

// Application includes
#include "small_displacement.h"

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

template<std::size_t TDim>
class SmallDisplacementShiftedBoundaryElement : public SmallDisplacement
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of SmallDisplacementShiftedBoundaryElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementShiftedBoundaryElement);

    typedef SmallDisplacement BaseType;

    static constexpr std::size_t NumNodes = TDim + 1;

    static constexpr std::size_t StrainSize = 3 * (TDim-1);

    static constexpr std::size_t BlockSize = TDim;

    static constexpr std::size_t LocalSize = NumNodes * TDim;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with geometry
    SmallDisplacementShiftedBoundaryElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /// Constructor with geometry and properties
    SmallDisplacementShiftedBoundaryElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~SmallDisplacementShiftedBoundaryElement() = default;

    /// Copy constructor.
    SmallDisplacementShiftedBoundaryElement(const SmallDisplacementShiftedBoundaryElement& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SmallDisplacementShiftedBoundaryElement& operator=(const SmallDisplacementShiftedBoundaryElement& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

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
    SmallDisplacementShiftedBoundaryElement() : SmallDisplacement()
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacement);
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
        const array_1d<double,TDim>& rUnitNormal,
        array_1d<double,TDim>& rCauchyTraction);

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
        const BoundedMatrix<double,StrainSize,LocalSize>& rB,
        const array_1d<double,TDim>& rUnitNormal,
        BoundedMatrix<double,TDim,LocalSize>& rAuxMat);

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
}; // Class SmallDisplacementShiftedBoundaryElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    SmallDisplacementShiftedBoundaryElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const SmallDisplacementShiftedBoundaryElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.
