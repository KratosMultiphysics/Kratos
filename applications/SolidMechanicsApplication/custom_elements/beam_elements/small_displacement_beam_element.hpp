//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SMALL_DISPLACEMENT_BEAM_ELEMENT_H_INCLUDED)
#define  KRATOS_SMALL_DISPLACEMENT_BEAM_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/beam_element.hpp"

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

/// Small displacements beam element for 2D and 3D space

/**
 * Implements a Small Displacement definition for structural analysis.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) SmallDisplacementBeamElement
    :public BeamElement
{
public:

    ///@name Type Definitions
    ///@{

    /// Counted pointer of SmallDisplacementBeamElement
    KRATOS_CLASS_POINTER_DEFINITION( SmallDisplacementBeamElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    SmallDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry);

    SmallDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    SmallDisplacementBeamElement(SmallDisplacementBeamElement const& rOther);

    /// Destructor.
    ~SmallDisplacementBeamElement() override;


    ///@}
    ///@name Operators
    ///@{

   /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;


    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Small Displacement Beam Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Small Displacement Beam Element #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      GetGeometry().PrintData(rOStream);
    }
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
    SmallDisplacementBeamElement() {};

    //constexpr const std::size_t& Dimension() const {return GetGeometry().WorkingSpaceDimension();}

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize Element General Variables
     */
    void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementDataType& rVariables,
                                     const unsigned int& rPointNumber) override;

    /**
     * Calculate Deformation Matrix
     */
    virtual void CalculateDeformationMatrix(Matrix& rB, const Vector& rN, const Matrix& rDN_DX);


    /**
     * Calculate Element Constitutive Matrix
     */
    void CalculateConstitutiveMatrix(ElementDataType& rVariables) override;


     /**
     * Calculate Element Stress Resultants and Couples
     */
    void CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber) override;

    /**
     * Calculation of the Rotation tensor
     */
    void CalculateTransformationMatrix(Matrix& RotationMatrix);


    /**
     * Transform Vector Variable form Spatial Frame to Global Frame
     */
    void MapLocalToGlobal(ElementDataType& rVariables, Matrix& rVariable) override;

    /**
     * Transform Vector Variable form Spatial Frame to Global Frame
     */
    void MapLocalToGlobal(ElementDataType& rVariables, VectorType& rVector) override;



    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                     ElementDataType& rVariables,
                                     double& rIntegrationWeight) override;


    /**
     * Calculation of the Material Stiffness Matrix
     */
    void CalculateLocalStiffnessMatrix(Matrix& LocalMatrix,ElementDataType& rVariables);


     /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
					       ElementDataType & rVariables,
					       double& rIntegrationWeight) override;


   ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

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
    ///@name Private  Access
    ///@{
    ///@}
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;


    // A private default constructor necessary for serialization


    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class SmallDisplacementBeamElement

} // namespace Kratos.
#endif //  KRATOS_SMALL_DISPLACEMENT_BEAM_ELEMENT_H_INCLUDED defined

