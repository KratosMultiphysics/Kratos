//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               April 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_V_P_ELEMENT_H_INCLUDED)
#define  KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_V_P_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/large_displacement_segregated_V_P_element.hpp"


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

/// Updated Lagrangian Segregated Solid Element for 3D and 2D geometries

/**
 * Implements a Large Displacement Lagrangian definition for V-P segregated analysis.
 * This works for linear Triangles and Tetrahedra (base class)
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) UpdatedLagrangianSegregatedVPElement : public LargeDisplacementSegregatedVPElement
{
public:

  ///@name Type Definitions
  ///@{
  ///Reference type definition for constitutive laws
  typedef ConstitutiveLaw ConstitutiveLawType;
  ///Pointer type for constitutive laws
  typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
  ///StressMeasure from constitutive laws
  typedef ConstitutiveLawType::StressMeasure StressMeasureType;
  ///Type definition for integration methods
  typedef GeometryData::IntegrationMethod IntegrationMethod;
  ///Type for size
  typedef GeometryData::SizeType SizeType;

  /// Counted pointer of UpdatedLagrangianSegregatedVPElement
  KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianSegregatedVPElement );
  ///@}

  ///@name Life Cycle
  ///@{

  /// Empty constructor needed for serialization
  UpdatedLagrangianSegregatedVPElement();

  /// Default constructors
  UpdatedLagrangianSegregatedVPElement(IndexType NewId, GeometryType::Pointer pGeometry);

  UpdatedLagrangianSegregatedVPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

  ///Copy constructor
  UpdatedLagrangianSegregatedVPElement(UpdatedLagrangianSegregatedVPElement const& rOther);


  /// Destructor.
  ~UpdatedLagrangianSegregatedVPElement() override;

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  UpdatedLagrangianSegregatedVPElement& operator=(UpdatedLagrangianSegregatedVPElement const& rOther);

  ///@}
  ///@name Operations
  ///@{
  /**
   * Returns the currently selected integration method
   * @return current integration method selected
   */
  /**
   * creates a new total lagrangian updated element pointer
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

  /**
   * clones the selected element variables, creating a new one
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

  //************* GETTING METHODS
  //SET
  /**
   * Set a double  Value on the Element Constitutive Law
   */
  void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  //GET:
  /**
   * Get on rVariable a double Value from the Element Constitutive Law
   */
  void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  //************* STARTING - ENDING  METHODS
  /**
   * Called to initialize the element.
   * Must be called before any calculation is done
   */
  void Initialize() override;

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
    buffer << "Updated Lagrangian Segregated V-P Solid Element #" << Id();
    return buffer.str();
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "Updated Lagrangian Segregated V-P Solid Element #" << Id();
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

  /**
   * Container for historical total elastic deformation measure F0 = dx/dX
   */
  std::vector< Matrix > mDeformationGradientF0;

  /**
   * Container for the total deformation gradient determinants
   */
  Vector mDeterminantF0;

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{

  /**
   * Finalize Element Internal Variables
   */
  void FinalizeStepVariables(ElementDataType & rVariables,
                             const double& rPointNumber ) override;

  /**
   * Initialize Element General Variables
   */
  void InitializeElementData(ElementDataType & rVariables,
                             const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Calculate Element Kinematics
   */
  void CalculateKinematics(ElementDataType& rVariables,
                           const double& rPointNumber) override;

  /**
   * Calculate Element Kinetics
   */
  void CalculateKinetics(ElementDataType& rVariables,
                         const double& rPointNumber) override;

  /**
   * Set Parameters for the Constitutive Law and Calculate Material Response
   */
  void CalculateMaterialResponse(ElementDataType& rVariables,
                                 ConstitutiveLaw::Parameters& rValues,
                                 const int & rPointNumber) override;

  /**
   * Calculation of the Pressure Stiffness Matrix. Kpp
   */
  void CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                          ElementDataType& rVariables) override;

  /**
   * Calculation of the Pressure Vector.
   */
  void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
                                     ElementDataType & rVariables) override;

  /**
   * Calculate stabilization factor
   */
  void CalculateStabilizationTau(ElementDataType& rVariables);

  /**
   * Get Faces with nodes in the freesurface
   */
  void GetFreeSurfaceFaces(std::vector<std::vector<SizeType> >& Faces);

  /**
   * Get Face normal
   */
  void GetFaceNormal(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, Vector& rNormal);

  /**
   * Get Face weight
   */
  void GetFaceWeight(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, double& rWeight, double& rNormalSize);

  /**
   * Get Face normal
   */
  void GetFaceNormal(const std::vector<SizeType>& rFace, Vector& rNormal);


  /**
   * Get the Historical Deformation Gradient to calculate after finalize the step
   */
  void GetHistoricalVariables( ElementDataType& rVariables,
                               const double& rPointNumber ) override;

  /**
   * Calculation of the Volume Change of the Element
   */
  double& CalculateVolumeChange(double& rVolumeChange, ElementDataType& rVariables) override;

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

  void save(Serializer& rSerializer) const override;

  void load(Serializer& rSerializer) override;


  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}

}; // Class UpdatedLagrangianSegregatedVPElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_V_P_ELEMENT_H_INCLUDED  defined
