//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

#if !defined(KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "compressible_potential_flow_application_variables.h"


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

class CompressiblePotentialFlowElement : public Element {
public:
    
    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double,TNumNodes> phis, distances;
        double rho;
        double vol;
        
        bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
    };



  ///@name Type Definitions
  ///@{

  typedef Element BaseType;

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of CompressiblePotentialFlowElement
  KRATOS_CLASS_POINTER_DEFINITION(CompressiblePotentialFlowElement);

  ///@}
  ///@name Life Cycle
  ///@{

  /**
   * Constructor.
   */
  CompressiblePotentialFlowElement(IndexType NewId = 0);

  /**
   * Constructor using an array of nodes
   */
  CompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes);

  /**
   * Constructor using Geometry
   */
  CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry);

  /**
   * Constructor using Properties
   */
  CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

  /**
   * Copy Constructor
   */
  CompressiblePotentialFlowElement(CompressiblePotentialFlowElement const& rOther);

  /**
   * Destructor
   */
  virtual ~CompressiblePotentialFlowElement();

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  CompressiblePotentialFlowElement & operator=(CompressiblePotentialFlowElement const& rOther);

  ///@}
  ///@name Operations
  ///@{

  /**
   * ELEMENTS inherited from this class have to implement next
   * Create and Clone methods: MANDATORY
   */

  /**
   * creates a new element pointer
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

  /**
   * creates a new element pointer
   * @param NewId: the ID of the new element
   * @param pGeom: the geometry to be employed
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

  /**
   * creates a new element pointer and clones the previous element data
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

  /**
   * this determines the elemental equation ID vector for all elemental
   * DOFs
   * @param rResult: the elemental equation ID vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

  /**
   * determines the elemental list of DOFs
   * @param ElementalDofList: the list of DOFs
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override;

  /**
   * ELEMENTS inherited from this class have to implement next
   * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
   * they can be managed internally with a private method to do the same calculations
   * only once: MANDATORY
   */

  /**
   * this is called during the assembling process in order
   * to calculate all elemental contributions to the global system
   * matrix and the right hand side
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rRightHandSideVector: the elemental right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateLocalSystem(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the elemental right hand side vector only
   * @param rRightHandSideVector: the elemental right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;


  /**
   * This method provides the place to perform checks on the completeness of the input
   * and the compatibility with the problem options as well as the contitutive laws selected
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   * this method is: MANDATORY
   */
  virtual int Check(const ProcessInfo& rCurrentProcessInfo) override;
  
  virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

  virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

            
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
  virtual std::string Info() const;

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const;

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const;

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
    void ComputeLHSGaussPointContribution(
                                       const double weight,
                                       bounded_matrix<double,3,3>& lhs,
                                       const ElementalData<3,2>& data);

    void ComputeRHSGaussPointContribution(
                                        const double weight,
                                        array_1d<double,3>& rhs,
                                        const ElementalData<3,2>& data);
    

    void GetValuesOnSplitElement(Vector& split_element_values );
    
    void ComputeSplitVolumes(array_1d<double,3>& distances, double& positive_vol, double& negative_vol);

    
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
  ///@name Serialization
  ///@{

  friend class Serializer;

  virtual void save(Serializer& rSerializer) const;
  virtual void load(Serializer& rSerializer);

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

}; // Class CompressiblePotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED  defined
