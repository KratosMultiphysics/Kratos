//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONTACT_DOMAIN_LM_2D_CONDITION_H_INCLUDED )
#define  KRATOS_CONTACT_DOMAIN_LM_2D_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_conditions/contact_domain_condition.hpp"

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


class ContactDomainLM2DCondition
    : public ContactDomainCondition
{
public:


    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    ///Tensor order 1 definition
    typedef array_1d<double, 3>    VectorType;
    ///NodeType
    typedef Node < 3 > NodeType;
    ///Geometry Type
    typedef Geometry<NodeType> GeometryType;
    ///Element Type
    typedef Element::ElementType ElementType;


    /// Counted pointer of ContactDomainLM2DCondition
    KRATOS_CLASS_POINTER_DEFINITION(ContactDomainLM2DCondition);

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    ContactDomainLM2DCondition() : ContactDomainCondition() {};

    ContactDomainLM2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    ContactDomainLM2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    ContactDomainLM2DCondition(ContactDomainLM2DCondition const& rOther);


    /// Destructor.
    virtual ~ContactDomainLM2DCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ContactDomainLM2DCondition& operator=(ContactDomainLM2DCondition const& rOther);


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    //************* GETTING METHODS


    //************* COMPUTING  METHODS



    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    //int Check(const ProcessInfo& rCurrentProcessInfo);

    //std::string Info() const;

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
    //      virtual String Info() const;

    /// Print information about this object.
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;
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

    /**
     * Calculation of the Contact Master Nodes and Mechanical variables
     */
    void SetMasterGeometry();


    /**
     * Calculate Tau stabilization or Penalty factor
     */
    virtual void CalculateContactFactor(ProcessInfo& rCurrentProcessInfo);
	

    /**
     * Calculation of the Contact Previous Gap
     */
    void CalculatePreviousGap();


    /**
     * Calculation of the Contact Multipliers or Penalty Factors
     */
    virtual void CalculateExplicitFactors(GeneralVariables& rVariables,
					  ProcessInfo& rCurrentProcessInfo);
    /**
     * Tangent Matrix construction methods:
     */
    void CalculateDomainShapeN(GeneralVariables& rVariables);


    /**
     * Calculate Integration Weight:
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);

    /**
     * Calculation of the Material Stiffness Matrix by components
     */
    virtual void CalcContactStiffness (double &Kcont,GeneralVariables& rVariables,
				       unsigned int& ndi,unsigned int& ndj,
				       unsigned int& idir,unsigned int& jdir);


    /**
     * Normal Force construction by components
     */
    virtual void CalculateNormalForce       (double &F,GeneralVariables& rVariables,
					     unsigned int& ndi,unsigned int& idir);

    /**
     * Tangent Stick Force construction by components
     */
    virtual void CalculateTangentStickForce (double &F,GeneralVariables& rVariables,
					     unsigned int& ndi,unsigned int& idir);
    /**
     * Tangent Slip Force construction by components
     */
    virtual void CalculateTangentSlipForce  (double &F,GeneralVariables& rVariables,
					     unsigned int& ndi,unsigned int& idir);

    ///@}
    ///@name Protected Operations
    ///@{

 
    /**
     * Util geometrical methods:
     */

    VectorType & CalculateCurrentTangent( VectorType &rTangent );


    void  CalcBaseDistances (BaseLengths& Base,VectorType& P1,VectorType& P2,VectorType& PS,VectorType& Normal);

    VectorType & ComputeFaceNormal(VectorType &Normal, VectorType& P1, VectorType &P2);

    VectorType & ComputeFaceTangent(VectorType &Tangent ,VectorType& P1, VectorType &P2);
   
    VectorType & ComputeFaceTangent(VectorType &Tangent ,VectorType& Normal);

    inline bool CheckFictiousContacts(GeneralVariables& rVariables);

    inline bool CalculatePosition(const double x0, const double y0,
				  const double x1, const double y1,
				  const double x2, const double y2,
				  const double xc, const double yc);

    inline bool CalculateObtuseAngle(const double x0, const double y0,
				     const double x1, const double y1,
				     const double xc, const double yc);


    inline double CalculateVol(const double x0, const double y0,
			       const double x1, const double y1,
			       const double x2, const double y2);
  


    void FSigmaP(GeneralVariables& rVariables, std::vector<Vector > &SigmaP, VectorType& AuxVector,unsigned int &ndi,unsigned int &ndj,unsigned int &ndk,unsigned int &ndr);

    void FSigmaPnd(GeneralVariables& rVariables, std::vector<Vector > &SigmaP, VectorType& AuxVector,unsigned int &ndi,unsigned int &ndj);



    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

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
    ///@name Serialization
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ContactDomainLM2DCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
    ContactDomainLM2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
    const ContactDomainLM2DCondition& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
///@}

} // namespace Kratos.
#endif // KRATOS_CONTACT_DOMAIN_LM_2D_CONDITION_H_INCLUDED  defined 
