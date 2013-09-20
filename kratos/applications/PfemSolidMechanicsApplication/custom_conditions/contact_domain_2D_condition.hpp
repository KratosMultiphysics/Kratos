//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONTACT_DOMAIN_2D_CONDITION_H_INCLUDED )
#define  KRATOS_CONTACT_DOMAIN_2D_CONDITION_H_INCLUDED


// System includes
#include <iostream>
#include <iomanip>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


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

/// Updated Total Lagrangian element for 2D and 3D geometries.

/**
 * Implements an Updated Lagrangian Condition based on the reference (or previous) configuration
 * This works for arbitrary geometries in 2D and 3D
 */

class ContactDomain2DCondition
    : public Condition
{
protected:

    /**
     * Parameters to be used in the Condition as they are. Direct interface to Parameters Struct
     */


    typedef struct
    {
        double L;    //base side lentgh
        double A;    //distance 2-3
        double B;    //distance 1-2

    private:
      
      friend class Serializer;
      
      void save(Serializer& rSerializer) const 
      {
	// rSerializer.save("BaseLength",L);
	// rSerializer.save("BaseA",A);
	// rSerializer.save("BaseB",B);
      }
      
      void load(Serializer& rSerializer) 
      {
	// rSerializer.load("BaseLength",L);
	// rSerializer.load("BaseA",A);
	// rSerializer.load("BaseB",B);
      }

    } BaseLengths;


    typedef struct
    {
        double Normal;        //Projection of the Traction Vector on the current normal
        double Tangent;       //Projection of the Traction Vector on the current tangent

    private:
      
      friend class Serializer;
      
      void save(Serializer& rSerializer) const 
      {
	// rSerializer.save("NormalTraction",Normal);
	// rSerializer.save("TangentTraction",Tangent);
      }
      
      void load(Serializer& rSerializer) 
      {
	// rSerializer.load("NormalTraction",Normal);
	// rSerializer.load("TangentTraction",Tangent);
      }


    } Tensil;


    typedef struct
    {

        double          gapN; //normal gap
        double          gapT; //tangential gap

        double          PreviousGapN;   //Configuration n-1 effective gap
        double          PreviousGapT;   //Configuration n-1 effective gap

        double          gapTsign; //sign or direction of the tangential gap

        double          active;  //active set on

        int             stick;

        double          LmN;            //Lagrange Multipliyer normal
        double          LmT;            //Lagrange Multipliyer tangential


        //variables of the contact element 2D
        Vector           Nn;            //Discrete variacion of the shape function  in the current normal direction
        Vector           Nt;            //Discrete variacion of the shape function  in the current tangent direction
        Vector          Nrn;           //Discrete variacion of the normal  in the reference normal direction

        std::vector<Vector >       Nsigma;
        std::vector<Vector >       Tsigma;

        //Geometric variables

        array_1d<double, 3>        CurNormal;    //Current Normal
        array_1d<double, 3>        CurTangent;   //Current Tangent
        array_1d<double, 3>        RefNormal;    //Reference Normal
        array_1d<double, 3>        RefTangent;   //Reference Tangent

        array_1d<double, 3>        PreNormal;    //Previous Normal
        array_1d<double, 3>        PreTangent;   //Previous Tangent

        std::vector<BaseLengths>   CurrentBase;    //Current Base Lengths variables
        std::vector<BaseLengths>   ReferenceBase;  //Reference Base Lengths variables

        //Resultant mechanical tractions
        array_1d<double, 3>        TractionVector;         //Traction Vector in the reference configuration
        Tensil  CurrentTensil;    //Traction Vectors Tangential and Normal

        //The stabilization parameter
        double Tau;

        //Penalty method or Lagrangian Multiplyers
        unsigned int            penalty;

        //friction:
        double          muCoefficient;
        unsigned int    friction_active;//friction active or inactive
        unsigned int    friction_on;

        //iteration
        int             iteration;


    private:
      
      friend class Serializer;
      
      void save(Serializer& rSerializer) const 
      {
	// rSerializer.save("NormalGap",gapN);
	// rSerializer.save("TangentGap",gapT);
	// rSerializer.save("PreviousNormalGap",PreviousGapN);
	// rSerializer.save("PreviousTangentGap",PreviousGapT);
	// rSerializer.save("TangentGapSign",gapTsign);
	// rSerializer.save("SetActive",active);
	// rSerializer.save("SetStick",stick);

	// rSerializer.save("NormalMultiplier",LmN);
	// rSerializer.save("TangentMultiplier",LmT);
	
	// rSerializer.save("NormalShapeFunc",Nn);
	// rSerializer.save("TangentShapeFunc",Nt);
	// rSerializer.save("RefNormalShapeFunc",Nrn);

	// rSerializer.save("NormalSigma",Nsigma);
	// rSerializer.save("TangentSigma",Tsigma);


	// rSerializer.save("CurrentNormal",CurNormal);
	// rSerializer.save("CurrentTangent",CurTangent);
	// rSerializer.save("ReferenceNormal",RefNormal);
	// rSerializer.save("ReferenceTangent",RefTangent);

	// rSerializer.save("PreviousNormal",PreNormal);
	// rSerializer.save("PreviousTangent",PreTangent);

	// rSerializer.save("CurrentBase",CurrentBase);
	// rSerializer.save("ReferenceBase",ReferenceBase);

	// rSerializer.save("TractionVector",TractionVector);

	// rSerializer.save("CurrentTensil",CurrentTensil);

	// rSerializer.save("StabilizationTau",Tau); 	
	// rSerializer.save("PenaltySelection",penalty); 

	// rSerializer.save("FrictionCoeffitient",muCoefficient); 
	// rSerializer.save("FrictionActive",friction_active); 
	// rSerializer.save("FrictionOn",friction_on); 
	
	// rSerializer.save("Iteration",iteration);
      }

      
      void load(Serializer& rSerializer) 
      {
	// rSerializer.load("NormalGap",gapN);
	// rSerializer.load("TangentGap",gapT);
	// rSerializer.load("PreviousNormalGap",PreviousGapN);
	// rSerializer.load("PreviousTangentGap",PreviousGapT);
	// rSerializer.load("TangentGapSign",gapTsign);
	// rSerializer.load("SetActive",active);
	// rSerializer.load("SetStick",stick);

	// rSerializer.load("NormalMultiplier",LmN);
	// rSerializer.load("TangentMultiplier",LmT);
	
	// rSerializer.load("NormalShapeFunc",Nn);
	// rSerializer.load("TangentShapeFunc",Nt);
	// rSerializer.load("RefNormalShapeFunc",Nrn);

	// rSerializer.load("NormalSigma",Nsigma);
	// rSerializer.load("TangentSigma",Tsigma);


	// rSerializer.load("CurrentNormal",CurNormal);
	// rSerializer.load("CurrentTangent",CurTangent);
	// rSerializer.load("ReferenceNormal",RefNormal);
	// rSerializer.load("ReferenceTangent",RefTangent);

	// rSerializer.load("PreviousNormal",PreNormal);
	// rSerializer.load("PreviousTangent",PreTangent);

	// rSerializer.load("CurrentBase",CurrentBase);
	// rSerializer.load("ReferenceBase",ReferenceBase);

	// rSerializer.load("TractionVector",TractionVector);

	// rSerializer.load("CurrentTensil",CurrentTensil);

	// rSerializer.load("StabilizationTau",Tau); 	
	// rSerializer.load("PenaltySelection",penalty); 

	// rSerializer.load("FrictionCoeffitient",muCoefficient); 
	// rSerializer.load("FrictionActive",friction_active); 
	// rSerializer.load("FrictionOn",friction_on); 
	
	// rSerializer.load("Iteration",iteration);
      }
      

    } ContactParameters;


    typedef struct
    {
        double  detF;
        double  detJ;
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;
        Matrix  F;
        Matrix  DN_DX;
        Matrix  ConstitutiveMatrix;

        ContactParameters Contact;

        std::vector<unsigned int> order;
        std::vector<unsigned int> nodes;
        std::vector<unsigned int> slaves;

   private:
      
      friend class Serializer;
      
      void save(Serializer& rSerializer) const 
      {
	// rSerializer.save("DetF",detF);
	// rSerializer.save("DetJ",detJ);
	// rSerializer.save("StrainVector",StrainVector);
	// rSerializer.save("StressVector",StressVector);
	// rSerializer.save("ShapeFunctions",N);
	// rSerializer.save("DeformationGradient",F);
	// rSerializer.save("ShapeFunctionsDerivatives",DN_DX);
	// rSerializer.save("ConstitutiveMatrix",ConstitutiveMatrix);

	// rSerializer.save("ContactParameters",Contact);
	// rSerializer.save("Order",order);
	// rSerializer.save("Nodes",nodes);
	// rSerializer.save("Slaves",slaves);
      }
      
      void load(Serializer& rSerializer) 
      {
	// rSerializer.load("DetF",detF);
	// rSerializer.load("DetJ",detJ);
	// rSerializer.load("StrainVector",StrainVector);
	// rSerializer.load("StressVector",StressVector);
	// rSerializer.load("ShapeFunctions",N);
	// rSerializer.load("DeformationGradient",F);
	// rSerializer.load("ShapeFunctionsDerivatives",DN_DX);
	// rSerializer.load("ConstitutiveMatrix",ConstitutiveMatrix);

	// rSerializer.load("ContactParameters",Contact);
	// rSerializer.load("Order",order);
	// rSerializer.load("Nodes",nodes);
	// rSerializer.load("Slaves",slaves);
      }

    } Standard;



public:


    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;


    /// Counted pointer of ContactDomain2DCondition
    KRATOS_CLASS_POINTER_DEFINITION(ContactDomain2DCondition);

    KRATOS_DEFINE_LOCAL_FLAG( ACTIVE );
    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    ContactDomain2DCondition() : Condition() {};

    ContactDomain2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    ContactDomain2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    ContactDomain2DCondition(ContactDomain2DCondition const& rOther);


    /// Destructor.
    virtual ~ContactDomain2DCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ContactDomain2DCondition& operator=(ContactDomain2DCondition const& rOther);


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

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod();

    /**
     * Sets on rConditionalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0);



    //on integration points:
    /**
     * Access for variables on Integration points.
     * This gives access to variables stored in the constitutive law on each integration point.
     * Specialisations of element.h (e.g. the TotalLagrangian) must specify the actual
     * interface to the constitutive law!
     * Note, that these functions expect a std::vector of values for the
     * specified variable type that contains a value for each integration point!
     * SetValueOnIntegrationPoints: set the values for given Variable.
     * GetValueOnIntegrationPoints: get the values for given Variable.
     */

    //SET
    /**
     * Set a Vector Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set a Matrix Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    //GET:
    /**
     * Set on rVariable a double Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set on rVariable a Vector Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set on rVariable a Matrix Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);



    //************* STARTING - ENDING  METHODS

    /**
     * Called to initialize the element.
     * Must be called before any calculation is done
     */
    void Initialize();


    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);



    //************* COMPUTING  METHODS

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side vector only
     * @param rLeftHandSideVector: the elemental left hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);


    //on integration points:
    /**
     * Calculate a double Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Vector Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Matrix Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo);


    void Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo);

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

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
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateConditionalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo,
                                            bool CalculateStiffnessMatrixFlag,
                                            bool CalculateResidualVectorFlag);
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
    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;

    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    /**
     * Variables stored in the element during the computation
     */

    GeometryType::Pointer mpMasterGeometry;
    Standard mVariables;

    ///@}
    ///@name Private Operators
    ///@{


    /**
     * Initialize Variables
     */
    void InitializeVariables ();


    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();


    /**
     * Initialize System Matrices
     */
    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
				  VectorType& rRightHandSideVector,
				  bool CalculateStiffnessMatrixFlag,
				  bool CalculateResidualVectorFlag );

    /**
     * Calculate Tau stabilization
     */
    void CalculateTauStab(ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculate Condition Kinematics
     */
    void CalculateKinematics(const double& rPointNumber,ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculation of the Deformation Gradient F
     */
    void CalculateDeformationGradient(const Matrix& rDN_DX,
                                      Matrix& rF);


    /**
     * Calculation of the Contact Master Nodes and Mechanical variables
     */
    void SetMasterGeometry();

    /**
     * Calculation of the Contact Previous Gap
     */
    void CalcPreviousGap();

    /**
     * Calculation of the Contact Multipliers
     */
    void CalcMultipliers(ProcessInfo& rCurrentProcessInfo);

    void CalcPenaltyParameters(ProcessInfo& rCurrentProcessInfo);


    /**
     * Util methods:
     */

  
    inline void  CalcBaseDistances (BaseLengths& Base,array_1d<double, 3>& P1,array_1d<double, 3>& P2,array_1d<double, 3>& PS,array_1d<double, 3>& Normal);

    inline array_1d<double, 3> & ComputeFaceNormal(array_1d<double, 3> &Normal, array_1d<double, 3>& P1, array_1d<double, 3> &P2);

    inline array_1d<double, 3> & ComputeFaceTangent(array_1d<double, 3> &Tangent ,array_1d<double, 3>& P1, array_1d<double, 3> &P2);

    inline array_1d<double, 3> & CalcCurrentTangent( array_1d<double, 3> &Tangent );

    inline array_1d<double, 3> & ComputeFaceTangent(array_1d<double, 3> &Tangent ,array_1d<double, 3>& Normal);

    inline bool CheckFictiousContacts();
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
  

    /**
     * Friction Parameters:
     */
    inline void CalcRelativeVelocity     (array_1d<double, 3 > & TangentVelocity);

    inline void CalcRelativeDisplacement (array_1d<double, 3 > & TangentDisplacement);

    void CalcFrictionCoefficient (const array_1d<double, 3 > & TangentVelocity);

//************************************************************************************
//************************************************************************************

    /**
     * Calculation of the Material Stiffness Matrix. Km = BT * D * B
     */
    void CalculateAndAddKm(MatrixType& rK,
                           double& rIntegrationWeight
                          );

    void CalculateAndAddPenaltyKm(MatrixType& rK,
                           double& rIntegrationWeight
                          );


    /**
     * Calculation of the Internal Forces Vector. Fi = B * sigma
     */
    void CalculateAndAddContactForces(VectorType& rRightHandSideVector,
				      double& rIntegrationWeight
                                      );

    void CalculateAndAddContactPenaltyForces(VectorType& rRightHandSideVector,
				      double& rIntegrationWeight
				     );



    /**
     * Tangent construction methods:
     */
    void CalcDomainShapeN();

    void FSigmaP(std::vector<Vector > &SigmaP, array_1d<double, 3>& AuxVector,unsigned int &ndi,unsigned int &ndj,unsigned int &ndk,unsigned int &ndr);

    void FSigmaPnd(std::vector<Vector > &SigmaP, array_1d<double, 3>& AuxVector,unsigned int &ndi,unsigned int &ndj);

    void CalcContactStiffness (double &Kcont,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir);

    void CalcContactPenaltyStiffness (double &Kcont,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir);
  

    /**
     * Force construction methods:
     */
    inline void CalcNormalForce       (double &F,unsigned int& ndi,unsigned int& idir);
    inline void CalcTangentStickForce (double &F,unsigned int& ndi,unsigned int& idir);
    inline void CalcTangentSlipForce  (double &F,unsigned int& ndi,unsigned int& idir);


     inline void CalcNormalPenaltyForce       (double &F,unsigned int& ndi,unsigned int& idir);
     inline void CalcTangentStickPenaltyForce (double &F,unsigned int& ndi,unsigned int& idir);
     inline void CalcTangentSlipPenaltyForce  (double &F,unsigned int& ndi,unsigned int& idir);
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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ContactDomain2DCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
    ContactDomain2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
    const ContactDomain2DCondition& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
///@}

} // namespace Kratos.
#endif // KRATOS_CONTACT_DOMAIN_2D__CONDITION_H_INCLUDED  defined 
