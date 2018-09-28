//
//   Project Name:        KratosMachiningApplication $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:        October 2013 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_THERMAL_CONTACT_DOMAIN_CONDITION_H_INCLUDED )
#define  KRATOS_THERMAL_CONTACT_DOMAIN_CONDITION_H_INCLUDED

// System includes
#include <iostream>
#include <iomanip>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"

#include "custom_utilities/contact_domain_utilities.hpp"

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

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) ThermalContactDomainCondition
    : public Condition
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

    ///NodeType
    typedef Node < 3 > NodeType;
    ///Geometry Type
    typedef Geometry<NodeType> GeometryType;
    ///Element Type
    typedef Element::ElementType ElementType;


    ///Tensor order 1 definition
    typedef ContactDomainUtilities::PointType               PointType;
    ///SurfaceVector
    typedef ContactDomainUtilities::SurfaceVector       SurfaceVector;
    ///SurfaceScalar
    typedef ContactDomainUtilities::SurfaceScalar       SurfaceScalar;
    ///BaseLengths
    typedef ContactDomainUtilities::BaseLengths           BaseLengths;


protected:

    /**
     * Parameters to be used in the Condition as they are. Direct interface to Parameters Struct
     */

    typedef struct
    {
        Flags           Options;      //calculation options

        double          ThermalGap;   //thermal gap


        //Geometric variables
	Vector           ProjectionsVector;

	SurfaceVector        CurrentSurface;      //Current Normal and Tangent
        SurfaceVector        ReferenceSurface;    //Reference Normal and Tangent

        std::vector<BaseLengths>   CurrentBase;    //Current Base Lengths variables
        std::vector<BaseLengths>   ReferenceBase;  //Reference Base Lengths variables


        //Axisymmetric
        double  CurrentRadius;
        double  ReferenceRadius;


        //friction:
        double          RelativeVelocityNorm;
	double          FrictionForceNorm;


    } GeneralVariables;


    typedef struct
    {

       //The stabilization parameter and penalty parameter
        double          StabilizationFactor;

        //Contact condition conectivities
        std::vector<unsigned int> nodes;
        std::vector<unsigned int> order;
        std::vector<unsigned int> slaves;

	//Pointer Variables
        GeometryType*         mpMasterGeometry;
        ElementType*          mpMasterElement;
        Condition*            mpMasterCondition;
        NodeType*             mpMasterNode;

	/**
         * sets the value of a specified pointer variable
	 */

	void SetMasterGeometry  (GeometryType& rGeometry){ mpMasterGeometry = &rGeometry; }
        void SetMasterElement   (ElementType& rElement){ mpMasterElement = &rElement; }
        void SetMasterCondition (ConditionType& rCondition){ mpMasterCondition = &rCondition; }
        void SetMasterNode      (NodeType& rNode){ mpMasterNode = &rNode; }

	/**
         * returns the value of a specified pointer variable
         */

        GeometryType& GetMasterGeometry()   { return (*mpMasterGeometry); }
	ElementType& GetMasterElement()     { return (*mpMasterElement); }
	ConditionType& GetMasterCondition() { return (*mpMasterCondition); }
	NodeType& GetMasterNode()           { return (*mpMasterNode); }


    } ContactVariables;



public:

    /// Counted pointer of ThermalContactDomainCondition
    KRATOS_CLASS_POINTER_DEFINITION(ThermalContactDomainCondition);

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    ThermalContactDomainCondition() : Condition() {};

    ThermalContactDomainCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    ThermalContactDomainCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    ThermalContactDomainCondition(ThermalContactDomainCondition const& rOther);


    /// Destructor.
    virtual ~ThermalContactDomainCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ThermalContactDomainCondition& operator=(ThermalContactDomainCondition const& rOther);


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new total lagrangian updated condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;


    /**
     * clones a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    //************* GETTING METHODS

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod() override;

    /**
     * Sets on rConditionalDofList the degrees of freedom of the considered condition geometry
     */
    void GetDofList(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rResult the ID's of the condition degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;



    //on integration points:
    /**
     * Access for variables on Integration points.
     * This gives access to variables stored in the constitutive law on each integration point.
     * Specialisations of condition.h (e.g. the TotalLagrangian) must specify the actual
     * interface to the constitutive law!
     * Note, that these functions expect a std::vector of values for the
     * specified variable type that contains a value for each integration point!
     * SetValueOnIntegrationPoints: set the values for given Variable.
     * GetValueOnIntegrationPoints: get the values for given Variable.
     */

    //SET
    /**
     * Set a double Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
    /**
     * Set a Vector Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Set a Matrix Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    //GET:
    /**
     * Set on rVariable a double Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Set on rVariable a Vector Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Set on rVariable a Matrix Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;



    //************* STARTING - ENDING  METHODS

    /**
     * Called to initialize the condition.
     * Must be called before any calculation is done
     */
    void Initialize() override;


    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;



    //************* COMPUTING  METHODS

    /**
     * this is called during the assembling process in order
     * to calculate all conditional contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the conditional left hand side matrix
     * @param rRightHandSideVector: the conditional right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the conditional right hand side vector only
     * @param rRightHandSideVector: the conditional right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the conditional left hand side vector only
     * @param rLeftHandSideVector: the conditional left hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;


    //on integration points:
    /**
     * Calculate a double Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate a Vector Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate a Matrix Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo) override;



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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Thermal Contact Domain Condition #" << Id();
        return buffer.str();

    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Thermal Contact Domain Condition #" << Id();
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
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;

    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;


    /**
     * Variables stored in the condition during the computation
     */
    ContactVariables mContactVariables;


    /**
     * Contact Domain Utilities
     */
    ContactDomainUtilities  mContactUtilities;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * Calculation of the Contact Master Nodes and Mechanical variables
     */
    virtual void SetMasterGeometry()
    {
      KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" );
    };

    /**
     * Calculate Thermal Conductivity
     */
    void CalculateHeatConductivity();


    /**
     * Calculates the conditional contributions
     */
    virtual void CalculateConditionalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo,
					    Flags& rCalculationFlags);



    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables,
				     ProcessInfo& rCurrentProcessInfo,
				     const unsigned int& rPointNumber)
    {
      KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" );
    };



    /**
     * Initialize Variables
     */
    void InitializeVariables ();



    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					  VectorType& rRightHandSideVector,
					  Flags& rCalculationFlags);


    /**
     * Calculate LHS
     */
    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
				    GeneralVariables& rVariables,
				    double& rIntegrationWeight);

    /**
     * Calculate RHS
     */
    virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                    GeneralVariables& rVariables,
                                    double& rIntegrationWeight);

    /**
     * Calculation of the Material Stiffness Matrix. Ktherm =  TauK * ([proj] * [proj]T)
     */
    virtual void CalculateAndAddThermalKm(MatrixType& rLeftHandSideMatrix,
					  GeneralVariables& rVariables,
					  double& rIntegrationWeight);

    /**
     * Calculation of the Internal Forces Vector. Fi = B * sigma
     */
    virtual void CalculateAndAddThermalContactForces(VectorType& rRightHandSideVector,
						     GeneralVariables& rVariables,
						     double& rIntegrationWeight);


    /**
     * Calculate Integration Weight:
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);


    /**
     * Force construction methods (frictional):
     */
    virtual void CalculateThermalFrictionForce(double &F, GeneralVariables& rVariables, unsigned int& ndi)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" );

	};


    /**
     * Force construction methods (thermal conduction):
     */
    virtual void CalculateThermalConductionForce(double &F, GeneralVariables& rVariables, unsigned int& ndi)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" );

	};


    /**
     * Calculate Relative Velocity:
     */
    void CalculateRelativeVelocity(GeneralVariables& rVariables, PointType& TangentVelocity, ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculate Relative Displacement:
     */
    void CalculateRelativeDisplacement(GeneralVariables& rVariables, PointType & TangentDisplacement, ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculate current tangent vector
     */
    virtual PointType & CalculateCurrentTangent(PointType &rTangent)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" );

	};

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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

}; // Class ThermalContactDomainCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_THERMAL_CONTACT_DOMAIN_CONDITION_H_INCLUDED  defined
