#if !defined(KRATOS_MESHLESS_PENALTY_CONTINUITY_CRACK_CONDITION_H_INCLUDED )
#define  KRATOS_MESHLESS_PENALTY_CONTINUITY_CRACK_CONDITION_H_INCLUDED



// System includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_conditions/meshless_base_coupling_condition.h"
#include "custom_conditions/meshless_base_condition.h"


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

/// Short class definition.
/** Detail class definition.
*/
class MeshlessPenaltyCouplingCrackCondition
    : public MeshlessBaseCouplingCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MeshlessPenaltyCouplingCrackCondition
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessPenaltyCouplingCrackCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshlessPenaltyCouplingCrackCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MeshlessPenaltyCouplingCrackCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MeshlessPenaltyCouplingCrackCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

	Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

	void Initialize();

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo) override;


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


    ///@}

private:
    ///@name Static Member Variables
    ///@{
	ConstitutiveLaw::Pointer mConstitutiveLawVector;
	//array_1d<double, 3> mg1_0_slave;
	//array_1d<double, 3> mg2_0_slave;
	//array_1d<double, 3> mg3_0_slave;
	//array_1d<double, 3> mg1_0_master;
	//array_1d<double, 3> mg2_0_master;
	//array_1d<double, 3> mg3_0_master;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    MeshlessPenaltyCouplingCrackCondition() : MeshlessBaseCouplingCondition()
    {
    }

  //  virtual void save(Serializer& rSerializer) const
  //  {
  //      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
		//rSerializer.save("g1_0_slave", mg1_0_slave);
		//rSerializer.save("g2_0_slave", mg2_0_slave);
		//rSerializer.save("g3_0_slave", mg3_0_slave);
		//rSerializer.save("g1_0_master", mg1_0_master);
		//rSerializer.save("g2_0_master", mg2_0_master);
		//rSerializer.save("g3_0_master", mg3_0_master);
  //  }

  //  virtual void load(Serializer& rSerializer)
  //  {
  //      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
		//rSerializer.load("g1_0_slave", mg1_0_slave);
		//rSerializer.load("g2_0_slave", mg2_0_slave);
		//rSerializer.load("g3_0_slave", mg3_0_slave);
		//rSerializer.load("g1_0_master", mg1_0_master);
		//rSerializer.load("g2_0_master", mg2_0_master);
		//rSerializer.load("g3_0_master", mg3_0_master);
  //  }
  //  
    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

	//void CaculateRotation(const Matrix &ShapeFunctionDerivatives,
	//	Vector &Phi_r, 
	//	Matrix &Phi_rs, 
	//	array_1d<double, 2> Phi,
	//	const array_1d<double, 2> &Tangents, 
	//	const bool Master);

	//void JacobianElement(const Matrix& DN_De,
	//	Matrix& Jacobian, const bool Master);

	//void MappingGeometricToParameterOneElement(const Matrix& DN_De,
	//	const array_1d<double, 2>& Tangents,
	//	double& JacobianMatrix);
    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //MeshlessPenaltyCouplingCrackCondition& operator=(const MeshlessPenaltyCouplingCrackCondition& rOther);

    /// Copy constructor.
    //MeshlessPenaltyCouplingCrackCondition(const MeshlessPenaltyCouplingCrackCondition& rOther);


    ///@}

}; // Class MeshlessPenaltyCouplingCrackCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    MeshlessPenaltyCouplingCrackCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshlessPenaltyCouplingCrackCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_PENALTY_CONTINUITY_CRACK_CONDITION_H_INCLUDED  defined 


