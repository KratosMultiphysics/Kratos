#if !defined(KRATOS_MESHLESS_BASE_COUPLING_CONDITION_H_INCLUDED )
#define  KRATOS_MESHLESS_BASE_COUPLING_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

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
class MeshlessBaseCouplingCondition
    : public MeshlessBaseCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MeshlessBaseCouplingCondition
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessBaseCouplingCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshlessBaseCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MeshlessBaseCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MeshlessBaseCouplingCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;
	
	void Initialize();

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
	array_1d<double, 3> mg1_0_slave;
	array_1d<double, 3> mg2_0_slave;
	array_1d<double, 3> mg3_0_slave;
	array_1d<double, 3> mg1_0_master;
	array_1d<double, 3> mg2_0_master;
	array_1d<double, 3> mg3_0_master;
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

	void CaculateRotationalShapeFunctions(Vector &Phi_r, Vector &Phi_r_Lambda, Matrix &Phi_rs, array_1d<double, 2> &Diff_Phi);

	void CaculateRotation(const Matrix &ShapeFunctionDerivatives,
		Vector &Phi_r,
		Matrix &Phi_rs,
		array_1d<double, 2> &Phi,
		array_1d<double, 3> &TrimTangent,
		const array_1d<double, 2> &Tangents,
		const bool Master);

	void JacobianElement(const Matrix& DN_De,
		Matrix& Jacobian, const bool Master);

	void MappingGeometricToParameterMasterElement(const Matrix& DN_De_Master,
		const array_1d<double, 2>& Tangents,
		double& JGeometricToParameter);

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
	///@name Serialization
	///@{
	friend class Serializer;

	// A private default constructor necessary for serialization
	MeshlessBaseCouplingCondition() : MeshlessBaseCondition()
	{
	}

	virtual void save(Serializer& rSerializer) const override
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
		rSerializer.save("g1_0_slave", mg1_0_slave);
		rSerializer.save("g2_0_slave", mg2_0_slave);
		rSerializer.save("g3_0_slave", mg3_0_slave);
		rSerializer.save("g1_0_master", mg1_0_master);
		rSerializer.save("g2_0_master", mg2_0_master);
		rSerializer.save("g3_0_master", mg3_0_master);
	}

	virtual void load(Serializer& rSerializer) override
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
		rSerializer.load("g1_0_slave", mg1_0_slave);
		rSerializer.load("g2_0_slave", mg2_0_slave);
		rSerializer.load("g3_0_slave", mg3_0_slave);
		rSerializer.load("g1_0_master", mg1_0_master);
		rSerializer.load("g2_0_master", mg2_0_master);
		rSerializer.load("g3_0_master", mg3_0_master);
	}
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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //MeshlessBaseCouplingCondition& operator=(const MeshlessBaseCouplingCondition& rOther);

    /// Copy constructor.
    //MeshlessBaseCouplingCondition(const MeshlessBaseCouplingCondition& rOther);


    ///@}

}; // Class MeshlessBaseCouplingCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    MeshlessBaseCouplingCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshlessBaseCouplingCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_LAGRANGE_COUPLING_CONDITION_H_INCLUDED  defined 


