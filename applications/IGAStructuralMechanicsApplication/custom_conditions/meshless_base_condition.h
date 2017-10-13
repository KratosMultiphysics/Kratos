#if !defined(KRATOS_MESHLESS_BASE_CONDITION_H_INCLUDED )
#define  KRATOS_MESHLESS_BASE_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


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
class MeshlessBaseCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MeshlessBaseCondition
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessBaseCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshlessBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MeshlessBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MeshlessBaseCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

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
	// A protected default constructor necessary for serialization
	MeshlessBaseCondition() : Condition()
	{
	}
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

	void Jacobian(const Matrix& DN_De, Matrix& Jacobian);

	void JacobianInitial(const Matrix& DN_De,
		Matrix& Jacobian);

	void MappingGeometricToParameter(const Matrix& DN_De,
		const array_1d<double, 2>& Tangents,
		double& JGeometricToParameter);


	void GetBasisVectors(
		const Matrix& DN_De,
		array_1d<double, 3>& g1,
		array_1d<double, 3>& g2,
		array_1d<double, 3>& g3);

	void GetInitialBasisVectors(
		const Matrix& DN_De,
		array_1d<double, 3>& g10,
		array_1d<double, 3>& g20,
		array_1d<double, 3>& g30);

	void CrossProduct(
		array_1d<double, 3>& cross,
		const array_1d<double, 3>& a,
		const array_1d<double, 3>& b);

	array_1d<double, 3> CrossProduct(
		const array_1d<double, 3>& a,
		const array_1d<double, 3>& b);

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
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization


    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }
    
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
    //MeshlessBaseCondition& operator=(const MeshlessBaseCondition& rOther);

    /// Copy constructor.
    //MeshlessBaseCondition(const MeshlessBaseCondition& rOther);


    ///@}

}; // Class MeshlessBaseCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    MeshlessBaseCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshlessBaseCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_BASE_CONDITION_H_INCLUDED  defined 


