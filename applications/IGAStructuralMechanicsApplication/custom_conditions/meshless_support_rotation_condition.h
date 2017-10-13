#if !defined(KRATOS_MESHLESS_SUPPORT_ROTATION_CONDITION_H_INCLUDED )
#define  KRATOS_MESHLESS_SUPPORT_ROTATION_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "solid_mechanics_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application_variables.h"
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
class MeshlessSupportRotationCondition
    : public MeshlessBaseCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MeshlessSupportRotationCondition
    KRATOS_CLASS_POINTER_DEFINITION(MeshlessSupportRotationCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshlessSupportRotationCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MeshlessSupportRotationCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MeshlessSupportRotationCondition();
	
    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

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
    ///@}
    ///@name Member Variables
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    MeshlessSupportRotationCondition() : MeshlessBaseCondition()
    {
    }

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
	void CalculateRotation(const Matrix &ShapeFunctionDerivatives,
		Vector &Phi_r, Matrix &Phi_rs, array_1d<double, 2>& Phi, 
		const array_1d<double, 2> &Tangents);
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
    //MeshlessSupportRotationCondition& operator=(const MeshlessSupportRotationCondition& rOther);

    /// Copy constructor.
    //MeshlessSupportRotationCondition(const MeshlessSupportRotationCondition& rOther);
    ///@}

}; // Class MeshlessSupportRotationCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    MeshlessSupportRotationCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshlessSupportRotationCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_SUPPORT_ROTATION_CONDITION_H_INCLUDED  defined 


