//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_STABILIZED_CFD_APPLICATION_H_INCLUDED )
#define  KRATOS_STABILIZED_CFD_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

// application includes
#include "custom_elements/dss.h"
#include "custom_elements/dss_fic.h"
#include "custom_elements/dss_fic_limited.h"
#include "custom_elements/dss_gls.h"
#include "custom_elements/dynss.h"
#include "custom_elements/dss_notau2.h"
#include "custom_elements/dynss_notau2.h"
#include "custom_elements/dss_ps.h"
#include "custom_conditions/dss_face.h"

namespace Kratos {

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
class KratosStabilizedCFDApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosStabilizedCFDApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosStabilizedCFDApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosStabilizedCFDApplication();

	/// Destructor.
	~KratosStabilizedCFDApplication() override {}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	void Register() override;



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
	std::string Info() const override {
		return "KratosStabilizedCFDApplication";
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const override {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	void PrintData(std::ostream& rOStream) const override {
  		KRATOS_WATCH("in my application");
  		KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
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

    /// Developement version of dynamic subscale elements
    const DSS<2> mDSS2D;
    const DSS<3> mDSS3D;
    const DSS<2> mDSS2D4N;
    const DSS<3> mDSS3D8N;

	// Base FIC implementation
    const DSS_FIC<2> mDSS2D_FIC;
    const DSS_FIC<3> mDSS3D_FIC;
    const DSS_FIC<2> mDSS2D4N_FIC;
    const DSS_FIC<3> mDSS3D8N_FIC;

	// FIC with dynamic limiter
    const DSS_FIC_LIMITED<2> mDSS2D_FIC_LIMITED;
    const DSS_FIC_LIMITED<3> mDSS3D_FIC_LIMITED;
    const DSS_FIC_LIMITED<2> mDSS2D4N_FIC_LIMITED;
    const DSS_FIC_LIMITED<3> mDSS3D8N_FIC_LIMITED;

	// Basic GLS
    const DSS_GLS<2> mDSS2D_GLS;
    const DSS_GLS<3> mDSS3D_GLS;
    const DSS_GLS<2> mDSS2D4N_GLS;
    const DSS_GLS<3> mDSS3D8N_GLS;

	// Dynamic subscales
    const DynSS<2> mDynSS2D;
    const DynSS<3> mDynSS3D;
    const DynSS<2> mDynSS2D4N;
    const DynSS<3> mDynSS3D8N;

	// Quasistatic subscales without pressure subscale term
    const DSS_notau2<2> mDSS2D_notau2;
    const DSS_notau2<3> mDSS3D_notau2;
    const DSS_notau2<2> mDSS2D4N_notau2;
    const DSS_notau2<3> mDSS3D8N_notau2;

	// Dynamic subscales without pressure subscale term
    const DYNSS_NOTAU2<2> mDYNSS2D_NOTAU2;
    const DYNSS_NOTAU2<3> mDYNSS3D_NOTAU2;
    const DYNSS_NOTAU2<2> mDYNSS2D4N_NOTAU2;
    const DYNSS_NOTAU2<3> mDYNSS3D8N_NOTAU2;

	// Proposed pressure subscale model
    const DSS_PS<2> mDSS2D_PS;
    const DSS_PS<3> mDSS3D_PS;
    const DSS_PS<2> mDSS2D4N_PS;
    const DSS_PS<3> mDSS3D8N_PS;

	// Boundary condition for elements with skew-symmetric convective term
    const DSSFace<2,2> mDSSFace2D;
    const DSSFace<3,3> mDSSFace3D;
    const DSSFace<3,4> mDSSFace3D4N;

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
	KratosStabilizedCFDApplication& operator=(KratosStabilizedCFDApplication const& rOther);

	/// Copy constructor.
	KratosStabilizedCFDApplication(KratosStabilizedCFDApplication const& rOther);


	///@}

}; // Class KratosStabilizedCFDApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_STABILIZED_CFD_APPLICATION_H_INCLUDED  defined
