//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2014-03-08 08:56:42 $
//
//

#if !defined(KRATOS_DERIVATIVE_RECOVERY_MESHING_TOOLS)
#define  KRATOS_DERIVATIVE_RECOVERY_MESHING_TOOLS

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

#include "includes/model_part.h"

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

/// This class performs simple meshing tools related to the derivative recovery algorithms.
/** @author  Guillermo Casas Gonzalez <gcasas@cimne.upc.edu>
*/

template<std::size_t TDim>
class KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecoveryMeshingTools
{
public:
///@name Type Definitions
///@{
    typedef ModelPart::ElementsContainerType::iterator ElementIteratorType;
    typedef ModelPart::NodesContainerType::iterator    NodeIteratorType;
    typedef Properties PropertiesType;

/// Pointer definition of DerivativeRecoveryMeshingTools
KRATOS_CLASS_POINTER_DEFINITION(DerivativeRecoveryMeshingTools);

///@}
///@name Life Cycle
///@{

/// Default constructor.

DerivativeRecoveryMeshingTools(){}

/// Destructor.
virtual ~DerivativeRecoveryMeshingTools(){}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{

//***************************************************************************************************************
//***************************************************************************************************************

void FillUpEdgesModelPartFromSimplicesModelPart(ModelPart& r_edges_model_part, ModelPart& r_tetra_model_part, std::string element_type);

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************


///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.
virtual std::string Info() const
{
    return "";
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const {}

///@}
///@name Friends
///@{
///@}

protected:

private:


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
DerivativeRecoveryMeshingTools& operator=(DerivativeRecoveryMeshingTools const& rOther);

///@}

}; // Class DerivativeRecoveryMeshingTools

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DerivativeRecoveryMeshingTools<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_DERIVATIVE_RECOVERY_MESHING_TOOLS  defined
