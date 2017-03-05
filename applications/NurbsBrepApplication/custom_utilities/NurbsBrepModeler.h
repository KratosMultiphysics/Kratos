//
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_NURBS_BREP_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_BREP_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <ctime>

// External includes 
//#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/kratos_parameters.h"

#include "BrepModelGeometryReader.h"
#include "BrepModel.h"

namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 

  // Variables definition 



  ///@} 
  ///@name Type Definitions
  ///@{ 
/**
 * Typedefs for search
 */
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
  class NurbsBrepModeler
  {
  public:
    ///@name Type Definitions
    ///@{
    typedef std::vector<BrepModel> BrepModelVector;
    
    /// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(NurbsBrepModeler);

    ///@}
    ///@name Life Cycle 
    ///@{ 
        
   
    /// Constructor.
        NurbsBrepModeler(Parameters& cad_geometry, ModelPart& model_part);

        //NurbsBrepModeler();
    /// Destructor.
        virtual ~NurbsBrepModeler();

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
    //       static const ApplicationCondition  msApplicationCondition; 
    ///@} 
        ///@name Member Variables
    ///@{ 

    Parameters& m_cad_geometry;
    ModelPart& m_model_part;
    //BrepModelGeometryReader m_brep_model_geometry_reader;
    BrepModelVector m_brep_model_vector;

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
    NurbsBrepModeler& operator=(NurbsBrepModeler const& rOther);

    /// Copy constructor.
    NurbsBrepModeler(NurbsBrepModeler const& rOther);


    ///@}    

  }; // Class NurbsBrepModeler 

  ///@} 


  ///@name Type Definitions       
  ///@{ 


  ///@} 
  ///@name Input and output 
  ///@{ 

  ///@} 


}  // namespace Kratos.

#endif // KRATOS_NURBS_BREP_MODELER_APPLICATION_H_INCLUDED  defined 


