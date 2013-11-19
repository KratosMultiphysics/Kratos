//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONTACT_DOMAIN_2D_MODELER_H_INCLUDED )
#define  KRATOS_CONTACT_DOMAIN_2D_MODELER_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_modelers/triangle_mesh_2D_modeler.hpp"

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
class ContactDomain2DModeler 
  : public TriangleMesh2DModeler
{
protected:

    struct ContactVariables
    {
      double  offset_factor;  
      double  penalty_stab;
      double  mu_static;
      double  mu_dynamic;
      
      unsigned int    friction_active;
      unsigned int    penalty_contact;
    };

public:

 
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( ContactDomain2DModeler );


    typedef Node<3>                                                       PointType;
    typedef Node<3>::Pointer                                       PointPointerType;
    typedef std::vector<PointPointerType>                        PointPointerVector;
    typedef ModelPart::PropertiesType                                PropertiesType;
    typedef ModelPart::PropertiesContainerType              PropertiesContainerType;
    typedef ModelPart::MeshType                                            MeshType;
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ContactDomain2DModeler() {} //

    /// Destructor.
    virtual ~ContactDomain2DModeler() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //*******************************************************************************************
    //*******************************************************************************************

    void TransferContactBoundaryData(ModelPart& rModelPart, bool initial);


    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateContactMesh (ModelPart& rModelPart,
			      Element   const& rReferenceElement,
			      Condition const& rReferenceCondition,
			      bool constrained = false,
			      double my_alpha  = 1.4,
			      double h_factor  = 0.5,
			      double my_offset = 1.0,
			      double penalty_stab = 0.01,
			      bool friction_active = false,
			      double mu_static = 0.3,
			      double mu_dynamic = 0.2,
			      bool penalty_contact = false);

       


    //*******************************************************************************************
    //*******************************************************************************************
    void GenerateContactDT(ModelPart& rModelPart,
			   TriangleMesh2DModeler::MeshingVariables& rMeshingVariables,
			   ContactVariables& rContactVariables);
    


    //*******************************************************************************************
    //*******************************************************************************************
    void GenerateContactCDT(ModelPart& rModelPart,
			    TriangleMesh2DModeler::MeshingVariables& rMeshingVariables,
			    ContactVariables& rContactVariables);



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
    ///@name Private Operators
    ///@{

    /// Assignment operator.
    ContactDomain2DModeler& operator=(ContactDomain2DModeler const& rOther);

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
    ///@name Unaccessible methods
    ///@{

    void SetTriangulateShrankNodes   (ModelPart& rModelPart,
				      ModelPart::NodesContainerType&  rBoundaryNodes,
				      TriangleMesh2DModeler::MeshingVariables& rMeshingVariables,
				      ContactVariables& rContactVariables,
				      struct triangulateio& in,
				      struct triangulateio& out);
    


    
  
    //Clear contact conditions in model_part
    void ClearContactConditions (ModelPart& rModelPart);

    //Set contact elements in model_part after the Delaunay Tesselation
    void SetContactConditions (ModelPart& rModelPart,
			       ModelPart::NodesContainerType& rBoundaryNodes,
			       TriangleMesh2DModeler::MeshingVariables& rMeshingVariables,
			       ContactVariables& rContactVariables,
			       struct triangulateio &out);
    ///@}

}; // Class ContactDomain2DModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      ContactDomain2DModeler& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const ContactDomain2DModeler& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_CONTACT_DOMAIN_2D_MODELER_H_INCLUDED  defined 


