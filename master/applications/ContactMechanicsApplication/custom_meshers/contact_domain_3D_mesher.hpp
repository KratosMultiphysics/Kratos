//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CONTACT_DOMAIN_3D_MESHER_H_INCLUDED )
#define  KRATOS_CONTACT_DOMAIN_3D_MESHER_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_meshers/tetrahedral_mesh_3D_mesher.hpp"

///VARIABLES used:
//Data:
//(props)
//StepData:
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)


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
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) ContactDomain3DMesher
  : public TetrahedralMesh3DMesher
{
protected:

    struct ContactVariables
    {
      double  OffsetFactor;
      double  PenaltyParameter;
      double  StabilityParameter;
      double  StaticFrictionCoefficient;
      double  DynamicFrictionCoefficient;

      unsigned int    FrictionFlag;
      unsigned int    PenaltyContactFlag;
    };

public:


    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( ContactDomain3DMesher );

    ///Tensor order 1 definition
    //typedef BoundedVector<double, 3>                       PointType;
    typedef array_1d<double, 3>                               PointType;

    typedef MesherUtilities::MeshingInfoParameters  InfoParametersType;
    typedef MesherUtilities::MeshingParameters   MeshingParametersType;
    typedef MesherUtilities::RefiningParameters   RefineParametersType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ContactDomain3DMesher() { std::cout<<"3D mesher for Contact"<<std::endl; } //

    /// Destructor.
    virtual ~ContactDomain3DMesher() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
    virtual std::string Info() const override
    {
	return "";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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
    ContactDomain3DMesher& operator=(ContactDomain3DMesher const& rOther);

    ///@}
    ///@name Private Operations
    ///@{


    //set nodes to a mesh
    void SetNodes(ModelPart& rModelPart,
		  MeshingParametersType& rMeshingVariables) override;


    //set faces in the triangulateio before the Delaunay Tesselation
    void SetFaces ( ModelPart &rModelPart,
		    MeshingParametersType & rMeshingVariables,
		    tetgenio &in ) override;

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}

}; // Class ContactDomain3DMesher

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      ContactDomain3DMesher& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const ContactDomain3DMesher& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_CONTACT_DOMAIN_3D_MESHER_H_INCLUDED  defined
