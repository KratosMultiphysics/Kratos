//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MODELER_UTILITIES_H_INCLUDED )
#define  KRATOS_MODELER_UTILITIES_H_INCLUDED

// External includes

// System includes

// Project includes
//#include "includes/kratos_flags.h"
#include "geometries/triangle_2d_3.h"
#include "custom_modelers/spatial_bounding_box.hpp"

#include "pfem_solid_mechanics_application.h"

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
class ModelerUtilities
{
public:
  
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( ModelerUtilities );


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
    ModelerUtilities() {} //

    /// Destructor.
    virtual ~ModelerUtilities() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
   
    void SetDomainLabels (ModelPart& rModelPart);


    void BuildTotalMesh (ModelPart& rModelPart);
    

    void CleanRemovedNodes (ModelPart& rModelPart,ModelPart::IndexType MeshId=0);


    //*******************************************************************************************
    //*******************************************************************************************

    bool CheckSubdomain(Geometry<Node<3> >& vertices);
    
    bool CheckInnerCentre(Geometry<Node<3> >& vertices);
    
    bool CheckOuterCentre(Geometry<Node<3> >& vertices,double &offset_factor);
         

    //*******************************************************************************************
    //*******************************************************************************************

      
    double FindBoundaryH (Node<3>& BoundaryPoint);

    //returns false if it should be removed
    bool ShrankAlphaShape(double alpha_param, Geometry<Node<3> >& vertices,double & offset_factor);

    //returns false if it should be removed
    bool AlphaShape(double alpha_param, Geometry<Node<3> >& vertices);
    
    void CheckParticles (ModelPart& rModelPart,ModelPart::IndexType MeshId=0);
    
    //*******************************************************************************************
    //*******************************************************************************************
   
    inline double CalculateSideLength(PointType& P1,PointType& P2)
    {
      return sqrt( (P1.X()-P2.X())*(P1.X()-P2.X()) + (P1.Y()-P2.Y())*(P1.Y()-P2.Y()) );
    };

    inline double CalculateCircRadius(Geometry< Node<3> >& rGeometry)
    {
      double L1 = CalculateSideLength (rGeometry[0],rGeometry[1]);
      double L2 = CalculateSideLength (rGeometry[1],rGeometry[2]);
      double L3 = CalculateSideLength (rGeometry[2],rGeometry[0]);
      
      
      double Area = rGeometry.Area();
      
      
      double  Rcrit = Area*2/(L1+L2+L3);
      
      return Rcrit;
      
    };
  

    inline double CalculateCircRadius(const double x0, const double y0,
				      const double x1, const double y1,
				      const double x2, const double y2,
				      double& Area)
    {
      
      double L1 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
      double L2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
      double L3 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));     
      
      Area = fabs( 0.5 * ( (x0*y1) - (x0*y2) - (x1*y0) + (x1*y2) + (x2*y0) - (x2*y1) ) );
      
      // std::cout<< " Area "<<Area<<" L1 "<<L1<<" L2 "<<L2<<" L3 "<<L3<<std::endl;
      
      double  Rcrit = Area*2/(L1+L2+L3);
      
      return Rcrit;
    }

    inline double CalculateAverageSideLength(const double x0, const double y0,
					     const double x1, const double y1,
					     const double x2, const double y2) 
    { 
      double length_0 = sqrt( x0*x0 + y0*y0 );
      double length_1 = sqrt( x1*x1 + y1*y1 );
      double length_2 = sqrt( x2*x2 + y2*y2 );
      
      return 0.5*( length_0 + length_1 + length_2 );
    };
  //*******************************************************************************************
  //*******************************************************************************************				   
    
      

    bool CheckConditionInBox(Condition::Pointer& pCond, SpatialBoundingBox& RefiningBox, ProcessInfo& CurrentProcessInfo);
    
    bool CheckElementInBox(Element::Pointer& pElem, SpatialBoundingBox& RefiningBox, ProcessInfo& CurrentProcessInfo);
    
    bool CheckVerticesInBox(Geometry<Node<3> >& rGeom, SpatialBoundingBox& RefiningBox, ProcessInfo& CurrentProcessInfo);
    

    Condition::Pointer FindMasterCondition(Condition::Pointer& pCond, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found);

    Condition::Pointer FindMasterCondition(Condition::Pointer& pCond, PointType& pSlaveNode, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found);
    
        

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
    ModelerUtilities& operator=(ModelerUtilities const& rOther);

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
     
    ///@}

}; // Class ModelerUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      ModelerUtilities& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const ModelerUtilities& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_MODELER_UTILITIES_H_INCLUDED  defined 

