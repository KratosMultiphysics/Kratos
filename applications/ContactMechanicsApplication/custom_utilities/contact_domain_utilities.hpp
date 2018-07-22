//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CONTACT_DOMAIN_UTILITIES_H_INCLUDED )
#define  KRATOS_CONTACT_DOMAIN_UTILITIES_H_INCLUDED

// External includes

// System includes

// Project includes
#include "utilities/math_utils.h"

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
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) ContactDomainUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ContactDomainUtilities
    KRATOS_CLASS_POINTER_DEFINITION( ContactDomainUtilities );

    /**
     * Flags related to the condition computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_FRICTION_FORCES );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_FRICTION_STIFFNESS );

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_NODAL_CONTACT_FORCES );

    ///Tensor order 1 definition
    //typedef BoundedVector<double, 3>            PointType;
    typedef array_1d<double, 3>                    PointType;


    typedef struct
    {
        double L;    //base side lentgh
        double A;    //distance 2-3
        double B;    //distance 1-2

    } BaseLengths;


    typedef struct
    {
        PointType Normal;        //normal direction
        PointType Tangent;       //tangent direction

    } SurfaceVector;


    typedef struct
    {
        double Normal;        //normal component
        double Tangent;       //tangent component

    } SurfaceScalar;


    typedef struct
    {
      double Covariant;       //covariant component
      double Contravariant;   //contravariant component

    } ScalarBaseType;


    typedef struct
    {
      Matrix    Metric;      //metric of the base
      PointType DirectionA;  //reference base direction a
      PointType DirectionB;  //reference base direction b

    } SurfaceBase;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ContactDomainUtilities() {} //

    /// Destructor.
    virtual ~ContactDomainUtilities() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    inline double CalculateVolume(const double x0, const double y0,
				  const double x1, const double y1,
				  const double x2, const double y2)
	{
		return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
	}


    //************************************************************************************
    //************************************************************************************

    inline bool CalculatePosition(const double x0, const double y0,
				  const double x1, const double y1,
				  const double x2, const double y2,
				  const double xc, const double yc)
	{
		double area = CalculateVolume(x0,y0,x1,y1,x2,y2);

		//std::cout<<" Area "<<area<<std::endl;

		if(area < 1e-15)
		{
			//KRATOS_THROW_ERROR( std::logic_error,"element with zero area found", "" )
			std::cout<<"element with zero area found: "<<area<<" position ("<<x0<<", "<<y0<<") ("<<x1<<", "<<y1<<") ("<<x2<<", "<<y2<<") "<<std::endl;
		}

		PointType N;

		N[0] = CalculateVolume(x1,y1,x2,y2,xc,yc)  / area;
		N[1] = CalculateVolume(x2,y2,x0,y0,xc,yc)  / area;
		N[2] = CalculateVolume(x0,y0,x1,y1,xc,yc)  / area;

		double tol = 1e-3;
		double upper_limit = 1.0+tol;
		double lower_limit = -tol;

		if(N[0] >= lower_limit && N[1] >= lower_limit && N[2] >= lower_limit && N[0] <= upper_limit && N[1] <= upper_limit && N[2] <= upper_limit) //if the xc yc is inside the triangle
			return true;

		return false;
	}

    //************************************************************************************
    //************************************************************************************

    inline bool CalculateObtuseAngle(const double x0, const double y0,
				     const double x1, const double y1,
				     const double xc, const double yc)
	{

		double side0 = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) ); //master side
		double side1 = sqrt( (x0-xc)*(x0-xc) + (y0-yc)*(y0-yc) );
		double side2 = sqrt( (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) );

		double cos_angle = 0;
		double aux       = (2*side1*side0);
		if(aux!=0)
			cos_angle = ((side1*side1) + (side0*side0) - (side2*side2)) / aux;

		if(cos_angle<(-0.1))
			return true;

		aux       = (2*side2*side0);
		if(aux!=0)
			cos_angle = ((side2*side2) + (side0*side0) - (side1*side1)) / aux;

		if(cos_angle<(-0.1))
			return true;

		return false;
	}


     void CalculateBaseArea (double& area,
			     double& a,
			     double& b,
			     double& c);

     void CalculateLineIntersection (double& a,
				     const PointType& P1,
				     const PointType& P2,
				     const PointType& V1,
				     const PointType& V2);


    void  CalculateEdgeDistances (std::vector<BaseLengths>& BaseVector,
				  const PointType& P1,
				  const PointType& P2,
				  const PointType& PS1,
				  const PointType& PS2,
				  const PointType& Normal);

    void  CalculateBaseDistances (std::vector<BaseLengths>& BaseVector,
				  const PointType& P1,
				  const PointType& P2,
				  const PointType& P3,
				  const PointType& PS,
				  const PointType& Normal);


    void  CalculateBaseDistances (BaseLengths& Base,
				  const PointType& P1,
				  const PointType& P2,
				  const PointType& PS,
				  const PointType& Normal);

    PointType & CalculateSurfaceNormal(PointType& Normal,
				       const PointType& P1,
				       const PointType& P2);


    PointType & CalculateFaceNormal(PointType& Normal,
				    const PointType& P1,
				    const PointType& P2);


    PointType & CalculateFaceTangent(PointType& Tangent,
				     const PointType& P1,
				     const PointType& P2);



    PointType & CalculateFaceTangent(PointType& Tangent,
				     PointType& Normal);


    void GetOppositeEdge(unsigned int& i, unsigned int& j, unsigned int& k, unsigned int& l);


    void BuildEdgeVector(std::vector<std::vector<std::vector<unsigned int> > >& rEdges);
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
    ContactDomainUtilities& operator=(ContactDomainUtilities const& rOther);

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

}; // Class ContactDomainUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      ContactDomainUtilities& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const ContactDomainUtilities& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_CONTACT_DOMAIN_UTILITIES_H_INCLUDED  defined
