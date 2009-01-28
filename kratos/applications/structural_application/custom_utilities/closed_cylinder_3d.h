/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
/* *********************************************************   
*          
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2009-01-15 08:29:58 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

#if !defined(KRATOS_CLOSED_CYLINDER_3D_INCLUDED )
#define  KRATOS_CLOSED_CYLINDER_3D_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"
#include <cmath>

//Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    class ClosedCylinder3D
    {
        public:
            typedef Dof<double> TDofType;
            typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType;
            typedef ModelPart::ElementsContainerType ElementsArrayType;
            typedef double* ContainerType;
            typedef Element::DofsVectorType DofsVectorType;
            typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
            typedef Geometry<Node<3> >::GeometryType GeometryType;
            typedef Geometry<Node<3> >::CoordinatesArrayType CoordinatesArrayType;
            typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
            
            /** 
             * Constructor.
             * @param M starting point of the cylinder, must be midpoint of the circle
             * @param e1 point on the circle above the midpoint
             * @param e2 point on the circle above the midpoint and orthogonal to e2
             * @param e3 point on the direction vector of the cylinder
             * e1, e2, e3 form a othogonal basis                  
             */
            ClosedCylinder3D( Point<3> M, Point<3> e1, Point<3> e2,  Point<3> e3)
            :mM(M),mE1(e1),mE2(e2),mE3(e3)
            {
                r1 = mE1-M;
                r2 = mE2-M;
                r3 = mE3-M;
				r1_normalized= r1/MathUtils<double>::Norm(r1);
                r2_normalized= r2/MathUtils<double>::Norm(r2);
				r3_normalized= r3/MathUtils<double>::Norm(r3);
                mLength= MathUtils<double>::Norm(r3);
                radius = MathUtils<double>::Norm(r1);
                std::cout << "Closed Cylinder created" << std::endl;
                KRATOS_WATCH( mM );
                KRATOS_WATCH( mE1 );
                KRATOS_WATCH( mE2 );
                KRATOS_WATCH( mE3 );
                KRATOS_WATCH( r1 );
                KRATOS_WATCH( r2 );
                KRATOS_WATCH( r3 );
                KRATOS_WATCH( radius );
                KRATOS_WATCH( mLength );
            }
            
			/**
             * Returns the starting point of the cylinder, is the midpoint of the circle
			 * @return mM
             */
            Point<3>& GetCenter()
            {
                return mM;
            }

            /**
             * Returns point on the circle above the midpoint
			 * @return mE1
             */
            array_1d<double,3>& GetPoint1()
            {
                return mE1;
            }
            
            /**
             * Returns point on the circle above the midpoint and orthogonal to e2
			 * @return mE2
             */
            array_1d<double,3>& GetPoint2()
            {
                return mE2;
            }
            
            /**
             * Returns point on the direction vector of the cylinder
			 * @return mE3
             */
            array_1d<double,3>& GetPoint3()
            {
                return mE3;
            }

            /**
             * Returns the direction vector r1
			 * @return r1
             */
            array_1d<double, 3> GetR1()
            {
                return r1;
            }
            
            /**
             * Returns the direction vector r2
			 * @return r2
             */
            array_1d<double, 3> GetR2()
            {
                return r2;
            }

            /**
             * Returns the direction vector r3
			 * @return r3
             */
            array_1d<double, 3> GetR3()
            {
                return r3;
            }

            /** 
             * Calculates a point on the cylinder with the parameters s and t.
             * @param s represents the line segment on the circle
             * @param t represents the distance on the direction vector of the cylinder from M
             * @param result a Point in 3D result(s,t) on the cylinder
             * @return result
             */
            Point<3>& GetPoint( Point<3>& result, double s , double t)
            {
                result = cos( s/radius )*r1 + sin( s/radius )*r2 + (mM+t*r3);
                return result;
            }

            /** 
             * Calculates the partial derivative after t of a point on the cylinder with the parameters s and t 
             * @param result represents the partial derivative after t
             * @return result
             */
            array_1d<double,3>& GetDerivative_t( array_1d<double,3>& result )
            {
                result = r3;
                return result;
            }

			/** 
             * Calculates the partial derivative after s of a point on the cylinder with the parameters s and t 
             * @param result represents the partial derivative after s
             * @return result
             */
            array_1d<double,3>& GetDerivative_s( array_1d<double,3>& result, double s )
            {
                result = (1.0/radius)*(-sin(s/radius)*r1+cos(s/radius)*r2);
                return result;
            }

            /** 
             * Calculates a point on the cylinder cap with the parameters s and t.
             * @param s represents the distance on the direction vector r1 from M
             * @param t represents the distance on the direction vector r2 from M
             * @param result a Point in 3D result(s,t) on the cylinder cap
             * @return result
             */
            Point<3>& GetPointOnCap( Point<3>& result, double s, double t)
            {
                result = r1*s + r2*t + (mM+r3);
                return result;
            }

			/** 
             * Calculates the partial derivative after t of a point on the cylinder cap with the parameters s and t 
             * @param result represents the partial derivative after t
             * @return result
             */
            array_1d<double,3>& GetDerivativeCap_t( array_1d<double,3>& result)
            {
                result = r2;
                return result;
            }

			/** 
             * Calculates the partial derivative after s of a point on the cylinder cap with the parameters s and t 
             * @param result represents the partial derivative after s
             * @return result
             */
            array_1d<double,3>& GetDerivativeCap_s( array_1d<double,3>& result)
            {
                result = r1;
                return result;
            }

            /**
             * Returns the length of the cylinder
			 * @return mLength
             */
            double GetLength()
            {
                return mLength;
            }

            /**
             * Returns the radius of the cylinder
			 * @return radius
             */
            double GetRadius()
            {
                return radius;
            }

            /**
             * Checks weather a point is inside or on the cylinder
			 * @param point that is checked
			 * @retval false if the point is not in or on the cylinder
			 * @retval true if the point is in or on the cylinder
             */
            bool IsInorOn(Point<3>& point)
            {
                Vector distance= point-mM;
                double t= MathUtils<double>::Dot3(distance,
                        r3_normalized);
                double s1= MathUtils<double>::Dot3(distance,
                        r1_normalized);
                double s2= MathUtils<double>::Dot3(distance,
                        r2_normalized);

                if(sqrt(s1*s1+s2*s2)-radius<=0.000001 && t>=-0.000001 && t-mLength<= 0.000001)
                    return true;
                else
                    return false;
            }

            /**
             * Checks weather a point is inside the cylinder
			 * @param point that is checked
			 * @retval false if the point is on or outside the cylinder
			 * @retval true if the point is inside the cylinder
             */
            bool IsIn(Point<3>& point)
            {
                Vector distance= point-mM;
                double t= MathUtils<double>::Dot3(distance,
                        r3_normalized);
                double s1= MathUtils<double>::Dot3(distance,
                        r1_normalized);
                double s2= MathUtils<double>::Dot3(distance,
                        r2_normalized);

                if(sqrt(s1*s1+s2*s2)-radius<-0.000001 && t>0.000001 && t-mLength<-0.000001)
                    return true;
                else
                    return false;
            }

            /**
             * Checks weather a point is inside a bigger cylinder
             * @param point that is checked
             * @retval false if the point is on or outside the cylinder
             * @retval true if the point is inside a bigger cylinder
             */
            bool IsInBigger(Point<3>& point)
            {
                Vector distance= point-mM;
                double t= MathUtils<double>::Dot3(distance,
                        r3_normalized);
                double s1= MathUtils<double>::Dot3(distance,
                        r1_normalized);
                double s2= MathUtils<double>::Dot3(distance,
                        r2_normalized);

                if(sqrt(s1*s1+s2*s2)-radius*2<0.000001 && t>-0.000001 && t-mLength*2<0.000001)
                    return true;
                else
                    return false;
            }

            /**
             * Checks weather a point is on the cylinder
			 * @param point that is checked
			 * @retval false if the point is inside or outside the cylinder
			 * @retval true if the point is on the cylinder
             */
            bool IsOn(Point<3>& point)
            {
                Vector distance= point-mM;
                double t= MathUtils<double>::Dot3(distance,
                        r3_normalized);
                double s1= MathUtils<double>::Dot3(distance,
                        r1_normalized);
                double s2= MathUtils<double>::Dot3(distance,
                        r2_normalized);

                double rad = sqrt(s1*s1+s2*s2);
                if((rad-radius<=0.000001 && rad-radius>=-0.000001) ||
                   (t-mLength<= 0.000001 && t-mLength>= -0.000001))
                    return true;
                else{
                    return false;
                }
            }

            /**
             * Checks weather a point is on the wall of the cylinder
			 * @param point that is checked
			 * @retval false if the point is not on the cylinder wall
			 * @retval true if the point is on the cylinder wall
             */
            bool IsOnWall(Point<3>& point)
            {
                Vector distance= point-mM;
                double t= MathUtils<double>::Dot3(distance,
                        r3_normalized);
                double s1= MathUtils<double>::Dot3(distance,
                        r1_normalized);
                double s2= MathUtils<double>::Dot3(distance,
                        r2_normalized);

                double rad = sqrt(s1*s1+s2*s2);
                if(rad-radius<=0.000001 && rad-radius>=-0.000001 && 
                 !(t-mLength<= 0.000001 && t-mLength>= -0.000001))
                    return true;
                else{
                    return false;
                }
            }

            /**
             * Checks weather a point is on the cap of the cylinder
			 * @param point that is checked
			 * @retval false if the point is not on the cylinder cap
			 * @retval true if the point is on the cylinder cap
             */
            bool IsOnCap(Point<3>& point)
            {
                Vector distance= point-mM;
                double t= MathUtils<double>::Dot3(distance,
                        r3_normalized);
                double s1= MathUtils<double>::Dot3(distance,
                        r1_normalized);
                double s2= MathUtils<double>::Dot3(distance,
                        r2_normalized);

                double rad = sqrt(s1*s1+s2*s2);
                if(t-mLength<= 0.000001 && t-mLength>= -0.000001 && 
                 !(rad-radius<=0.000001 && rad-radius>=-0.000001))
                    return true;
                else{
                    return false;
                }
            }

			/**
             * Returns to a given point the closest point on the plane spanned by the cylinder cap 
			 * @param point for which the closest point on the cap is searched
			 * @return closestpoint
             */
			Point<3> ClosestPointOnCap(Point<3>& point){
				Vector distance= point-mM;
				double t= MathUtils<double>::Dot3(distance,
					   r3_normalized);

				Point<3> closestpoint;

				closestpoint = point + (mLength-t)*r3_normalized;

				return closestpoint;
			}

			/**
             * Returns the closest point on the cylinder wall to a given point if the projected point will be on the closed cylinder
			 * @param point for which the closest point on the wall is searched
			 * @return closestpoint
             */
            Point<3> ClosestPointOnWall(Point<3>& point){
				Vector distance= point-mM;
				double t= MathUtils<double>::Dot3(distance,
					   r3_normalized);
				double s1= MathUtils<double>::Dot3(distance,
					   r1_normalized);
				double s2= MathUtils<double>::Dot3(distance,
					   r2_normalized);

				Point<3> closestpoint;

                if (t <= mLength + 0.000001 && t >= -0.000001){
                //point is above cylinder wall
                    Vector helpvector = s1 * r1_normalized + s2 * r2_normalized;
                    helpvector = helpvector / MathUtils<double>::Norm(helpvector);
                    closestpoint = mM + radius * (helpvector) + t*r3_normalized;
                }
                else
				   closestpoint = point;
				
				return closestpoint;
			}

            /**
             * Returns the closest point on the cylinder to a given point 
			 * @param point for which the closest point is searched
			 * @return closestpoint
             */
			Point<3> ClosestPoint(Point<3>& point){
                Vector distance= point-mM;
                double t= MathUtils<double>::Dot3(distance,
                        r3_normalized);
                double s1= MathUtils<double>::Dot3(distance,
                        r1_normalized);
                double s2= MathUtils<double>::Dot3(distance,
                        r2_normalized);

                Point<3> closestpoint;
				double rad = sqrt(s1*s1+s2*s2);

				if (mLength < t){
                //outside cylidner cap				
                    if (radius < rad){
                    //outside cylinder wall and cap
                        Vector helpvector = s1 * r1_normalized + s2 * r2_normalized;
                        helpvector = helpvector / MathUtils<double>::Norm(helpvector);
                        closestpoint = mM + radius * (helpvector) + mLength*r3_normalized;
                    }
					else
                    //only ouside cylinder cap
                        closestpoint = point + (mLength-t)*r3_normalized;
				}
				else if (radius < rad){
				//only outside cylinder wall
					Vector helpvector = s1 * r1_normalized + s2 * r2_normalized;
					helpvector = helpvector / MathUtils<double>::Norm(helpvector);
		    		closestpoint = mM + radius * (helpvector) + t*r3_normalized;
				}
				else{
				//node is inside cylidner
					if(fabs(mLength-t) < fabs(radius-rad)){
					//node is closer to cap than on wall
						closestpoint = point + (mLength-t)*r3_normalized;
					}
                	else{
					//node is closer on wall than on cap
						Vector helpvector = s1 * r1_normalized + s2 * r2_normalized;
						helpvector = helpvector / MathUtils<double>::Norm(helpvector);
		    			closestpoint = mM + radius * (helpvector) + t*r3_normalized;
					}
				}
				return closestpoint;
			}

        private:
            
            Point<3> mM;
            Point<3> mE1;
            Point<3> mE2;
            Point<3> mE3;
            array_1d<double, 3> r1;
            array_1d<double, 3> r2;
            array_1d<double, 3> r3;
			Vector r1_normalized;
			Vector r2_normalized;
			Vector r3_normalized;
            double radius;
            double mLength;
    };//Class ClosedCylinder3D
}//namespace Kratos.

#endif /* KRATOS_CIRCLE_3D_INCLUDED  defined */
