/*
==============================================================================
KratosMultiScaleApplication
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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-12-30 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(KRATOS_GID_EXTENDED_GAUSS_POINT_CONTAINER_H_INCLUDED)
#define  KRATOS_GID_EXTENDED_GAUSS_POINT_CONTAINER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>

// External includes
#include "gidpost/source/gidpost.h"

// Project includes
#include "includes/define.h"
#include "includes/gid_io.h"
#include "geometries/geometry_data.h"

#include "multiscale_application.h"

#define tet10_a 0.108103018168070
#define tet10_b 0.445948490915965
#define tet10_c 0.816847572980459

namespace Kratos
{
/**
 * Type definitions
 */
typedef ModelPart::ElementsContainerType ElementsArrayType;
typedef ModelPart::NodesContainerType NodesArrayType;
typedef ModelPart::ConditionsContainerType ConditionsArrayType;
typedef GeometryData::IntegrationMethod IntegrationMethodType;
typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;

/**
 * Auxiliary class to store gauss point containers and perform result printing
 * on gauss points
 */
class GidExtendedGaussPointsContainer
{
public:
    ///Constructor
    GidExtendedGaussPointsContainer( const char * gp_title, KratosGeometryFamily geometryFamily,
                             GiD_ElementType gid_element_type,
                             int number_of_integration_points,
                             std::vector<int> index_container )
        :
        mGPTitle(gp_title),mKratosElementFamily(geometryFamily),
        mGidElementFamily(gid_element_type), mSize(number_of_integration_points),
        mIndexContainer(index_container) {}

    bool AddElement( const ModelPart::ElementsContainerType::iterator pElemIt )
    {
        KRATOS_TRY
        if( pElemIt->GetGeometry().GetGeometryFamily() == mKratosElementFamily
                && pElemIt->GetGeometry().IntegrationPoints(
                    pElemIt->GetIntegrationMethod() ).size() == mSize )
        {
            mMeshElements.push_back( *(pElemIt.base() ) );
            return true;
        }
        else return false;
        KRATOS_CATCH("")
    }

    bool AddCondition( const ModelPart::ConditionsContainerType::iterator pCondIt )
    {
        KRATOS_TRY
        if( pCondIt->GetGeometry().GetGeometryFamily() == mKratosElementFamily
                && pCondIt->GetGeometry().IntegrationPoints().size() == mSize )
        {
            mMeshConditions.push_back( *(pCondIt.base() ) );
            return true;
        }
        else return false;
        KRATOS_CATCH("")
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<double> rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);

            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );

            std::vector<double> ValuesOnIntPoint(mSize);

            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); it != mMeshElements.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), ValuesOnIntPoint[index] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), ValuesOnIntPoint[index] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<array_1d<double,3> > rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
			WriteGaussPoints(ResultFile);
			
            GiD_fBeginResult( ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                              GiD_Vector, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
							 
            std::vector<array_1d<double,3> > ValuesOnIntPoint(mSize,ZeroVector(3));
			
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); it != mMeshElements.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[index][0], ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[index][0], ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<array_1d<double,6> > rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);

            GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), ( char*)("Kratos"),
                             SolutionTag, GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );

            std::vector<array_1d<double, 6> > ValuesOnIntPoint(mSize);

            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); it != mMeshElements.end(); ++it )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        //KRATOS_WATCH(it->Id())
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
                                               ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2],
                                               ValuesOnIntPoint[index][3], ValuesOnIntPoint[index][4],
                                               ValuesOnIntPoint[index][5] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
                                               ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2],
                                               ValuesOnIntPoint[index][3], ValuesOnIntPoint[index][4],
                                               ValuesOnIntPoint[index][5] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<Vector> rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);

            GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Vector, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );

            std::vector<Vector> ValuesOnIntPoint(mSize);

            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); it != mMeshElements.end(); ++it )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if( ValuesOnIntPoint[0].size() == 3 )
								GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
												  ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if( ValuesOnIntPoint[0].size() == 3 )
                                GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
                                                 ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<Matrix> rVariable, ModelPart& r_model_part,
                               double SolutionTag, int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);

            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"),
                             SolutionTag, GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );

            std::vector<Matrix> ValuesOnIntPoint(mSize);

            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); it != mMeshElements.end(); ++it )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if(ValuesOnIntPoint[index].size1() ==3 && ValuesOnIntPoint[index].size2() ==3)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                                   ValuesOnIntPoint[index](0,2) );
							else if(ValuesOnIntPoint[index].size1() ==2 && ValuesOnIntPoint[index].size2() ==2)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), 0.0,
                                                   ValuesOnIntPoint[index](0,1), 0.0, 0.0);
                            else if(ValuesOnIntPoint[index].size1() ==1 && ValuesOnIntPoint[index].size2() ==3)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), 0.0,
                                                   ValuesOnIntPoint[index](0,2), 0.0, 0.0);
                            else if(ValuesOnIntPoint[index].size1() ==1 && ValuesOnIntPoint[index].size2() ==4)
 								GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), 0.0, 0.0);
							else if(ValuesOnIntPoint[index].size1() ==1 && ValuesOnIntPoint[index].size2() ==6)
								GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
												   ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                                   ValuesOnIntPoint[index](0,5) );

                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if(ValuesOnIntPoint[index].size1() ==3 && ValuesOnIntPoint[index].size2() ==3)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                                   ValuesOnIntPoint[index](0,2) );
                            else if(ValuesOnIntPoint[index].size1() ==1 && ValuesOnIntPoint[index].size2() ==6)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                                   ValuesOnIntPoint[index](0,5) );
                            else if(ValuesOnIntPoint[index].size1() ==1 && ValuesOnIntPoint[index].size2() ==3)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), 0.0,
                                                   ValuesOnIntPoint[index](0,2), 0.0, 0.0);
							else if(ValuesOnIntPoint[index].size1() ==1 && ValuesOnIntPoint[index].size2() ==4)
 								GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), 0.0, 0.0);
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    void Reset()
    {
        mMeshElements.clear();
        mMeshConditions.clear();
    }

protected:

    void WriteGaussPoints(GiD_FILE MeshFile)
    {
        //setting up gauss points
        if( mGidElementFamily == GiD_Tetrahedra && mSize == 4 )
        {
            GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Tetrahedra, NULL, 4, 0, 0 );
            GiD_fWriteGaussPoint3D( MeshFile, 0.58541020,0.13819660,0.13819660 );
            GiD_fWriteGaussPoint3D( MeshFile, 0.13819660,0.58541020,0.13819660 );
            GiD_fWriteGaussPoint3D( MeshFile, 0.13819660,0.13819660,0.58541020 );
            GiD_fWriteGaussPoint3D( MeshFile, 0.13819660,0.13819660,0.13819660 );
            GiD_fEndGaussPoint(MeshFile);
        }
        if( mGidElementFamily == GiD_Tetrahedra && mSize == 5 )
        {
            GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Tetrahedra, NULL, 4, 0, 0 );
            GiD_fWriteGaussPoint3D( MeshFile, 1.0/6.0, 1.0/6.0, 1.0/6.0 );
            GiD_fWriteGaussPoint3D( MeshFile, 1.0/2.0, 1.0/6.0, 1.0/6.0 );
            GiD_fWriteGaussPoint3D( MeshFile, 1.0/6.0, 1.0/2.0, 1.0/6.0 );
            GiD_fWriteGaussPoint3D( MeshFile, 1.0/6.0, 1.0/6.0, 1.0/2.0 );
            GiD_fEndGaussPoint(MeshFile);
        }
        else if( mGidElementFamily == GiD_Tetrahedra && mSize == 10 )
        {
            GiD_fBeginGaussPoint(MeshFile, "tet10_element_gp", GiD_Tetrahedra, NULL, 10, 0, 0);
            GiD_fWriteGaussPoint3D(MeshFile, tet10_a,  tet10_a,  tet10_a );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_c,  tet10_a,  tet10_a );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_a,  tet10_c,  tet10_a );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_a,  tet10_a,  tet10_c );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_b,  tet10_a,  tet10_a );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_b,  tet10_b,  tet10_a );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_a,  tet10_b,  tet10_a );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_a,  tet10_a,  tet10_b );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_b,  tet10_a,  tet10_b );
            GiD_fWriteGaussPoint3D(MeshFile, tet10_a,  tet10_b,  tet10_b );
            GiD_fEndGaussPoint(MeshFile);
        }
        else if( mGidElementFamily == GiD_Tetrahedra && mSize == 11 )
        {
            GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Tetrahedra, NULL, 4, 0, 0 );
            //GiD_WriteGaussPoint3D( 1.0/4.0, 1.0/4.0, 1.0/4.0 );
            GiD_fWriteGaussPoint3D( MeshFile, 1.0/14.0, 1.0/14.0, 1.0/14.0 );
            GiD_fWriteGaussPoint3D( MeshFile, 11.0/14.0, 1.0/14.0, 1.0/14.0 );
            GiD_fWriteGaussPoint3D( MeshFile, 1.0/14.0, 11.0/14.0, 1.0/14.0 );
            GiD_fWriteGaussPoint3D( MeshFile, 1.0/14.0, 1.0/14.0, 11.0/14.0 );
            //GiD_WriteGaussPoint3D( (1.0+std::sqrt(5.0/14.0))/4.0, (1.0-std::sqrt(5.0/14.0))/4.0, (1.0-std::sqrt(5.0/14.0))/4.0 );
            //GiD_WriteGaussPoint3D( (1.0+std::sqrt(5.0/14.0))/4.0,(1.0+std::sqrt(5.0/14.0))/4.0, (1.0-std::sqrt(5.0/14.0))/4.0 );
            //GiD_WriteGaussPoint3D( (1.0-std::sqrt(5.0/14.0))/4.0,(1.0+std::sqrt(5.0/14.0))/4.0, (1.0-std::sqrt(5.0/14.0))/4.0 );
            //GiD_WriteGaussPoint3D( (1.0-std::sqrt(5.0/14.0))/4.0,(1.0-std::sqrt(5.0/14.0))/4.0, (1.0+std::sqrt(5.0/14.0))/4.0 );
            //GiD_WriteGaussPoint3D( (1.0+std::sqrt(5.0/14.0))/4.0,(1.0-std::sqrt(5.0/14.0))/4.0, (1.0+std::sqrt(5.0/14.0))/4.0 );
            //GiD_WriteGaussPoint3D( (1.0-std::sqrt(5.0/14.0))/4.0,(1.0+std::sqrt(5.0/14.0))/4.0, (1.0+std::sqrt(5.0/14.0))/4.0 );
            GiD_fEndGaussPoint(MeshFile);
        }
        else if ( mGidElementFamily == GiD_Triangle && mSize == 3 )
        {
            GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Triangle, NULL, 3, 0, 0 );
            GiD_fWriteGaussPoint2D( MeshFile, 1.0/6.0, 1.0/6.0 );
            GiD_fWriteGaussPoint2D( MeshFile, 2.0/3.0, 1.0/6.0 );
            GiD_fWriteGaussPoint2D( MeshFile, 1.0/6.0, 2.0/3.0 );
            GiD_fEndGaussPoint(MeshFile);
        }
        else
        {
            GiD_fBeginGaussPoint(MeshFile, mGPTitle, mGidElementFamily, NULL,
                                mSize, 0, 1);
            GiD_fEndGaussPoint(MeshFile);
        }
    }
    

    ///member variables
    const char * mGPTitle;
    KratosGeometryFamily mKratosElementFamily;
    GiD_ElementType mGidElementFamily;
    unsigned int mSize;
    std::vector<int> mIndexContainer;
    ModelPart::ElementsContainerType mMeshElements;
    ModelPart::ConditionsContainerType mMeshConditions;

};

}

#endif // KRATOS_GID_EXTENDED_GAUSS_POINT_CONTAINER_H_INCLUDED defined 
