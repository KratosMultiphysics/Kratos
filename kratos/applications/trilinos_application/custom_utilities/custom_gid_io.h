/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2009-01-13 16:51:36 $
//   Revision:            $Revision: 1.3 $
//
//

#if !defined(KRATOS_GID_TRILINOS_GAUSS_POINT_CONTAINER_H_INCLUDED)
#define  KRATOS_GID_TRILINOS_GAUSS_POINT_CONTAINER_H_INCLUDED

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
#include "trilinos_application.h"

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
class TrilinosGidGaussPointsContainer : public GidGaussPointsContainer
{
public:
    typedef GidGaussPointsContainer BaseType;
    ///Constructor
    TrilinosGidGaussPointsContainer( const char* gp_title, KratosGeometryFamily geometryFamily,
                                     GiD_ElementType gid_element_type,
                                     int number_of_integration_points,
                                     std::vector<int> index_container )
        :BaseType( gp_title, geometryFamily, gid_element_type, number_of_integration_points,
                   index_container) {}

    virtual void PrintResults( GiD_FILE ResultFile, Variable<double> rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int parameter_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);
            if(rVariable == SATURATION)
            {
                GiD_fBeginRangeTable(ResultFile, "saturation");
                GiD_fWriteMinRange( ResultFile, 0.0001, "unsaturated" );
                GiD_fWriteRange( ResultFile, 0.0001, 0.9999, "partially saturated" );
                GiD_fWriteMaxRange(ResultFile, 0.9999, "saturated");
                GiD_fEndRangeTable(ResultFile );
                GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), "Kratos", SolutionTag,
                                 GiD_Scalar, GiD_OnGaussPoints, mGPTitle,
                                 "saturation", 0, NULL );
            }
            else
                GiD_fBeginResult( ResultFile,  (char *)(rVariable.Name()).c_str(), "Kratos", SolutionTag,
                                 GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<double> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         r_model_part.GetProcessInfo() );
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
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), ValuesOnIntPoint[index] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile );
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<array_1d<double, 3> > rVariable, ModelPart& r_model_part, double
                               SolutionTag, unsigned int parameter_index )
    {
        if( mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);
            GiD_fBeginResult( ResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                             SolutionTag, GiD_Vector,
                             GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<array_1d<double,3> > ValuesOnIntPoint(mSize);


            for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); it++ )
            {
                if( ! it->GetValue( IS_INACTIVE ) )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    KRATOS_WATCH( ValuesOnIntPoint[0] );
                    for(unsigned int i=0; i<mIndexContainer.size(); i++)
                    {
                        GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[i][0],
                                         ValuesOnIntPoint[i][1], ValuesOnIntPoint[i][2] );
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<Vector> rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int parameter_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);

            if( rVariable == INSITU_STRESS )
                GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), "Kratos", SolutionTag,
                                 GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            else if( (rVariable == MATERIAL_PARAMETERS) || (rVariable == INTERNAL_VARIABLES) )
            {
                std::stringstream param_index;
                param_index << parameter_index;
                GiD_fBeginResult( ResultFile, (char *)(rVariable.Name() + param_index.str() ).c_str(),
                                 "Kratos", SolutionTag, GiD_Scalar, GiD_OnGaussPoints,
                                 mGPTitle, NULL, 0, NULL );
            }
            else
                GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), "Kratos", SolutionTag,
                                 GiD_Vector, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<Vector> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if( rVariable == INSITU_STRESS )
                            {
                                if(ValuesOnIntPoint[i].size() ==6 )
                                    GiD_fWrite3DMatrix(ResultFile, it->Id(), ValuesOnIntPoint[index](0),
                                                       ValuesOnIntPoint[index](1), ValuesOnIntPoint[index](2),
                                                       ValuesOnIntPoint[index](3), ValuesOnIntPoint[index](4),
                                                       ValuesOnIntPoint[index](5) );
                            }
                            else if( (rVariable == MATERIAL_PARAMETERS)
                                     || (rVariable == INTERNAL_VARIABLES) )
                            {
                                double value = 0.0;
                                if( ValuesOnIntPoint[index].size() > parameter_index )
                                    value = ValuesOnIntPoint[index][parameter_index];
                                GiD_fWriteScalar( ResultFile, it->Id(), value );
                            }
                            else if( ValuesOnIntPoint[0].size() == 3 )
                                GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
                                                 ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if( rVariable == INSITU_STRESS )
                            {
                                if(ValuesOnIntPoint[i].size() ==6 )
                                    GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0),
                                                       ValuesOnIntPoint[index](1), ValuesOnIntPoint[index](2),
                                                       ValuesOnIntPoint[index](3), ValuesOnIntPoint[index](4),
                                                       ValuesOnIntPoint[index](5) );
                            }
                            else if( (rVariable == MATERIAL_PARAMETERS)
                                     || (rVariable == INTERNAL_VARIABLES) )
                            {
                                double value = 0.0;
                                if( ValuesOnIntPoint[index].size() > parameter_index )
                                    value = ValuesOnIntPoint[index][parameter_index];
                                GiD_fWriteScalar(ResultFile,  it->Id(), value );
                            }
                            else if( ValuesOnIntPoint[0].size() == 3 )
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
                               double SolutionTag, unsigned int parameter_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile);
            GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), "Kratos", SolutionTag,
                             GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<Matrix> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if(ValuesOnIntPoint[index].size1() ==3
                                    && ValuesOnIntPoint[index].size2() ==3)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                                   ValuesOnIntPoint[index](0,2) );
                            if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==6)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                                   ValuesOnIntPoint[index](0,5) );
                            if(ValuesOnIntPoint[index].size1() == 1
                                    && ValuesOnIntPoint[index].size2() == 3)
                                GiD_fWrite2DMatrix( ResultFile,  it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2) );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); it++ )
                {
                    if( ! it->GetValue( IS_INACTIVE ) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         r_model_part.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++)
                        {
                            int index = mIndexContainer[i];
                            if(ValuesOnIntPoint[index].size1() ==3
                                    && ValuesOnIntPoint[index].size2() ==3)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                                   ValuesOnIntPoint[index](0,2) );
                            if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==6)
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                                   ValuesOnIntPoint[index](0,5) );
                            if(ValuesOnIntPoint[index].size1() == 1
                                    && ValuesOnIntPoint[index].size2() == 3)
                                GiD_fWrite2DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2) );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

};//class TrilinosGidGaussPointsContainer

}// namespace Kratos.

#endif // KRATOS_GID_TRILINOS_GAUSS_POINT_CONTAINER_H_INCLUDED defined 
