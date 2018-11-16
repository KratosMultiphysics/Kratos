//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Pooyan Dadvand
//

#if !defined(KRATOS_GID_GAUSS_POINT_CONTAINER_H_INCLUDED)
#define  KRATOS_GID_GAUSS_POINT_CONTAINER_H_INCLUDED

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
#include "geometries/geometry_data.h"

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
class GidGaussPointsContainer
{
public:
    ///Constructor
    GidGaussPointsContainer( const char * gp_title, KratosGeometryFamily geometryFamily,
                             GiD_ElementType gid_element_type,
                             int number_of_integration_points,
                             std::vector<int> index_container )
        :
        mGPTitle(gp_title),mKratosElementFamily(geometryFamily),
        mGidElementFamily(gid_element_type), mSize(number_of_integration_points),
        mIndexContainer(index_container) {}

    ///Destructor
  virtual ~GidGaussPointsContainer(){};

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
        if( pCondIt->GetGeometry().GetGeometryFamily() == mKratosElementFamily &&
            pCondIt->GetGeometry().IntegrationPoints(pCondIt->GetIntegrationMethod()).size() == mSize )
        {
            mMeshConditions.push_back( *(pCondIt.base() ) );
            return true;
        }
        else return false;
        KRATOS_CATCH("")
    }


//            virtual void PrintResults( Variable<array_1d<double,3> > rVariable, ModelPart& rModelPart,
//                                        double SolutionTag, unsigned int ValueIndex )

    virtual void PrintFlagsResults(
        GiD_FILE ResultFile,
        Kratos::Flags rFlag,
        std::string rFlagName,
        ModelPart& rModelPart,
        double SolutionTag
        )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 ) {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult(ResultFile,  (char *)(rFlagName).c_str(), (char *)("Kratos"), SolutionTag, GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); it++ ) {
                    const double double_flag = static_cast<double>(it->Is(rFlag));
                    for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                        GiD_fWriteScalar( ResultFile, it->Id(), double_flag);
                    }
                }
            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    const double double_flag = static_cast<double>(it->Is(rFlag));
                    for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                        GiD_fWriteScalar( ResultFile, it->Id(), double_flag );
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<bool> rVariable, ModelPart& rModelPart,
                               double SolutionTag, unsigned int ValueIndex )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 ) {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );

            std::vector<bool> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), static_cast<double>(ValuesOnIntPoint[index]) );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), static_cast<double>(ValuesOnIntPoint[index]) );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<int> rVariable, ModelPart& rModelPart,
                               double SolutionTag, unsigned int ValueIndex )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 ) {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<int> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), static_cast<double>(ValuesOnIntPoint[index]) );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), static_cast<double>(ValuesOnIntPoint[index]) );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<double> rVariable, ModelPart& rModelPart,
                               double SolutionTag, unsigned int ValueIndex )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<double> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), ValuesOnIntPoint[index] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWriteScalar( ResultFile, it->Id(), ValuesOnIntPoint[index] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<array_1d<double,3> > rVariable, ModelPart& rModelPart,
                               double SolutionTag, unsigned int ValueIndex )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 ) {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult( ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Vector, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<array_1d<double,3> > ValuesOnIntPoint(mSize,ZeroVector(3));
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            if( ValuesOnIntPoint[0].size() == 3 )
                                GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
                                                 ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
                                             ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<array_1d<double,6> > rVariable, ModelPart& rModelPart,
                               double SolutionTag, unsigned int ValueIndex )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), ( char*)("Kratos"),
                             SolutionTag, GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<array_1d<double, 6> > ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); ++it ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) )
                    {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index][0],
                                               ValuesOnIntPoint[index][1], ValuesOnIntPoint[index][2],
                                               ValuesOnIntPoint[index][3], ValuesOnIntPoint[index][4],
                                               ValuesOnIntPoint[index][5] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
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


    virtual void PrintResults( GiD_FILE ResultFile, Variable<Vector> rVariable, ModelPart& rModelPart,
                               double SolutionTag, unsigned int ValueIndex )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );

            std::vector<Vector> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); ++it ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            const auto& values = ValuesOnIntPoint[index];
                            if (values.size() ==3 )
                                GiD_fWrite2DMatrix(ResultFile, it->Id(), values[0], values[1], values[2]);
                            else if (values.size() == 6 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), values[0], values[1], values[2],
                                    values[3], values[4], values[5] );
                        }
                    }
                }
            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            const auto& values = ValuesOnIntPoint[index];
                            if (values.size() ==3 )
                                GiD_fWrite2DMatrix(ResultFile, it->Id(), values[0], values[1], values[2]);
                            else if (values.size() == 6 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), values[0], values[1], values[2],
                                    values[3], values[4], values[5] );
                        }
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, Variable<Matrix> rVariable, ModelPart& rModelPart,
                               double SolutionTag, int ValueIndex )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 ) {
            //WriteGaussPoints(ResultFile);
            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"),
                             SolutionTag, GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<Matrix> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 ) {
                for( auto it = mMeshElements.begin(); it != mMeshElements.end(); ++it ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );
                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            if(ValuesOnIntPoint[index].size1() ==3
                                    && ValuesOnIntPoint[index].size2() ==3) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                                   ValuesOnIntPoint[index](0,2) );
                            }
                            else if(ValuesOnIntPoint[index].size1() ==2
                                    && ValuesOnIntPoint[index].size2() ==2) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), 0.0,
                                                   ValuesOnIntPoint[index](0,1), 0.0, 0.0);
                            }
                            else if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==3) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), 0.0,
                                                   ValuesOnIntPoint[index](0,2), 0.0, 0.0);
                            }
                            else if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==4) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), 0.0, 0.0);
                            }
                            else if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==6) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                                   ValuesOnIntPoint[index](0,5) );
                            }
                        }
                    }
                }

                // Resize first matrix to (0,0) for test below
                ValuesOnIntPoint[0].resize(0, 0, false);

            }
            if( mMeshConditions.size() != 0 ) {
                for( auto it = mMeshConditions.begin(); it != mMeshConditions.end(); it++ ) {
                    if( !(it->IsDefined(ACTIVE)) || it->Is(ACTIVE) ) {
                        it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                         rModelPart.GetProcessInfo() );

                        if (ValuesOnIntPoint[0].size1() == 0 && ValuesOnIntPoint[0].size2() == 0) {
                            // If we aren't getting any results, break
                            break;
                        }

                        for(unsigned int i=0; i<mIndexContainer.size(); i++) {
                            int index = mIndexContainer[i];
                            if(ValuesOnIntPoint[index].size1() ==3
                                    && ValuesOnIntPoint[index].size2() ==3) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](1,1), ValuesOnIntPoint[index](2,2),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](1,2),
                                                   ValuesOnIntPoint[index](0,2) );
                            }
                            else if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==6) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), ValuesOnIntPoint[index](0,4),
                                                   ValuesOnIntPoint[index](0,5) );
                            }
                            else if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==3) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), 0.0,
                                                   ValuesOnIntPoint[index](0,2), 0.0, 0.0);
                            }
                            else if(ValuesOnIntPoint[index].size1() ==1
                                    && ValuesOnIntPoint[index].size2() ==4) {
                                GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[index](0,0),
                                                   ValuesOnIntPoint[index](0,1), ValuesOnIntPoint[index](0,2),
                                                   ValuesOnIntPoint[index](0,3), 0.0, 0.0);
                            }
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

public:
    void WriteGaussPoints(GiD_FILE MeshFile)
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
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
            else if( mGidElementFamily == GiD_Quadrilateral && mSize == 4 )
            {
                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Quadrilateral, NULL, 4, 0, 0 );
                GiD_fWriteGaussPoint2D( MeshFile, - 1.00/std::sqrt(3.0), - 1.00/std::sqrt(3.0) );
                GiD_fWriteGaussPoint2D( MeshFile,   1.00/std::sqrt(3.0), - 1.00/std::sqrt(3.0) );
                GiD_fWriteGaussPoint2D( MeshFile,   1.00/std::sqrt(3.0),   1.00/std::sqrt(3.0) );
                GiD_fWriteGaussPoint2D( MeshFile, - 1.00/std::sqrt(3.0),   1.00/std::sqrt(3.0) );
                GiD_fEndGaussPoint(MeshFile);
            }
            else if( mGidElementFamily == GiD_Quadrilateral && mSize == 9 )
            {
                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Quadrilateral, NULL, 9, 0, 0 );
                GiD_fWriteGaussPoint2D( MeshFile, -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00));
                GiD_fWriteGaussPoint2D( MeshFile,  0.00 , -std::sqrt(3.00/5.00) );
                GiD_fWriteGaussPoint2D( MeshFile,  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00) );
                GiD_fWriteGaussPoint2D( MeshFile, -std::sqrt(3.00/5.00), 0.00 );
                GiD_fWriteGaussPoint2D( MeshFile,   0.00 , 0.00 );
                GiD_fWriteGaussPoint2D( MeshFile,  std::sqrt(3.00/5.00), 0.00);
                GiD_fWriteGaussPoint2D( MeshFile, -std::sqrt(3.00/5.00), std::sqrt(3.00/5.00) );
                GiD_fWriteGaussPoint2D( MeshFile,  0.00, std::sqrt(3.00/5.00) );
                GiD_fWriteGaussPoint2D( MeshFile,  std::sqrt(3.00/5.00), std::sqrt(3.00/5.00) );
                GiD_fEndGaussPoint(MeshFile);
            }
            else if( mGidElementFamily == GiD_Tetrahedra && mSize == 5 )
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
            /* START: Adding manually the custom prism */
            else if ( mGidElementFamily == GiD_Prism  && mSize > 1)
            {
                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Prism, NULL, 6, 0, 0 );

                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 - std::sqrt(1.00 / 3.00)));
                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5 * (1.0 - std::sqrt(1.00 / 3.00)));
                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 - std::sqrt(1.00 / 3.00)));

                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 + std::sqrt(1.00 / 3.00)));
                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5 * (1.0 + std::sqrt(1.00 / 3.00)));
                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 + std::sqrt(1.00 / 3.00)));

                GiD_fEndGaussPoint(MeshFile);
            }
//            else if ( mGidElementFamily == GiD_Prism  && mSize == 2 )
//            {
//                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Prism, NULL, 6, 0, 0 );

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 - std::sqrt(1.00 / 3.00)));
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5 * (1.0 - std::sqrt(1.00 / 3.00)));
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 - std::sqrt(1.00 / 3.00)));

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 + std::sqrt(1.00 / 3.00)));
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5 * (1.0 + std::sqrt(1.00 / 3.00)));
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5 * (1.0 + std::sqrt(1.00 / 3.00)));

//                GiD_fEndGaussPoint(MeshFile);
//            }
//            else if ( mGidElementFamily == GiD_Prism  && mSize == 3 )
//            {
//                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Prism, NULL, 9, 0, 0 );

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.1127016653792583114821);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.1127016653792583114821);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.1127016653792583114821);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.8872983346207416885180);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.8872983346207416885180);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.8872983346207416885180);

//                GiD_fEndGaussPoint(MeshFile);
//            }
//            else if ( mGidElementFamily == GiD_Prism  && mSize == 5 )
//            {
//                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Prism, NULL, 15, 0, 0 );

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.0469100770306680036012);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.0469100770306680036012);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.0469100770306680036012);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.2307653449471584544819);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.2307653449471584544819);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.2307653449471584544819);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.7692346550528415455182);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.7692346550528415455182);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.7692346550528415455182);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.9530899229693319963988);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.9530899229693319963988);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.9530899229693319963988);

//                GiD_fEndGaussPoint(MeshFile);
//            }
//            else if ( mGidElementFamily == GiD_Prism  && mSize == 7 )
//            {
//                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Prism, NULL, 21, 0, 0 );

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.0254460438286207377369);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.0254460438286207377369);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.0254460438286207377369);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.1292344072003027800681);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.1292344072003027800681);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.1292344072003027800681);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.2970774243113014165467);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.2970774243113014165467);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.2970774243113014165467);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.7029225756886985834533);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.7029225756886985834533);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.7029225756886985834533);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.8707655927996972199320);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.8707655927996972199320);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.8707655927996972199320);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.9745539561713792622631);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.9745539561713792622631);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.9745539561713792622631);

//                GiD_fEndGaussPoint(MeshFile);
//            }
//            else if ( mGidElementFamily == GiD_Prism  && mSize == 11 )
//            {
//                GiD_fBeginGaussPoint( MeshFile, mGPTitle, GiD_Prism, NULL, 33, 0, 0 );

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.0108856709269715035981);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.0108856709269715035981);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.0108856709269715035981);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.0564687001159523504624);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.0564687001159523504624);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.0564687001159523504624);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.1349239972129753379533);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.1349239972129753379533);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.1349239972129753379533);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.2404519353965940920372);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.2404519353965940920372);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.2404519353965940920372);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.3652284220238275138343);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.3652284220238275138343);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.3652284220238275138343);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.5000000000000000000000);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.6347715779761724861657);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.6347715779761724861657);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.6347715779761724861657);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.7595480646034059079628);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.7595480646034059079628);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.7595480646034059079628);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.8650760027870246620467);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.8650760027870246620467);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.8650760027870246620467);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.9435312998840476495376);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.9435312998840476495376);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.9435312998840476495376);

//                GiD_fWriteGaussPoint3D( MeshFile, 1.00 / 3.00 , 2.00 / 3.00, 0.9891143290730284964020);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 1.00 / 3.00, 0.9891143290730284964020);
//                GiD_fWriteGaussPoint3D( MeshFile, 2.00 / 3.00 , 2.00 / 3.00, 0.9891143290730284964020);

//                GiD_fEndGaussPoint(MeshFile);
//            }
            /* END: Adding manually the custom prism */
            else if ( mGidElementFamily == GiD_Point ||  mGidElementFamily == GiD_Sphere ||  mGidElementFamily == GiD_Circle )
            {
                //Gid does not accept gauss points on Points, Circles or Spheres! (october 18th 2014)
            }
            else
            {
                GiD_fBeginGaussPoint(MeshFile, mGPTitle, mGidElementFamily, NULL, mSize, 0, 1);
                GiD_fEndGaussPoint(MeshFile);
                /* NOTE: for linear elements, GiD gauss point coordinates are equispaced, visualization coordinates
                 * are wrong. Is there a GiD_fWriteGaussPoint2D equivalent for lines? Not according to documentation!
                 * gidpost documentation does not list one,
                 * GiD 12 explicitly states in the post format documentation that line coordinates cannot be given.
                 */
            }
        }
    }


protected:


    ///member variables
    const char * mGPTitle;
    KratosGeometryFamily mKratosElementFamily;
    GiD_ElementType mGidElementFamily;
    unsigned int mSize;
    std::vector<int> mIndexContainer;
    ModelPart::ElementsContainerType mMeshElements;
    ModelPart::ConditionsContainerType mMeshConditions;
};//class GidGaussPointsContainer
}// namespace Kratos.

#endif // KRATOS_GID_GAUSS_POINT_CONTAINER_H_INCLUDED defined

