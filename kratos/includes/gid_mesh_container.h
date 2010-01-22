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

#if !defined(KRATOS_GID_MESH_CONTAINER_H_INCLUDED)
#define  KRATOS_GID_MESH_CONTAINER_H_INCLUDED

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
     * Auxiliary class to store meshes of different element types and to
     * write these meshes to an output file
     */
    class GidMeshContainer
    {
        public:
            ///Constructor
            GidMeshContainer( GeometryData::KratosGeometryType geometryType, 
                              GiD_ElementType elementType, const char* mesh_title )
            :mGeometryType(geometryType), mGidElementType(elementType), mMeshTitle(mesh_title){}
            
            bool AddElement( const ModelPart::ElementsContainerType::iterator pElemIt )
            {
                KRATOS_TRY
                if( pElemIt->GetGeometry().GetGeometryType() == mGeometryType )
                {
                    mMeshElements.push_back( *(pElemIt.base() ) );
                    Geometry<Node<3> >&geom = pElemIt->GetGeometry();
                    for( Element::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
                    {
                        mMeshNodes.push_back( *(it.base() ) );
                    }
                    return true;
                }
                else
                    return false;
                KRATOS_CATCH("")
            }
            
            bool AddCondition(const ModelPart::ConditionsContainerType::iterator pCondIt)
            {
                KRATOS_TRY
                        if( pCondIt->GetGeometry().GetGeometryType() == mGeometryType )
                {
                    mMeshConditions.push_back( *(pCondIt.base() ) );
                    Geometry<Node<3> >&geom = pCondIt->GetGeometry();
//                     for( int i=0; i<geom.size(); i++ )
//                     {
//                         mMeshNodes.push_back( geom.pGetPoint(i) );
//                     }
                    for( Condition::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
                    {
                        mMeshNodes.push_back( *(it.base() ) );
                    }
                    return true;
                }
                else
                    return false;
                KRATOS_CATCH("")
            }
            
            void FinalizeMeshCreation()
            {
                if( mMeshElements.size() != 0 )
                {
                    mMeshNodes.Unique();
                }
                if( mMeshConditions.size() != 0 )
                {
                    mMeshNodes.Unique();
                }
            }
            
            void WriteMesh(bool deformed)
            {
                KRATOS_TRY
                if( mMeshElements.size() != 0 )
                {
                    if( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        std::cout << "writing a 2D mesh" << std::endl;
                        GiD_BeginMesh( mMeshTitle, GiD_2D, mGidElementType,
                                       mMeshElements.begin()->GetGeometry().size() );
                    }
                    else if( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        std::cout << "writing a 3D mesh" << std::endl;
                        GiD_BeginMesh( mMeshTitle, GiD_3D, mGidElementType,
                                       mMeshElements.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_ERROR(std::logic_error,"check working space dimension of model","");
                    //printing nodes
                    GiD_BeginCoordinates();
                    for( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin(); 
                         it != mMeshNodes.end(); ++it )
                    {
                        if( deformed )
                            GiD_WriteCoordinates( (it)->Id(), (it)->X(),
                                (it)->Y(), (it)->Z());
                        else
                            GiD_WriteCoordinates( (it)->Id(), (it)->X0(),
                                (it)->Y0(), (it)->Z0());
                    }
                    GiD_EndCoordinates();
                    //printing elements
                    GiD_BeginElements();
                    int* nodes_id = new int[mMeshElements.begin()->GetGeometry().size()+1];
                    for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin(); 
                         it != mMeshElements.end(); ++it )
                    {
                        for( unsigned int i=0; i<(it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[i].Id();
                        
                        //workaround: reordering node ids for Hexahedra20 elements
                        if( mGeometryType == GeometryData::Kratos_Hexahedra3D20 )
                        {
                            nodes_id[12] = (it)->GetGeometry()[16].Id();
                            nodes_id[13] = (it)->GetGeometry()[17].Id();
                            nodes_id[14] = (it)->GetGeometry()[18].Id();
                            nodes_id[15] = (it)->GetGeometry()[19].Id();
                            nodes_id[16] = (it)->GetGeometry()[12].Id();
                            nodes_id[17] = (it)->GetGeometry()[13].Id();
                            nodes_id[18] = (it)->GetGeometry()[14].Id();
                            nodes_id[19] = (it)->GetGeometry()[15].Id();
                        }
                        if( mGeometryType == GeometryData::Kratos_Line2D3 
                            || mGeometryType == GeometryData::Kratos_Line3D3 )
                        {
                            nodes_id[0] = (it)->GetGeometry()[0].Id();
                            nodes_id[1] = (it)->GetGeometry()[2].Id();
                            nodes_id[2] = (it)->GetGeometry()[1].Id();
                        }
                        nodes_id[(it)->GetGeometry().size()]= (it)->GetProperties().Id();
                        if( it->Has( IS_INACTIVE ) )
                        {
                            if( ! it->GetValue( IS_INACTIVE ) )
                            {
                                KRATOS_WATCH( it->GetValue(IS_INACTIVE) );
                                GiD_WriteElementMat((it)->Id(), nodes_id);
                            }
                        }
                        else
                        {
                            KRATOS_WATCH( "IS_INACTIVE not defined" );
                            GiD_WriteElementMat((it)->Id(), nodes_id);
                        }
                    }
                    delete [] nodes_id;
                    GiD_EndElements();
                    GiD_EndMesh();
                }
                
                if( mMeshConditions.size() != 0 )
                {
                    if( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        std::cout << "writing a 2D mesh" << std::endl;
                        GiD_BeginMesh( mMeshTitle, GiD_2D, mGidElementType,
                                       mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else if( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        std::cout << "writing a 3D mesh" << std::endl;
                        GiD_BeginMesh( mMeshTitle, GiD_3D, mGidElementType,
                                       mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_ERROR(std::logic_error,"check working space dimension of model","");
                    //printing nodes
                    GiD_BeginCoordinates();
                    for( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin(); 
                         it != mMeshNodes.end(); ++it )
                    {
                        if( deformed )
                            GiD_WriteCoordinates( (it)->Id(), (it)->X(),
                                (it)->Y(), (it)->Z());
                        else
                            GiD_WriteCoordinates( (it)->Id(), (it)->X0(),
                                (it)->Y0(), (it)->Z0());
                    }
                    GiD_EndCoordinates();
                    //printing elements
                    GiD_BeginElements();
                    int* nodes_id = new int[mMeshConditions.begin()->GetGeometry().size()+1];
                    for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(	); 
                         it != mMeshConditions.end(); ++it )
                    {
                        
                        for( unsigned int i=0; i<(it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[i].Id();
                        //workaround: reordering node ids for Hexahedra20 elements
                        if( mGeometryType == GeometryData::Kratos_Hexahedra3D20 )
                        {
                            nodes_id[12] = (it)->GetGeometry()[16].Id();
                            nodes_id[13] = (it)->GetGeometry()[17].Id();
                            nodes_id[14] = (it)->GetGeometry()[18].Id();
                            nodes_id[15] = (it)->GetGeometry()[19].Id();
                            nodes_id[16] = (it)->GetGeometry()[12].Id();
                            nodes_id[17] = (it)->GetGeometry()[13].Id();
                            nodes_id[18] = (it)->GetGeometry()[14].Id();
                            nodes_id[19] = (it)->GetGeometry()[15].Id();
                        }
                        nodes_id[(it)->GetGeometry().size()]= (it)->GetProperties().Id();
                        if( it->Has( IS_INACTIVE ) )
                        {
                            if( ! it->GetValue( IS_INACTIVE ) )
                            {
                                KRATOS_WATCH( it->GetValue(IS_INACTIVE) );
                                GiD_WriteElementMat((it)->Id(), nodes_id);
                            }
                        }
                        else
                        {
                            KRATOS_WATCH( "IS_INACTIVE not defined" );
                            GiD_WriteElementMat((it)->Id(), nodes_id);
                        }

                    }
                    delete [] nodes_id;
                    GiD_EndElements();
                    GiD_EndMesh();
                }
                KRATOS_CATCH("")
            }
            
            void Reset()
            {
                mMeshNodes.clear();
                mMeshElements.clear();
                mMeshConditions.clear();
            }
            
        protected:
            ///member variables
            GeometryData::KratosGeometryType mGeometryType;
            GiD_ElementType mGidElementType;
            ModelPart::NodesContainerType mMeshNodes;
            ModelPart::ElementsContainerType mMeshElements;
            ModelPart::ConditionsContainerType mMeshConditions;
            const char* mMeshTitle;
    };//class GidMeshContainer
}// namespace Kratos.

#endif // KRATOS_GID_MESH_CONTAINER_H_INCLUDED defined 
