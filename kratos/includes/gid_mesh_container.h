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

#pragma once
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
#include "includes/deprecated_variables.h"


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
    GidMeshContainer ( GeometryData::KratosGeometryType geometryType,
                       GiD_ElementType elementType, const char* mesh_title )
        :mGeometryType (geometryType), mGidElementType (elementType), mMeshTitle (mesh_title) {}
    bool AddElement ( const ModelPart::ElementConstantIterator pElemIt )
    {
        KRATOS_TRY
        if ( pElemIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
            mMeshElements.push_back ( * (pElemIt.base() ) );
            Geometry<Node >&geom = pElemIt->GetGeometry();
            for ( Element::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
            {
                mMeshNodes.push_back ( * (it.base() ) );
            }
            return true;
        }
        else
            return false;
        KRATOS_CATCH ("")
    }
    bool AddCondition (const ModelPart::ConditionConstantIterator pCondIt)
    {
        KRATOS_TRY
        if ( pCondIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
            mMeshConditions.push_back ( * (pCondIt.base() ) );
            Geometry<Node >&geom = pCondIt->GetGeometry();
            for ( Condition::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
            {
                mMeshNodes.push_back ( * (it.base() ) );
            }
            return true;
        }
        else
            return false;
        KRATOS_CATCH ("")
    }
    bool AddGeometry (const ModelPart::GeometryConstantIterator pGeomIt)
    {
        KRATOS_TRY
        if ( pGeomIt->GetGeometryType() == mGeometryType )
        {
            mMeshGeometries.push_back ( * (pGeomIt.base() ) );
            Geometry<Node> const&geom = *pGeomIt;
            for ( Geometry<Node>::const_iterator it = geom.begin(); it != geom.end(); it++)
            {
                mMeshNodes.push_back ( * (it.base() ) );
            }
            return true;
        }
        else
            return false;
        KRATOS_CATCH ("")
    }
    void FinalizeMeshCreation()
    {
        if ( mMeshElements.size() != 0 )
        {
            mMeshNodes.Unique();
        }
        if ( mMeshConditions.size() != 0 )
        {
            mMeshNodes.Unique();
        }
        if ( mMeshGeometries.size() != 0 )
        {
            mMeshNodes.Unique();
        }
    }
    void WriteMesh (GiD_FILE MeshFile, bool deformed)
    {
        KRATOS_TRY

        bool nodes_written = false;
        if ( mMeshElements.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
            {
                const int prop_id = it->HasProperties() ? it->GetProperties().Id() : 0;
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> elements_per_layer (max_id+1,0);
            //KRATOS_WATCH(max_id);

            //fill layer list
            for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
            {
                const int prop_id = it->HasProperties() ? (it)->GetProperties().Id() : 0;
                elements_per_layer[prop_id] += 1;
            }
            //std::cout << "start printing elements" <<std::endl;
            for (unsigned int current_layer = 0; current_layer < elements_per_layer.size(); current_layer++)
            {
                if (elements_per_layer[current_layer] > 0)
                {
                    //create an appropriate name
                    std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
                    current_layer_name << mMeshTitle << "_" << current_layer ;
                    if ( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        //std::cout << " -print element 2D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, mGidElementType,mMeshElements.begin()->GetGeometry().size() );
                    }
                    else if ( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        //std::cout << " -print element 3D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, mGidElementType,mMeshElements.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_THROW_ERROR (std::logic_error,"check working space dimension of model","");
                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
                        {
                            if ( deformed )
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(),
                                                       (it)->Y(), (it)->Z() );
                            else
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(),
                                                       (it)->Y0(), (it)->Z0() );
                        }
                        GiD_fEndCoordinates(MeshFile);

                        nodes_written = true;
                    }
                    //printing elements
                    GiD_fBeginElements(MeshFile);
                    int* nodes_id = new int[mMeshElements.begin()->GetGeometry().size() + 1];
                    for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                            it != mMeshElements.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry() [i].Id();

                        if (mGeometryType == GeometryData::KratosGeometryType::Kratos_Line2D3 || mGeometryType == GeometryData::KratosGeometryType::Kratos_Line3D3)
                        {
                            nodes_id[0] = (it)->GetGeometry() [0].Id();
                            nodes_id[1] = (it)->GetGeometry() [1].Id();
                            nodes_id[2] = (it)->GetGeometry() [2].Id();
                        }
                        const unsigned int elem_layer = it->HasProperties() ? it->GetProperties().Id() : 0;
                        nodes_id[ (it)->GetGeometry().size()]= elem_layer + 1;

                        if (it->IsActive())
                            if (elem_layer == current_layer)
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id);

                    }
                    delete [] nodes_id;
                    GiD_fEndElements(MeshFile);
                    GiD_fEndMesh(MeshFile);
                }
            }
            //std::cout << "end printing elements" <<std::endl;
        }
        if ( mMeshConditions.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
            {
                const int prop_id = it->HasProperties() ? (it)->GetProperties().Id() : 0;
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> conditions_per_layer (max_id+1,0);
            //fill layer list
            for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
            {
                const int prop_id = it->HasProperties() ? it->GetProperties().Id() : 0;
                conditions_per_layer[prop_id] += 1;
            }
            //std::cout << "start printing conditions" <<std::endl;
            for (unsigned int current_layer = 0; current_layer < conditions_per_layer.size(); current_layer++)
            {
                if (conditions_per_layer[current_layer] > 0)
                {
                    std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
                    current_layer_name << mMeshTitle << "_" << current_layer ;

                    if ( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        //std::cout << " -print condition 2D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, mGidElementType,
                                        mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else if ( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        //std::cout << " -print condition 3D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, mGidElementType,
                                        mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_THROW_ERROR (std::logic_error,"check working space dimension of model","");
                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
                        {
                            if ( deformed )
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(),
                                                       (it)->Y(), (it)->Z() );
                            else
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(),
                                                       (it)->Y0(), (it)->Z0() );
                        }
                        GiD_fEndCoordinates(MeshFile);
                        nodes_written = true;
                    }
                    else  //printing these headers avoids the exception raised by GidPost. There's an Assert stopping the program if we don't put this here.
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        GiD_fEndCoordinates(MeshFile);
                    }
                    //printing elements
                    GiD_fBeginElements(MeshFile);
                    int* nodes_id = new int[mMeshConditions.begin()->GetGeometry().size() + 1];
                    for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(  );
                            it != mMeshConditions.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry() [i].Id();
                        //workaround: reordering node ids for Hexahedra20 elements
                        if (mGeometryType == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
                        {
                            nodes_id[12] = (it)->GetGeometry() [16].Id();
                            nodes_id[13] = (it)->GetGeometry() [17].Id();
                            nodes_id[14] = (it)->GetGeometry() [18].Id();
                            nodes_id[15] = (it)->GetGeometry() [19].Id();
                            nodes_id[16] = (it)->GetGeometry() [12].Id();
                            nodes_id[17] = (it)->GetGeometry() [13].Id();
                            nodes_id[18] = (it)->GetGeometry() [14].Id();
                            nodes_id[19] = (it)->GetGeometry() [15].Id();
                        }
                        const unsigned int cond_layer = it->HasProperties() ? it->GetProperties().Id() : 0;
                        nodes_id[(it)->GetGeometry().size()]= cond_layer + 1;

                        if (it->IsActive())
                            if (cond_layer == current_layer)
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id);
                    }
                    delete [] nodes_id;
                    GiD_fEndElements(MeshFile);
                    GiD_fEndMesh(MeshFile);
                }
            }
            //std::cout << "end printing conditions" <<std::endl;
        }
        if ( mMeshGeometries.size() != 0 )
        {
            std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
            current_layer_name << mMeshTitle << "_0" ;
            if ( mMeshGeometries.front()->WorkingSpaceDimension() == 2 )
            {
                GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, mGidElementType, mMeshGeometries.front()->size() );
            }
            else if ( mMeshGeometries.front()->WorkingSpaceDimension() == 3 )
            {
                GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, mGidElementType, mMeshGeometries.front()->size() );
            }
            else
                KRATOS_THROW_ERROR (std::logic_error,"check working space dimension of model","");
            //printing nodes
            if(nodes_written == false)
            {
                GiD_fBeginCoordinates(MeshFile);
                for ( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin(); it != mMeshNodes.end(); ++it )
                {
                    if ( deformed )
                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(), (it)->Y(), (it)->Z() );
                    else
                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(), (it)->Y0(), (it)->Z0() );
                }
                GiD_fEndCoordinates(MeshFile);
                nodes_written = true;
            }
            else
            {
                GiD_fBeginCoordinates(MeshFile);
                GiD_fEndCoordinates(MeshFile);
            }
            //printing geometry connectivity
            GiD_fBeginElements(MeshFile);
            int* nodes_id = new int[mMeshGeometries.front()->size() + 1];
            for ( std::vector<ModelPart::GeometryType::Pointer>::iterator it = mMeshGeometries.begin(); it != mMeshGeometries.end(); ++it )
            {
                ModelPart::GeometryType const& geom = **it;
                for ( unsigned int i = 0; i < geom.size(); i++ )
                    nodes_id[i] = geom[i].Id();
                nodes_id[ geom.size() ] = 1;
                GiD_fWriteElementMat ( MeshFile, geom.Id(), nodes_id );
            }
            delete [] nodes_id;
            GiD_fEndElements(MeshFile);
            GiD_fEndMesh(MeshFile);
        }
        KRATOS_CATCH ("")
    }
    void Reset()
    {
        mMeshNodes.clear();
        mMeshElements.clear();
        mMeshConditions.clear();
        mMeshGeometries.clear();
    }
    ModelPart::NodesContainerType GetMeshNodes()
    {
        return mMeshNodes;
    }
protected:
    ///member variables
    GeometryData::KratosGeometryType mGeometryType;
    GiD_ElementType mGidElementType;
    ModelPart::NodesContainerType mMeshNodes;
    ModelPart::ElementsContainerType mMeshElements;
    ModelPart::ConditionsContainerType mMeshConditions;
    std::vector<ModelPart::GeometryType::Pointer> mMeshGeometries;
    const char* mMeshTitle;
};//class GidMeshContainer
}// namespace Kratos.
