//
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include "utilities/openmp_utils.h"
// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/entity_erase_process.h"

#include "utilities/binbased_fast_point_locator.h"

#include "u_qualityMetrics.h"
#include "Math3D.h"
#include "u_Types.h"
#include "u_TetraFunctions.h"
#include "u_ShowMetrics.h"
#include "u_ParallelFunctions.h"
#include "u_MeshLoaders.h"
#include "u_elementCluster.h"
#include "u_ProcessTime.h"
#include "u_TetGenInterface.h"

#ifdef _OPENMP
#include <omp.h>
#endif


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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
class TetrahedraReconnectUtility
{
public:
    ///@name Type Definitions
    ///@{
    int maxNumThreads ;
    int blockSize ;
    bool debugMode;
    /// Pointer definition of TetrahedraReconnectUtility
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedraReconnectUtility);

    ///@}
    ///@name Life Cycle
    ///@{
    // inner Mesh
    TVolumeMesh *m ;

    ModelPart &refMP ;

    /// Default constructor.
    TetrahedraReconnectUtility(ModelPart& r_model_part):
        refMP(r_model_part)
    {
        std::cout << "Creating mesh" << "\n";
        m = new TVolumeMesh();
        // Convert to inner format
        innerConvertFromKratos(r_model_part , m );
        //refMP = r_model_part;

        maxNumThreads = 0;
        blockSize = 2048;

    }

    /// Destructor.
    virtual ~TetrahedraReconnectUtility() = default;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@{

    ///@brief function innerConvertFromKratos
    /// This function converts from Kratos format, to the inner structure
    ///@param r_model_part the input mesh
    ///@param m the output mesh
    void innerConvertFromKratos(ModelPart& r_model_part, TVolumeMesh *m)
    {
        if (debugMode) std::cout << "Reading nodes"<< "\n";
        //reorder node Ids consecutively
        unsigned int id=1;
        for (ModelPart::NodesContainerType::iterator i_node = r_model_part.NodesBegin() ; i_node != r_model_part.NodesEnd() ; i_node++)
            i_node->SetId(id++);

        //loop on nodes
        for (ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            float4 fPos ;
            fPos.x = it->X();
            fPos.y = it->Y();
            fPos.z = it->Z();
            auto v = new TVertex(fPos);
            v->setID( it->Id() );
            v->blockID = true;
            m->vertexes->Add(v);
        }
        if (debugMode) std::cout << "Reading elements"<< "\n";

        for (ModelPart::ElementsContainerType::iterator el_it=r_model_part.ElementsBegin(); el_it!=r_model_part.ElementsEnd(); el_it++)
        {
            Geometry< Node >& geom = el_it->GetGeometry();
            if (geom.size() != 4)
            {
                std::cout << "Invalid element size" <<  el_it->Id();
                continue;
            }
            TVertex *v0 = m->findVertexById(geom[0].Id());
            if (  v0 == nullptr )
            {
                std::cout << "Invalid element reference" <<  el_it->Id();
                continue;
            }
            TVertex *v1 = m->findVertexById(geom[1].Id());
            if (  v1 == nullptr )
            {
                std::cout << "Invalid element reference" <<  el_it->Id();
                continue;
            }
            TVertex *v2 = m->findVertexById(geom[2].Id());
            if (  v2 == nullptr )
            {
                std::cout << "Invalid element reference" <<  el_it->Id();
                continue;
            }
            TVertex *v3 = m->findVertexById(geom[3].Id());
            if (  v3 == nullptr )
            {
                std::cout << "Invalid element reference" <<  el_it->Id();
                continue;
            }

            auto t = new TTetra(nullptr, v0,v1,v2,v3);
            m->elements->Add(t);
        }

        if (debugMode) std::cout << " Number of vertexes read :"<< m->vertexes->Count() <<"\n";
        if (debugMode) std::cout << " Number of elements read :"<< m->elements->Count() << "\n";
        m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
        if (debugMode) std::cout << " Number of faces read :"<< m->fFaces->Count() << "\n";

    }

    ///@brief function innerConvertToKratos
    /// This function converts back to Kratos format
    ///@param m the input mesh
    ///@param mrModelPart the output mesh
    ///@param removeFreeVertexes Free vertexes are not added to the new structure
    void innerConvertToKratos(ModelPart& mrModelPart , TVolumeMesh *m, bool removeFreeVertexes)
    {
        if (debugMode) std::cout << "-------------Generating for Kratos----------------" << "\n";
        m->vertexes->Sort(sortByID);
        Element::Pointer pReferenceElement = *(mrModelPart.Elements().begin()).base();

        // Mark elements to delete
        // 0 Not remove
        // 1 Remove
        for (int i=0; i< m->vertexes->Count(); i++)
        {
            if (m->vertexes->structure[i]->elementsList->Count() == 0)
                m->vertexes->elementAt(i)->flag = 1;
            else
                m->vertexes->elementAt(i)->flag = 0;
        }
        if (debugMode) std::cout << " Nodes : Mark elements to remove : " << "\n";
        bool shownMessage = false;
        /* Control de Indices
        int indexv = 0;
        TMeshLoader* ml2 = new TVMWLoader();
        std::string s("d:/MeshKratos_innerConvertion.vwm");
        ml2->save( s, m);
        delete ml2;
        */
        // Create new Nodes
        for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.Nodes().begin() ; i_node != mrModelPart.Nodes().end() ; i_node++)
        {
            TVertex* v = m->findVertexById(i_node->Id());
            if (v == nullptr)
            {
                if (!shownMessage)
                    std::cout << "Error at id "<<i_node->Id() << "\n";
                shownMessage = true;
                continue;
            }
            i_node->Set(TO_ERASE ,v->flag == 1);
        }
        if (debugMode) std::cout << " Generating Elements for Kratos " << "\n";
        //generate new nodes
        mrModelPart.Elements().clear();
        //add preserved elements to the kratos
        Properties::Pointer properties = mrModelPart.GetMesh().pGetProperties(1);
        (mrModelPart.Elements()).reserve(m->elements->Count());
        bool messageShown = false;

        for (int i=0; i< m->elements->Count() ; i++)
        {

            TTetra* t = (TTetra*)(m->elements->elementAt(i));
            Node::Pointer v0 = mrModelPart.pGetNode(t->vertexes[0]->getID());
            Node::Pointer v1 = mrModelPart.pGetNode(t->vertexes[1]->getID());
            Node::Pointer v2 = mrModelPart.pGetNode(t->vertexes[2]->getID());
            Node::Pointer v3 = mrModelPart.pGetNode(t->vertexes[3]->getID());
            if ((v0.get() == nullptr) || (v1.get()==nullptr) || (v2.get() == nullptr) || (v3.get() == nullptr))
            
            {
                if (!messageShown)
                    std::cout << "Invalid vertex access " << t->getID() <<"\n";
                messageShown = true;
                continue ;
            }
            Tetrahedra3D4<Node > geom(
                v0,v1,v2,v3
            );

            Element::Pointer p_element = pReferenceElement->Create(i+1, geom, properties);
            (mrModelPart.Elements()).push_back(p_element);
            //KRATOS_WATCH(*p_element);
        }
        if (debugMode) std::cout << "Generation OK " << "\n";
        mrModelPart.Elements().Sort();

        if (removeFreeVertexes )
        {
            std::cout << "Removing free vertexes " << "\n";
            (EntitiesEraseProcess<Node>(mrModelPart)).Execute();
        }
        if (debugMode) std::cout << "-----------------Generation Finished OK-------------------" << "\n";

        //  TVolumeMesh *testM = new TVolumeMesh();
        //  innerConvertFromKratos( mrModelPart , testM);

    }

    bool EvaluateQuality()
    {
        auto qt = new TetQuality(m) ;
        qt->refresh();
        qt->print();
        int numNegElements = qt->nonPositive;
        delete qt;
        return numNegElements == 0;
    }

    /**
    * This function performs the meshing optimization by Cluster reconnection.
    */
    void TestRemovingElements()
    {
        for (int i=0 ; i<100 ; i++)
        {
            m->elementsToRemove->Add( m->elements->elementAt(i));
        }
        m->updateRefs();
    }

    ///@brief function updateNodesPositions
    /// Update only nodes positions without regenerating the structure
    ///@param r_model_part the input mesh
    void updateNodesPositions(ModelPart& r_model_part)
    {
        std::cout << "Updating nodes"<< "\n";
        //loop on nodes
        for (ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            float4 fPos ;
            fPos.x = it->X();
            fPos.y = it->Y();
            fPos.z = it->Z();
            fPos.w = 1.0f;

            int id = it->Id();
            TVertex *v = m->findVertexById(id);
            if (v == nullptr)
                v->fPos = fPos;
        }
    }
    ///@brief function OptimizeQuality
    ///@param r_model_part the input mesh
    ///@param iterations amount of iterations to optimize
    ///@param processByNode boolean that indicates if processing is done by Node
    ///@param processByFace boolean that indicates if processing is done by Face
    ///@param processByEdge boolean that indicates if processing is done by Edge
    ///@param saveToFile boolean to use as Debug Mode and see intermediate generated meshes
    ///@param removeFreeVertexes boolean Removes vertexes that do not have elements referencing them
    ///@param evaluateInParallel boolean Activate/Deactivate parallel processing mode
    ///@param reinsertNodes boolean Activate/Deactivate to try re inserting removed nodes

    void setMaxNumThreads(int mxTh)
    {
        maxNumThreads = mxTh;
    }

    void setBlockSize(int bs)
    {
        blockSize = bs;
    }

    bool isaValidMesh()
    {
        auto qt = new TetQuality(m) ;
        qt->refresh();
        int numNegElements = qt->nonPositive;
        delete qt;
        return numNegElements == 0;
    }

    void InterpolateAndAddNewNodes(ModelPart& rModelPart, TVolumeMesh* m, BinBasedFastPointLocator<3>& element_finder)
    {
        unsigned int n_points_before_refinement = rModelPart.Nodes().size();

        //if the refinement was performed, we need to add it to the model part.
        if (m->vertexes->Count()>(int)n_points_before_refinement)
        {
            //definitions for spatial search
//             typedef Node PointType;
//             typedef Node ::Pointer PointTypePointer;
            array_1d<double, 4 > N;
            const int max_results = 10000;
            BinBasedFastPointLocator<3>::ResultContainerType results(max_results);

            Node::DofsContainerType& reference_dofs = (rModelPart.NodesBegin())->GetDofs();

            int step_data_size = rModelPart.GetNodalSolutionStepDataSize();
            std::cout <<"...Adding Nodes :"<<  m->vertexes->Count() -  n_points_before_refinement<<"\n";
            //TODO: parallelize this loop
            for (int i = n_points_before_refinement; i<m->vertexes->Count(); i++)
            {
                int id=i+1;

                TVertex *v = m->vertexes->elementAt(i);
                double x= (float)(v->fPos.x);
                double y= (float)(v->fPos.y);
                double z= (float)(v->fPos.z);

                Node::Pointer pnode = rModelPart.CreateNewNode(id,x,y,z);

                //putting the new node also in an auxiliary list
                //KRATOS_WATCH("adding nodes to list")
                //list_of_new_nodes.push_back( pnode );

                //std::cout << "new node id = " << pnode->Id() << std::endl;
                //generating the dofs
                for (Node::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
                {
                    Node::DofType &rDof = **iii;
                    Node::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

                    (p_new_dof)->FreeDof();
                }

                //do interpolation
                auto result_begin = results.begin();
                Element::Pointer pelement;

                bool is_found = element_finder.FindPointOnMesh(pnode->Coordinates(), N, pelement, result_begin, max_results);

                if (is_found == true)
                {
                    Geometry<Node >& geom = pelement->GetGeometry();

                    Interpolate( geom, N, step_data_size, pnode);
                }


            }
        }
        //	std::cout << "During refinement we added " << tet.numberofpoints-n_points_before_refinement<< "nodes " <<std::endl;
    }

    void setVertexExpectedSize(ModelPart& rModelPart, TVolumeMesh* m )
    {
        for (int el = 0; el< m->elements->Count(); el++)
        {

            //calculate the prescribed h
            /*                double prescribed_h = (nodes_begin + out.tetrahedronlist[old_base]-1)->FastGetSolutionStepValue(NODAL_H);
            prescribed_h += (nodes_begin + out.tetrahedronlist[old_base+1]-1)->FastGetSolutionStepValue(NODAL_H);
            prescribed_h += (nodes_begin + out.tetrahedronlist[old_base+2]-1)->FastGetSolutionStepValue(NODAL_H);
            prescribed_h += (nodes_begin + out.tetrahedronlist[old_base+3]-1)->FastGetSolutionStepValue(NODAL_H);
            prescribed_h *= 0.25;*/
            TTetra *t = (TTetra*)(m->elements->elementAt(el));
            ModelPart::NodesContainerType& rNodes  = rModelPart.Nodes();
            ModelPart::NodesContainerType::iterator it1 = (rNodes).find( t->vertexes[0]->getID());
            ModelPart::NodesContainerType::iterator it2 = (rNodes).find( t->vertexes[1]->getID());
            ModelPart::NodesContainerType::iterator it3 = (rNodes).find( t->vertexes[2]->getID());
            ModelPart::NodesContainerType::iterator it4 = (rNodes).find( t->vertexes[3]->getID());

            if ( it1 == rModelPart.Nodes().end() )
                KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it1->Id());
            if ( it2 == rModelPart.Nodes().end() )
                KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it2->Id());
            if ( it3 == rModelPart.Nodes().end() )
                KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it3->Id());
            if ( it4 == rModelPart.Nodes().end() )
                KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it4->Id());

            Node::Pointer pn1 =  *it1.base();
            Node::Pointer pn2 =  *it2.base();
            Node::Pointer pn3 =  *it3.base();
            Node::Pointer pn4 =  *it4.base();

            t->vertexes[0]->expectedSize = (pn1)->FastGetSolutionStepValue(NODAL_H);
            t->vertexes[1]->expectedSize = (pn2)->FastGetSolutionStepValue(NODAL_H);
            t->vertexes[2]->expectedSize = (pn3)->FastGetSolutionStepValue(NODAL_H);
            t->vertexes[3]->expectedSize = (pn4)->FastGetSolutionStepValue(NODAL_H);

            //if h is the height of a perfect tetrahedra, the edge size is edge = sqrt(3/2) h
            //filling in the list of "IDEAL" tetrahedron volumes=1/12 * (edge)^3 * sqrt(2)~0.11785* h^3=
            //0.2165063509*h^3
//					double prescribed_h = (t->vertexes[0]->expectedSize + t->vertexes[1]->expectedSize +
//										  t->vertexes[2]->expectedSize  +t->vertexes[3]->expectedSize) * 0.25;
//
            //out.tetrahedronvolumelist[counter] = 0.217*prescribed_h*prescribed_h*prescribed_h;


        }
    }
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    void Interpolate( Geometry<Node >& geom, const array_1d<double,4>& N,
                      unsigned int step_data_size,
                      Node::Pointer pnode)
    {
        unsigned int buffer_size = pnode->GetBufferSize();


        for (unsigned int step = 0; step<buffer_size; step++)
        {

            //getting the data of the solution step
            double* step_data = (pnode)->SolutionStepData().Data(step);


            double* node0_data = geom[0].SolutionStepData().Data(step);
            double* node1_data = geom[1].SolutionStepData().Data(step);
            double* node2_data = geom[2].SolutionStepData().Data(step);
            double* node3_data = geom[3].SolutionStepData().Data(step);

            //copying this data in the position of the vector we are interested in
            for (unsigned int j= 0; j<step_data_size; j++)
            {

                step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j] + N[3]*node3_data[j];


            }
        }
        pnode->FastGetSolutionStepValue(IS_BOUNDARY)=0.0;
    }

    void OptimizeQuality(ModelPart& r_model_part,int simIter, int iterations ,
                         bool processByNode, bool processByFace, bool processByEdge,
                         bool saveToFile, bool removeFreeVertexes ,
                         bool evaluateInParallel , bool reinsertNodes , bool debugMode, int minAngle)
    {
        this->debugMode = debugMode;
        m->vertexes->Sort(sortByID);


        // Save the mesh as generated from Kratos
        if (saveToFile)
        {
            TMeshLoader* ml2 = new TVMWLoader();
            std::string s("out_MeshFromKratos.vwm");
            ml2->save( s, m);
            delete ml2;
        }
        if (debugMode && saveToFile)
        {
            EvaluateQuality();
            setVertexExpectedSize(r_model_part , m);
            std::cout <<"...Start Optimization..." <<"\n";
            auto  ml2 = new TVMWLoader();
            std::string s("");
            s = "../out_MeshFromKratos" + intToString(simIter)+".vwm";

            ml2->save( s.data(), m);
            delete ml2;

            return ;
        }

        if (debugMode) std::cout <<"...Start Optimization..." <<"\n";
        // maxNumThreads works as a FLAG
        //    when equals to 0, take "max available threads"
        //    in other case, take "set num threads"

        if ( maxNumThreads == 0)
        {
#ifdef _OPENMP
            OpenMPUtils::SetNumThreads(omp_get_max_threads());
#endif
        }
        else
            OpenMPUtils::SetNumThreads(maxNumThreads);

        if (debugMode)
        {
            std::cout <<"Number of active threads"<< OpenMPUtils::GetNumThreads() <<"\n";
            std::cout <<"Debug mode is Active" <<"\n";
            startTimers();
        }
        else
        {
            stopTimers();
        }
        std::cout <<"...Trying to allocate memory\n";
        preparePool();
        std::cout <<"...Allocate memory OK\n";

        for (int iter = 0 ; iter< iterations ; iter ++)
        {
            //ParallelEvaluateClusterByNode(m,vrelaxQuality);
            if (processByNode)
            {

                if (evaluateInParallel )
                {
                    if (debugMode) std::cout <<"...Parallel optimizing by Node. Iteration : "<< iter <<"\n";
                    ParallelEvaluateClusterByNode((TVolumeMesh*)(m),vrelaxQuality,minAngle);
                }
                else
                {
                    if (debugMode)  std::cout <<"...Optimizing by Node. Iteration : "<< iter <<"\n";
                    evaluateClusterByNode( (TVolumeMesh*)(m),minAngle,vrelaxQuality);
                }
                if (debugMode)
                    m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
            }

            if (processByFace)
            {
                if (evaluateInParallel )
                {
                    if (debugMode)  std::cout <<"...Parallel optimizing by Face. Iteration : "<< iter <<"\n";
                    ParallelEvaluateClusterByFace((TVolumeMesh*)(m),vrelaxQuality,minAngle);
                    if (debugMode)  std::cout <<"...End. Iteration : "<< iter <<"\n";
                }
                else
                {
                    if (debugMode)  std::cout <<"...Optimizing by Face. Iteration : "<< iter <<"\n";
                    evaluateClusterByFace(m,minAngle,vrelaxQuality);
                    if (debugMode)  std::cout <<"...End. Iteration : "<< iter <<"\n";
                }
                if (debugMode)
                    m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
            }

            if (processByEdge)
            {

                if (evaluateInParallel )
                {
                    if (debugMode)  std::cout <<"...Parallel optimizing by Edge. Iteration : "<< iter <<"\n";
                    ParallelEvaluateClusterByEdge((TVolumeMesh*)(m),vrelaxQuality,minAngle);
                    if (debugMode)  std::cout <<"...End. Iteration : "<< iter <<"\n";
                }
                else
                {
                    if (debugMode)  std::cout <<"...Optimizing by Edge. Iteration : "<< iter <<"\n";
                    evaluateClusterByEdge( (TVolumeMesh*)(m),minAngle,vrelaxQuality);
                }
                if (debugMode)
                    m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
            }

            if (debugMode)
            {
                std :: cout<< "Number of faces:" << m->fFaces->Count() << "\n";
                EvaluateQuality();
                m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
                m->validate(true);
                std :: cout<< "Number of faces:" << m->fFaces->Count() << "\n";
                showProcessTime();
            }
            // Save the mesh as generated from Kratos
            if (saveToFile)
            {
                TMeshLoader* ml2 = new TVMWLoader();
                std::string s("");
                s = "out_MeshFromKratos" + intToString(iter)+".vwm";
                // BUG Linux/Windows
                ml2->save(s , m);
                delete ml2;
            }
        }

        if (saveToFile)
        {
            TMeshLoader* ml2 = new TVMWLoader();
            std::string s("out_MeshC.vwm");
            ml2->save(s, m);
            delete ml2;
        }

        if (reinsertNodes)
        {
            if (debugMode)  std::cout <<"...Trying to reinsert nodes..." <<"\n";
            tryToReinsertNodes();
        }

        // Sino pudo mejorar la malla
        /*
        TVolumeMesh* m2 = (TVolumeMesh*)(GenerateMesh(m));
        m= m2;
        */
        ///------------------------------
        //construct spatial structure with an auxiliary model part

        std::cout <<"...Interpolating and adding new Nodes..." <<"\n";
        BinBasedFastPointLocator<3> element_finder(r_model_part);
        element_finder.UpdateSearchDatabase();
        InterpolateAndAddNewNodes(r_model_part, m,element_finder);

    }
    ///@brief function tryToReinsertNodes
    /// Reinsert removed nodes into the structure
    void tryToReinsertNodes()
    {
        int vToR =m->vertexesToRemove->Count();
        int ri = vertexTetraReInsertion(m , m->vertexesToRemove);
        m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
        std :: cout<< " Reinsert vertexes " << ri << " of " <<  vToR <<"\n";
        EvaluateQuality();
        m->validate(true);
        std :: cout<< "........................................"<<"\n";
    }
    ///@brief function FinalizeOptimization
    /// Destroy the structure
    void FinalizeOptimization(bool removeFreeVertexes )
    {
        if (m == nullptr ) return ;
        if (debugMode)  std::cout <<"...Output to Kratos Format" <<"\n";
        // Get back in Kratos
        innerConvertToKratos(refMP , m , removeFreeVertexes);
        delete m;
        m = nullptr;
        std::cout <<"...Trying to release memory\n";
        clearPool();
        std::cout <<"...Release memory OK\n";

    }


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
        std::stringstream buffer;
        buffer << "TetrahedraReconnectUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TetrahedraReconnectUtility";
    }

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TetrahedraReconnectUtility& operator=(TetrahedraReconnectUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    TetrahedraReconnectUtility(TetrahedraReconnectUtility const& rOther);


    ///@}

}; // Class TetrahedraReconnectUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TetrahedraReconnectUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TetrahedraReconnectUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED  defined 


