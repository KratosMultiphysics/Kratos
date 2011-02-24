//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_TETGEN_VOLUME_MESHER_H_INCLUDED )
#define  KRATOS_TETGEN_VOLUME_MESHER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 

#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "utilities/geometry_utilities.h"
#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"


namespace Kratos
{
    ///@addtogroup MeshingApplication
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

    /// This class performs the (constrained) Delaunay meshing of the internal volume given a surface discretized in terms of triangular conditions

    /** This class performs the meshing of the internal volume given a surface discretized in terms of triangular conditions.
     * Meshing is performed using the Tetgen mesher by Hang Si
     */
    class TetgenVolumeMesher
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of TetgenVolumeMesher
        KRATOS_CLASS_POINTER_DEFINITION(TetgenVolumeMesher);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        TetgenVolumeMesher(ModelPart& rmodel_part)
        :mrModelPart(rmodel_part)
        {
        }

        /// Destructor.

        virtual ~TetgenVolumeMesher()
        {
        }

        /**
         * This function is designed to add volumes, to carve out cavities inside a domain.
         * The user is expected to manually prescribe the position of the cavities.
         * As many "Holes" as needed can be prescribed"
         * @param x - x coordinate of the cavity
         * @param y - y coordinate of the cavity
         * @param z - z coordinate of the cavity
         */
        void AddHole(double x, double y, double z)
        {
            array_1d<double, 3 >  hole_coordinates;
            hole_coordinates[0] = x;
            hole_coordinates[1] = y;
            hole_coordinates[2] = z;
            mholes.push_back(hole_coordinates);
        }

        /**
         * This function performs the actual meshing. The conditions used to define the
         * surface will be mantained in the final mesh.
         * The user is expected to prescribe the command line settings to be used accordin to the Tetgen manual.
         * WARNING: please note that YY is needed between the parameters
         * @param settings - string containing the command line settings
         * @param ElementName - name of the type of elements to be created \n
         *        please note that the element type given should be registered in the Kratos
         */
        void GenerateMesh(std::string settings, std::string ElementName)
        {
            tetgenio in, out;

            if(mrModelPart.Conditions().size() == 0)
                KRATOS_ERROR(std::logic_error,"model part does not have faces","")
            if(mrModelPart.Elements().size() != 0)
                KRATOS_ERROR(std::logic_error,"model part already has volume elements","")

            //reorder node Ids consecutively
            ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();
            for (unsigned int i = 0; i < mrModelPart.Nodes().size(); i++)
            {
                (nodes_begin + i)->SetId(i + 1);
            }

            // All indices start from 1.
            in.initialize();
            in.firstnumber = 1;

            in.numberofpoints = mrModelPart.Nodes().size();
            in.pointlist = new REAL[in.numberofpoints * 3];

            //give the coordinates to the mesher
            
            for (unsigned int i = 0; i < mrModelPart.Nodes().size(); i++)
            {
                int base = i * 3;
                in.pointlist[base] = (nodes_begin + i)->X();
                in.pointlist[base + 1] = (nodes_begin + i)->Y();
                in.pointlist[base + 2] = (nodes_begin + i)->Z();
            }

            //prepare the list of faces
            in.numberoffacets = mrModelPart.Conditions().size();
            in.facetlist = new tetgenio::facet[in.numberoffacets];
            in.facetmarkerlist = new int[in.numberoffacets];

            //give the surface connectivity to the mesher
            ModelPart::ConditionsContainerType::iterator condition_begin = mrModelPart.ConditionsBegin();
            for (unsigned int i = 0; i < mrModelPart.Conditions().size(); i++)
            {
                tetgenio::facet* f = &in.facetlist[i];
                f->numberofpolygons = 1;
                f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
                f->numberofholes = 0;
                f->holelist = NULL;
                tetgenio::polygon* p = &f->polygonlist[0];
                p->numberofvertices = 3;
                p->vertexlist = new int[p->numberofvertices];
                p->vertexlist[0] = (condition_begin + i)->GetGeometry()[0].Id();
                p->vertexlist[1] = (condition_begin + i)->GetGeometry()[1].Id();
                p->vertexlist[2] = (condition_begin + i)->GetGeometry()[2].Id();

                in.facetmarkerlist[i] = static_cast<int>( (condition_begin + i)->GetValue(FLAG_VARIABLE) );
            }

            //give the hole list to the mesher
            in.numberofholes = mholes.size();
            in.holelist = new double[in.numberofholes*3];
            for (unsigned int i = 0; i < mholes.size(); i++)
            {
                int base = i * 3;
                in.holelist[base] =     mholes[i][0];
                in.holelist[base + 1] = mholes[i][1];
                in.holelist[base + 2] = mholes[i][2];
            }


//            char tetgen_options[] = settings.c_str() ; //"pqYY";
//                        //perform meshing with tetgen
//            tetrahedralize(tetgen_options, &in, &out);

            char *tmp=new char[settings.size()+1];
            tmp[settings.size()]=0;
            memcpy(tmp,settings.c_str(),settings.size());

            tetrahedralize(tmp, &in, &out);

            delete [] tmp;




            //generate new nodes
            Node < 3 > ::DofsContainerType& reference_dofs = (mrModelPart.NodesBegin())->GetDofs();
            for(int i=in.numberofpoints; i< out.numberofpoints; i++)
            {
                int base = i * 3;
                double x = out.pointlist[base];
                double y = out.pointlist[base+1];
                double z = out.pointlist[base+2];

                Node<3>::Pointer p_new_node = Node<3>::Pointer(new Node<3>(i+1, x, y, z));

                // Giving model part's variables list to the node
                p_new_node->SetSolutionStepVariablesList(&(mrModelPart.GetNodalSolutionStepVariablesList()));
                p_new_node->SetBufferSize(mrModelPart.NodesBegin()->GetBufferSize());
                for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                {
                    Node < 3 > ::DofType& rDof = *iii;
                    Node < 3 > ::DofType::Pointer p_new_dof = p_new_node->pAddDof(rDof);
                    (p_new_dof)->FreeDof(); //the variables are left as free for the internal node
                }

                mrModelPart.Nodes().push_back(p_new_node);
            }

            //generate new Elements
            Element const& rReferenceElement = KratosComponents<Element>::Get(ElementName);
            nodes_begin = mrModelPart.NodesBegin();
            Properties::Pointer properties = mrModelPart.GetMesh().pGetProperties(0);

            for(int i=0; i!=out.numberoftetrahedra; i++)
            {
                int base=i*4;

                 Tetrahedra3D4<Node<3> > geom(
						*( (nodes_begin +  out.tetrahedronlist[base]-1).base() 	),
						*( (nodes_begin +  out.tetrahedronlist[base+1]-1).base() 	),
						*( (nodes_begin +  out.tetrahedronlist[base+2]-1).base() 	),
						*( (nodes_begin +  out.tetrahedronlist[base+3]-1).base() 	)
						);

                 Element::Pointer p_element = rReferenceElement.Create(i+1, geom, properties);
	        (mrModelPart.Elements()).push_back(p_element);
            }

            in.deinitialize();
            out.deinitialize();
            in.initialize(); //better deinitialize and initialize again
            out.initialize();


        }

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

        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "TetgenVolumeMesher";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "TetgenVolumeMesher";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
        {
        }


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
        std::vector< array_1d<double, 3 > > mholes;
        ModelPart& mrModelPart;


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

        TetgenVolumeMesher & operator=(TetgenVolumeMesher const& rOther)
        {
            return *this;
        }

        /// Copy constructor.

        TetgenVolumeMesher(TetgenVolumeMesher const& rOther):
        mrModelPart(rOther.mrModelPart)
        {
        }


        ///@}

    }; // Class TetgenVolumeMesher

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function

    inline std::istream & operator >>(std::istream& rIStream,
            TetgenVolumeMesher& rThis)
    {
        return rIStream;
    }

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const TetgenVolumeMesher& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TETGEN_VOLUME_MESHER_H_INCLUDED  defined


