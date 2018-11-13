//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:        LMonforte $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      	 March 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined( KRATOS_CPT_REFINE_MESH_ELEMENTS_ON_SIZE_PROCESS_H_INCLUDED )
#define KRATOS_CPT_REFINE_MESH_ELEMENTS_ON_SIZE_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"

///VARIABLES used:
//Data:     
//StepData: NODAL_H, CONTACT_FORCE
//Flags:    (checked) TO_REFOME, BOUNDARY, NEW_ENTITY
//          (set)     
//          (modified)  
//          (reset)   
//(set):=(set in this process)

namespace Kratos
{

    ///@name Kratos Classes
    ///@{

    /// Refine Mesh Elements Process 2D and 3D
    /** The process labels the elements to be refined in the mesher
        it applies a size constraint to elements that must be refined.

    */
    class CptRefineMeshElementsOnSizeProcess
      : public Process
    {
        public:
         ///@name Type Definitions
         ///@{

         /// Pointer definition of Process
         KRATOS_CLASS_POINTER_DEFINITION( CptRefineMeshElementsOnSizeProcess );

         typedef ModelPart::ConditionType         ConditionType;
         typedef ModelPart::PropertiesType       PropertiesType;
         typedef ConditionType::GeometryType       GeometryType;

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         CptRefineMeshElementsOnSizeProcess(ModelPart& rModelPart,
               ModelerUtilities::MeshingParameters& rRemeshingParameters,
               int EchoLevel) 
            : mrModelPart(rModelPart),
            mrRemesh(rRemeshingParameters)
      {
         mEchoLevel = EchoLevel;
      }


         /// Destructor.
         virtual ~CptRefineMeshElementsOnSizeProcess() {}


         ///@}
         ///@name Operators
         ///@{

         /// This operator is provided to call the process as a function and simply calls the Execute method.
         void operator()()
         {
            Execute();
         }


         ///@}
         ///@name Operations
         ///@{


        /// Execute method is used to execute the Process algorithms.
        virtual void Execute()
        {
            KRATOS_TRY

            if( mEchoLevel > 0 ){
               std::cout<<" [ LMV SELECT ELEMENTS TO REFINE : "<<std::endl;
               //std::cout<<"   refine selection "<<std::endl;
            }

            //***SIZES :::: parameters do define the tolerance in mesh size: 
            double size_for_inside_elements   = 0.75 * mrRemesh.Refine->CriticalRadius;
            double size_for_boundary_elements = 1.50 * mrRemesh.Refine->CriticalRadius; 

            double nodal_h_refining_factor     = 0.75;
            double nodal_h_non_refining_factor = 2.00;

            ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();

            // 1. try to look for the cone position, or something similar
            
            /* loop over contacting nodes in order to find the 
             *      tip position yCone & xConeMin (smallest x and y)
             *      shaft radius xConeMax - xConeMin
             * 
             * define 2 distances of 2.0 and 4.0 times the cone radius
             */
             
            double yCone = 10.0;
            double xConeMax = 0.0;
            double xConeMin = 1000.0;
            int counter = 0;
            
            ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();
            for (unsigned int i = 0; i < mrModelPart.Nodes().size(); i++) 
            {
                if ( ( nodes_begin + i)->SolutionStepsDataHas( CONTACT_FORCE) ) {
                    array_1d<double, 3 > & ContactForce = (nodes_begin + i)->FastGetSolutionStepValue(CONTACT_FORCE) ;
                    if ( norm_2(ContactForce) > 0) {
                        double y = (nodes_begin + i)->Y();
                        double x = (nodes_begin + i)->X();
                        counter += 1;

                        if ( y < yCone)
                            yCone = y; 
                        if ( x < xConeMin)
                            xConeMin = x; 
                        if ( x > xConeMax )
                            xConeMax = x;
                    }
                }
            }
            
            double Radi = xConeMax - xConeMin; 
            double xCone = xConeMin + 0.5* Radi ;

            const double FirstDistance = 3.0 * Radi;
            const double SecondDistance = 6.0 * Radi;

            std::cout << " position of the cone: yCone = " << yCone << std::endl;


            // 2. Try to set the size depending on the relative position with the cone
            int id = 0; 
            if ( (mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ELEMENTS)
                    && mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ADD_NODES) ) )
            {
                std::cout<<"CHECK REFINE"<<std::endl;

                ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();

                unsigned int nds = (*element_begin).GetGeometry().size();

                ModelerUtilities::MeshContainer& InMesh = mrRemesh.InMesh;

                InMesh.CreateElementList(mrRemesh.Info->NumberOfElements, nds); //number of preserved elements
                InMesh.CreateElementSizeList(mrRemesh.Info->NumberOfElements);

                int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();

                int* InElementList        = mrRemesh.InMesh.GetElementList();
                double* InElementSizeList = mrRemesh.InMesh.GetElementSizeList();

                int* OutElementList       = mrRemesh.OutMesh.GetElementList();

                ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();


                //SET THE REFINED ELEMENTS AND THE AREA (NODAL_H)
                //*********************************************************************

                for(int el = 0; el< OutNumberOfElements; el++)
                {
                    if(mrRemesh.PreservedElements[el])
                    {
                        unsigned int    count_dissipative = 0;
                        /*double prescribed_h      = 0;
                        bool   dissipative       = false;  //dissipative means reference threshold is overwhelmed
                        bool   refine_size       = false;

                        unsigned int    count_boundary_inserted = 0;
                        unsigned int    count_boundary = 0;
                        unsigned int    count_contact_boundary = 0;*/ // LMV. not using it

                        Geometry<Node<3> > vertices;
                        double elementX = 0;
                        double elementY = 0;
                        for(unsigned int pn=0; pn<nds; pn++)
                        {

                            InElementList[id*nds+pn]= OutElementList[el*nds+pn];

                            vertices.push_back(*(nodes_begin + OutElementList[el*nds+pn]-1).base());

                            elementX += (nodes_begin+OutElementList[el*nds+pn]-1)->X();
                            elementY += (nodes_begin+OutElementList[el*nds+pn]-1)->Y();

                            if ( (nodes_begin+OutElementList[el*nds+pn]-1)->Is(TO_REFINE) )
                                count_dissipative+=1;
                        }

                        elementX /= 3.0;
                        elementY /= 3.0;
                        elementX -= xCone; 
                        elementY -= yCone;

                        double distance = elementX*elementX + elementY*elementY;
                        distance = std::sqrt(distance);

                        bool refine_candidate = true;

                        double element_size = 0;
                        double element_radius = mModelerUtilities.CalculateElementRadius(vertices, element_size);

                        double trianglearea;
                        double MyAlpha = 20.0;
                        if (distance < FirstDistance) {
                            MyAlpha = 0.5;
                            if (element_radius < 0.95 * size_for_inside_elements)
                                MyAlpha = 1.5;
                        } else if ( distance < SecondDistance) {
                            MyAlpha = 0.7;
                            if ( element_radius < 1.3* size_for_inside_elements )
                                MyAlpha = 2.0;
                        }
                        trianglearea = element_size * MyAlpha;
                        InElementSizeList[id] = trianglearea;
                        id += 1;

                        //********* PLASTIC POWER ENERGY REFINEMENT CRITERION (A)
                        if(count_dissipative>=nds-1){
                            mrRemesh.Info->CriticalElements += 1;
                        }
                    }
                }
            }
            else
            {
                // 
                std::cout<<"CHECK REFINE: FALSE"<<std::endl;
                
                ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();

                unsigned int nds = (*element_begin).GetGeometry().size();

                ModelerUtilities::MeshContainer& InMesh = mrRemesh.InMesh;

                InMesh.CreateElementList(mrRemesh.Info->NumberOfElements, nds); //number of preserved elements
                InMesh.CreateElementSizeList(mrRemesh.Info->NumberOfElements);

                int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();

                int* InElementList        = mrRemesh.InMesh.GetElementList();
                double* InElementSizeList = mrRemesh.InMesh.GetElementSizeList();

                int* OutElementList       = mrRemesh.OutMesh.GetElementList();

                ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();	  

                for(int el = 0; el< OutNumberOfElements; el++)
                {
                    if(mrRemesh.PreservedElements[el])
                    {
                        Geometry<Node<3> > vertices;
                        for(unsigned int pn=0; pn<nds; pn++)
                        {
                            InElementList[id*nds+pn]= OutElementList[el*nds+pn];
                            vertices.push_back(*(nodes_begin + OutElementList[el*nds+pn]-1).base());
                        }

                        double element_size = 0;
                        mModelerUtilities.CalculateElementRadius(vertices, element_size);

                        InElementSizeList[id] = nodal_h_non_refining_factor * element_size;

                        id++;
                    }
                }
            }
            
            if( mEchoLevel > 0 )
            {
                std::cout<<"   Visited Elements: "<<id<<std::endl;
                std::cout<<"   LMV SELECT ELEMENTS TO REFINE ]; "<<std::endl;
            }

            KRATOS_CATCH(" ")
        }


         /// this function is designed for being called at the beginning of the computations
         /// right after reading the model and the groups
         virtual void ExecuteInitialize()
         {
         }

         /// this function is designed for being execute once before the solution loop but after all of the
         /// solvers where built
         virtual void ExecuteBeforeSolutionLoop()
         {
         }

         /// this function will be executed at every time step BEFORE performing the solve phase
         virtual void ExecuteInitializeSolutionStep()
         {	
         }

         /// this function will be executed at every time step AFTER performing the solve phase
         virtual void ExecuteFinalizeSolutionStep()
         {
         }

         /// this function will be executed at every time step BEFORE  writing the output
         virtual void ExecuteBeforeOutputStep()
         {
         }

         /// this function will be executed at every time step AFTER writing the output
         virtual void ExecuteAfterOutputStep()
         {
         }

         /// this function is designed for being called at the end of the computations
         /// right after reading the model and the groups
         virtual void ExecuteFinalize()
         {
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
            return "CptRefineMeshElementsOnSizeProcess";
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const
         {
            rOStream << "CptRefineMeshElementsOnSizeProcess";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const
         {
         }


         ///@}
         ///@name Friends
         ///@{

         ///@}


      private:
         ///@name Static Member Variables
         ///@{

         ///@}
         ///@name Static Member Variables
         ///@{
         ModelPart& mrModelPart;

         ModelerUtilities::MeshingParameters& mrRemesh;

         ModelerUtilities mModelerUtilities;  

         int mEchoLevel;

         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}



         /// Assignment operator.
         CptRefineMeshElementsOnSizeProcess& operator=(CptRefineMeshElementsOnSizeProcess const& rOther);


         /// this function is a private function


         /// Copy constructor.
         //Process(Process const& rOther);


         ///@}

    }; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
      CptRefineMeshElementsOnSizeProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
      const CptRefineMeshElementsOnSizeProcess& rThis)
{
   rThis.PrintInfo(rOStream);
   rOStream << std::endl;
   rThis.PrintData(rOStream);

   return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_CPT_REFINE_MESH_ELEMENTS_ON_SIZE_PROCESS_H_INCLUDED defined 

