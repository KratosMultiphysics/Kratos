//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined( KRATOS_SELECT_FLUID_ELEMENTS_MESHER_PROCESS_H_INCLUDED )
#define KRATOS_SELECT_FLUID_ELEMENTS_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "custom_utilities/mesher_utilities.hpp"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:
//StepData: NODAL_H, CONTACT_FORCE
//Flags:    (checked) TO_ERASE, BOUNDARY, NEW_ENTITY
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
class SelectFluidElementsMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( SelectFluidElementsMesherProcess );

  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  SelectFluidElementsMesherProcess(ModelPart& rModelPart,
                                     MesherUtilities::MeshingParameters& rRemeshingParameters,
                                     int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~SelectFluidElementsMesherProcess() {}


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
  void Execute() override
  {
    KRATOS_TRY

    if( mEchoLevel > 1 ){
      std::cout<<" [ SELECT MESH ELEMENTS in PfemFluid: ("<<mrRemesh.OutMesh.GetNumberOfElements()<<") "<<std::endl;
      std::cout<<"MODEL PART InNumberOfElements "<<mrRemesh.InMesh.GetNumberOfElements()<<std::endl;
      std::cout<<"MODEL PART InNumberOfPoints "<<mrRemesh.InMesh.GetNumberOfPoints()<<std::endl;
      std::cout<<"MODEL PART OutNumberOfElements "<<mrRemesh.OutMesh.GetNumberOfElements()<<std::endl;
      std::cout<<"MODEL PART OutNumberOfPoints "<<mrRemesh.OutMesh.GetNumberOfPoints()<<std::endl;
    }

    int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
    mrRemesh.PreservedElements.clear();
    mrRemesh.PreservedElements.resize(OutNumberOfElements);
    std::fill( mrRemesh.PreservedElements.begin(), mrRemesh.PreservedElements.end(), 0 );
    mrRemesh.MeshElementsSelectedFlag = true;


    mrRemesh.Info->NumberOfElements=0;

    const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
    double currentTime = rCurrentProcessInfo[TIME];
    double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    bool firstMesh=false;
    if(currentTime<2*timeInterval){
      firstMesh=true;
    }

    bool box_side_element = false;
    bool wrong_added_node = false;

    int number_of_slivers = 0;

    if(mrRemesh.ExecutionOptions.IsNot(MesherUtilities::SELECT_TESSELLATION_ELEMENTS))
    {
      for(int el=0; el<OutNumberOfElements; ++el)
      {
        mrRemesh.PreservedElements[el]=1;
        mrRemesh.Info->NumberOfElements+=1;
      }
    }
    else
    {
      if( mEchoLevel > 1 )
        std::cout<<"   Start Element Selection "<<OutNumberOfElements<<std::endl;

      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();

      const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

      unsigned int nds = 3;  //linear triangle
      if(dimension==3) //linear tetrahedron
        nds = 4;

      int* OutElementList = mrRemesh.OutMesh.GetElementList();

      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

      int el = 0;
      int number = 0;
      // std::cout<<" num nodes "<<rNodes.size()<<std::endl;
      // std::cout<<" NodalPreIdsSize "<<mrRemesh.NodalPreIds.size()<<std::endl;

      //#pragma omp parallel for reduction(+:number) private(el)
      for(el=0; el<OutNumberOfElements; ++el)
      {
        Geometry<Node<3> > vertices;

        unsigned int  numfreesurf =0;
        unsigned int  numboundary =0;
        unsigned int  numrigid =0;
        unsigned int  numsolid =0;
        unsigned int  numfluid =0;
        unsigned int  numinlet =0;
        unsigned int  numisolated =0;
        // unsigned int  numinsertednodes =0;
        std::vector<double > normVelocityP;
        normVelocityP.resize(nds);
        unsigned int  checkedNodes =0;
        box_side_element = false;


        for(unsigned int pn=0; pn<nds; ++pn)
        {
          //set vertices
          if(mrRemesh.NodalPreIds[OutElementList[el*nds+pn]]<0){
            if(mrRemesh.Options.IsNot(MesherUtilities::CONTACT_SEARCH))
              std::cout<<" ERROR: something is wrong: nodal id < 0 "<<std::endl;
            box_side_element = true;
            break;
          }

          if(OutElementList[el*nds+pn]<=0)
            std::cout<<" ERROR: something is wrong: nodal id < 0 "<<el<<std::endl;


          if( (unsigned int)OutElementList[el*nds+pn] > mrRemesh.NodalPreIds.size() ){
            wrong_added_node = true;
            std::cout<<" ERROR: something is wrong: node out of bounds "<<std::endl;
            break;
          }

          //vertices.push_back( *((rNodes).find( OutElementList[el*nds+pn] ).base() ) );
          vertices.push_back(rNodes(OutElementList[el*nds+pn]));

          //check flags on nodes
          if(vertices.back().Is(ISOLATED)){
            numisolated++;
          }

          if(vertices.back().Is(BOUNDARY)){
            numboundary++;
          }

          if(vertices.back().Is(RIGID)){
            numrigid++;
          }

          if(vertices.back().Is(SOLID)){
            numsolid++;
          }

          if( vertices.back().Is(FREE_SURFACE) ){
            numfreesurf++;
          }

          if(vertices.back().Is(FLUID)){
            numfluid++;
          }

          if(vertices.back().Is(FREE_SURFACE) || vertices.back().Is(ISOLATED)){
            const array_1d<double,3> &velocityP0=vertices.back().FastGetSolutionStepValue(VELOCITY,0);
            normVelocityP[pn]=norm_2(velocityP0);
            checkedNodes++;
          }

          if(vertices.back().Is(INLET)){
            numinlet++;
          }
        }


        if(box_side_element || wrong_added_node){
          std::cout<<" ,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Box_Side_Element "<<std::endl;
          continue;
        }


        double Alpha =  mrRemesh.AlphaParameter; //*nds;

        if(dimension==2){

          if( (numisolated+numfreesurf)==nds && firstMesh==false ){
            if(checkedNodes==nds){
              const double maxValue=1.5;
              const double minValue=1.0/maxValue;
              if(normVelocityP[0]/normVelocityP[1]>maxValue || normVelocityP[0]/normVelocityP[1]<minValue ||
                 normVelocityP[0]/normVelocityP[2]>maxValue || normVelocityP[0]/normVelocityP[2]<minValue ||
                 normVelocityP[1]/normVelocityP[2]>maxValue || normVelocityP[1]/normVelocityP[2]<minValue){
                Alpha*=0;
              }
            }else{
              KRATOS_INFO( "ATTENTION!!! CHECKED NODES= " ) <<checkedNodes<<" != "<<nds<<" [numfreesurf:"<<numfreesurf<<" numisolated:"<<numisolated<<std::endl;
              Alpha*=0;
            }
          }



          if(numrigid==0 && numsolid==0 && numfreesurf==0 && numisolated==0){
            Alpha*=1.75;
          }
          else{

            if( numboundary==nds && (numsolid !=0 || numrigid!=0) ){
              Alpha*=0.95;
            }
            else{
              Alpha*=1.04;
            }
          }

        }
        else if(dimension==3){

           if( (numisolated+numfreesurf)==nds ){
            if(checkedNodes==nds){
              const double maxValue=1.5;
              const double minValue=1.0/maxValue;
              if(normVelocityP[0]/normVelocityP[1]<minValue || normVelocityP[0]/normVelocityP[2]<minValue || normVelocityP[0]/normVelocityP[3]<minValue ||
                 normVelocityP[0]/normVelocityP[1]>maxValue || normVelocityP[0]/normVelocityP[2]<maxValue || normVelocityP[0]/normVelocityP[3]>maxValue ||
                 normVelocityP[1]/normVelocityP[2]<minValue || normVelocityP[1]/normVelocityP[3]<minValue ||
                 normVelocityP[1]/normVelocityP[2]>maxValue || normVelocityP[1]/normVelocityP[3]<maxValue ||
                 normVelocityP[2]/normVelocityP[3]<minValue ||
                 normVelocityP[2]/normVelocityP[3]>maxValue){
                Alpha*=0;
              }
            }else{
              KRATOS_INFO( "ATTENTION!!! CHECKED NODES= " ) <<checkedNodes<<" != "<<nds<<" [numfreesurf:"<<numfreesurf<<" numisolated:"<<numisolated<<std::endl;
              Alpha*=0;
            }

          }

          if(numrigid==0 && numfreesurf==0 && numisolated==0){
            Alpha*=1.75;
          }else{
            Alpha*=1.125;
          }

          if(numrigid==nds || (numrigid+numsolid)>=nds){
            Alpha*=0.95;
          }

        }
        if(firstMesh==true){
          Alpha*=1.15;
        }

        // Alpha*=1.175;

        bool accepted=false;

        MesherUtilities MesherUtils;

        if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
        {
          accepted=MesherUtils.ShrankAlphaShape(Alpha,vertices,mrRemesh.OffsetFactor,dimension);
        }
        else
        {
          double MeanMeshSize=mrRemesh.Refine->CriticalRadius;
          accepted=MesherUtils.AlphaShape(Alpha,vertices,dimension,MeanMeshSize);
        }

        //3.1.-
        bool self_contact = false;
        if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
          self_contact = MesherUtils.CheckSubdomain(vertices);

        //4.- to control that the element is inside of the domain boundaries
        if(accepted)
        {
          if(mrRemesh.Options.Is(MesherUtilities::CONTACT_SEARCH))
          {
            accepted=MesherUtils.CheckOuterCentre(vertices,mrRemesh.OffsetFactor, self_contact);
          }
        }

        if(numrigid==nds || (numrigid+numsolid)>=nds){
          accepted=false;
          //   if(isolatedWallElement==true || (dimension==2 && countIsolatedWallNodes==0)){
          //   accepted=false;
          // }
        }

        //5.- to control that the element has a good shape
        if(accepted && (numfreesurf>0 || numboundary == nds || numboundary-(numrigid+numsolid) > 0))
        {

          if(dimension==3 && nds==4){
            Geometry<Node<3> >* tetrahedron = new Tetrahedra3D4<Node<3> > (vertices);
            double Volume = tetrahedron->Volume();
            double CriticalVolume=0.01*mrRemesh.Refine->MeanVolume;
            //std::cout<<"CriticalVolume "<<Volume<<" MeanVolume "<<std::endl;
            if(Volume<CriticalVolume){
              KRATOS_INFO("SLIVER")<<" Volume="<<Volume<<" VS Critical Volume="<<CriticalVolume<<std::endl;
              accepted = false;
              number_of_slivers++;
            }
            delete tetrahedron;
          }

        }

        if(accepted)
        {
          //std::cout<<" Element ACCEPTED after cheking Center "<<number<<std::endl;
          number+=1;
          mrRemesh.PreservedElements[el] = number;
        }

      }

      mrRemesh.Info->NumberOfElements=number;

    }

    if( mEchoLevel > 1 ){
      std::cout<<"Number of Preserved Fluid Elements "<<mrRemesh.Info->NumberOfElements<<" (slivers detected: "<<number_of_slivers<<") "<<std::endl;
      std::cout<<"TOTAL removed nodes "<<mrRemesh.Info->RemovedNodes<<std::endl;
    }
    if(mrRemesh.ExecutionOptions.IsNot(MesherUtilities::KEEP_ISOLATED_NODES)){


      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
      const unsigned int nds = (*element_begin).GetGeometry().size();

      int* OutElementList = mrRemesh.OutMesh.GetElementList();

      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

      //check engaged nodes
      for(int el=0; el<OutNumberOfElements; ++el)
      {
        if( mrRemesh.PreservedElements[el] ){
          for(unsigned int pn=0; pn<nds; ++pn)
          {
            //set vertices
            rNodes[OutElementList[el*nds+pn]].Set(BLOCKED);

          }

        }

      }

      int count_release = 0;
      for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; ++i_node)
      {

        if( i_node->IsNot(BLOCKED)  ){
          if(!(i_node->Is(FREE_SURFACE) || i_node->Is(RIGID))){
            i_node->Set(TO_ERASE);
            if( mEchoLevel > 0 )
              std::cout<<" NODE "<<i_node->Id()<<" RELEASE "<<std::endl;
            if( i_node->Is(BOUNDARY) )
              std::cout<<" ERROR: node "<<i_node->Id()<<" IS BOUNDARY RELEASE "<<std::endl;
          }
        }
        if( i_node->Is(TO_ERASE)  ){
          count_release++;
        }

        i_node->Reset(BLOCKED);
      }

      if( mEchoLevel > 0 )
        std::cout<<"   fluid NUMBER OF RELEASED NODES "<<count_release<<std::endl;

    }
    else{

      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

      for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; ++i_node)
      {
        i_node->Reset(BLOCKED);
      }

    }


    mrRemesh.InputInitializedFlag = false;
    // mMesherUtilities.SetNodes(mrModelPart,mrRemesh);
    mMesherUtilities.SetNodes(mrModelPart,mrRemesh);
    mrRemesh.InputInitializedFlag = true;

    if( mEchoLevel > 1 ){
      std::cout<<"   Generated_Elements :"<<OutNumberOfElements<<std::endl;
      std::cout<<"   Passed_AlphaShape  :"<<mrRemesh.Info->NumberOfElements<<std::endl;
      std::cout<<"   SELECT MESH ELEMENTS ]; "<<std::endl;
    }

    KRATOS_CATCH( "" )

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
  std::string Info() const override
  {
    return "SelectFluidElementsMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "SelectFluidElementsMesherProcess";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
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

  MesherUtilities::MeshingParameters& mrRemesh;

  MesherUtilities mMesherUtilities;

  int mEchoLevel;

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
  SelectFluidElementsMesherProcess& operator=(SelectFluidElementsMesherProcess const& rOther);


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
                                  SelectFluidElementsMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SelectFluidElementsMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SELECT_FLUID_ELEMENTS_MESHER_PROCESS_H_INCLUDED defined
