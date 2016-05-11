//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_VOLUME_CALCULATION_PROCESS_H_INCLUDED )
#define  KRATOS_VOLUME_CALCULATION_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"

///VARIABLES used:
//Data:     
//StepData: 
//Flags:    (checked) 
//          (set)     
//          (modified)  
//          (reset)   


namespace Kratos
{

///@name Kratos Classes
///@{

class ModelVolumeCalculationProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( ModelVolumeCalculationProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModelVolumeCalculationProcess(ModelPart& rModelPart, bool Axisymmetric, int EchoLevel) 
      : mrModelPart(rModelPart)
    {
      mAxisymmetric = Axisymmetric;
      mEchoLevel = EchoLevel;
    }

    ModelVolumeCalculationProcess(ModelPart& rModelPart) 
      : mrModelPart(rModelPart)
    {
      mAxisymmetric = false;
    } 

    /// Destructor.
    virtual ~ModelVolumeCalculationProcess() {}


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
    virtual void Execute() {}

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

      KRATOS_TRY

      mModelVolume = ComputeModelVolume(mMeshVolume);
      
      if( mEchoLevel > 0 )
	std::cout<<"  [ Model Volume : "<<mModelVolume<<" ]"<<std::endl;


      // unsigned int start=0;
      // unsigned int NumberOfMeshes=mrModelPart.NumberOfMeshes();
      // if(NumberOfMeshes>1) 
      // 	start=1;
      // for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      // 	{
      // 	  std::cout<<"   [ Mesh Volume ["<<MeshId<<"]: "<<mMeshVolume[MeshId]<<" ]"<<std::endl;
      // 	}

      KRATOS_CATCH( "" )
	
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
      KRATOS_TRY	 

      std::vector<double> MeshVolume;
      double ModelVolume = ComputeModelVolume(MeshVolume);
      
      if( mEchoLevel > 0 )
	std::cout<<"  [ Model Volume : "<<ModelVolume<<" ] [ Step increment : "<<ModelVolume-mModelVolume<<" ] "<<std::endl;

      // unsigned int start=0;
      // unsigned int NumberOfMeshes=mrModelPart.NumberOfMeshes();
      // if(NumberOfMeshes>1) 
      // 	start=1;
      // for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      // 	{
      // 	  std::cout<<"   [ Mesh Volume ["<<MeshId<<"]: "<<MeshVolume[MeshId]<<" ] [ Step increment : "<<MeshVolume[MeshId]-mMeshVolume[MeshId]<<" ] "<<std::endl;
      // 	}

      KRATOS_CATCH( "" )      
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
        return "ModelVolumeCalculationProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ModelVolumeCalculationProcess";
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
    ModelPart&  mrModelPart;
 
    bool mAxisymmetric;
 
    double mModelVolume;
  
    std::vector<double> mMeshVolume;
  
    int mEchoLevel;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ModelVolumeCalculationProcess& operator=(ModelVolumeCalculationProcess const& rOther);


    /// this function is a private function

    double ComputeModelVolume(std::vector<double> & rMeshVolume)
    {

      KRATOS_TRY

      unsigned int start=0;
      unsigned int NumberOfMeshes=mrModelPart.NumberOfMeshes();
      if(NumberOfMeshes>1) 
	start=1;

      
      rMeshVolume.resize(NumberOfMeshes+start);
      std::fill( rMeshVolume.begin(), rMeshVolume.end(), 0.0 );

      double ModelVolume = 0;

      if( mAxisymmetric ){

	//std::cout<<"  AXISYMMETRIC MODEL "<<std::endl;

	double radius = 0;
	double two_pi = 6.28318530717958647693;

	//By the way: set meshes options from bools
	for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	  {
	    for(ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin(MeshId); ie != mrModelPart.ElementsEnd(MeshId); ie++)
	      {
		
		const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

		if( dimension > 2 )
		  std::cout<<" Axisymmetric problem with dimension: "<<dimension<<std::endl;

		radius = 0;
		for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
		  radius += ie->GetGeometry()[i].X();

		radius/=double(ie->GetGeometry().size());
		
		rMeshVolume[MeshId] += ie->GetGeometry().Area() * two_pi * radius ;


	      }
	  
	    ModelVolume += rMeshVolume[MeshId];
	  }

      }
      else{

	//By the way: set meshes options from bools
	for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	  {
	    for(ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin(MeshId); ie != mrModelPart.ElementsEnd(MeshId); ie++)
	      {
		const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();
		if( dimension == 2){
		  rMeshVolume[MeshId] += ie->GetGeometry().Area() *  ie->GetProperties()[THICKNESS];
		}
		else if( dimension == 3){
		  rMeshVolume[MeshId] += ie->GetGeometry().Volume();
		}
		else{
		  //do nothing.
		}

	      }
	  
	    ModelVolume += rMeshVolume[MeshId];

	  }

      }

      return ModelVolume;

      KRATOS_CATCH( "" )

    }


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
                                  ModelVolumeCalculationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ModelVolumeCalculationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MODEL_VOLUME_CALCULATION_PROCESS_H_INCLUDED  defined 


