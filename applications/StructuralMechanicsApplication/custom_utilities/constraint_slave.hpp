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
//   Last Modified by:    $Author: Aditya Ghantasala $
//   Date:                $Date: 14-03-2017 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(CONSTRAINT_SLAVE_H_INCLUDED)
#define CONSTRAINT_SLAVE_H_INCLUDED

#include<vector>
#include<map>
#include <iostream>

namespace Kratos
{

    /** \brief Quaternion
	* A simple class that implements the main features of quaternion algebra
	*/
	class MpcData : public Flags
	{

	public:
		
		typedef Dof<double>::Pointer DofPointerType;
		typedef Dof<double> DofType;

		///@name Life Cycle
		///@{

		/**
		Creates a Quaternion from its coefficients.
		@param id : The ID of the DOF which is considered as a slave
		*/
		MpcData(){
			this->totalNumberOfMasterDOFs = 0;
		}
        /// Destructor.
        virtual ~MpcData(){};
		
		///@}
		
	public:
		
		///@name Operators
		///@{
		
		
		///@}
		
	public:

		///@name Access
		///@{

		
		/**
		Adds a master and a corresponding weight to the current slave
		@param slaveDOF the DOF number which is to be added
		*/
		void AddSlave(DofType slaveDOF){
			
            if(this->slaveDofMasterDataMap.count(slaveDOF) > 0){
                
            }else{
                slaveDofVector.push_back(slaveDOF);
                masterData p;
                this->slaveDofMasterDataMap[slaveDOF] = p;                
            }
		}
		/**
		Adds a slave to the data structure
		@param slaveDOF the slave DOF to which master DOF is to be added
		@param masterDOF the DOF of the master DOF which is to be added
		@param masterWeight the weight of the masterDOF
		*/
		void AddMaster(DofType slaveDOF, DofType masterDOF, double masterWeight){
			this->slaveDofMasterDataMap[slaveDOF].masterDOFs.push_back(masterDOF);
			this->slaveDofMasterDataMap[slaveDOF].masterWeights.push_back(masterWeight);
		}
		
		/**
		Get the MasterDOFs vector for this slave
		@return MasterDOFs vector for this slave
		*/
		const std::vector<DofType> GetMasterDOFs(DofType slaveDOF){
			return this->slaveDofMasterDataMap[slaveDOF].masterDOFs;
		}

		/**
		Get the weights for corresponding masters of a slave
		@return weights for corresponding masters of a slave
		*/
		const std::vector<double> GetMasterWeightsForSlave(DofType slaveDOF){
			return this->slaveDofMasterDataMap[slaveDOF].masterWeights;
		}

		/**
		Get the Total number of MasterDOFs in the current Node
		@return Total number of MasterDOFs in the current Node
		 */
		const int GetTotalNumbeOfMasterDOFs(){
			this->totalNumberOfMasterDOFs = 0;
			for( const auto it : this->slaveDofMasterDataMap ){
				this->totalNumberOfMasterDOFs += (it.second).masterDOFs.size();
			}
			return this->totalNumberOfMasterDOFs;
		}


		/**
		Get the SlaveDOFs vector for this Node
		@return SlaveDOFs vector for this Node
		*/
		const std::vector<DofType>& GetSlaveDOFs(){
			return this->slaveDofVector;
		}

		///@

		///@name Static Operations
		///@{
		/**
		 * Returns the string containing a detailed description of this object.
		 * @return the string with informations
		 */
		virtual void GetInfo()const {
			std::cout << std::endl;
			std::cout << "===============================================================" << std::endl;
			std::cout << "This is adfa a MpcData Object" << std::endl;
 			std::cout << "Total number of slaves :: "<<(slaveDofVector).size()<< std::endl;
			for( const auto it : this->slaveDofMasterDataMap ){
				std::cout<<"SlaveDOF :: "<<(it.first).EquationId()<<std::endl;
				for(unsigned int i=0; i<(it.second).masterDOFs.size(); ++i ){
					std::cout<<"\t MasterDOFs :: "<< ((it.second).masterDOFs[i]).EquationId()<<std::endl;
					std::cout<<"\t MasterWeights ::"<< (it.second).masterWeights[i]<< std::endl;
				}
			}
			std::cout << "===============================================================" << std::endl;
			std::cout << std::endl;
		}
		
		///@}
		virtual void PrintInfo(std::ostream& rOStream) const {
			rOStream <<" MpcData object "<<std::endl;
		}

	private:

		///@name Member Variables
		///@{
		std::vector<DofType> slaveDofVector;
		unsigned int totalNumberOfMasterDOFs;

		struct masterData{
			std::vector<DofType> masterDOFs;
			std::vector<double> masterWeights;
		};

		std::map<DofType, masterData> slaveDofMasterDataMap;
		///@}

		///@name Serialization
		///@{
		friend class Serializer;

		    virtual void save(Serializer& rSerializer) const
		    {
		        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );

		    }

		    virtual void load(Serializer& rSerializer)
		    {
		        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
		    }

		///@}

		
	};

	///@name Input/Output funcitons
	///@{

	inline std::istream & operator >> (std::istream & rIStream, MpcData & rThis);

	inline std::ostream & operator << (std::ostream & rOStream, const MpcData & rThis)
	{
	    return rOStream;
	}

	///@}

        
} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
