//Author: M.A. Celigueta (CIMNE, maceli@cimne.upc.edu)
// Date: January 2017


// System includes
#include <string>
#include <iostream> 

// Project includes
#include "bonding_spheric_continuum_particle.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{             
	void BondingSphericContinuumParticle::UpdateContinuumNeighboursVector(ProcessInfo& r_process_info) {
            KRATOS_TRY
            if(mNeighbourElements.size() == mContinuumInitialNeighborsSize) return;
            
            const unsigned int initial_number_of_cont_neighbours = mContinuumInitialNeighborsSize;
            
            std::vector<SphericParticle*> cont_ini_neighbour_elems;
            std::vector<int> cont_ini_ids;
            std::vector<double> cont_ini_deltas;
            std::vector<int> cont_ini_failure_ids;  
            std::vector<SphericParticle*> discont_ini_neighbour_elems; 
            std::vector<int> discont_ini_ids;
            std::vector<double> discont_ini_deltas;
            std::vector<SphericParticle*> discont_neighbour_elems;
                        
            //Assuming that the neighbours are ordered (first initial continuum, then initial discontinuum, then other discontinuum)
            for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
                cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                cont_ini_ids.push_back(mIniNeighbourIds[i]);
                cont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                cont_ini_failure_ids.push_back(mIniNeighbourFailureId[i]);  
            }            

            //// Discontinuum and initial
            for (unsigned int i = mContinuumInitialNeighborsSize; i < mInitialNeighborsSize; i++) {
                if(mNeighbourElements[i] == NULL) {
                    discont_ini_ids.push_back(mIniNeighbourIds[i]);
                    discont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    discont_ini_neighbour_elems.push_back(mNeighbourElements[i]);   
                    continue;
                }
                
                if(NeighbourIsToBeBonded(mNeighbourElements[i]->Id())){
                    cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                    cont_ini_ids.push_back(mIniNeighbourIds[i]);
                    cont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    cont_ini_failure_ids.push_back(0);                        
                } else {
                    discont_ini_ids.push_back(mIniNeighbourIds[i]);
                    discont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    discont_ini_neighbour_elems.push_back(mNeighbourElements[i]);                    
                }            
            }            
            
            //The rest (discontinuum not initial)
            for (unsigned int i = mInitialNeighborsSize; i < mNeighbourElements.size(); i++) {               
                if(NeighbourIsToBeBonded(mNeighbourElements[i]->Id())){
                    cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                    cont_ini_ids.push_back(mNeighbourElements[i]->Id());
                    cont_ini_deltas.push_back(0.0);
                    cont_ini_failure_ids.push_back(0);                        
                } else {
                    discont_neighbour_elems.push_back(mNeighbourElements[i]);                    
                }            
            }

            //Adding everything up in single vectors mNeighbourElements, mIniNeighbourIds, mIniNeighbourDelta  
            mContinuumInitialNeighborsSize = cont_ini_neighbour_elems.size();
            mInitialNeighborsSize = mContinuumInitialNeighborsSize + discont_ini_neighbour_elems.size();
            
            mIniNeighbourIds.resize(mInitialNeighborsSize);
            mIniNeighbourDelta.resize(mInitialNeighborsSize);
            mIniNeighbourFailureId.resize(mContinuumInitialNeighborsSize);
            
            //Continuum 
            for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
                mNeighbourElements[i] = cont_ini_neighbour_elems[i];
                mIniNeighbourIds[i]   = cont_ini_ids[i];
                mIniNeighbourDelta[i] = cont_ini_deltas[i]; 
                mIniNeighbourFailureId[i] = cont_ini_failure_ids[i];
            }

            mContinuumConstitutiveLawArray.resize(mContinuumInitialNeighborsSize);

            for (unsigned int i = initial_number_of_cont_neighbours; i < mContinuumInitialNeighborsSize; i++) {
                DEMContinuumConstitutiveLaw::Pointer NewContinuumConstitutiveLaw = GetProperties()[DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER]-> Clone();
                mContinuumConstitutiveLawArray[i] = NewContinuumConstitutiveLaw;
                mContinuumConstitutiveLawArray[i]->Initialize();
                
                mBondElements.push_back(NULL);
            }                      

            KRATOS_CATCH("")
	};        
        
        bool BondingSphericContinuumParticle::NeighbourIsToBeBonded(const int neighbour_id) {
            for(int i=0; i<(int)mIdsOfElementsToBeBonded.size(); i++) {
                if(mIdsOfElementsToBeBonded[i] == neighbour_id) {
                    return true;
                }
            }
            return false;
        }
        
        void BondingSphericContinuumParticle::ComputeForceWithNeighbourFinalOperations(){
            SphericContinuumParticle::ComputeForceWithNeighbourFinalOperations();
            //TODO: ADD SOMETHING HERE like 'mIdsOfElementsToBeBonded.push_back(id_of_neighbour)'
            
             for(int i=mContinuumInitialNeighborsSize; i<(int)mNeighbourElements.size(); i++) {
                  //if (> ) { 
                mIdsOfElementsToBeBonded.push_back(i);
                //mIdsOfElementsToBeBonded.push_back(mNeighbourElements[i]);

                 // }
             }
        }
                       
}  // namespace Kratos.
