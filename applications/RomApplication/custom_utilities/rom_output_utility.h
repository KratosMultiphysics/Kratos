// KRATOS ┏━━━┓
//        ┃┏━┓┃
//        ┃┗━┛┣━━┳┓┏┓
//        ┃┏┓┏┫┏┓┃┗┛┃
//        ┃┃┃┗┫┗┛┃┃┃┃
//        ┗┛┗━┻━━┻┻┻┛Application
//
//  License:         BSD License
//
//  Main authors:       Sebastian  Ares de Parga
//                      Raul Bravo
//

#if !defined( ROM_OUTPUT_UTILITY_H_INCLUDED )
#define  ROM_OUTPUT_UTILITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"

/* Application includes */
#include "rom_application_variables.h"

namespace Kratos
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;
    typedef Variable<double>             DoubleVarType;

    // This utility returns the Snapshot of a given time-step for Structural and Fluid Mechanics simulations.
    class RomOutputUtility
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(RomOutputUtility);

        RomOutputUtility(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ): mpModelPart(rModelPart){
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "main_model_part": "Structure"
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        //Assigning the process' variables
        mModelPartName = ThisParameters["main_model_part"].GetString();
    }

        ~RomOutputUtility()= default;

        void Start(){
            const auto &r_process_info = mpModelPart.GetProcessInfo();
            mDimension = r_process_info[DOMAIN_SIZE];
            int number_of_nodes = mpModelPart.NumberOfNodes();
            if (mModelPartName=="Structure"){
                mVectorDimension = mDimension*number_of_nodes;
            }
            else if (mModelPartName=="FluidModelPart"){
                mVectorDimension = (mDimension+1)*number_of_nodes;
            }
        }

        Vector GetSnapshotOfDisplacements(){
            Vector snapshot_vector =  ZeroVector(mVectorDimension);
            #pragma omp parallel for 
            for(auto& node : mpModelPart.Nodes()){
                Vector aux = ZeroVector(mDimension);
                aux = node.FastGetSolutionStepValue(DISPLACEMENT);
                int index_node = node.Id()-1;
                int index_dof = index_node*mDimension;
                for(int i=0;i<mDimension;i++){
                    snapshot_vector(index_dof+i) = aux[i]; 
                }
            }
            return snapshot_vector;
        }

        Vector GetSnapshotOfVelocityAndPressure(){
            Vector snapshot_vector =  ZeroVector(mVectorDimension);
            #pragma omp parallel for
            for(auto& node : mpModelPart.Nodes()){
                Vector aux = ZeroVector(mDimension);
                aux = node.FastGetSolutionStepValue(VELOCITY);
                int index_node = node.Id()-1;
                int index_dof = index_node*(mDimension+1);
                for(int i=0;i<mDimension;i++){
                    snapshot_vector(index_dof+i) = aux[i]; 
                }
                snapshot_vector(index_dof+mDimension) = node.FastGetSolutionStepValue(PRESSURE);
            }
            return snapshot_vector;
        }

        protected:
            ModelPart& mpModelPart;
            int mVectorDimension,mDimension;
            std::string mModelPartName;
    };



} // namespace Kratos



#endif // ROM_OUTPUT_UTILITY_H_INCLUDED  defined