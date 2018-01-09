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

/* *********************************************************
*
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-11-11 14:03:41 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_TRILINOS_STATIC_SCHEME_VARIABLE_PROPERTY )
#define  KRATOS_TRILINOS_STATIC_SCHEME_VARIABLE_PROPERTY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"
#include "Epetra_Import.h"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
// #include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/c2c_variables.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.

This class provides the implementation of the basic tasks that are needed by the solution strategy.
It is intended to be the place for tailoring the solution strategies to problem specific tasks.

Detail class definition.

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


*/
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme : public TrilinosResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme);

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme():
        TrilinosResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>(),
        mImporterIsInitialized(false)
    {}

    /** Destructor.
    */
    virtual ~TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */
//     /**
//    Performing the Initialize of the solution.
//    */
//    void Initialize(
//    ModelPart& r_model_part
//    ) {
//        //Initialize variables
//          ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
//          ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
//          const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
//          const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
//          const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
//          const Variable<double>& rTransferCoef = my_settings->GetTransferCoefficientVariable();
//          // const double amb_temp = CurrentProcessInfo[AMBIENT_TEMPERATURE];


//        ModelPart::TableType rDensityVar_table = r_model_part.GetTable(1);
//        ModelPart::TableType C_table = r_model_part.GetTable(2);
//        ModelPart::TableType F_table = r_model_part.GetTable(3);
//        ModelPart::TableType DF_DT_table = r_model_part.GetTable(4);
//        ModelPart::TableType rDiffusionVar_table = r_model_part.GetTable(5);
//        ModelPart::TableType HTC_table = r_model_part.GetTable(7);

//        double density_var = r_model_part.GetProcessInfo()[DENSITY];


//// 	      for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
//          #pragma omp parallel for
//          for (int k = 0; k< static_cast<int> (r_model_part.Nodes().size()); k++)
//          {
//        ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin() + k;
//              const double unknown_val = ind->FastGetSolutionStepValue(rUnknownVar);
//              const double dist = ind->FastGetSolutionStepValue(DISTANCE);

//              /*double htc_var = ind->FastGetSolutionStepValue(rTransferCoef);
//              double rho = rDensityVar_table.GetValue(unknown_val);
//              double cc =C_table.GetValue(unknown_val);
//              ind->FastGetSolutionStepValue(rTransferCoef) = htc_var/(rho*cc);	*/
//              double htc_var = HTC_table.GetValue(unknown_val);
//              ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;

//            if(dist < 0){
//             // double density_var = rDensityVar_table.GetValue(unknown_val);
//              double specific_heat_var =C_table.GetValue(unknown_val);
//              double solid_fraction_var = F_table.GetValue(unknown_val);
//              double solid_fraction_rate_var = DF_DT_table.GetValue(unknown_val);
//              double conductvity_var = rDiffusionVar_table.GetValue(unknown_val);
//              //double htc_var = HTC_table.GetValue(unknown_val);


//              ind->FastGetSolutionStepValue(rDensityVar) = density_var;
//              ind->FastGetSolutionStepValue(rDensityVar,1) = density_var;

//              ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_var;
//              ind->FastGetSolutionStepValue(SPECIFIC_HEAT,1) = specific_heat_var;

//              ind->FastGetSolutionStepValue(SOLIDFRACTION) = solid_fraction_var;
//              ind->FastGetSolutionStepValue(SOLIDFRACTION,1) = solid_fraction_var;

//              ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = solid_fraction_rate_var;
//              ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE,1) = solid_fraction_rate_var;

//              ind->FastGetSolutionStepValue(rDiffusionVar) = conductvity_var;
//              ind->FastGetSolutionStepValue(rDiffusionVar,1) = conductvity_var;

//             // ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;

//            }
//            else
//            {
//              ind->FastGetSolutionStepValue(rDensityVar) = 1.0;
//              ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = 1000.0;
//              ind->FastGetSolutionStepValue(SOLIDFRACTION) = 1.0;
//              ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = 0.0;
//              ind->FastGetSolutionStepValue(rDiffusionVar) = 1.0;

//              ind->FastGetSolutionStepValue(rTransferCoef) = 1.0;
//             // ind->FastGetSolutionStepValue(rUnknownVar) = amb_temp;
//            }

//          }
////               mSchemeIsInitialized = true;
//      }
    void Initialize(
        ModelPart& r_model_part
    )
    {
        //Initialize variables
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
        const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
        // const Variable<double>& rTransferCoef = my_settings->GetTransferCoefficientVariable();
        // const double amb_temp = CurrentProcessInfo[AMBIENT_TEMPERATURE];


        ModelPart::TableType rDensityVar_table = r_model_part.GetTable(1);
        ModelPart::TableType C_table = r_model_part.GetTable(2);
        ModelPart::TableType F_table = r_model_part.GetTable(3);
        ModelPart::TableType DF_DT_table = r_model_part.GetTable(4);
        ModelPart::TableType rDiffusionVar_table = r_model_part.GetTable(5);
        //ModelPart::TableType HTC_table = r_model_part.GetTable(7);

        double density_var = r_model_part.GetProcessInfo()[DENSITY];
        const double latent_heat = r_model_part.GetProcessInfo()[LATENT_HEAT];
        const unsigned int buffer_size = r_model_part.GetBufferSize();

// 	      for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (r_model_part.Nodes().size()); k++)
        {
            ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin() + k;
            const double unknown_val = ind->FastGetSolutionStepValue(rUnknownVar);
            const double dist = ind->FastGetSolutionStepValue(DISTANCE);

            /*double htc_var = ind->FastGetSolutionStepValue(rTransferCoef);
            double rho = rDensityVar_table.GetValue(unknown_val);
            double cc =C_table.GetValue(unknown_val);
            ind->FastGetSolutionStepValue(rTransferCoef) = htc_var/(rho*cc);	*/
            double specific_heat_var =C_table.GetValue(unknown_val);
            //double htc_var = HTC_table.GetValue(unknown_val);
            //ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;
            double conductivity_var = rDiffusionVar_table.GetValue(unknown_val);

             if(dist <= 0)
             {
                // double density_var = rDensityVar_table.GetValue(unknown_val);
                // double specific_heat_var =C_table.GetValue(unknown_val);
                double solid_fraction_var = F_table.GetValue(unknown_val);
                double solid_fraction_rate_var = DF_DT_table.GetValue(unknown_val);

                //double htc_var = HTC_table.GetValue(unknown_val);


                ind->FastGetSolutionStepValue(rDensityVar) = density_var;
                ind->FastGetSolutionStepValue(rDensityVar,1) = density_var;

                ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_var;
                ind->FastGetSolutionStepValue(SPECIFIC_HEAT,1) = specific_heat_var;

                ind->FastGetSolutionStepValue(SOLIDFRACTION) = solid_fraction_var;
                ind->FastGetSolutionStepValue(SOLIDFRACTION,1) = solid_fraction_var;

                ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = solid_fraction_rate_var;
                ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE,1) = solid_fraction_rate_var;

                ind->FastGetSolutionStepValue(rDiffusionVar) = conductivity_var;
                ind->FastGetSolutionStepValue(rDiffusionVar,1) = conductivity_var;

                //ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;


                //assign an initial value to the enthalpy
                const double initial_enthalpy = specific_heat_var*unknown_val + (1.0-ind->FastGetSolutionStepValue(SOLIDFRACTION))*latent_heat;
                for(unsigned int i=0; i<buffer_size; i++)
                    ind->FastGetSolutionStepValue(ENTHALPY,i) = initial_enthalpy;


             }
             else
             {
                 const double specific_heat_air = 1000.0;
                 ind->FastGetSolutionStepValue(rDensityVar) = 1.0;
                 ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_air;
                 ind->FastGetSolutionStepValue(SOLIDFRACTION) = 0.0;
                 ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = 0.0;
                 ind->FastGetSolutionStepValue(rDiffusionVar) = conductivity_var*1000.0; //0.05

                 //ind->FastGetSolutionStepValue(rTransferCoef) = 0.0; ///(density_var*specific_heat_var);
                 // ind->FastGetSolutionStepValue(rUnknownVar) = amb_temp;

                 //assign an initial value to the enthalpy
                 for(unsigned int i=0; i<buffer_size; i++)
                     ind->FastGetSolutionStepValue(ENTHALPY,i) = specific_heat_air*unknown_val;
              }

        }
//               mSchemeIsInitialized = true;
    }



    //***************************************************************************
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);

        //update variables based on resolved unknowns
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
        const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
        //const Variable<double>& rTransferCoef = my_settings->GetTransferCoefficientVariable();
        //const double amb_temp = CurrentProcessInfo[AMBIENT_TEMPERATURE];

        ModelPart::TableType rDensityVar_table = r_model_part.GetTable(1);
        ModelPart::TableType C_table = r_model_part.GetTable(2);
        ModelPart::TableType F_table = r_model_part.GetTable(3);
        ModelPart::TableType DF_DT_table = r_model_part.GetTable(4);
        ModelPart::TableType rDiffusionVar_table = r_model_part.GetTable(5);
        ModelPart::TableType HTC_table = r_model_part.GetTable(7);

        double density_var = r_model_part.GetProcessInfo()[DENSITY];

        //const double latent_heat = r_model_part.GetProcessInfo()[LATENT_HEAT];

        //for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (r_model_part.Nodes().size()); k++)
        {
            ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin() + k;
            const double unknown_val = ind->FastGetSolutionStepValue(rUnknownVar);
            const double dist = ind->FastGetSolutionStepValue(DISTANCE);

            double specific_heat_var =C_table.GetValue(unknown_val);
            //double htc_var = HTC_table.GetValue(unknown_val);
            //ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;

             if(dist <= 0)
             {
                //double density_var = rDensityVar_table.GetValue(unknown_val);
                double solid_fraction_var = F_table.GetValue(unknown_val);
                double solid_fraction_rate_var = DF_DT_table.GetValue(unknown_val);
                double conductvity_var = rDiffusionVar_table.GetValue(unknown_val);

                ind->FastGetSolutionStepValue(rDensityVar) = density_var;
                ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_var;
                ind->FastGetSolutionStepValue(SOLIDFRACTION) = solid_fraction_var;
                ind->GetValue(SOLIDFRACTION) = solid_fraction_var; //also save in database without history
                ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = solid_fraction_rate_var;
                ind->FastGetSolutionStepValue(rDiffusionVar) = conductvity_var;
                //ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;

                //here compute the ENTHALPY = int(c dT,0,T) + L(T)
                //which should be computed incrementally as   Hn+1 = Hn + 1/2*(Tn+1 - Tn)*(cn+1 - cn) + L(Tn+1) - L(Tn)
//                const double Delta_T = unknown_val - ind->GetValue(rUnknownVar);
//                const double delta_solid_fraction = 0.0; //(1-Sn+1) - (1-Sn)
                const double delta_enthalpy = 0.0; //Delta_T*specific_heat_var + delta_solid_fraction*latent_heat;
                ind->FastGetSolutionStepValue(ENTHALPY) = /*ind->FastGetSolutionStepValue(ENTHALPY,1) +*/ delta_enthalpy;
                ind->GetValue(ENTHALPY) = ind->FastGetSolutionStepValue(ENTHALPY);
                ind->GetValue(SOLIDFRACTION) = solid_fraction_var;
                ind->GetValue(SOLIDFRACTION_RATE) = solid_fraction_rate_var;
            }
            else
            {
                const double conductivity_var = rDiffusionVar_table.GetValue(unknown_val);
                const double specific_heat_air = 1000.0;
                ind->FastGetSolutionStepValue(rDensityVar) = 1.0; //density_var;  //1.0;
                ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_air;
                ind->FastGetSolutionStepValue(SOLIDFRACTION) = 0.0;
                ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = 0.0;
                ind->FastGetSolutionStepValue(rDiffusionVar) = conductivity_var*1000.0;  //0.05
                //ind->FastGetSolutionStepValue(rTransferCoef) = 0.0; //htc_var/(50.0); //*density_var*specific_heat_air);

                ind->FastGetSolutionStepValue(ENTHALPY) = specific_heat_air*unknown_val; // * (ind->FastGetSolutionStepValue(rUnknownVar)) ;
                ind->GetValue(ENTHALPY) = specific_heat_air*unknown_val;
                ind->GetValue(SOLIDFRACTION) = 0.0;
                ind->GetValue(SOLIDFRACTION_RATE) = 0.0;
            }

        }

        KRATOS_CATCH("")
    }

    /**
    Performing the update of the solution.
    */
    //***************************************************************************
    virtual void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY;

        if (!DofImporterIsInitialized())
            this->InitializeDofImporter(rDofSet,Dx);

        int system_size = TSparseSpace::Size1(A);

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp( mpDofImporter->TargetMap() );

        //importing in the new temp vector the values
        int ierr = temp.Import(Dx,*mpDofImporter,Insert);
        if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");

        double* temp_values; //DO NOT make delete of this one!!
        temp.ExtractView( &temp_values );

        Dx.Comm().Barrier();


// ModelPart::NodesContainerType::iterator node_it = r_model_part.Nodes().find(2756);
// std::cout << A.Comm().MyPID() << " node 2756 " << node_it->FastGetSolutionStepValue(PARTITION_INDEX) << " disp_x id " << node_it->pGetDof(DISPLACEMENT_X)->EquationId() << std::endl;

// std::cout << "rank=" << A.Comm().MyPID() << "dof with id 117 "<< rDofSet.find(117) << std::endl;

        //performing the update
        typename DofsArrayType::iterator dof_begin = rDofSet.begin();
        for(unsigned int iii=0; iii<rDofSet.size(); iii++)
        {
            int global_id = (dof_begin+iii)->EquationId();
            if(global_id < system_size)
            {
                double aaa = temp[mpDofImporter->TargetMap().LID(global_id)];
                /*		if(global_id == 117) std::cout << "rank = " << b.Comm().MyPID() << " global id" << global_id << "local num" << iii << " map num " << dof_update_map.LID(global_id) << " value= " << aaa << std::endl;*/
                (dof_begin+iii)->GetSolutionStepValue() += aaa;
            }
        }

//                        delete [] temp_values;  //deleting this is WRONG! do not do it!!


       AddiotionalTablePropertyUpdate(r_model_part);

        KRATOS_CATCH("")
    }




    /*@} */
    /**@name Operations */
    /*@{ */


    virtual void Clear()
    {
        mpDofImporter.reset();
        mImporterIsInitialized = false;
    }

    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    bool DofImporterIsInitialized()
    {
        return mImporterIsInitialized;
    }

    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    virtual void InitializeDofImporter(DofsArrayType& rDofSet,
                                       TSystemVectorType& Dx)
    {
        int system_size = TSparseSpace::Size(Dx);
        int number_of_dofs = rDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        //filling the array with the global ids
        int counter = 0;
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            int id = i_dof->EquationId();
            if( id < system_size )
            {
                index_array[counter] = id;
                counter += 1;
            }
        }

        std::sort(index_array.begin(),index_array.end());
        std::vector<int>::iterator NewEnd = std::unique(index_array.begin(),index_array.end());
        index_array.resize(NewEnd-index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        Dx.Comm().SumAll(&tot_update_dofs,&check_size,1);
        if ( (check_size < system_size) &&  (Dx.Comm().MyPID() == 0) )
        {
            std::stringstream Msg;
            Msg << "Dof count is not correct. There are less dofs then expected." << std::endl;
            Msg << "Expected number of active dofs = " << system_size << " dofs found = " << check_size << std::endl;
            KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"")
        }

        //defining a map as needed
        Epetra_Map dof_update_map(-1,index_array.size(), &(*(index_array.begin())),0,Dx.Comm() );

        //defining the importer class
        boost::shared_ptr<Epetra_Import> pDofImporter( new Epetra_Import(dof_update_map,Dx.Map()) );
        mpDofImporter.swap(pDofImporter);

        mImporterIsInitialized = true;
    }


    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /// Get pointer Epetra_Import instance that can be used to import values from Dx to the owner of each Dof.
    /**
     * @note Important: always check that the Importer is initialized before calling using
     * DofImporterIsInitialized or initialize it with InitializeDofImporter.
     * @return Importer
     */
    boost::shared_ptr<Epetra_Import> pGetImporter()
    {
        return mpDofImporter;
    }

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    bool mImporterIsInitialized;

    boost::shared_ptr<Epetra_Import> mpDofImporter;

    /*@} */
    /**@name Private Operators*/
    /*@{ */
   void  AddiotionalTablePropertyUpdate(ModelPart& r_model_part){
       KRATOS_TRY;

       //update variables based on resolved unknowns
       ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
       ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
       const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
       //const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
       //const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
       //const Variable<double>& rTransferCoef = my_settings->GetTransferCoefficientVariable();
       //const double amb_temp = CurrentProcessInfo[AMBIENT_TEMPERATURE];

       ModelPart::TableType rDensityVar_table = r_model_part.GetTable(1);
       ModelPart::TableType C_table = r_model_part.GetTable(2);
       ModelPart::TableType F_table = r_model_part.GetTable(3);
       ModelPart::TableType DF_DT_table = r_model_part.GetTable(4);
       ModelPart::TableType rDiffusionVar_table = r_model_part.GetTable(5);
//         ModelPart::TableType HTC_table = r_model_part.GetTable(7);

       //double density_var = r_model_part.GetProcessInfo()[DENSITY];

       const double latent_heat = r_model_part.GetProcessInfo()[LATENT_HEAT];


       //for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
       #pragma omp parallel for
       for (int k = 0; k< static_cast<int> (r_model_part.Nodes().size()); k++)
       {
           ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin() + k;
           const double unknown_val = ind->FastGetSolutionStepValue(rUnknownVar);
           const double dist = ind->FastGetSolutionStepValue(DISTANCE);

//            double specific_heat_var =C_table.GetValue(unknown_val);
//             double htc_var = HTC_table.GetValue(unknown_val);
//             ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;

            if(dist <= 0)
            {
               //double density_var = rDensityVar_table.GetValue(unknown_val);
               double solid_fraction_var = F_table.GetValue(unknown_val);
               double solid_fraction_rate_var = DF_DT_table.GetValue(unknown_val);
//                double conductvity_var = rDiffusionVar_table.GetValue(unknown_val);

//                ind->FastGetSolutionStepValue(rDensityVar) = density_var;
//                ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_var;
               ind->FastGetSolutionStepValue(SOLIDFRACTION) = solid_fraction_var;
               ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = solid_fraction_rate_var;
//                ind->FastGetSolutionStepValue(rDiffusionVar) = conductvity_var;
//                 ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;

               //here compute the ENTHALPY = int(c dT,0,T) + L(T)
               //which should be computed incrementally as   Hn+1 = Hn + 1/2*(Tn+1 - Tn)*(cn+1 - cn) + L(Tn+1) - L(Tn)
               const double Delta_T = unknown_val - ind->GetValue(rUnknownVar);
               const double avg_c = ind->FastGetSolutionStepValue(SPECIFIC_HEAT);
               const double delta_solid_fraction = ind->GetValue(SOLIDFRACTION) - ind->FastGetSolutionStepValue(SOLIDFRACTION); //(1-Sn+1) - (1-Sn)
               const double delta_enthalpy = Delta_T*avg_c + delta_solid_fraction*latent_heat;
               ind->FastGetSolutionStepValue(ENTHALPY) = /*ind->GetValue(ENTHALPY) +*/ delta_enthalpy;
           }
           else
           {
               const double specific_heat_air = ind->FastGetSolutionStepValue(SPECIFIC_HEAT); //1000.0;
//                 ind->FastGetSolutionStepValue(rDensityVar) = 1.0;
//                 ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_air;
//                 ind->FastGetSolutionStepValue(SOLIDFRACTION) = 0.0;
//                 ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = 0.0;
//                 ind->FastGetSolutionStepValue(rDiffusionVar) = 1.0;
//                 ind->FastGetSolutionStepValue(rTransferCoef) = ind->FastGetSolutionStepValue(rTransferCoef)/(density_var*specific_heat_air);

               ind->FastGetSolutionStepValue(ENTHALPY) = specific_heat_air*unknown_val; // * (ind->FastGetSolutionStepValue(rUnknownVar)) ;

           }

       }

       KRATOS_CATCH("")

//       //update variables based on resolved unknowns
//        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
//        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
//        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
//        const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
//        const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
//        const Variable<double>& rTransferCoef = my_settings->GetTransferCoefficientVariable();
//        //const double amb_temp = CurrentProcessInfo[AMBIENT_TEMPERATURE];

//       ModelPart::TableType rDensityVar_table = r_model_part.GetTable(1);
//       ModelPart::TableType C_table = r_model_part.GetTable(2);
//       ModelPart::TableType F_table = r_model_part.GetTable(3);
//       ModelPart::TableType DF_DT_table = r_model_part.GetTable(4);
//       ModelPart::TableType rDiffusionVar_table = r_model_part.GetTable(5);
//       ModelPart::TableType HTC_table = r_model_part.GetTable(7);

//       double density_var = r_model_part.GetProcessInfo()[DENSITY];

//    //for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
//     #pragma omp parallel for
//    for (int k = 0; k< static_cast<int> (r_model_part.Nodes().size()); k++)
//    {
//        ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin() + k;
//        const double unknown_val = ind->FastGetSolutionStepValue(rUnknownVar);
//        const double dist = ind->FastGetSolutionStepValue(DISTANCE);

//        double htc_var = HTC_table.GetValue(unknown_val);
//        ind->FastGetSolutionStepValue(rTransferCoef) = htc_var;

//        if(dist < 0){
//        //double density_var = rDensityVar_table.GetValue(unknown_val);
//        double specific_heat_var =C_table.GetValue(unknown_val);
//        double solid_fraction_var = F_table.GetValue(unknown_val);
//        double solid_fraction_rate_var = DF_DT_table.GetValue(unknown_val);
//        double conductvity_var = rDiffusionVar_table.GetValue(unknown_val);

//        ind->FastGetSolutionStepValue(rDensityVar) = density_var;
//        ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = specific_heat_var;
//        ind->FastGetSolutionStepValue(SOLIDFRACTION) = solid_fraction_var;
//        ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = solid_fraction_rate_var;
//        ind->FastGetSolutionStepValue(rDiffusionVar) = conductvity_var;
//        }
//        else
//        {
//        ind->FastGetSolutionStepValue(rDensityVar) = 1.0;
//        ind->FastGetSolutionStepValue(SPECIFIC_HEAT) = 1000.0;
//        ind->FastGetSolutionStepValue(SOLIDFRACTION) = 1.0;
//        ind->FastGetSolutionStepValue(SOLIDFRACTION_RATE) = 0.0;
//        ind->FastGetSolutionStepValue(rDiffusionVar) = 1.0;
//        ind->FastGetSolutionStepValue(rTransferCoef) = 1.0;
//        //ind->FastGetSolutionStepValue(rUnknownVar) = amb_temp;
//        }

//    }
   }

    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_STANDARD_STATIC_SCHEME  defined */
