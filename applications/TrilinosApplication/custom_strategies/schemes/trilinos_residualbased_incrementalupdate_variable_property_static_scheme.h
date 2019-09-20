//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_TRILINOS_STATIC_SCHEME_VARIABLE_PROPERTY )
#define  KRATOS_TRILINOS_STATIC_SCHEME_VARIABLE_PROPERTY


/* System includes */

/* External includes */
#include "Epetra_Import.h"


/* Project includes */
#include "includes/c2c_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

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
class TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
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
        ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>()
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
    ) override
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
    ) override
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
    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY;

        BaseType::Update(r_model_part,rDofSet,A,Dx,b);

        AddiotionalTablePropertyUpdate(r_model_part);

        KRATOS_CATCH("")
    }




    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


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

    /*@} */
    /**@name Protected  Access */
    /*@{ */


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
