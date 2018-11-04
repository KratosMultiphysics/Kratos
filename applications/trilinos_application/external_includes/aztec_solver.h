//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_AZTEC_SOLVER_H_INCLUDED )
#define  KRATOS_AZTEC_SOLVER_H_INCLUDED

// External includes
#include "string.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack.h"
#include "Ifpack_ConfigDefs.h"



namespace Kratos
{
enum AztecScalingType {NoScaling,LeftScaling,SymmetricScaling};

template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AztecSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of AztecSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION( AztecSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    AztecSolver(Parameters settings)
    {
        Parameters default_settings( R"(
        {
        "solver_type": "AztecSolver",
        "tolerance" : 1.0e-6,
        "max_iteration" : 200,
        "preconditioner_type" : "None",
        "overlap_level":1,
        "gmres_krylov_space_dimension" : 100,
        "scaling":false,
        "verbosity":0,
        "trilinos_aztec_parameter_list": {},
        "trilinos_preconditioner_parameter_list": {}
        }  )" );

        settings.ValidateAndAssignDefaults(default_settings);


        //settings for the AZTEC solver
        mtol = settings["tolerance"].GetDouble();
        mmax_iter = settings["max_iteration"].GetDouble();

        //IFpack settings
        moverlap_level = settings["overlap_level"].GetInt();

        //scaling settings
        if (settings["scaling"].GetBool() == false)
            mscaling_type = NoScaling;
        else
            mscaling_type = LeftScaling;

        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mparameter_list
        maztec_parameter_list = Teuchos::ParameterList();

        if(settings["verbosity"].GetInt() == 0)
            maztec_parameter_list.set("AZ_output", "AZ_none");
        else
            maztec_parameter_list.set("AZ_output", settings["verbosity"].GetInt());

        //choose the solver type
        if(settings["solver_type"].GetString() == "CGSolver")
        {
            maztec_parameter_list.set("AZ_solver", "AZ_cg");
        }
        else if(settings["solver_type"].GetString() == "BICGSTABSolver")
        {
            maztec_parameter_list.set("AZ_solver", "AZ_bicgstab");
        }
        else if(settings["solver_type"].GetString() == "GMRESSolver")
        {
            maztec_parameter_list.set("AZ_solver", "AZ_gmres");
            maztec_parameter_list.set("AZ_kspace", settings["gmres_krylov_space_dimension"].GetInt());
        }
        else if(settings["solver_type"].GetString() == "AztecSolver")
        {
            //do nothing here. Leave full control to the user through the "trilinos_aztec_parameter_list"
        }
        else
        {
            KRATOS_ERROR << " the solver type specified : " << settings["solver_type"].GetString()  << " is not supported";
        }


        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        for(auto it = settings["trilinos_aztec_parameter_list"].begin(); it != settings["trilinos_aztec_parameter_list"].end(); it++)
        {
            if(it->IsString()) maztec_parameter_list.set(it.name(), it->GetString());
            else if(it->IsInt()) maztec_parameter_list.set(it.name(), it->GetInt());
            else if(it->IsBool()) maztec_parameter_list.set(it.name(), it->GetBool());
            else if(it->IsDouble()) maztec_parameter_list.set(it.name(), it->GetDouble());
        }

        mpreconditioner_parameter_list = Teuchos::ParameterList();
        if(settings["preconditioner_type"].GetString() == "DiagonalPreconditioner")
        {
            mIFPreconditionerType = "None";
        }
        else if(settings["preconditioner_type"].GetString() == "ILU0")
        {
            mIFPreconditionerType = "ILU";
            mpreconditioner_parameter_list.set("fact: drop tolerance", 1e-9);
            mpreconditioner_parameter_list.set("fact: level-of-fill", 1);
        }
        else if(settings["preconditioner_type"].GetString() == "ILUT")
        {
            mIFPreconditionerType = "ILU";
            mpreconditioner_parameter_list.set("fact: drop tolerance", 1e-9);
            mpreconditioner_parameter_list.set("fact: level-of-fill", 10);
        }
        else if(settings["preconditioner_type"].GetString() == "ICC")
        {
            mIFPreconditionerType = "ICC";
            mpreconditioner_parameter_list.set("fact: drop tolerance", 1e-9);
            mpreconditioner_parameter_list.set("fact: level-of-fill", 10);
        }
        else if(settings["preconditioner_type"].GetString() == "AmesosPreconditioner")
        {
            mIFPreconditionerType = "Amesos";
            mpreconditioner_parameter_list.set("amesos: solver type", "Amesos_Klu");
        }
        else if(settings["preconditioner_type"].GetString() == "None")
        {
            mIFPreconditionerType = "AZ_none";
        }
        else
            KRATOS_ERROR << "wrong choice for preconditioner_type. Selction was :" << settings["preconditioner_type"].GetString() << " Available choices are: None,ILU0,ILUT,ICC,AmesosPreconditioner";

        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        for(auto it = settings["trilinos_preconditioner_parameter_list"].begin(); it != settings["trilinos_preconditioner_parameter_list"].end(); it++)
        {
            if(it->IsString()) mpreconditioner_parameter_list.set(it.name(), it->GetString());
            else if(it->IsInt()) mpreconditioner_parameter_list.set(it.name(), it->GetInt());
            else if(it->IsBool()) mpreconditioner_parameter_list.set(it.name(), it->GetBool());
            else if(it->IsDouble()) mpreconditioner_parameter_list.set(it.name(), it->GetDouble());
        }


    }

    /**
     * Default constructor
     */
    AztecSolver(Teuchos::ParameterList& aztec_parameter_list, std::string IFPreconditionerType, Teuchos::ParameterList& preconditioner_parameter_list, double tol, int nit_max, int overlap_level)
    {
        //settings for the AZTEC solver
        maztec_parameter_list = aztec_parameter_list;
        mtol = tol;
        mmax_iter = nit_max;

        //IFpack settings
        mIFPreconditionerType = IFPreconditionerType;
        mpreconditioner_parameter_list = preconditioner_parameter_list;
        moverlap_level = overlap_level;

        mscaling_type = LeftScaling;

        /*			if(overlap_level == 0)
        				KRATOS_THROW_ERROR(std::logic_error,"the overlap level for the Aztecsolver with IFPackshould be greater than 0","");*/
    }

    /**
     * Destructor
     */
    virtual ~AztecSolver() {}

    //function to set the scaling typedef
    void SetScalingType(AztecScalingType scaling_type)
    {
        mscaling_type = scaling_type;
    }

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        KRATOS_TRY
        rA.Comm().Barrier();

        Epetra_LinearProblem AztecProblem(&rA,&rX,&rB);


        //perform GS1 scaling if required
        if(mscaling_type == SymmetricScaling)
        {
            KRATOS_THROW_ERROR(std::logic_error,"somethign wrong with the scaling to be further teststed","")
            Epetra_Vector scaling_vect(rA.RowMap());
            rA.InvColSums(scaling_vect);

            int MyLength = scaling_vect.MyLength();
            for( int i=0 ; i<MyLength ; ++i ) scaling_vect[i] = sqrt(scaling_vect[i]);

            AztecProblem.LeftScale(scaling_vect);
            AztecProblem.RightScale(scaling_vect);
        }
        else if (mscaling_type == LeftScaling)
        {
            Epetra_Vector scaling_vect(rA.RowMap());
            rA.InvColSums(scaling_vect);

            AztecProblem.LeftScale(scaling_vect);

        }

        AztecOO aztec_solver(AztecProblem);
        aztec_solver.SetParameters(maztec_parameter_list);

        //here we verify if we want a preconditioner
        if( mIFPreconditionerType!=std::string("AZ_none") )
        {

            //ifpack preconditioner type
            Ifpack Factory;

            std::string PrecType = mIFPreconditionerType;
            Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &rA, moverlap_level);
            assert(Prec != 0);

            IFPACK_CHK_ERR(Prec->SetParameters(mpreconditioner_parameter_list));
            IFPACK_CHK_ERR(Prec->Initialize());
            IFPACK_CHK_ERR(Prec->Compute());

            // HERE WE SET THE IFPACK PRECONDITIONER
            aztec_solver.SetPrecOperator(Prec);

            //and ... here we solve
            aztec_solver.Iterate(mmax_iter,mtol);

            delete Prec;
        }
        else
        {
            aztec_solver.Iterate(mmax_iter,mtol);
        }

// 		for( int i=0 ; i<(rX).MyLength() ; ++i )
// 		{
// 		     (&rX)[i] = (&rX)[i]/scaling_vect[i] ;
// 		}


        rA.Comm().Barrier();

        return true;
        KRATOS_CATCH("");
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {

        return false;
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
//                rOStream << "Aztec solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    //aztec solver settings
    Teuchos::ParameterList maztec_parameter_list;
    double mtol;
    int mmax_iter;
    AztecScalingType mscaling_type;

    std::string mIFPreconditionerType;

    Teuchos::ParameterList mpreconditioner_parameter_list;
    int moverlap_level;



    /**
     * Assignment operator.
     */
    AztecSolver& operator=(const AztecSolver& Other);

    /**
     * Copy constructor.
     */
    AztecSolver(const AztecSolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, AztecSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AztecSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_AZTEC_SOLVER_H_INCLUDED  defined
