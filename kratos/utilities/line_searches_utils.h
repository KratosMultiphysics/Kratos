//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//
//


#if !defined(KRATOS_LINE_SEARCHES_UTILS)
#define  KRATOS_LINE_SEARCHES_UTILS


/* System includes */


/* External includes */



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

#include <cmath>

namespace Kratos
{
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class LineSearchesUtils

{

public:

    typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    LineSearchesUtils(
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TSchemeType::Pointer pScheme,
        unsigned  int MaxLineSearchIterations   = 20,
        double tolls  = 0.10,    // energy tolerance factor on LineSearch (0.8 is ok)
        double amp    = 1.618,   // maximum amplification factor
        double etmxa  = 5,       // maximum allowed step length
        double etmna  = 0.1,     // minimum allowed step length
        bool MoveMeshFlag       = true,
        bool ApplyLineSearches  = true
    )
        :rmodel_part(model_part)

    {
        KRATOS_TRY
        SetParametersLineSearches(MaxLineSearchIterations,
                                  tolls,
                                  amp,
                                  etmxa,
                                  etmna,
                                  MoveMeshFlag,
                                  ApplyLineSearches);

        //saving the scheme
        mpScheme = pScheme;

        //saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        //setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                 new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(mpLinearSolver)
                             );





        KRATOS_CATCH("")
    }


    /** Destructor.
    */
    virtual ~LineSearchesUtils () {}


protected:

    // parameters of line searches
    unsigned int  mMaxLineSearchIterations;
    double mtolls;
    double mamp;
    double meta;
    double metmxa;
    double metmna;
    bool   mMoveMeshFlag;
    bool   mApplyLineSearches;

    typename TSchemeType::Pointer mpScheme;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    typename TLinearSolver::Pointer mpLinearSolver;

    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    ModelPart& rmodel_part;

    // methods of line searches
    void SetParametersLineSearches(unsigned int& MaxLineSearchIterations,
                                   double& tolls,
                                   double& amp,
                                   double& etmxa,
                                   double& etmna,
                                   bool& MoveMeshFlag,
                                   bool& ApplyLineSearches)
    {
        KRATOS_TRY

        mMaxLineSearchIterations = MaxLineSearchIterations;
        mtolls  = mtolls;
        mamp    = mamp;
        metmxa  = etmxa;
        metmna  = etmna;
        mMoveMeshFlag      = MoveMeshFlag;
        mApplyLineSearches = ApplyLineSearches;
        meta= 1.00;


        KRATOS_CATCH("")
    }




    bool LineSearches( TSystemVectorType&  X_old,
                       TSystemVectorType&  Delta_p,
                       TSystemVectorType&  mDx,
                       TSystemVectorType&  mb,
                       TSystemMatrixType&  mA
                     )
    {

        KRATOS_TRY
        double seta               = 0.00;
        double so                 = 0.00;
        unsigned int ils          = 0;
        unsigned int ico          = 0;
        unsigned iteration_number = 1;
        Vector prodr(mMaxLineSearchIterations+1,false);
        Vector lseta(mMaxLineSearchIterations+2,false);
        lseta(1) =  1.00;
        prodr(0) =  1.00;
        meta     =  1.00;

        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();


        //TSystemMatrixType& mA        = *(mpA);
        //TSystemVectorType& mDx       = *(mpDx);
        //TSystemVectorType& mb        = *(mpb);
        //TSystemVectorType& X_old     = *(pX_old);
        //TSystemVectorType& Delta_p   = *(pDelta_p);

        TSparseSpace::SetToZero(mb);
        pBuilderAndSolver->BuildRHS(pScheme,rmodel_part,mb);
        so = TSparseSpace::Dot(Delta_p,mb); // mb inicial

        KRATOS_WATCH(so);

        while (iteration_number <= mMaxLineSearchIterations)
        {

            pScheme->Update(rmodel_part,rDofSet,mA,mDx,mb);
            //if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
            pScheme->FinalizeNonLinIteration(rmodel_part,mA,mDx,mb);

            std::cout<<"Line_Search_Iteration_Number:"<<(iteration_number)<<std::endl;
            TSparseSpace::SetToZero(mb);
            pBuilderAndSolver->BuildRHS(pScheme,rmodel_part,mb); //mb para calculat eta

            TSparseSpace::Copy(Delta_p, mDx);

            seta = TSparseSpace::Dot(Delta_p,mb); //Delta_p =  Const
            //KRATOS_WATCH(mb);
            //KRATOS_WATCH(Delta_p);
            KRATOS_WATCH(seta);

            // Posible restriccion en el denominador
            ils = iteration_number -1;
            prodr(ils+1) = (seta/so);
            KRATOS_WATCH(prodr(ils+1))
            if (fabs(prodr(ils+1)) < mtolls)
            {
                return true;
                break;
            }
            else
            {
                //Search devuelve lseta(ils+2)
                SetDatabaseToValue(rDofSet, X_old);
                //if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
                Search(ils,prodr,lseta,ico);
                std::cout<<"Eta:"<< lseta(ils+2)<<std::endl;
                meta = lseta(ils+2);
                ////KRATOS_WATCH(mDx);
                TSparseSpace::Assign(mDx,meta, mDx); //mDx =lseta(iteration_number+1)*mDx
                ////KRATOS_WATCH(mDx);
            }

            if (ico==2)
            {
                //ilfail = 2;
                std::cout<<"********************************"<<std::endl;
                std::cout<<"******Line Search Trouble*******"<<std::endl;
                std::cout<<"********************************"<<std::endl;
                return false;
                break;
            }

            if (iteration_number>= mMaxLineSearchIterations)
            {
                this->MaxIterationsExceeded();
                return false;
                break;
            }


            iteration_number++;
        }
        KRATOS_CATCH("")
    }


    void SetDatabaseToValue(
        DofsArrayType& rDofSet,
        const TSystemVectorType& X_old
    )
    {
        KRATOS_TRY

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() = X_old[i_dof->EquationId()];
            }
        }
        KRATOS_CATCH("")
    }


    // Permite escribir los desplazamientos antiguos en el X_old
    void BackupDatabase(
        DofsArrayType const& rDofSet,
        TSystemVectorType& X_old
    )
    {
        KRATOS_TRY

        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                X_old[i_dof->EquationId()] = i_dof->GetSolutionStepValue();
            }
        }
        KRATOS_CATCH("")
    }

//realiza Line-Search local para encontrar la longitud del paso en eta(ils+2)
//eta: historia de longitudes de pasos en la iteracion de L-S (eta(1)=0.,eta(2)=1.)
//prodr: productos r
//amp: factor de amplificacion maximo para el paso de carga en caso de extrapolacion
//etmxa y etmna: longitudes de pasos maximos y minimos permitidos
//ico: 1 si en el paso anterior ha sido necesario usar etmxa o etmna, sino es 0
// al salir ico es 1 si en esta iteracion se ha usado etmxa o etmna y 2 si se han
//usado por segunda vez.

//obtener ineg=No del S-L actual que este mas cerca del origen y que su prodr sea negativo
//tambien busca el S-L maximo en etmaxp
// si no hay existe un eta con prodr<0 ineg=999

    void Search(unsigned int& ils,Vector& prodr, Vector& eta, unsigned int& ico)
    {
        KRATOS_TRY

        double etaint      = 0.00;
        double etaneg      = 1e5;
        double etmxt       = 0.00;
        double etmaxp      = 0.00;
        double etaalt      = 0.00;
        double etaext      = 0.00;
        unsigned int ipos  = 0;
        unsigned int ineg  = 999;

        for (unsigned int i=0; i<=ils+1; i++)
        {
            if (eta(i)>etmaxp)
            {
                etmaxp=eta(i);
            }
            if (prodr(i)<0 && eta(i)<=etaneg)
            {
                etaneg=eta(i);
                ineg=i;
            }
        }

        //Decide si interpola o extrapola
        if (ineg !=999)
        {
            //Interpola
            //entontrar ipos=No del S-L con prodr positivo mas cercano a ineg pero con S-L mas peque�o
            ipos = 0;
            for( unsigned int i=0; i<=ils+1; i++)
            {
                if (prodr(i)>=0 && eta(i)<=eta(ineg) && eta(i)>=eta(ipos))
                {
                    ipos=i;
                }
            }
            //interpola para encontrar S-L etaint
            etaint=prodr(ineg)*eta(ipos)-prodr(ipos)*eta(ineg);
            if ((prodr(ineg)-prodr(ipos))==0)
            {
                ico=2;
                return;
            }
            else
            {
                etaint=etaint/(prodr(ineg)-prodr(ipos));
            }
            //alternativamente encuentra etaalt para asegurar un cambio razonable
            etaalt = eta(ipos) + 0.2*(eta(ineg)-eta(ipos));
            //toma el maximo entre etaint y etaalt
            etaint=std::max(etaint,etaalt);
            //la longitud de paso minimo
            if (etaint<metmna)
            {
                etaint=metmna;
                if (ico==1)
                {
                    ico=2;
                }
                else if (ico==0)
                {
                    ico=1;
                }
            }

            eta(ils+2)=etaint;
        }

        else if (ineg==999)
        {

            //extrapola
            //asigna temporalmente la longitud de paso maxima
            etmxt = mamp*etmaxp;
            if(etmxt>metmxa)
            {
                etmxt=metmxa;
            }
            //extrapola la longitud de paso entre el actual y el anterior
            etaext=prodr(ils+1)*eta(ils)-prodr(ils)*eta(ils+1);
            if ((prodr(ils+1)-prodr(ils))==0)
            {
                ico=2;
                return;
            }
            else
            {
                etaext=etaext/(prodr(ils+1)-prodr(ils));
            }

            eta(ils+2)=etaext;
            //se acepta eta si esta dentro de los limites
            if (etaext<=0 || etaext> etmxt)
            {
                eta(ils+2) = etmxt;
            }
            if (eta(ils+2)==metmxa && ico==1)
            {
                ico=2;
                return;
            }
            if (eta(ils+2)==metmxa)
            {
                ico=1;
            }
        }
        KRATOS_CATCH("")
    }

    void MaxIterationsExceeded()
    {
        std::cout << "*****************************************************************" << std::endl;
        std::cout << "******* ATTENTION: Max Iterations Line Searches Exceeded ********" << std::endl;
        std::cout << "*****************************************************************" << std::endl;
    }


    void MoveMesh()
    {
        KRATOS_TRY

        for(ModelPart::NodeIterator i = rmodel_part.NodesBegin() ;
                i != rmodel_part.NodesEnd() ; ++i)
        {
            array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
            (i)->X() = (i)->X0() + disp[0];
            (i)->Y() = (i)->Y0() + disp[1];
            (i)->Z() = (i)->Z0() + disp[2];
        }
        KRATOS_CATCH("")
    }

}; // end class lineseaches
} // namespace Kratos
#endif









