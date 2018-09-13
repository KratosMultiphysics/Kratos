// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_LINE_SEARCHES_UTILITY)
#define  KRATOS_LINE_SEARCHES_UTILITY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

//#include <cmath>

namespace Kratos
{
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class LineSearchesUtility

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

    /** Counted pointer of ClassName */
    /*KRATOS_CLASS_POINTER_DEFINITION(LineSearchesUtility);

    typedef Dof<TDataType> TDofType;
    typedef PointerVectorSet<TDofType, IdentityFunction<TDofType> > DofsArrayType;
    typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
    typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;
    */

    /************************************* CONSTRUCTOR *********************************/
    /***********************************************************************************/

    LineSearchesUtility(
        //typename TSchemeType::Pointer pScheme,
        //typename TLinearSolver::Pointer pNewLinearSolver,
        unsigned  int& MaxLineSearchIterations,
        double& tolls,         // Energy tolerance factor on LineSearch (0.8 is ok)
        double& amp,           // Maximum amplification factor
        double& etmxa,         // Maximum allowed step length
        double& etmna,         // Minimum allowed step length
        bool& MoveMeshFlag,
        bool& ApplyLineSearches
    )
    //:rmodel_part(model_part)

    {
        KRATOS_TRY;
        SetParametersLineSearches(MaxLineSearchIterations,
                                  tolls,
                                  amp,
                                  etmxa,
                                  etmna,
                                  MoveMeshFlag,
                                  ApplyLineSearches
                                 );

        KRATOS_CATCH("");
    }

    /************************************* DESTRUCTOR **********************************/
    /***********************************************************************************/

    virtual ~LineSearchesUtility () {}

    /***********************************************************************************/
    /***********************************************************************************/

protected:

    // Parameters of line searches
    unsigned int  mMaxLineSearchIterations;
    double mtolls;
    double mamp;
    double meta;
    double metmxa;
    double metmna;
    bool   mMoveMeshFlag;
    bool   mApplyLineSearches;

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Methods of line searches. It set the parameters of the line search
    * @param MaxLineSearchIterations: Maximum number of line searchs
    * @param tolls:
    * @param amp:
    * @param etmxa:
    * @param etmna:
    * @param MoveMeshFlag: Boolean ("flag") that determines if the mesh is moved
    * @param ApplyLineSearches: Boolean that determines if the line search is applied
    */

    void SetParametersLineSearches(
            unsigned int& MaxLineSearchIterations,
            double& tolls,
            double& amp,
            double& etmxa,
            double& etmna,
            bool& MoveMeshFlag,
            bool& ApplyLineSearches
            )
    {
        KRATOS_TRY;

        mMaxLineSearchIterations = MaxLineSearchIterations;
        mtolls  = tolls;
        mamp    = amp;
        metmxa  = etmxa;
        metmna  = etmna;
        mMoveMeshFlag      = MoveMeshFlag;
        mApplyLineSearches = ApplyLineSearches;
        meta = 1.00;

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It computes the line searches
    * @param rmodel_part: Model part of the problem
    * @param pScheme: The scheme used to solve the problem
    * @param pBuilderAndSolver: The Builder and Solver used to solve the problem
    * @param rDofSet: DOF that are set (BC, free or fix)
    * @param X_old: The previous value of the "variables"
    * @param Delta_p:
    * @param mDx: The incremenent of "displacements" (variables of the problem)
    * @param mb: The RHS of the system
    * @param mA: The LHS of the system
    */

    bool LineSearches(
            ModelPart& rmodel_part,
            typename TSchemeType::Pointer& pScheme,
            typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
            DofsArrayType& rDofSet,
            TSystemVectorType&  X_old,
            TSystemVectorType&  Delta_p,
            TSystemVectorType&  mDx,
            TSystemVectorType&  mb,
            TSystemMatrixType&  mA
            )
    {

        KRATOS_TRY;
        double seta               = 0.00;
        double so                 = 0.00;
        unsigned int ils          = 0;
        unsigned int ico          = 0;
        unsigned int iteration_number = 1;
        TSystemVectorType prodr(mMaxLineSearchIterations+1);
        TSystemVectorType lseta(mMaxLineSearchIterations+2);
        TSparseSpace::SetToZero(prodr);
        TSparseSpace::SetToZero(lseta);

        lseta(0) =  0.00;
        lseta(1) =  1.00;
        prodr(0) =  1.00;
        meta     = lseta(1);

        so = TSparseSpace::Dot(Delta_p,mb); // Initial mb
        //if (so < 0.00)
        {
            while (iteration_number <= mMaxLineSearchIterations)
            {
                rDofSet = pBuilderAndSolver->GetDofSet();
                pScheme->Update(rmodel_part,rDofSet,mA,mDx,mb);
                BackupDatabase(rDofSet,mDx);
                //KRATOS_WATCH(mDx)
                if(mMoveMeshFlag == true) MoveMesh(rmodel_part);
                TSparseSpace::SetToZero(mb);
                pBuilderAndSolver->BuildRHS(pScheme,rmodel_part,mb); //mb para calculat eta
                TSparseSpace::Copy(Delta_p, mDx);
                seta = TSparseSpace::Dot(mDx,mb); //Delta_p =  Const

                ils = iteration_number - 1;
                prodr(ils+1) = (seta/so);
                std::cout<< "Line_Search_Iteration:" << iteration_number << " Reaches Toler. = " <<  prodr(ils + 1) << "  Required Toler. = " << mtolls<<std::endl;

                if (fabs(prodr(ils+1)) < mtolls)
                {
                    return true;
                }
                else
                {
                    // Search gives lseta(ils+2)
                    SetDatabaseToValue(rDofSet, X_old);

                    if(mMoveMeshFlag == true)
                    {
                        MoveMesh(rmodel_part);
                    }

                    Search(ils, prodr, lseta, ico);

                    if (ico==2)
                    {
                        //ilfail = 2;
                        meta = 1.00;
                        std::cout<<"****** Line Search Trouble. No Convergence Was Reached *******"<<std::endl;
                        return false;
                    }

                    meta = lseta(ils+2);
                    TSparseSpace::Assign(mDx,meta, mDx); //mDx =lseta(iteration_number+1)*mDx
                }

                if (iteration_number >= mMaxLineSearchIterations)
                {
                    std::cout << "*****************************************************************" << std::endl;
                    std::cout << "******* ATTENTION: Max Iterations Line Searches Exceeded ********" << std::endl;
                    std::cout << "*****************************************************************" << std::endl;
                    return false;
                    //break;
                }

                iteration_number++;
            }
        }

        return true;
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It set the values in the database in the values that are free
    * @param rDofSet: DOF that are set (BC, free or fix)
    * @param X_old: The previous value of the "variables"
    */

    void SetDatabaseToValue(
            DofsArrayType& rDofSet,
            const TSystemVectorType& X_old
            )
    {
        KRATOS_TRY;

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() = X_old[i_dof->EquationId()];
            }
        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It allows to write the old displacements in X_old
    * @param rDofSet: DOF that are set (BC, free or fix)
    * @param X_old: The previous value of the "variables"
    */

    void BackupDatabase(
            DofsArrayType const& rDofSet,
            TSystemVectorType& X_old
            )
    {
        KRATOS_TRY;

        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                X_old[i_dof->EquationId()] = i_dof->GetSolutionStepValue();
            }
        }
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It realises a local Line-Search to look for the length in the step eta(ils+2)
    * eta: History of the length of the steps in the iteration of the S-L(eta(1)=0.,eta(2)=1.)
    * amp: Amplification factor for the load step in case of  extrapolation
    * etmxa and etmna: mMaximal and minimal step lengths allowed
    * Obtains ineg=false of the current S-L which is closer to the origin and with a negative prodr.
    * It looks for the S-L maximal in etmaxp. There is a eta with prodr < 0 and ineg = 999
    * @param ils
    * @param prodr: r products
    * @param lseta:
    * @param ico: 1 when in the previous iteration etmxa or etmna have been employed, otherwise is 0.
    * Once ended if etmxa or etmna have been employed the value of ico es 1 if it is the first time, and 2 if it is the second time
    */

    void Search(
            unsigned int& ils,
            TSystemVectorType& prodr,
            TSystemVectorType& lseta,
            unsigned int& ico
            )
    {
        KRATOS_TRY;

        double etaint      = 0.00;
        double etaneg      = 1e5;
        double etmxt       = 0.00;
        double etmaxp      = 0.00;
        double etaalt      = 0.00;
        double etaext      = 0.00;
        unsigned int ipos  = 0;
        unsigned int ineg  = 999;
        //KRATOS_WATCH(ils);

        for (unsigned int i = 0; i <= ils+1; i++)
        {
            if (lseta(i) > etmaxp)
            {
                etmaxp = lseta(i);
            }
            if (prodr(i) < 0 && lseta(i) <= etaneg)
            {
                etaneg = lseta(i);
                ineg = i;
            }
        }

        // It decided if interpole or extrapole
        if (ineg != 999)
        {
            // Interpole
            // Find ipos = false of the S-L with prodr positive closer to ineg, but with a smaller S-L
            ipos = 0;
            for( unsigned int i = 0; i<=ils+1; i++)
            {
                if (prodr(i) >= 0 && lseta(i) <= lseta(ineg) && lseta(i) >= lseta(ipos))
                {
                    ipos = i;
                }
            }

            // It interpoles to find a S-L etaint
            etaint = prodr(ineg) * lseta(ipos) - prodr(ipos) * lseta(ineg);
            if ((prodr(ineg)-prodr(ipos)) == 0)
            {
                ico = 2;
                std::cout << " Warning: Division By Zero. Line Searches Failed" << std::endl;
                return;
            }
            else
            {
                etaint = etaint/(prodr(ineg)-prodr(ipos));
            }

            // Alternatively finds etaalt fo ensure a feasible change
            etaalt = lseta(ipos) + 0.2*(lseta(ineg)-lseta(ipos));

            // It tooks the maximal between etaint and etaalt
            etaint=std::max(etaint,etaalt);

            // Length of the minimal step
            if (etaint < metmna)
            {
                etaint = metmna;
                if (ico == 0)
                {
                    ico = 1;
                }
                else if (ico == 1)
                {
                    ico = 2;
                    std::cout << "Min Step Length Reached Twice" << std::endl;
                }
            }

            lseta(ils + 2) = etaint;
            meta = lseta(ils+2);
            return;
        }
        else
        {
            // Extrapole
            // It assigns temporally the length of the maximal step
            etmxt = mamp * etmaxp;
            if(etmxt > metmxa)
            {
                etmxt = metmxa;
            }

            // Extrapole the length between the current and previous step
            etaext=prodr(ils + 1) * lseta(ils) - prodr(ils) * lseta(ils + 1);
            if ((prodr(ils + 1) - prodr(ils)) == 0)
            {
                ico = 2;
                std::cout << " Warning: Division By Zero. Line Searches Filed" << std::endl;
                return;
            }
            else
            {
                etaext = etaext/(prodr(ils + 1) - prodr(ils));
            }

            lseta(ils + 2) = etaext;

            // It is accepted that is between the limits
            if (etaext <= 0 || etaext > etmxt)
            {
                lseta(ils + 2) = etmxt;
                meta = lseta(ils + 2);
                return;
            }
            if (lseta(ils + 2) == metmxa && ico == 1)
            {
                ico = 2;
                std::cout << "Max Step Length Again" << std::endl;
                return;
            }
            if (lseta(ils + 2) == metmxa)
            {
                ico = 1;
            }

            return;
        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It moves the nodes of the mesh when is invoqued
    * @param rmodel_part: The model part
    */

    void MoveMesh(ModelPart& rmodel_part)
    {
        KRATOS_TRY;

        for(ModelPart::NodeIterator i = rmodel_part.NodesBegin() ; i != rmodel_part.NodesEnd() ; i++)
        {
            array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
            (i)->X() = (i)->X0() + disp[0];
            (i)->Y() = (i)->Y0() + disp[1];
            (i)->Z() = (i)->Z0() + disp[2];
        }

        KRATOS_CATCH("");
    }

}; // end class lineseaches
} // namespace Kratos
#endif









