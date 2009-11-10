/* 
 * File:   AMGpreconditioner.h
 * Author: isaac
 *
 * Created on 1 / octubre / 2009, 09:39
 */

#ifndef _KRATOS_AMGPRECONDITIONER_H
#define	_KRATOS_AMGPRECONDITIONER_H

#include "AMGpreconditioner.h"
#include "includes/ublas_interface.h"
#include <cstdio>
#include <cstdlib>


class Kratos_AMGpreconditioner : public AMGpreconditioner {
public:
	Kratos_AMGpreconditioner() {};
    	Kratos_AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, const Kratos::Vector _preSweeps, const Kratos::Vector _postSweeps) : AMGpreconditioner(_W, _numLevelsRoh, _assumeZerosForEachStep, _numMaxHierarchyLevels, _minimumSizeAllowed)
	{
		 preSww = new size_t [_numMaxHierarchyLevels];
		 postSww = new size_t [_numMaxHierarchyLevels];


		//size_t pp [_numMaxHierarchyLevels];
		//size_t qq [_numMaxHierarchyLevels];
		for(size_t i=0; i < _numMaxHierarchyLevels; i++){
			preSww[i] = (size_t)_preSweeps[i];
			postSww[i] = (size_t)_postSweeps[i];
		}
		//printf("\n\nEstem a KRATOS_AMGPReconditioner, i posem ara la W = %f, i despres els sweeps %u, i %u, i %u, i %u\n", _W, preSww[0], preSww[1], preSww[2], preSww[3]);

		setPreSweeps( preSww);
		setPostSweeps( postSww);
		
		/*printf("PRINTING from KratosAMGPreconditioner variable values:\n W = %f, Roh = %u, Zeros = %s, HierarchyLevels = %u, minimumSize = %u, firstPre = %u, secondPre = %u\n", _W, _numLevelsRoh, (_assumeZerosForEachStep)?"true":"false", _numMaxHierarchyLevels, _minimumSizeAllowed, preSww[0], preSww[1]);*/

		//exit(0);
	};
    	virtual ~Kratos_AMGpreconditioner(){
		delete[] preSww;
		delete[] postSww;
    	};

private:
	size_t* preSww;
	size_t* postSww;
};

#endif	/* _AMGPRECONDITIONER_H */

