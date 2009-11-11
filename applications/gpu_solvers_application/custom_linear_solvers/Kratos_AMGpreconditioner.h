/*
==============================================================================
KratosGPUApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2009
Pooyan Dadvand, Riccardo Rossi, Isaac Gallego, Farshid Mossaiby 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
isaac.gallego.pla@gmail.com
mossaiby@yahoo.com
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

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

#ifndef _KRATOS_AMGPRECONDITIONER_H
#define	_KRATOS_AMGPRECONDITIONER_H

#include "AMGpreconditioner.h"
#include "includes/ublas_interface.h"
#include <cstdio>
#include <cstdlib>


class Kratos_AMGpreconditioner : public AMGpreconditioner {
public:
	Kratos_AMGpreconditioner() {};
    	Kratos_AMGpreconditioner(double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, const Kratos::Vector _preSweeps, const Kratos::Vector _postSweeps, bool isPreconditioner) : AMGpreconditioner(_W, _numLevelsRoh, _assumeZerosForEachStep, _numMaxHierarchyLevels, _minimumSizeAllowed, isPreconditioner)
	{
		 preSww = new size_t [_numMaxHierarchyLevels];
		 postSww = new size_t [_numMaxHierarchyLevels];

		for(size_t i=0; i < _numMaxHierarchyLevels; i++){
			preSww[i] = (size_t)_preSweeps[i];
			postSww[i] = (size_t)_postSweeps[i];
		}

		setPreSweeps( preSww);
		setPostSweeps( postSww);
		
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

