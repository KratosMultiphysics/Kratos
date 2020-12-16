/*

	 _____  __ _____________ _______  ______ ___________
	/     \|  |  \____ \__  \\_  __ \/  ___// __ \_  __ \
   |  Y Y  \  |  /  |_> > __ \|  | \/\___ \\  ___/|  | \/
   |__|_|  /____/|   __(____  /__|  /____  >\___  >__|
		 \/      |__|       \/           \/     \/
   Copyright (C) 2004 - 2020 Ingo Berg

	Redistribution and use in source and binary forms, with or without modification, are permitted
	provided that the following conditions are met:

	  * Redistributions of source code must retain the above copyright notice, this list of
		conditions and the following disclaimer.
	  * Redistributions in binary form must reproduce the above copyright notice, this list of
		conditions and the following disclaimer in the documentation and/or other materials provided
		with the distribution.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
	IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
	FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
	IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
	OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef MU_PARSER_CALLBACK_H
#define MU_PARSER_CALLBACK_H

#include "muParserDef.h"

/** \file
	\brief Definition of the parser callback class.
*/

namespace mu
{

	/** \brief Encapsulation of prototypes for a numerical parser function.

		Encapsulates the prototyp for numerical parser functions. The class
		stores the number of arguments for parser functions as well
		as additional flags indication the function is non optimizeable.
		The pointer to the callback function pointer is stored as void*
		and needs to be casted according to the argument count.
		Negative argument counts indicate a parser function with a variable number
		of arguments.
	*/
	class API_EXPORT_CXX ParserCallback final
	{
	public:
		ParserCallback(fun_type0  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type1  a_pFun, bool a_bAllowOpti, int a_iPrec = -1, ECmdCode a_iCode = cmFUNC);
		ParserCallback(fun_type2  a_pFun, bool a_bAllowOpti, int a_iPrec, EOprtAssociativity a_eAssociativity);
		ParserCallback(fun_type2  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type3  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type4  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type5  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type6  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type7  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type8  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type9  a_pFun, bool a_bAllowOpti);
		ParserCallback(fun_type10 a_pFun, bool a_bAllowOpti);

		ParserCallback(bulkfun_type0  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type1  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type2  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type3  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type4  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type5  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type6  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type7  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type8  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type9  a_pFun, bool a_bAllowOpti);
		ParserCallback(bulkfun_type10 a_pFun, bool a_bAllowOpti);

		ParserCallback(multfun_type a_pFun, bool a_bAllowOpti);
		ParserCallback(strfun_type1 a_pFun, bool a_bAllowOpti);
		ParserCallback(strfun_type2 a_pFun, bool a_bAllowOpti);
		ParserCallback(strfun_type3 a_pFun, bool a_bAllowOpti);
		ParserCallback(strfun_type4 a_pFun, bool a_bAllowOpti);
		ParserCallback(strfun_type5 a_pFun, bool a_bAllowOpti);
		ParserCallback();
		ParserCallback(const ParserCallback& a_Fun);
    ParserCallback & operator=(const ParserCallback& a_Fun);

		ParserCallback* Clone() const;

		bool  IsOptimizable() const;
		void* GetAddr() const;
		ECmdCode  GetCode() const;
		ETypeCode GetType() const;
		int GetPri()  const;
		EOprtAssociativity GetAssociativity() const;
		int GetArgc() const;

	private:
		void* m_pFun;                   ///< Pointer to the callback function, casted to void

		/** \brief Number of numeric function arguments

			This number is negative for functions with variable number of arguments. in this cases
			they represent the actual number of arguments found.
		*/
		int   m_iArgc;
		int   m_iPri;                   ///< Valid only for binary and infix operators; Operator precedence.
		EOprtAssociativity m_eOprtAsct; ///< Operator associativity; Valid only for binary operators 
		ECmdCode  m_iCode;
		ETypeCode m_iType;
		bool  m_bAllowOpti;             ///< Flag indication optimizeability 
	};


	/** \brief Container for Callback objects. */
	typedef std::map<string_type, ParserCallback> funmap_type;

} // namespace mu

#endif

