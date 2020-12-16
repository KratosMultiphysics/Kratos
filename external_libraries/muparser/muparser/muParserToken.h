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

#ifndef MU_PARSER_TOKEN_H
#define MU_PARSER_TOKEN_H

#include <string>
#include <stack>
#include <vector>
#include <memory>

#if defined(_MSC_VER)
	#pragma warning(push)
	#pragma warning(disable : 26812) 
#endif

#include "muParserError.h"
#include "muParserCallback.h"

/** \file
	\brief This file contains the parser token definition.
*/

namespace mu
{
	/** \brief Encapsulation of the data for a single formula token.

		Formula token implementation. Part of the Math Parser Package.
		Formula tokens can be either one of the following:
		<ul>
			<li>value</li>
			<li>variable</li>
			<li>function with numerical arguments</li>
			<li>functions with a string as argument</li>
			<li>prefix operators</li>
			<li>infix operators</li>
			<li>binary operator</li>
		</ul>
	*/
	template<typename TBase, typename TString>
	class ParserToken final
	{
	private:

		ECmdCode  m_iCode;  ///< Type of the token; The token type is a constant of type #ECmdCode.
		ETypeCode m_iType;
		void* m_pTok;		///< Stores Token pointer; not applicable for all tokens
		int  m_iIdx;		///< An otional index to an external buffer storing the token data
		TString m_strTok;   ///< Token string
		TString m_strVal;   ///< Value for string variables
		value_type m_fVal;  ///< the value 
		std::unique_ptr<ParserCallback> m_pCallback;

	public:

		/** \brief Constructor (default).

			Sets token to an neutral state of type cmUNKNOWN.
			\throw nothrow
			\sa ECmdCode
		*/
		ParserToken()
			:m_iCode(cmUNKNOWN)
			, m_iType(tpVOID)
			, m_pTok(0)
			, m_iIdx(-1)
			, m_strTok()
			, m_strVal()
			, m_fVal(0)
			, m_pCallback()
		{}

		//------------------------------------------------------------------------------
		/** \brief Create token from another one.

			Implemented by calling Assign(...)
			\throw nothrow
			\post m_iType==cmUNKNOWN
			\sa #Assign
		*/
		ParserToken(const ParserToken& a_Tok)
		{
			Assign(a_Tok);
		}

		
		/** \brief Assignment operator.

			Copy token state from another token and return this.
			Implemented by calling Assign(...).
			\throw nothrow
		*/
		ParserToken& operator=(const ParserToken& a_Tok)
		{
			Assign(a_Tok);
			return *this;
		}


		/** \brief Copy token information from argument.

			\throw nothrow
		*/
		void Assign(const ParserToken& a_Tok)
		{
			m_iCode = a_Tok.m_iCode;
			m_pTok = a_Tok.m_pTok;
			m_strTok = a_Tok.m_strTok;
			m_iIdx = a_Tok.m_iIdx;
			m_strVal = a_Tok.m_strVal;
			m_iType = a_Tok.m_iType;
			m_fVal = a_Tok.m_fVal;
			// create new callback object if a_Tok has one 
			m_pCallback.reset(a_Tok.m_pCallback.get() ? a_Tok.m_pCallback->Clone() : 0);
		}

		//------------------------------------------------------------------------------
		/** \brief Assign a token type.

		  Token may not be of type value, variable or function. Those have separate set functions.

		  \pre [assert] a_iType!=cmVAR
		  \pre [assert] a_iType!=cmVAL
		  \pre [assert] a_iType!=cmFUNC
		  \post m_fVal = 0
		  \post m_pTok = 0
		*/
		ParserToken& Set(ECmdCode a_iType, const TString& a_strTok = TString())
		{
			// The following types can't be set this way, they have special Set functions
			MUP_ASSERT(a_iType != cmVAR);
			MUP_ASSERT(a_iType != cmVAL);
			MUP_ASSERT(a_iType != cmFUNC);

			m_iCode = a_iType;
			m_iType = tpVOID;
			m_pTok = 0;
			m_strTok = a_strTok;
			m_iIdx = -1;

			return *this;
		}

		//------------------------------------------------------------------------------
		/** \brief Set Callback type. */
		ParserToken& Set(const ParserCallback& a_pCallback, const TString& a_sTok)
		{
			MUP_ASSERT(a_pCallback.GetAddr());

			m_iCode = a_pCallback.GetCode();
			m_iType = tpVOID;
			m_strTok = a_sTok;
			m_pCallback.reset(new ParserCallback(a_pCallback));

			m_pTok = 0;
			m_iIdx = -1;

			return *this;
		}

		//------------------------------------------------------------------------------
		/** \brief Make this token a value token.

			Member variables not necessary for value tokens will be invalidated.
			\throw nothrow
		*/
		ParserToken& SetVal(TBase a_fVal, const TString& a_strTok = TString())
		{
			m_iCode = cmVAL;
			m_iType = tpDBL;
			m_fVal = a_fVal;
			m_strTok = a_strTok;
			m_iIdx = -1;

			m_pTok = 0;
			m_pCallback.reset(0);

			return *this;
		}

		//------------------------------------------------------------------------------
		/** \brief make this token a variable token.

			Member variables not necessary for variable tokens will be invalidated.
			\throw nothrow
		*/
		ParserToken& SetVar(TBase* a_pVar, const TString& a_strTok)
		{
			m_iCode = cmVAR;
			m_iType = tpDBL;
			m_strTok = a_strTok;
			m_iIdx = -1;
			m_pTok = (void*)a_pVar;
			m_pCallback.reset(0);
			return *this;
		}

		//------------------------------------------------------------------------------
		/** \brief Make this token a variable token.

			Member variables not necessary for variable tokens will be invalidated.
			\throw nothrow
		*/
		ParserToken& SetString(const TString& a_strTok, std::size_t a_iSize)
		{
			m_iCode = cmSTRING;
			m_iType = tpSTR;
			m_strTok = a_strTok;
			m_iIdx = static_cast<int>(a_iSize);

			m_pTok = 0;
			m_pCallback.reset(0);
			return *this;
		}

		//------------------------------------------------------------------------------
		/** \brief Set an index associated with the token related data.

			In cmSTRFUNC - This is the index to a string table in the main parser.
			\param a_iIdx The index the string function result will take in the bytecode parser.
			\throw exception_type if #a_iIdx<0 or #m_iType!=cmSTRING
		*/
		void SetIdx(int a_iIdx)
		{
			if (m_iCode != cmSTRING || a_iIdx < 0)
				throw ParserError(ecINTERNAL_ERROR);

			m_iIdx = a_iIdx;
		}

		//------------------------------------------------------------------------------
		/** \brief Return Index associated with the token related data.

			In cmSTRFUNC - This is the index to a string table in the main parser.

			\throw exception_type if #m_iIdx<0 or #m_iType!=cmSTRING
			\return The index the result will take in the Bytecode calculatin array (#m_iIdx).
		*/
		int GetIdx() const
		{
			if (m_iIdx < 0 || m_iCode != cmSTRING)
				throw ParserError(ecINTERNAL_ERROR);

			return m_iIdx;
		}

		//------------------------------------------------------------------------------
		/** \brief Return the token type.

			\return #m_iType
			\throw nothrow
		*/
		ECmdCode GetCode() const
		{
			if (m_pCallback.get())
			{
				return m_pCallback->GetCode();
			}
			else
			{
				return m_iCode;
			}
		}

		//------------------------------------------------------------------------------
		ETypeCode GetType() const
		{
			if (m_pCallback.get())
			{
				return m_pCallback->GetType();
			}
			else
			{
				return m_iType;
			}
		}

		//------------------------------------------------------------------------------
		int GetPri() const
		{
			if (!m_pCallback.get())
				throw ParserError(ecINTERNAL_ERROR);

			if (m_pCallback->GetCode() != cmOPRT_BIN && m_pCallback->GetCode() != cmOPRT_INFIX)
				throw ParserError(ecINTERNAL_ERROR);

			return m_pCallback->GetPri();
		}

		//------------------------------------------------------------------------------
		EOprtAssociativity GetAssociativity() const
		{
			if (m_pCallback.get() == nullptr || m_pCallback->GetCode() != cmOPRT_BIN)
				throw ParserError(ecINTERNAL_ERROR);

			return m_pCallback->GetAssociativity();
		}

		//------------------------------------------------------------------------------
		/** \brief Return the address of the callback function assoziated with
				   function and operator tokens.

			\return The pointer stored in #m_pTok.
			\throw exception_type if token type is non of:
				   <ul>
					 <li>cmFUNC</li>
					 <li>cmSTRFUNC</li>
					 <li>cmPOSTOP</li>
					 <li>cmINFIXOP</li>
					 <li>cmOPRT_BIN</li>
				   </ul>
			\sa ECmdCode
		*/
		generic_fun_type GetFuncAddr() const
		{
			return (m_pCallback.get()) ? (generic_fun_type)m_pCallback->GetAddr() : 0;
		}

		//------------------------------------------------------------------------------
		/** \biref Get value of the token.

			Only applicable to variable and value tokens.
			\throw exception_type if token is no value/variable token.
		*/
		TBase GetVal() const
		{
			switch (m_iCode)
			{
			case cmVAL:  return m_fVal;
			case cmVAR:  return *((TBase*)m_pTok);
			default:     throw ParserError(ecVAL_EXPECTED);
			}
		}

		//------------------------------------------------------------------------------
		/** \brief Get address of a variable token.

		  Valid only if m_iType==CmdVar.
		  \throw exception_type if token is no variable token.
		*/
		TBase* GetVar() const
		{
			if (m_iCode != cmVAR)
				throw ParserError(ecINTERNAL_ERROR);

			return (TBase*)m_pTok;
		}

		//------------------------------------------------------------------------------
		/** \brief Return the number of function arguments.

		  Valid only if m_iType==CmdFUNC.
		*/
		int GetArgCount() const
		{
			MUP_ASSERT(m_pCallback.get());

			if (!m_pCallback->GetAddr())
				throw ParserError(ecINTERNAL_ERROR);

			return m_pCallback->GetArgc();
		}

		//------------------------------------------------------------------------------
		/** \brief Return the token identifier.

			If #m_iType is cmSTRING the token identifier is the value of the string argument
			for a string function.
			\return #m_strTok
			\throw nothrow
			\sa m_strTok
		*/
		const TString& GetAsString() const
		{
			return m_strTok;
		}
	};
} // namespace mu

#if defined(_MSC_VER)
	#pragma warning(pop)
#endif

#endif
