//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Pooyan Dadvand
//
//


// System includes


// External includes


// Project includes
#include "includes/code_location.h"
#include "includes/checks.h"


namespace Kratos
{
    CodeLocation::CodeLocation() :
		mFileName("Unknown"), mFunctionName("Unknown"), mLineNumber(-1) {}

	CodeLocation::CodeLocation(std::string const& FileName, std::string const& FunctionName, std::size_t LineNumber) :
		mFileName(FileName), mFunctionName(FunctionName), mLineNumber(LineNumber) {}

	CodeLocation::CodeLocation(CodeLocation const & Other) :
		mFileName(Other.mFileName), mFunctionName(Other.mFunctionName), mLineNumber(Other.mLineNumber) {}


	const std::string& CodeLocation::GetFileName() const	{
		return mFileName;
	}

	const std::string& CodeLocation::GetFunctionName()	const {
		return mFunctionName;
	}

	int CodeLocation::GetLineNumber() const {
		return mLineNumber;
	}

	std::string CodeLocation::CleanFileName() const
	{
		std::string clean_file_name(mFileName);
		ReplaceAll(clean_file_name, "\\", "/");

		std::size_t kratos_root_position = clean_file_name.rfind("/application/");
		if (kratos_root_position != std::string::npos)
			clean_file_name.erase(0, kratos_root_position);


		if (kratos_root_position == std::string::npos)
			kratos_root_position = clean_file_name.rfind("/kratos/");
		if (kratos_root_position != std::string::npos)
			clean_file_name.erase(0, kratos_root_position + 1 );
		return clean_file_name;
	}

	std::string CodeLocation::CleanFunctionName() const
	{
		std::string clean_function_name(mFunctionName);

		// Applying filters. Note that The sequence is important
		// NOTE: Please add new fliters at the end of the list!!
		RemoveNamespace(clean_function_name, "Kratos");
		RemoveNamespace(clean_function_name, "std");
		ReduceTemplateArgumentsToFirstN(clean_function_name, "ublas::vector", 1);
		ReduceTemplateArgumentsToFirstN(clean_function_name, "ublas::matrix", 1);
		ReduceTemplateArgumentsToFirstN(clean_function_name, "iterators::indirect_iterator", 1);
		ReduceTemplateArgumentsToFirstN(clean_function_name, "PointerVectorSet", 1);
		ReduceTemplateArgumentsToFirstN(clean_function_name, "basic_string", 1);
		ReplaceAll(clean_function_name, "__int64", "int");
		ReplaceAll(clean_function_name, "basic_string<char,...>", "string");
		ReduceTemplateArgumentsToFirstN(clean_function_name, "compressed_matrix", 0);
		ReplaceAll(clean_function_name, "ublas::vector<double,...>", "Vector");
		ReplaceAll(clean_function_name, "ublas::DenseMatrix<double,...>", "Matrix");
		ReduceTemplateArgumentsToFirstN(clean_function_name, "ResidualBasedBlockBuilderAndSolver", 1);
		ReduceTemplateArgumentsToFirstN(clean_function_name, "ResidualBasedLinearStrategy", 1);
		ReplaceAll(clean_function_name, "Dof<double>", "Dof");
		ReplaceAll(clean_function_name, "Node<3, Dof >", "Node");



		return clean_function_name;
	}

	void CodeLocation::RemoveNamespace(std::string& FunctionName, const std::string& Namespace)
	{
		std::string namespace_string = Namespace + "::";
		ReplaceAll(FunctionName, namespace_string, "");

	}

	void CodeLocation::ReduceTemplateArgumentsToFirstN(std::string& FunctionName, const std::string& TemplateName, std::size_t NumberOfArgumentsToKeep)
	{
		std::size_t start_position = 0;
		while ((start_position = FunctionName.find(TemplateName, start_position)) != std::string::npos)	{
			start_position += TemplateName.size();
			std::size_t template_position = GetNextPositionSkippingWhiteSpaces(FunctionName, start_position);
			//if (FunctionName[template_position] == '<') // It is not a template
			//	continue;
			std::string::iterator arguments_iterator = FunctionName.begin() + template_position;
			std::size_t number_of_open_pranthesis = 0;
			std::size_t number_of_open_templates = 1;
			arguments_iterator++;
			std::size_t number_of_arguments = 0;
			if (*arguments_iterator != '>')
				number_of_arguments = 1;
			std::size_t replace_position = std::string::npos;
			if ((number_of_arguments > NumberOfArgumentsToKeep))
				replace_position = arguments_iterator - FunctionName.begin();
			while (arguments_iterator != FunctionName.end() && number_of_open_templates > 0)
			{
				if (*arguments_iterator == '<')
					number_of_open_templates++;
				else if (*arguments_iterator == '>')
					number_of_open_templates--;
				else if (*arguments_iterator == '(')
					number_of_open_pranthesis++;
				else if (*arguments_iterator == ')')
					number_of_open_pranthesis--;
				else if (*arguments_iterator == ',')
				{
					if (number_of_open_pranthesis == 0 && number_of_open_templates == 1)
						number_of_arguments++;
					//else
					//	std::cout << "number_of_open_pranthesis " << number_of_open_pranthesis << "number_of_open_templates" << number_of_open_templates << std::endl;
					if ((number_of_arguments > NumberOfArgumentsToKeep) && (replace_position == std::string::npos))
						replace_position = arguments_iterator - FunctionName.begin() + 1;
				}
				//std::cout << *arguments_iterator;
				arguments_iterator++;

			}
			if (replace_position != std::string::npos)
			{
				std::size_t replace_size = arguments_iterator - FunctionName.begin() - 1 - replace_position;
				FunctionName.replace(replace_position, replace_size, "...");

			}
			//std::cout << number_of_arguments << " arguments " << std::endl;
		}
	}



	void CodeLocation::ReplaceAll(std::string& ThisString, const std::string& FromString, const std::string& ToString)
	{
		std::size_t start_position = 0;
		while ((start_position = ThisString.find(FromString, start_position)) != std::string::npos)
		{
			ThisString.replace(start_position, FromString.length(), ToString);
			start_position += ToString.length(); // ...
		}

	}

	std::size_t CodeLocation::GetNextPositionSkippingWhiteSpaces(std::string const& ThisString, std::size_t Position)
	{
		char c = ThisString[Position];
		while (IsWhiteSpace(c))
			c = ThisString[++Position];
		return Position;
	}

	bool CodeLocation::IsWhiteSpace(char C)
	{
		return ((C == ' ') || (C == '\t') || (C == '\r') || (C == '\n'));
	}

	std::ostream & operator <<(std::ostream& rOStream,
		const CodeLocation& Location)
	{
		rOStream << Location.CleanFileName() << ":" << Location.GetLineNumber() << ":" << Location.CleanFunctionName();
		return rOStream;
	}

}  // namespace Kratos.
