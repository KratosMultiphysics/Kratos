//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// Project includes
#include "includes/exception.hpp"

namespace CoSimIO {
namespace Internals {

Exception::Exception(const std::string& rWhat)
    : std::exception(),
        mMessage(rWhat),
        mCallStack()
{
    update_what();
}

Exception::Exception(const std::string& rWhat, const CodeLocation& rLocation)
    : Exception(rWhat)
{
    add_to_call_stack(rLocation);
}

Exception::Exception(const Exception& Other)
    : std::exception(Other),
      mWhat(Other.mWhat),
      mMessage(Other.mMessage),
      mCallStack(Other.mCallStack)
{
}

const char* Exception::what() const noexcept
{
    return mWhat.c_str();
}

const std::string& Exception::message() const
{
    return mMessage;
}

Exception& Exception::operator << (std::ostream& (*pf)(std::ostream&))
{
    std::stringstream buffer;
    pf(buffer);

    append_message(buffer.str());

    return *this;
}

Exception& Exception::operator << (const char* pString)
{
    append_message(pString);
    return *this;
}

Exception& Exception::operator << (const CodeLocation& rLocation)
{
    add_to_call_stack(rLocation);
    return *this;
}

void Exception::append_message(const std::string& rMessage)
{
    mMessage.append(rMessage);
    update_what();
}

void Exception::add_to_call_stack(const CodeLocation& rLocation)
{
    mCallStack.push_back(rLocation);
    update_what();
}

void Exception::update_what(){
    std::stringstream buffer;
    buffer << message() << "\n";
    if (mCallStack.empty()) {
        buffer << "in Unknown Location";
    } else {
        buffer << "in " << mCallStack[0] << "\n";
        for (auto i = mCallStack.begin()+1; i != mCallStack.end(); ++i) {
            buffer << "   " << *i << "\n";
        }
    }
    mWhat = buffer.str();
}

} // namespace Internals
} // namespace CoSimIO
