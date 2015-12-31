// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef KRATOS_EXPORT_API_H
#define KRATOS_EXPORT_API_H

#undef KRATOS_API_EXPORT
#undef KRATOS_API_IMPORT
#ifdef _WIN32
  #define KRATOS_API_EXPORT __declspec(dllexport)
  #define KRATOS_API_IMPORT __declspec(dllimport)
#else
  #define KRATOS_API_EXPORT __attribute__((visibility("default")))
  #define KRATOS_API_IMPORT __attribute__((visibility("default")))
#endif

// This fixes MSVC not expanding __VA_ARGS__ as defined in the C99 standard
#define KRATOS_EXPAND(A) A

// This expands the API call to either import or export based on the
// number of the arguments in API() call.
#define KRATOS_API_CALL(x,T1,T2,T3,...) T3
#define KRATOS_API(...) \
  KRATOS_EXPAND(KRATOS_API_CALL(,##__VA_ARGS__,KRATOS_API_EXPORT,KRATOS_API_IMPORT))
#define KRATOS_NO_EXPORT(...)

#endif
