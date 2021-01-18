/*
!=========================================================================================
!Copyright (c) 2009-2019, The Regents of the University of Massachusetts, Amherst.
!Developed by E. Polizzi
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification, 
!are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list of conditions 
!   and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
!   and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
!    products derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
!ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!==========================================================================================
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// External declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//PFEAST
extern void pfeastinit_(int *fpm,MPI_Fint* COMM_WORLD,int *nL3);
extern void pfeast_distribution_type_(int *N,int *isa,int *jsa,MPI_Fint* L3_COMM_WORLD,int *distype);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FEAST interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//PFEAST

void pfeastinit(int *fpm,MPI_Fint* COMM_WORLD,int *nL3) {
     pfeastinit_(fpm,COMM_WORLD,nL3);
}
void PFEASTINIT(int *fpm,MPI_Fint* COMM_WORLD,int *nL3) {
     pfeastinit_(fpm,COMM_WORLD,nL3);
}
void pfeast_distribution_type(int *N,int *isa,int *jsa,MPI_Fint* L3_COMM_WORLD,int *distype) {
     pfeast_distribution_type_(N,isa,jsa,L3_COMM_WORLD,distype);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

















