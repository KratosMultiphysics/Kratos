/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Template for use with dynamic code generation of a
    fixedValue fvPatchField.

    - without state

SourceFiles
    fixedValueFvPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueFvPatchTemplateVectorField_H
#define fixedValueFvPatchTemplateVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated FixedValueFvPatch
\*---------------------------------------------------------------------------*/

class swirlFixedValueFvPatchVectorField
:
    public fixedValueFvPatchField<vector>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("swirl");


    // Constructors

        //- Construct from patch and internal field
        swirlFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swirlFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        swirlFixedValueFvPatchVectorField
        (
            const swirlFixedValueFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swirlFixedValueFvPatchVectorField
        (
            const swirlFixedValueFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField >
            (
                new swirlFixedValueFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        swirlFixedValueFvPatchVectorField
        (
            const swirlFixedValueFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField >
            (
                new swirlFixedValueFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~swirlFixedValueFvPatchVectorField();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

