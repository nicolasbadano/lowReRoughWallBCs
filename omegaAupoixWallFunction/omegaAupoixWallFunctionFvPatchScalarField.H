/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

Class
    Foam::omegaAupoixWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    Provides a rough wall low-reynolds boundary condition for the omega field
    of k-omega SST model based on:

    Aupoix, B. (2015). Roughness Corrections for the k-ω Shear Stress Transport
    Model: Status and Proposals. Journal of Fluids Engineering, Vol. 137, 2.
    DOI: 10.1115/1.4028122

    Use together with:
      - kAupoixLowReWallFunction for k
      - calculated for nut

Usage
    \table
        Property     | Description             | Required    | Default value
        colebrook    | use colebrook fit       | no          | true
        Cmu          | model coefficient       | no          | 0.09
        E            | model coefficient       | no          | 9.8
        Ceps2        | model coefficient       | no          | 1.9
        Ks           |                         | yes
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            omegaAupoixWallFunction;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    omegaAupoixWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaAupoixWallFunctionFvPatchScalarField_H
#define omegaAupoixWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class omegaAupoixWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaAupoixWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        Switch colebrook_;

        //- Ks coefficient
        scalar Ks_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- Ceps2 coefficient
        scalar Ceps2_;

        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the Y+ at the edge of the laminar sublayer
        scalar yPlusLam(const scalar kappa, const scalar E);


public:

    //- Runtime type information
    TypeName("omegaAupoixWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaAupoixWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaAupoixWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given omegaAupoixWallFunctionFvPatchScalarField
        //  onto a new patch
        omegaAupoixWallFunctionFvPatchScalarField
        (
            const omegaAupoixWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        omegaAupoixWallFunctionFvPatchScalarField
        (
            const omegaAupoixWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaAupoixWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        omegaAupoixWallFunctionFvPatchScalarField
        (
            const omegaAupoixWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaAupoixWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes);


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
