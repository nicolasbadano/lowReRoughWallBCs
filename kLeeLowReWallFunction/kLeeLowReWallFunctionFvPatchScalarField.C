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

\*---------------------------------------------------------------------------*/

#include "kLeeLowReWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kLeeLowReWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


scalar kLeeLowReWallFunctionFvPatchScalarField::yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kLeeLowReWallFunctionFvPatchScalarField::kLeeLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Ks_(0.1),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    Ceps2_(1.9),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


kLeeLowReWallFunctionFvPatchScalarField::kLeeLowReWallFunctionFvPatchScalarField
(
    const kLeeLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Ks_(ptf.Ks_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    Ceps2_(ptf.Ceps2_),
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


kLeeLowReWallFunctionFvPatchScalarField::kLeeLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Ks_(dict.lookupOrDefault<scalar>("Ks", 0.1)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    Ceps2_(dict.lookupOrDefault<scalar>("Ceps2", 1.9)),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


kLeeLowReWallFunctionFvPatchScalarField::kLeeLowReWallFunctionFvPatchScalarField
(
    const kLeeLowReWallFunctionFvPatchScalarField& kwfpsf
)
:
    fixedValueFvPatchField<scalar>(kwfpsf),
    Ks_(kwfpsf.Ks_),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_)
{
    checkType();
}


kLeeLowReWallFunctionFvPatchScalarField::kLeeLowReWallFunctionFvPatchScalarField
(
    const kLeeLowReWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(kwfpsf, iF),
    Ks_(kwfpsf.Ks_),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kLeeLowReWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradU(mag(Uw.snGrad()));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const scalar beta5 = pow(0.075, 0.5);

    scalarField& kw = *this;
    forAll(kw, facei)
    {
        scalar ut = sqrt( (nuw[facei]+nutw[facei]) * magGradU[facei] );
        scalar KsPlus = Ks_*ut/nuw[facei];
        kw[facei] = min(1.,KsPlus/90.)*ut*ut/beta5;
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void kLeeLowReWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void kLeeLowReWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("Ks") << Ks_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ceps2") << Ceps2_ << token::END_STATEMENT << nl;
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    kLeeLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
