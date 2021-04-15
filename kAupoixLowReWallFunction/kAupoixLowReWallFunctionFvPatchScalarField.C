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

#include "kAupoixLowReWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kAupoixLowReWallFunctionFvPatchScalarField::checkType()
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


scalar kAupoixLowReWallFunctionFvPatchScalarField::yPlusLam
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

kAupoixLowReWallFunctionFvPatchScalarField::kAupoixLowReWallFunctionFvPatchScalarField
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
    yPlusLam_(yPlusLam(kappa_, E_)),
    colebrook_(true)
{
    checkType();
}


kAupoixLowReWallFunctionFvPatchScalarField::kAupoixLowReWallFunctionFvPatchScalarField
(
    const kAupoixLowReWallFunctionFvPatchScalarField& ptf,
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
    yPlusLam_(ptf.yPlusLam_),
    colebrook_(ptf.colebrook_)
{
    checkType();
}


kAupoixLowReWallFunctionFvPatchScalarField::kAupoixLowReWallFunctionFvPatchScalarField
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
    yPlusLam_(yPlusLam(kappa_, E_)),
    colebrook_(dict.lookupOrDefault<Switch>("colebrook", true))
{
    checkType();
}


kAupoixLowReWallFunctionFvPatchScalarField::kAupoixLowReWallFunctionFvPatchScalarField
(
    const kAupoixLowReWallFunctionFvPatchScalarField& kwfpsf
)
:
    fixedValueFvPatchField<scalar>(kwfpsf),
    Ks_(kwfpsf.Ks_),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_),
    colebrook_(kwfpsf.colebrook_)
{
    checkType();
}


kAupoixLowReWallFunctionFvPatchScalarField::kAupoixLowReWallFunctionFvPatchScalarField
(
    const kAupoixLowReWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(kwfpsf, iF),
    Ks_(kwfpsf.Ks_),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_),
    colebrook_(kwfpsf.colebrook_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kAupoixLowReWallFunctionFvPatchScalarField::updateCoeffs()
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
    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradU(mag(Uw.snGrad()));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const scalar Cmu25 = pow025(Cmu_);

    scalarField& kw = *this;

    // Set k wall values
    forAll(kw, facei)
    {
        scalar ut = sqrt( (nuw[facei]+nutw[facei]) * magGradU[facei] );
        scalar KsPlus = Ks_*ut/nuw[facei];
        if (colebrook_)
        {
            //  Colebrook
            kw[facei] = sqr(ut) / sqrt(Cmu_) * tanh( (log(KsPlus/30)/log(10.0) + 1 - tanh(KsPlus/125))*tanh(KsPlus/125));
        } else {
            //  Nikuradse
            kw[facei] = sqr(ut) / sqrt(Cmu_) * tanh( (log(KsPlus/30)/log(8.0) + 0.5 * (1 - tanh(KsPlus/100)))*tanh(KsPlus/75));
        }
    }

    // Limit kw to avoid failure of the turbulence model due to division by kw
    kw = max(kw, SMALL);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void kAupoixLowReWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void kAupoixLowReWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("Ks") << Ks_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ceps2") << Ceps2_ << token::END_STATEMENT << nl;
    os.writeKeyword("colebrook") << colebrook_ << token::END_STATEMENT << nl;
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    kAupoixLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
