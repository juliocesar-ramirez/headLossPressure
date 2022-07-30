/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "headLossPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::headLossPressureFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::headLossPressureFvPatchScalarField::
headLossPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    scalarData_(0.0),
    data_(Zero),
    fieldData_(p.size(), Zero),
    timeVsData_(),
    wordData_("wordDefault"),
    labelData_(-1),
    boolData_(false)
{
}


Foam::headLossPressureFvPatchScalarField::
headLossPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    scalarData_(dict.lookup<scalar>("scalarData")),
    data_(dict.lookup<scalar>("data")),
    fieldData_("fieldData", dict, p.size()),
    timeVsData_(Function1<scalar>::New("timeVsData", dict)),
    wordData_(dict.lookupOrDefault<word>("wordName", "wordDefault")),
    labelData_(-1),
    boolData_(false)
{


    fixedValueFvPatchScalarField::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    */
}


Foam::headLossPressureFvPatchScalarField::
headLossPressureFvPatchScalarField
(
    const headLossPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(mapper(ptf.fieldData_)),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


Foam::headLossPressureFvPatchScalarField::
headLossPressureFvPatchScalarField
(
    const headLossPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::headLossPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(fieldData_, fieldData_);
}


void Foam::headLossPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const headLossPressureFvPatchScalarField& tiptf =
        refCast<const headLossPressureFvPatchScalarField>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::headLossPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fixedValueFvPatchScalarField::operator==
    (
        data_
      + fieldData_
      + scalarData_*timeVsData_->value(t())
    );


    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::headLossPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "scalarData", scalarData_);
    writeEntry(os, "data", data_);
    writeEntry(os, "fieldData", fieldData_);
    writeEntry(os, timeVsData_());
    writeEntry(os, "wordData", wordData_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        headLossPressureFvPatchScalarField
    );
}

// ************************************************************************* //
