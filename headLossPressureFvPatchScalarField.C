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

#include "DimensionedScalarField.H"
#include "UList.H"
#include "doubleScalar.H"
#include "fvPatchField.H"
#include "headLossPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "quaternion.H"
#include "scalar.H"
#include "surfaceFieldsFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "volFieldsFwd.H"
#include <cmath>
#include <math.h>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::headLossPressureFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}

Foam::scalar Foam::headLossPressureFvPatchScalarField::f(scalar di,scalar epsilon, scalar U) {
    scalar f1=Zero;

    const dictionary &transportProperties =
      db().lookupObject<dictionary>("transportProperties");
  const dimensionedScalar nuDimScal("nu", dimViscosity, transportProperties);
  const scalar &nu = nuDimScal.value();

  scalar reynolds = di * U / nu;
  if (reynolds == 0) {
      f1=0;
  }
  else if(reynolds<2300){
      f1=64/reynolds;
  }
  else {
      scalar fInicial=sqr(1/(-1.8*log10((6.11/reynolds)+pow((epsilon/di)/3.7,1.11))));
      f1=sqr(1/(-2*log10((epsilon/di)/3.7+2.51/(reynolds*sqr(fInicial)))));
  }
  return f1;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::headLossPressureFvPatchScalarField::headLossPressureFvPatchScalarField(
    const fvPatch &p, const DimensionedField<scalar, volMesh> &iF)
    : fixedValueFvPatchScalarField(p, iF), scalarData_(0.0), d_(0.0),g_(0.0), H_(0.0),
      data_(Zero), fieldData_(p.size(), Zero), timeVsData_(),
      wordData_("wordDefault"), phiName_("phi"),UName_("U"), labelData_(-1),
      boolData_(false), minorLossFactor_() {}

Foam::headLossPressureFvPatchScalarField::headLossPressureFvPatchScalarField(
    const fvPatch &p, const DimensionedField<scalar, volMesh> &iF,
    const dictionary &dict)
    : fixedValueFvPatchScalarField(p, iF),
      scalarData_(dict.lookup<scalar>("scalarData")),
      d_(dict.lookup<scalar>("d")),g_(dict.lookup<scalar>("g")),
      H_(dict.lookup<scalar>("H")),
      data_(dict.lookup<scalar>("data")),
      fieldData_("fieldData", dict, p.size()),
      timeVsData_(Function1<scalar>::New("timeVsData", dict)),
      wordData_(dict.lookupOrDefault<word>("wordName", "wordDefault")),
      phiName_(dict.lookupOrDefault<word>("phi", "phi")),
      UName_(dict.lookupOrDefault<word>("U", "U")),
      labelData_(-1), boolData_(false),minorLossFactor_(dict.lookup("minorLossFactor")) {

  fixedValueFvPatchScalarField::evaluate();

  /*
  // Initialise with the value entry if evaluation is not possible
  fvPatchScalarField::operator=
  (
      scalarField("value", dict, p.size())
  );
  */
}

Foam::headLossPressureFvPatchScalarField::headLossPressureFvPatchScalarField(
    const headLossPressureFvPatchScalarField &ptf, const fvPatch &p,
    const DimensionedField<scalar, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchScalarField(ptf, p, iF, mapper),
      scalarData_(ptf.scalarData_), d_(ptf.d_),g_(ptf.g_),H_(ptf.H_), data_(ptf.data_),
      fieldData_(mapper(ptf.fieldData_)), timeVsData_(ptf.timeVsData_, false),
      wordData_(ptf.wordData_), phiName_(ptf.phiName_),UName_(ptf.UName_),labelData_(-1), boolData_(ptf.boolData_),minorLossFactor_(ptf.minorLossFactor_) {}

Foam::headLossPressureFvPatchScalarField::headLossPressureFvPatchScalarField(
    const headLossPressureFvPatchScalarField &ptf,
    const DimensionedField<scalar, volMesh> &iF)
    : fixedValueFvPatchScalarField(ptf, iF), scalarData_(ptf.scalarData_),
      d_(ptf.d_),g_(ptf.g_),H_(ptf.H_) ,data_(ptf.data_), fieldData_(ptf.fieldData_),
      timeVsData_(ptf.timeVsData_, false), wordData_(ptf.wordData_), phiName_(ptf.phiName_),UName_(ptf.UName_),
      labelData_(-1), boolData_(ptf.boolData_),minorLossFactor_(ptf.minorLossFactor_) {}

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

    const vectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    scalar sumPhi = gSum(patch().lookupPatchField<surfaceScalarField,  scalar>( phiName_));
    const scalar totalArea=gSum(patch().magSf());
    scalar Uavg=mag(sumPhi)/totalArea;

    // perdidas menores
    scalar dpMinor=0;
    scalar k;
    scalar dMinor;
    for (int i=0; i<minorLossFactor_.size(); i++) {
        dMinor=minorLossFactor_[i].first().component(0);
        k=minorLossFactor_[i].first().component(1);
        dpMinor+=k*0.5*sqrt(Uavg*sqrt(d_*dMinor));
    }

    //perdidas por friction
    scalar dpFriction=0;
    scalar dFriccion,epsilon,L; 
    for (int i=0; i<frictionLossFactor_.size(); i++) {
        dFriccion = frictionLossFactor_[i].first().component(0);
         epsilon = frictionLossFactor_[i].first().component(1);
         L = frictionLossFactor_[i].first().component(2);

         dpFriction=f(dFriccion,  epsilon,Uavg)*L/dFriccion*0.5*sqrt(Uavg*sqr(d_/dFriccion))+dpFriction;

    }

    fixedValueFvPatchScalarField::operator==(patm_ + g_ * H_ - dpFriction -
                                             dpMinor-0.5*magSqr(Up));

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::headLossPressureFvPatchScalarField::write(Ostream &os) const {
        fvPatchScalarField::write(os);
  writeEntry(os, "scalarData", scalarData_);
  writeEntry(os, "d", d_);
  writeEntry(os, "g", g_);
  writeEntry(os, "H", H_);
  writeEntry(os, "data", data_);
  writeEntry(os, "fieldData", fieldData_);
  writeEntry(os, timeVsData_());
  writeEntry(os, "wordData", wordData_);
  writeEntry(os, "phi", phiName_);
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
