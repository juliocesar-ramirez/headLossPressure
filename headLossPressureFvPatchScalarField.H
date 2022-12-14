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

Class
    Foam::headLossPressureFvPatchScalarField

Description
    This boundary condition provides a headLossPressure condition,
    calculated as:

        \f[
            Q = Q_{0} + Q_{p} + s*Q_{t}
        \f]

    where
    \vartable
        s      | single scalar value [units]
        Q_{0}  | single scalar value [units]
        Q_{p}  | scalar field across patch [units]
        Q_{t}  | scalar function of time [units]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd? | Default
        scalarData | single scalar value          | yes    |
        data       | single scalar value          | yes    |
        fieldData  | scalar field across patch    | yes    |
        timeVsData | scalar function of time      | yes    |
        wordData   | word, eg name of data object | no     | wordDefault
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       headLossPressure;
        scalarData -1;
        data       1;
        fieldData  uniform 3;
        timeVsData table (
                             (0 0)
                             (1 2)
                          );
        wordName   anotherName;
        value      uniform 4; // optional initial value
    }
    \endverbatim

SourceFiles
    headLossPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef BC_H
#define BC_H

#include "Function1.H"
#include "fixedValueFvPatchFields.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class headLossPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class headLossPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Single valued scalar quantity, e.g. a coefficient
        scalar scalarData_;

        scalar d_;
        scalar g_;
        scalar H_;

        //- Single valued Type quantity, e.g. reference pressure pRefValue_
        //  Other options include vector, tensor
        scalar data_;
        //- Field of Types, typically defined across patch faces
        //  e.g. total pressure p0_.  Other options include vectorField
        scalarField fieldData_;

        scalarField patm_;

        //- Type specified as a function of time for time-varying BCs
        autoPtr<Function1<scalar>> timeVsData_;

        //- Word entry, e.g. pName_ for name of the pressure field on database
        word wordData_;

        word phiName_;

        word UName_;
        //- Label, e.g. patch index, current time index
        label labelData_;

        //- Boolean for true/false, e.g. specify if flow rate is volumetric_
        bool boolData_;

        typedef Tuple2<vector, word> indexedVector;
        const List<indexedVector> minorLossFactor_;

        typedef Tuple2<vector, word> indexedVector2;
        const List<indexedVector2> frictionLossFactor_;

        // Private Member Functions

        //- Return current time
        scalar t() const;
        scalar f(scalar di, scalar epsilon, scalar U);

      public:
        //- Runtime type information
        TypeName("headLossPressure");

        // Constructors

        //- Construct from patch and internal field
        headLossPressureFvPatchScalarField(
            const fvPatch &, const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        headLossPressureFvPatchScalarField(
            const fvPatch &, const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        headLossPressureFvPatchScalarField(
            const headLossPressureFvPatchScalarField &, const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Disallow copy without setting internal field reference
        headLossPressureFvPatchScalarField(
            const headLossPressureFvPatchScalarField &) = delete;

        //- Copy constructor setting internal field reference
        headLossPressureFvPatchScalarField(
            const headLossPressureFvPatchScalarField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField>
        clone(const DimensionedField<scalar, volMesh> &iF) const {
          return tmp<fvPatchScalarField>(
              new headLossPressureFvPatchScalarField(*this, iF));
        }

    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
