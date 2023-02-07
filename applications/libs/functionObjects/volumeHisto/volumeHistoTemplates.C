/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "volumeHisto.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::volumeHisto::foundObject
(
    const word& name,
    const bool verbose
) const
{
    if (fvMeshFunctionObject::foundObject<Type>(name))
    {
        return true;
    }
    else
    {
        if (debug || verbose)
        {
            Warning
                << "    functionObjects::" << type() << " " << this->name()
                << " cannot find required object " << name << " of type "
                << Type::typeName << endl;
        }

        return false;
    }
}

template<class Type>
void Foam::functionObjects::volumeHisto::integrateFields
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);
        
        Type volIntegral = pTraits<Type>::zero;
        scalar volume = 0;

        scalar interval = maxVal_ - minVal_;
        scalar binX = interval/static_cast<double>(numBins_);
        
        Info<<"Total interval: "<<interval<<nl;
        Info<<"Bin wigth: "<<binX<<nl;


        scalarList histo(numBins_, 0.0);
            // volume
            forAll (field, cellI)
            {
            }

        forAll (field, cellI)
        {
            vector pos=mesh_.C()[cellI];
            if( cellInsideTheBox(pos) )
            {
                scalar magField = Foam::mag(field[cellI]);
                int bin = static_cast<int>( (magField - minVal_)/binX );
                if(bin>=numBins_) bin = numBins_ - 1;
                histo[bin] += mesh_.V()[cellI];
                volume += mesh_.V()[cellI];
            }
        }

        forAll(histo, hI)
        {
            histo[hI] = histo[hI]/volume;
        }

        Info<<"Volume:  "<<volume<<endl;
        Info<<"Histogram of "<<fieldName<<":  "<<histo<<endl;
    }
}



// ************************************************************************* //
