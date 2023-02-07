/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "tomo2cfd.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "volFields.H"

#include "memInfo.H"
#include "OFstreamMod.H"
#include "interpolation.H"
#include "interpolationCellPoint.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

#define NUMBER_OF_COLUMNS 10

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(tomo2cfd, 0);
    addToRunTimeSelectionTable(functionObject, tomo2cfd, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::functionObjects::tomo2cfd::printInfo() const
{
    Info<<"***************************************************************"<<nl;
    Info<<"Dictionary parameters"<<nl;
    Info<< "field names: " << fieldNames_ << nl;

    Info<< "displacement: " << displacement_ <<nl;

    Info<< "dimensions: " << dimensions_ <<nl;
    Info<<"***************************************************************"<<nl;
}

bool Foam::functionObjects::tomo2cfd::cellInsideTheBox(point& pos) const
{
    bool res = true;

    for(int i=0; i<3; i++)
    {
        res = res  &&  ( pos[i]>=minPoint_[i] && pos[i]<=maxPoint_[i] );
    }

    return res;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::tomo2cfd::tomo2cfd
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_)
{
    if (debug)
    {
        Info<< "functionObjects::tomo2cfd Constructor"<<nl;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::tomo2cfd::~tomo2cfd()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::tomo2cfd::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }
    
    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("tomo2cfd::read")
              << "There is no fields parameter in tomo2cfd dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<vector>("displacement", displacement_) ){
        SeriousErrorIn("::read")
              << "There is no displacement parameter in tomo2cfd dictionary"
              << exit(FatalError);
    }

    
    if( !dict.readIfPresent< List<int> >("dimensions", dimensions_) ){
        SeriousErrorIn("tomo2cfd::read")
              << "There is no minVal parameter in tomo2cfd dictionary"
              << exit(FatalError);
    }
    
    fieldSet_.read(dict);
    
    Info<<"END READ"<<nl<<nl;
    
    printInfo();
    
    Info<<nl<<nl;
    
    return true;
}


bool Foam::functionObjects::tomo2cfd::execute()
{
    if(debug)
    {
        Info<<"Calculate histograms of fields"<<nl;
        printInfo();
    }
    
    //for (const word& fieldName : fieldSet_.selection())
    for (const word& fieldName : fieldNames_)
    {
        integrateFields<scalar>(fieldName);
        integrateFields<vector>(fieldName);
        integrateFields<sphericalTensor>(fieldName);
        integrateFields<symmTensor>(fieldName);
        integrateFields<tensor>(fieldName);
    }
    
    return true;
}

bool Foam::functionObjects::tomo2cfd::write()
{
    Info<<"Currently there is no write into a file"<<nl;
    return true;
}


// ************************************************************************* //
