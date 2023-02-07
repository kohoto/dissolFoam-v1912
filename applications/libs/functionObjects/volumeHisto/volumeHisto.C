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

#include "volumeHisto.H"
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
    defineTypeNameAndDebug(volumeHisto, 0);
    addToRunTimeSelectionTable(functionObject, volumeHisto, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::functionObjects::volumeHisto::printInfo() const
{
    Info<<"***************************************************************"<<nl;
    Info<<"Dictionary parameters"<<nl;
    Info<< "field names: " << fieldNames_ << nl;

    Info<< "min point: " << minPoint_ <<nl;
    Info<< "max point: " << maxPoint_ <<nl;

    Info<< "minVal: " << minVal_ <<nl;
    Info<< "maxVal: " << maxVal_ <<nl;
    Info<< "number of bins: " << numBins_ <<nl;
    Info<<"***************************************************************"<<nl;
}

bool Foam::functionObjects::volumeHisto::cellInsideTheBox(point& pos) const
{
    bool res = true;

    for(int i=0; i<3; i++)
    {
        res = res  &&  ( pos[i]>=minPoint_[i] && pos[i]<=maxPoint_[i] );
    }

    return res;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeHisto::volumeHisto
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
        Info<< "functionObjects::volumeHisto Constructor"<<nl;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeHisto::~volumeHisto()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumeHisto::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }
    
    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("volumeHisto::read")
              << "There is no fields parameter in volumeHisto dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<point>("minPoint", minPoint_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no minPoint parameter in volumeIntegrate dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<point>("maxPoint", maxPoint_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no maxPoint parameter in volumeIntegrate dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<scalar>("minVal", minVal_) ){
        SeriousErrorIn("volumeHisto::read")
              << "There is no minVal parameter in volumeHisto dictionary"
              << exit(FatalError);
    }
    if( !dict.readIfPresent<scalar>("maxVal", maxVal_) ){
        SeriousErrorIn("volumeHisto::read")
              << "There is no maxVal parameter in volumeHisto dictionary"
              << exit(FatalError);
    }
    if( !dict.readIfPresent<int>("nBins", numBins_) ){
        SeriousErrorIn("volumeHisto::read")
              << "There is no nBins parameter in volumeHisto dictionary"
              << exit(FatalError);
    }
    
    if(minVal_>maxVal_)
        SeriousErrorIn("volumeHisto::read")
              << "min value is larger than max value"<< exit(FatalError);
    
    fieldSet_.read(dict);
    
    Info<<"END READ"<<nl<<nl;
    
    printInfo();
    
    Info<<nl<<nl;
    
    return true;
}


bool Foam::functionObjects::volumeHisto::execute()
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

bool Foam::functionObjects::volumeHisto::write()
{
    Info<<"Currently there is no write into a file"<<nl;
    return true;
}


// ************************************************************************* //
