/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::tomo2cfd

Group
    grpFieldFunctionObjects

Description
    Calculates the tomo2cfd, the curl of the velocity.

    The field is stored on the mesh database so that it can be retrieved
    and used for other applications.

Usage
    \verbatim
    tomo2cfd1
    {
        type        tomo2cfd;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: tomo2cfd    | yes         |
        U            | Name of velocity field  | no          | U
        result       | Name of Courant number field | no    | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    tomo2cfd.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_tomo2cfd_H
#define functionObjects_tomo2cfd_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "point.H"
#include "scalar.H"
#include "scalarField.H"
#include "interpolation.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class tomo2cfd Declaration
\*---------------------------------------------------------------------------*/

class tomo2cfd
:
    public fvMeshFunctionObject
{
protected:
    // Protected member data
        wordList fieldNames_;

        vector displacement_;

        //- Fields to integrate
        volFieldSelection fieldSet_;

        List<int>  dimensions_;

    // Protected Member Functions

        void printInfo() const;
        bool cellInsideTheBox(point& pos) const;
        
        template<class Type>
        bool foundObject(const word& name, const bool verbose = true) const;
    
        //- Calculate the field integral
        template<class Type>
        void interpolateFields
        (
            const word& fieldName
        );
        
private:
        //- Disallow default bitwise copy construct
        tomo2cfd(const tomo2cfd&);
        
        //- Disallow default bitwise assignment
        void operator=(const tomo2cfd&);

public:

    //- Runtime type information
    TypeName("tomo2cfd");

    // Constructors

        //- Construct from Time and dictionary
        tomo2cfd
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~tomo2cfd();
    
    // Member Functions

        //- Read the fieldExpression data
        virtual bool read(const dictionary&);

        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();

        //- Clear the result field from the objectRegistry
        //virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "tomo2cfdTemplates.C"
#endif

#endif

// ************************************************************************* //
