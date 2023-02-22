/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    getSurfCoords

Group
    ?

Description
    Utility to get rough surface coordinates.

    More specifically it:


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "coupledPatchInterpolation.H"

/*
 #######################################################################################
 *    Main program body
 #######################################################################################
*/

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H" // -> runTime
  #include "createDynamicFvMesh.H" // -> mesh.
  
    // Loop over boundary patches
  word solubleWall("solubleWall");
  // const pointField& points = mesh.points();
  
  forAll(mesh.boundary(), patch)
  {
    const word& patchName = mesh.boundary()[patch].name();            // Boundary patch name
    Info<< patchName << endl;
    if (patchName == solubleWall) {
      forAll(mesh.boundary()[patch], facei)
      {
        // const label& point = boundaryMesh[patch].meshPoints()[pointi];  // Node index
        vector pointX(0,0,0);
        const label& faceID = mesh.boundary()[patch].start() + facei;
        
        forAll (mesh.faces()[faceID], nodei) // since the same points are used for many facies, nodes data will have duplicates.
        {
          const label& nodeID = mesh.faces()[faceID][nodei];
          Info<< mesh.points()[nodeID] << endl;
        }
      }
    }
    else
    {
      forAll(mesh.boundary()[patch], facei)
      {
        // const label& point = boundaryMesh[patch].meshPoints()[pointi];  // Node index
        const label& faceID = mesh.boundary()[patch].start() + facei;
        
        forAll (mesh.faces()[faceID], nodei) // since the same points are used for many facies, nodes data will have duplicates.
        {
          const label& nodeID = mesh.faces()[faceID][nodei];
          Info << nodeID << endl;
        }
      }
    }
  }
  return 0;
}


// ************************************************************************* //
