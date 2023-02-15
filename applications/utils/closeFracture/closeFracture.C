/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
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

Application
  closeFracture

Description
  Read GSLIB output and assign to the mesh.

  Limitations now are the next:
    - a geometry should include one or two flat plates

Usage
  - surfRoughGen

Needs dictionary
  system/closeFractureDict

\*---------------------------------------------------------------------------*/

#include <complex>
#include <vector>
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "coupledPatchInterpolation.H"
#include "pointPatchField.H"

#include "normalMotionSlipPointPatchVectorField.H"
#include "fixedValuePointPatchField.H"

#include "velocityMotionSolver.H"

#include "motionSolver.H"

/*
 #######################################################################################
 *    Main program body
 #######################################################################################
*/

class closeFracture
{
protected:
  // Protected data
  int seed;
  int M;
  int N;
  double majLen;
  double minLen;
  double rgh;
  double dHurst; // Fractal dimension D = 3 - dHurst
  double cutLen;
  double maxDisp;
  word wayToApply;

public:
  // Constructor
  closeFracture(
      int seed_,
      int M_,
      int N_,
      double majLen_,
      double minLen_,
      double rgh_,
      double dHurst,
      double cutLen_,
      double maxDisp_,
      word wayToApply_)
      : seed(seed_),
        M(M_),
        N(N_),
        majLen(majLen_),
        minLen(minLen_),
        rgh(rgh_),
        dHurst(dHurst),
        cutLen(cutLen_),
        maxDisp(maxDisp_),
        wayToApply(wayToApply_)
  {
  }

  void getSurfaceDisplacement( // (mesh, faceDisp, patchID, majDir, minDir)
      dynamicFvMesh &mesh,     // in?
      scalarField &wd,         // inout: zeros before
      label &patchID,          // in
      int majDir,              // in
      int minDir               // in
  )
  {
    scalarField sFn(M * N, 0.0);
    fftDisp(sFn); // calculations are here
    // sFn: displacement
    scalarField sFp(M * N, 0.0);

    sFp = sFn; // sFn: displacement

    pointField pointFace = mesh.boundaryMesh()[patchID].faceCentres();
    pointField normlFace = mesh.boundaryMesh()[patchID].faceNormals();
    scalar maxMaj = max(pointFace.component(majDir));
    scalar maxLat = max(pointFace.component(minDir));
    scalar minMaj = min(pointFace.component(majDir));
    scalar minLat = min(pointFace.component(minDir));
    double Llat = maxLat - minLat;
    double Lmaj = maxMaj - minMaj;

    int mveDir = 3 - (minDir + majDir);

    forAll(pointFace, i)
    {
      // what is sign?
      scalar sign = normlFace[i].component(mveDir) /
                    mag(normlFace[i].component(mveDir));

      // get ind
      scalar curMaj = pointFace[i].component(majDir) - minMaj;
      scalar curLat = pointFace[i].component(minDir) - minLat;
      int curm = std::floor(curMaj / Lmaj * (M - 1));
      int curn = std::floor(curLat / Llat * (N - 1));
      int ind = curn + N * curm;

      if (wayToApply == "synchronous") // my surface will be always synchronous
        wd[i] = sign * sFn[ind];
    }

    Info << "Displacement loaded" << nl;
  }

private:
  // converts indexes
  label index(label m, label n) { return n + N * m; }

  scalar power(double ksq)
  {
    if (ksq == 0)
      return 0; // <rad^2> ~ 1/ksq^(1+H)
    if (ksq > 1)
      return 0; // cutoff wavelength = cutLen
    scalar p = Foam::pow(ksq, -(dHurst + 1));
    //    p *= Foam::exp(-ksq);
    return std::sqrt(p);
  }

  void fftDisp(scalarField &disp) // compute displacement
  {
    unsigned int MN = M * N;
    std::vector<std::complex<double>> f, F;
    f.resize(MN);
    F.resize(MN);

    Random rnd(seed); // seed is global var in this solver
    scalar TwoPi = constant::mathematical::twoPi;

    Info << "Displacement calc starts...." << nl;
    /*
     *   --- ---
     *  | 1 | 2 |
     *   --- ---
     *  | 3 | 4 |
     *   --- ---
     */
    // calculating 1 and 4
    for (int m = 0; m < M / 2 + 1; ++m)
    {
      for (int n = 0; n < N / 2 + 1; ++n)
      {
        scalar p = TwoPi * rnd.sample01<scalar>();
        scalar rad;
        double majk, mink, ksq;
        majk = m * cutLen / majLen;
        mink = n * cutLen / minLen;
        ksq = majk * majk + mink * mink;
        rad = power(ksq) * rnd.GaussNormal<scalar>();

        f[index(m, n)] =
            rad * std::complex<double>(Foam::cos(p), Foam::sin(p));
        f[index(((M - m) % M), (N - n) % N)] =
            rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }
    f[index(M / 2, 0)].imag(0.0);
    f[index(0, N / 2)].imag(0.0);
    f[index(M / 2, N / 2)].imag(0.0);

    for (int m = 1; m < M / 2; ++m)
    {
      for (int n = 1; n < N / 2; ++n)
      {
        scalar p = TwoPi * rnd.sample01<scalar>();
        scalar rad;
        double majk, mink, ksq;
        majk = TwoPi * m / majLen;
        mink = TwoPi * n / minLen;
        ksq = majk * majk + mink * mink;
        rad = power(ksq) * rnd.GaussNormal<scalar>();

        f[index(m, N - n)] =
            rad * std::complex<double>(Foam::cos(p), Foam::sin(p));
        f[index(M - m, n)] =
            rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }

    scalarField sF(MN);
    forAll(sF, ii)
    {
      sF[ii] = F[ii].real();
    }

    scalarField sF2 = sqr(sF);
    scalar avSF = average(sF);
    scalar avSF2 = average(sF2);

    scalar factor = rgh / Foam::sqrt(mag(avSF2 - sqr(avSF)));

    sF *= factor;

    disp = sF;
  }

}; // end of class closeFracture

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createDynamicFvMesh.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // reading dictionary closeFractureDict
  IOdictionary closeFractureDict // generate dictionary data class object
      (
          IOobject(
              "closeFractureDict",
              runTime.system(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE));

  word wayToApply;
  if (!closeFractureDict.readIfPresent<word>("apply", wayToApply))
  {
    SeriousErrorIn("main")
        << "There is no `synchronous` parameter in dictionary"
        << exit(FatalError);
  }
  word patchName;
  if (!closeFractureDict.readIfPresent<word>("patchName", patchName))
  {
    SeriousErrorIn("main")
        << "There is no `patchName` parameter in dictionary"
        << exit(FatalError);
  }

  int majDir;
  if (!closeFractureDict.readIfPresent<int>("majDir", majDir))
  {
    SeriousErrorIn("main")
        << "There is no `majDir` parameter in dictionary"
        << exit(FatalError);
  }
  int minDir;
  if (!closeFractureDict.readIfPresent<int>("minDir", minDir))
  {
    SeriousErrorIn("main")
        << "There is no `minDir` parameter in dictionary"
        << exit(FatalError);
  }
  int majLen;
  if (!closeFractureDict.readIfPresent<int>("majLen", majLen))
  {
    SeriousErrorIn("main")
        << "There is no `majLen` parameter in dictionary"
        << exit(FatalError);
  }
  int minLen;
  if (!closeFractureDict.readIfPresent<int>("minLen", minLen))
  {
    SeriousErrorIn("main")
        << "There is no `minLen` parameter in dictionary"
        << exit(FatalError);
  }
  int M;
  if (!closeFractureDict.readIfPresent<int>("majNum", M))
  {
    SeriousErrorIn("main")
        << "There is no `majNum` parameter in dictionary"
        << exit(FatalError);
  }
  int N;
  if (!closeFractureDict.readIfPresent<int>("minNum", N))
  {
    SeriousErrorIn("main")
        << "There is no `minNum` parameter in dictionary"
        << exit(FatalError);
  }

  int seed;
  if (!closeFractureDict.readIfPresent<int>("seed", seed))
  {
    SeriousErrorIn("main")
        << "There is no `seed` parameter in dictionary"
        << exit(FatalError);
  }
  scalar rgh;
  if (!closeFractureDict.readIfPresent<scalar>("roughness", rgh))
  {
    SeriousErrorIn("main")
        << "There is no `roughness` parameter in dictionary"
        << exit(FatalError);
  }
  double dHurst;
  if (!closeFractureDict.readIfPresent<double>("dHurst", dHurst))
  {
    SeriousErrorIn("main")
        << "There is no `dHurst` parameter in dictionary"
        << exit(FatalError);
  }
  double cutLen;
  if (!closeFractureDict.readIfPresent<double>("cutLen", cutLen))
  {
    SeriousErrorIn("main")
        << "There is no `cutLen` parameter in dictionary"
        << exit(FatalError);
  }
  double maxDisp;
  if (!closeFractureDict.readIfPresent<double>("maxDisp", maxDisp))
  {
    SeriousErrorIn("main")
        << "There is no `maxDisp` parameter in dictionary"
        << exit(FatalError);
  }

  Info << "patch:         " << patchName << endl;

  label patchID = mesh.boundaryMesh().findPatchID(patchName);
  const pointField &boundaryPoints = mesh.boundaryMesh()[patchID].localPoints();
  int i2 = 0;
  for (int m = 0; m < M; ++m)
  {
    for (int n = 0; n < N; ++n) // for same xs.
    {
      i2 += 1;
    }
    Info << boundaryPoints[i2][2] << endl; // this has to be the same as the first one in closureFracture
  }
  /////?????????????????????????????????/
  
  double cpuTime = runTime.elapsedCpuTime();

  ////////////////////////// Repeat the procedure after this for the mirrored plane. /////////////////////////////

  // Get patch ID for moving boundaries
  //label patchID = mesh.boundaryMesh().findPatchID(patchName);

  if (patchID == -1)
  {
    SeriousErrorIn("main")
        << "patch " << patchName << " is missing" << exit(FatalError);
  }

  label patchID_mirrored = mesh.boundaryMesh().findPatchID(patchName + "_mirrored");

  coupledPatchInterpolation patchInterpolator(mesh.boundaryMesh()[patchID_mirrored], mesh);

  // get coordinates of moving surfaces
  //const pointField &boundaryPoints = mesh.boundaryMesh()[patchID].localPoints(); // this is the coordinates of points of {patchID}.
  const pointField &boundaryPoints_mirrored = mesh.boundaryMesh()[patchID_mirrored].localPoints();

  // get apatures by {y+} - {y-}
  vectorField pointDispWall(boundaryPoints.size(), vector::zero);
  vectorField pointNface = mesh.boundaryMesh()[patchID].faceNormals();
  vector motionN(0.0, 0.0, -1.0);

  scalarField faceDisp(pointNface.size(), 0.0);
  scalarField pointDisp = patchInterpolator.faceToPointInterpolate(faceDisp); // I think I can make pointDisp with zeros, not by using faceDisp, but I'm just keeping it as is

  scalarField aperture(N, 0.0);
  double maxAperture = 0.0; // negative, mag is max
  int i = 0;
  int j = 0;
  for (int m = 0; m < M; ++m)
  {
    for (int n = 0; n < N; ++n) // for same xs.
    {
      aperture[n] = boundaryPoints[i][2] - boundaryPoints_mirrored[i][2];
      i += 1;
    }
    maxAperture = max(aperture); // since aperture will be all negative, max is the min in magnitude
    for (int n = 0; n < N; ++n) // for same xs.
    {
      pointDisp[j] = maxAperture;
      j += 1;
    }
  }

  Info << "Maximum and minimum point displacements  "
       << max(pointDisp) << "  " << min(pointDisp) << endl; // went to until this point

  // TODO: !!! Run DICPCG somewhere in between here
  forAll(pointDispWall, i) pointDispWall[i] = pointDisp[i] * motionN;
  //forAll(motionN, i) Info << motionN[i] << endl;

  pointVectorField &pointVelocity = const_cast<pointVectorField &>(
      mesh.objectRegistry::lookupObject<pointVectorField>("pointMotionU"));

  pointVelocity.boundaryFieldRef()[patchID] == pointDispWall; // update rough surface // move only the main surface

  mesh.update(); // if mesh is updated, return true // this line is not working!!!!!!

  cpuTime = runTime.elapsedCpuTime() - cpuTime;
  // !!! End Run DICPCG somewhere in between here
  Info << nl << "Time statistics:" << nl;

  // int wlNP = mesh.boundaryMesh()[patchID].nPoints();
  Info << "Total number of points:                   " << mesh.nPoints() << nl;
  // Info << "Number of points on the walls:            " << wlNP << nl;
  Info << "Running time:                             " << cpuTime << nl << endl;

  Info << "Overwriting points in current time directory." << nl;
  runTime.writeNow();

  ////////////////////////// END Repeat the procedure after this for the mirrored plane. /////////////////////////////

  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
