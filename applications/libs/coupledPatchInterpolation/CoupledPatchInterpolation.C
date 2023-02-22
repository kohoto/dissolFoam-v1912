/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "CoupledPatchInterpolation.H"
#include "faceList.H"
#include "demandDrivenData.H"

#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Patch>
const scalarListList&
CoupledPatchInterpolation<Patch>::faceToPointWeights() const
{
    if (!faceToPointWeightsPtr_)
    {
        makeFaceToPointWeights();
    }

    return *faceToPointWeightsPtr_;
}

template<class Patch>
const scalarList&
CoupledPatchInterpolation<Patch>::faceToPointSumWeights() const
{
    if (!faceToPointSumWeightsPtr_)
    {
        makeFaceToPointWeights();
    }

    return *faceToPointSumWeightsPtr_;
}

template<class Patch>
void CoupledPatchInterpolation<Patch>::makeFaceToPointWeights() const
{
    if (faceToPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "CoupledPatchInterpolation<Patch>::makeFaceToPointWeights() const"
        )   << "Face-to-edge weights already calculated"
            << abort(FatalError);
    }

    const pointField& points = patch_.localPoints();
    const List<typename Patch::FaceType>& faces = patch_.localFaces();

    faceToPointWeightsPtr_ = new scalarListList(points.size());
    scalarListList& weights = *faceToPointWeightsPtr_;

    // get reference to addressing
    const labelListList& pointFaces = patch_.pointFaces();

    faceToPointSumWeightsPtr_ = new scalarList( pointFaces.size() );
    scalarList& sumw = *faceToPointSumWeightsPtr_;

    forAll(pointFaces, pointi)
    {
        const labelList& curFaces = pointFaces[pointi];

        scalarList& pw = weights[pointi];
        pw.setSize(curFaces.size());

        scalar sumwloc = 0.0;

        forAll(curFaces, facei)
        {
          pw[facei] = 1.0/mag(faces[curFaces[facei]].centre(points) - points[pointi]);
          sumwloc += pw[facei];
        }
        sumw[pointi] = sumwloc;
    }
    
    // synchronization over coupled boundaries
    syncTools::syncPointList( mesh_,
                              patch_.meshPoints(),
                              sumw, 
                              plusEqOp<scalar>(),
                              0.0);
    
    // !! Important faceToPointWeightsPtr_ is not normalized !!
}


template<class Patch>
const scalarList&
CoupledPatchInterpolation<Patch>::faceToEdgeWeights() const
{
    if (!faceToEdgeWeightsPtr_)
    {
        makeFaceToEdgeWeights();
    }

    return *faceToEdgeWeightsPtr_;
}


template<class Patch>
void CoupledPatchInterpolation<Patch>::makeFaceToEdgeWeights() const
{
    if (faceToEdgeWeightsPtr_)
    {
        FatalErrorIn
        (
            "CoupledPatchInterpolation<Patch>::makeFaceToEdgeWeights() const"
        )   << "Face-to-edge weights already calculated"
            << abort(FatalError);
    }

    const pointField& points = patch_.localPoints();
    const List<typename Patch::FaceType>& faces = patch_.localFaces();
    const edgeList& edges = patch_.edges();
    const labelListList& edgeFaces = patch_.edgeFaces();

    faceToEdgeWeightsPtr_ = new scalarList(patch_.nInternalEdges());
    scalarList& weights = *faceToEdgeWeightsPtr_;

    forAll(weights, edgei)
    {
        vector P = faces[edgeFaces[edgei][0]].centre(points);
        vector N = faces[edgeFaces[edgei][1]].centre(points);
        vector S = points[edges[edgei].start()];
        vector e = edges[edgei].vec(points);

        scalar alpha =
            -(((N - P)^(S - P))&((N - P)^e))/(((N - P)^e )&((N - P)^e));

        vector E = S + alpha*e;

        weights[edgei] = mag(N - E)/(mag(N - E) + mag(E - P));
    }
}


template<class Patch>
void CoupledPatchInterpolation<Patch>::clearWeights()
{
  deleteDemandDrivenData(faceToPointWeightsPtr_);
  deleteDemandDrivenData(faceToPointSumWeightsPtr_);
  deleteDemandDrivenData(faceToEdgeWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Patch>
CoupledPatchInterpolation<Patch>::CoupledPatchInterpolation(const Patch& p, const fvMesh& mesh)
:
    patch_(p),
    mesh_(mesh),
    faceToPointWeightsPtr_(NULL),
    faceToPointSumWeightsPtr_(NULL),
    faceToEdgeWeightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class Patch>
CoupledPatchInterpolation<Patch>::~CoupledPatchInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Patch>
template<class Type> // template program
tmp<Field<Type> > CoupledPatchInterpolation<Patch>::faceToPointInterpolate // used in dissolCase
(
    const Field<Type>& ff //cellMotionU for solubleWall
) const
{
    // Check size of the given field
    if (ff.size() != patch_.size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > CoupledPatchInterpolation::"
            "faceToPointInterpolate(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << patch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }
    // Info << "CoupledPatchInterpolation<Patch>::faceToPointInterpolate" << endl;
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
//            patch_.nPoints(), pTraits<Type>::zero
            patch_.nPoints(), Zero
        )
    );

    Field<Type>& result = tresult.ref();

    List<Type> pointValue;
    pointValue.setSize( patch_.nPoints() );
    
    const labelListList& pointFaces = patch_.pointFaces();
    const scalarListList& weights = faceToPointWeights();
    //Type tvalue(0.0, 0.0, 0.0); // If I write like this, this is just taking the last value. So I'm not multiplying vector.
    forAll(pointFaces, pointi)
    {
        const labelList& curFaces = pointFaces[pointi];
        const scalarList& w = weights[pointi];
        
        // Type tvalue = pTraits<Type>::zero; // original
        // Vector<Type> tvalue(0.0, 0.0, 0.0);
        Type tvalue = pTraits<Type>::zero;
        Type fvalue = pTraits<Type>::zero;

        
        // Type a = ff[curFaces[0]].x();
        // Info << "a = " <<  a << endl;
        forAll(curFaces, facei) // facei:0-3
        {
          // Type a = ff[curFaces[facei]];
          //if (isA<vector>(a)) {Info << a.x() << nl;}
          tvalue += w[facei] * ff[curFaces[facei]]; // w: scalar. ff: vector [3 x 1], ff.size() = 19040
          
        //   label n=0;
        //   // vector hello(0.0, 1.0, 2.0);
        //   for(auto&& z : tvalue)
        //   {
        //     if(n==2) continue;
        //     z=0;
        //     n++;
        //   }
          // if (isA<vector>(tvalue)) tvalue.x() = tvalue.y() = scalar(0);
          //tvalue = cmptMultiply(tvalue, fvalue);
        }
        /* put in assign values

        */
        AssignZeros(tvalue);

        pointValue[pointi] = tvalue;
        // Info << "tvalue in main = " << tvalue << nl;
        // Info << "pointValue[pointi] = " <<  pointValue[pointi] << endl; // i noticed that if I set as [pointi, 1], then it shows the same array everytime.
        // If I change it to [i, pointi], then the value changed
        // If I change it to [1, pointi, 1] then all values are the same
        // If I change it to [1, 1, pointi] then the value changed
    }

    // Info << pTraits<Type>::dim <<nl; // 3
    // Info << pTraits<Type>::typeName <<nl; // vector
    Info << "Using CoupledPatchInterpolation.C - Tohoko 3" << nl;    
    syncTools::syncPointList( mesh_, patch_.meshPoints(), pointValue, plusEqOp<Type>(), pTraits<Type>::zero);
    
    // normalization
    const scalarList& sumw = faceToPointSumWeights();
    forAll(pointFaces, pointi)
    {
      result[pointi] = pointValue[pointi] / sumw[pointi]; // pointValue is vector
      // result[pointi].x() = 0.0; // doesn't work
    }

    return tresult;
}

template<class Patch>
void CoupledPatchInterpolation<Patch>::AssignZeros
(
    vector& tvalue
) const
{
    tvalue.x() = 0.0;
    tvalue.y() = 0.0;
    // Info << "tvalue in AssignZero = " << tvalue << nl;
}

template<class Patch>
void CoupledPatchInterpolation<Patch>:: AssignZeros
(
    double& tvalue
) const
{
    // Info << "AssignZeros for scalar" << nl; 
}

template<class Patch>
template<class Type>
tmp<Field<Type> > CoupledPatchInterpolation<Patch>::faceToPointInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = faceToPointInterpolate(tff());
    tff.clear();
    return tint;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > CoupledPatchInterpolation<Patch>::pointToFaceInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != patch_.nPoints())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > CoupledPatchInterpolation::"
            "pointToFaceInterpolate(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << patch_.nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            patch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    const List<typename Patch::FaceType>& localFaces = patch_.localFaces();

    forAll(result, facei)
    {
        const labelList& curPoints = localFaces[facei];

        forAll(curPoints, pointi)
        {
            result[facei] += pf[curPoints[pointi]];
        }

        result[facei] /= curPoints.size();
    }

    return tresult;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > CoupledPatchInterpolation<Patch>::pointToFaceInterpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = pointToFaceInterpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > CoupledPatchInterpolation<Patch>::faceToEdgeInterpolate
(
    const Field<Type>& pf
) const
{
    // Check size of the given field
    if (pf.size() != patch_.size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > CoupledPatchInterpolation::"
            "faceToEdgeInterpolate(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << patch_.size() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>(patch_.nEdges(), pTraits<Type>::zero)
    );

    Field<Type>& result = tresult();

    const edgeList& edges = patch_.edges();
    const labelListList& edgeFaces = patch_.edgeFaces();

    const scalarList& weights = faceToEdgeWeights();

    for (label edgei = 0; edgei < patch_.nInternalEdges(); edgei++)
    {
        result[edgei] =
            weights[edgei]*pf[edgeFaces[edgei][0]]
          + (1.0 - weights[edgei])*pf[edgeFaces[edgei][1]];
    }

    for (label edgei = patch_.nInternalEdges(); edgei < edges.size(); edgei++)
    {
        result[edgei] = pf[edgeFaces[edgei][0]];
    }

    return tresult;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > CoupledPatchInterpolation<Patch>::faceToEdgeInterpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = faceToEdgeInterpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Patch>
bool CoupledPatchInterpolation<Patch>::movePoints()
{
    clearWeights();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
