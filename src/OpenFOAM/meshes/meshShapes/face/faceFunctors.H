#include "vector.H"

namespace Foam
{

struct faceCentreFunctor
{
    const label* labels;
    const point* points;

    faceCentreFunctor
    (
        const label* _labels,
        const point* _points
    ):
        labels(_labels),
        points(_points)
    {}

    __host__ __device__
    vector operator()(const faceData& face) const
    {
	// Calculate the centre by breaking the faceData into triangles and
        // area-weighted averaging their centres

        const label start = face.start();
        const label nPoints = face.size();

        // If the faceData is a triangle, do a direct calculation
        if (nPoints == 3)
        {
            return
                (1.0/3.0)
               *(
                   points[labels[start]]
                 + points[labels[start+1]]
                 + points[labels[start+2]]
                );
        }


        point centrePoint(0,0,0);
        for (label pI=0; pI<nPoints; ++pI)
        {
            centrePoint += points[labels[pI+start]];
        }
        centrePoint /= nPoints;

        scalar sumA = 0;
        vector sumAc(0,0,0);

        for (label pI=0; pI<nPoints; ++pI)
        {
            const point& nextPoint = points[labels[((pI + 1) % nPoints)+start]];

            // Calculate 3*triangle centre
            const vector ttc
            (
                points[labels[pI+start]]
              + nextPoint
              + centrePoint
            );

            // Calculate 2*triangle area
            const scalar ta = Foam::mag
            (
                (points[labels[pI+start]] - centrePoint)
              ^ (nextPoint - centrePoint)
            );

            sumA += ta;
            sumAc += ta*ttc;
        }

        if (sumA > VSMALL)
        {
            return sumAc/(3.0*sumA);
        }
        else
        {
            return centrePoint;
        }
    }
};


struct faceSweptVolFunctor
{
    const label* labels;
    const point* oldPoints;
    const point* newPoints;
    const faceCentreFunctor oldCentre;
    const faceCentreFunctor newCentre;

    faceSweptVolFunctor
    (
        const label* _labels,
        const point* _oldPoints,
        const point* _newPoints
    ):
        labels(_labels),
        oldPoints(_oldPoints),
        newPoints(_newPoints),
        oldCentre(_labels,_oldPoints),
        newCentre(_labels,_newPoints)
    {}

    __host__ __device__
    scalar operator()(const faceData& face)
    {
        const label start = face.start();
        const label nPoints = face.size();

        scalar sv = 0;

        // Calculate the swept volume by breaking the faceData into triangles and
        // summing their swept volumes.
        // Changed to deal with small concavity by using a central decomposition

        point centreOldPoint = oldCentre(face);
        point centreNewPoint = newCentre(face);

        for (register label pi=0; pi<nPoints-1; ++pi)
        {
            // Note: for best accuracy, centre point always comes last
            sv += triPointRef
            (
                centreOldPoint,
                oldPoints[labels[pi+start]],
                oldPoints[labels[pi+start + 1]]
            ).sweptVol
            (
                triPointRef
                (
                    centreNewPoint,
                    newPoints[labels[pi+start]],
                    newPoints[labels[pi+start + 1]]
                )
            );
        }

        sv += triPointRef
        (
            centreOldPoint,
            oldPoints[labels[start+nPoints-1]],
            oldPoints[labels[start]]
        ).sweptVol
        (
            triPointRef
            (
                centreNewPoint,
                newPoints[labels[start+nPoints-1]],
                newPoints[labels[start]]
            )
        );

        return sv;
    }
};

}
