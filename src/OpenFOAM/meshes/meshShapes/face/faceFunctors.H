#include "vector.H"
#include "intersection.H"

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

struct faceRayFunctor
{
    const point* meshPoints;
    const intersection::algorithm alg;
    const intersection::direction dir;

    faceRayFunctor
    (
        const point* _meshPoints,
        const intersection::algorithm _alg,
        const intersection::direction _dir
    ):
        meshPoints(_meshPoints),
        alg(_alg),
        dir(_dir)
    {}

    __host__ __device__
    pointHit operator()(const faceData& tgtFace, const point& p, const vector& n)
    {
        const label start = tgtFace.start();
        const label nPoints = tgtFace.size();

        // If the face is a triangle, do a direct calculation
        if (nPoints == 3)
        {
            return triPointRef
            (
                meshPoints[start],
                meshPoints[start + 1],
                meshPoints[start + 2]
            ).ray(p, n, alg, dir);
        }

        point ctr(0, 0, 0);
        for (label i = 0; i < tgtFace.size(); i++)
        {
            ctr += meshPoints[start + i];
        }
        ctr /= tgtFace.size();

        scalar nearestHitDist = GREAT;

        scalar nearestMissDist = GREAT;
        bool eligible = false;

        // Initialize to miss, distance = GREAT
        pointHit nearest(p);

        point nextPoint = ctr;

        for (label pI = 0; pI < nPoints; pI++)
        {
            nextPoint = meshPoints[start + (pI + 1) % nPoints];

            // Note: for best accuracy, centre point always comes last
            //
            pointHit curHit = triPointRef
            (
                meshPoints[start + pI],
                nextPoint,
                ctr
            ).ray(p, n, alg, dir);

            if (curHit.hit())
            {
                if (Foam::mag(curHit.distance()) < Foam::mag(nearestHitDist))
                {
                    nearestHitDist = curHit.distance();
                    nearest.setHit();
                    nearest.setPoint(curHit.hitPoint());
                }
            }
            else if (!nearest.hit())
            {
                // Miss and no hit yet. Update miss statistics.
                if (curHit.eligibleMiss())
                {
                    eligible = true;

                    // Miss distance is the distance between the plane intersection
                    // point and the nearest point of the triangle
                    scalar missDist =
                        Foam::mag
                        (
                            p + curHit.distance()*n
                            - curHit.missPoint()
                        );

                    if (missDist < nearestMissDist)
                    {
                        nearestMissDist = missDist;
                        nearest.setDistance(curHit.distance());
                        nearest.setPoint(curHit.missPoint());
                    }
                }
            }
        }

        if (nearest.hit())
        {
            nearest.setDistance(nearestHitDist);
        }
        else
        {
            // Haven't hit a single face triangle
            nearest.setMiss(eligible);
        }

        return nearest;
    }
};

struct faceNormalFunctor
{
    const point* meshPoints;

    faceNormalFunctor
    (
        const point* _meshPoints
    ):
        meshPoints(_meshPoints)
    {}

    __host__ __device__
    vector operator()(const faceData& face)
    {
        const label start = face.start();
        const label nPoints = face.size();

        if (nPoints == 3)
        {
            return triPointRef
            (
                meshPoints[start],
                meshPoints[start + 1],
                meshPoints[start + 2]
            ).normal();
        }

        label pI;

        point centrePoint(0,0,0);
        for (pI = 0; pI < nPoints; ++pI)
        {
            centrePoint += meshPoints[start + pI];
        }
        centrePoint /= nPoints;

        vector n(0,0,0);

        point nextPoint = centrePoint;

        for (pI = 0; pI < nPoints; ++pI)
        {
            nextPoint = meshPoints[start + (pI + 1) % nPoints];

            n += triPointRef
            (
                meshPoints[start + pI],
                nextPoint,
                centrePoint
            ).normal();
        }

        return n;
    }
};
}
