#include "backgroundMesh.H"

namespace Bashyal
{
    void backgroundBlock::intersectClosedSurface(const faceList &faces, const pointField &points, point insidePoint)
    {

        Foam::HashTable<label, point> pointMap;
        Foam::HashTable<label, face> faceMap;

        pointField outputPoints;
        face outputFace;

        // List<List<point>> blockFaceHitPoints;
        // HashTable<pointFaceHit, point> blockHitMap;
        // Loop through each face in cubeAggregate's face list
        for (const face &faceI : faces)
        {
            this->intersectSurfaceFace(faceI, points, outputFace, outputPoints, insidePoint);

        }

        for (int i = 0; i < faces_.size(), i++)
        {
            List<point> blockFaceHitPointsI = blockFaceHitPoints[i];
            this->generateBlockFace(blockFaceHitPointsI, blockHitMap, outputFace, outputPoints, insidePoint);
        }
    }

    void backgroundBlock::intersectBlockFace(List<point> &blockFaceHitPoints, HashTable<pointFaceHit, point> &blockHitMap, const face &faceI, const pointField &points)
    {

    }

    void backgroundBlock::intersectSurfaceFace(const face &faceI, const pointField &points, face &outputFace, pointField &outputPoints, point targetPoint)
    {
        int countVerticesInside = this->countVerticesInsideBlock(faceI, points);
        if (countVerticesInside == faceI.size())
        {
            face reversedFace = faceI.reverseFace();
            outputPoints = reversedFace.points(points);
            outputFace.setSize(reversedFace.size());
            for (label i = 0; i < outputFace.size(); ++i)
            {
                // Set the vertex index directly, assuming ordered sequentially
                outputFace[i] = i;
            }
            return;
        }

        HashTable<pointFaceHit, point> hitMap;

        List<point> pts;

        // Loop through each face in the backgroundBlock
        for (const face &blockFace : this->getFaces())
        {

            // Call findIntersections between faceI and the current blockFace
            this->findIntersections(
                faceI, blockFace,          // Faces to intersect
                points, this->getPoints(), // Points for each face
                pts,                       // Store intersections
                hitMap);
        }

        if (pts.size() == 0)
        {
            return;
        }

        this->generateFace(pts, hitMap, outputFace, outputPoints, targetPoint);

        // Append points to pointMap
        // Append face to faceMap

        // Process intersections if needed, or pass them back to caller if needed
        // (Currently intersections is only collected within this function scope)
    }

    void backgroundBlock::findIntersections(
        const face &face1,
        const face &face2,
        const UList<point> &points1,
        const UList<point> &points2,
        List<point> &pts,
        HashTable<pointFaceHit, point> &hitMap)
    {
        // Use shootSurfaceRays method to shoot rays from face1 to face2 and vice versa
        shootSurfaceRays(face1, points1, face2, points2, pts, hitMap);
        shootBlockRays(face2, points2, face1, points1, pts, hitMap);
    }

    void backgroundBlock::shootSurfaceRays(
        const face &srcFace,
        const UList<point> &srcPoints,
        const face &tgtFace,
        const UList<point> &tgtPoints,
        List<point> &pts,
        HashTable<pointFaceHit, point> &hitMap)
    {
        int countHits = 0;
        for (label i = 0; i < srcFace.size(); ++i)
        {
            // Start point of the ray (vertex of source face)
            point startPoint = srcPoints[srcFace[i]];

            if (this->contains(startPoint))
            {
                if (this->isPointInsideSurface(startPoint))
                {
                    if (!hitMap.found(startPoint))
                    {
                        pointFaceHit intersection1 = pointFaceHit(startPoint);
                        intersection1.setInsideFlag();
                        pts.append(startPoint);
                        hitMap.insert(startPoint, intersection1);
                    }
                }
            }
            // Compute the edge vector (to the next vertex)
            vector edgeVec = srcFace.edge(i, srcPoints);

            // Calculate the edge distance as the magnitude of the edge vector
            scalar edgeDistance = mag(edgeVec);

            // Perform ray intersection with the target face
            pointHit hit = tgtFace.ray(
                startPoint,            // Start point of ray
                edgeVec,               // Direction vector (edge)
                tgtPoints,             // Points of the target face
                intersection::FULL_RAY // Use full-ray intersection
            );

            // If the ray hits the face, check the distance
            if (hit.hit())
            {
                // Calculate distance from the startPoint to the hit point
                scalar hitDistance = hit.distance();

                // Only include the hit if the distance is less than the edge distance
                if (hitDistance >= 0 && hitDistance <= edgeDistance)
                {
                    if (!hitMap.found(hit.point()))
                    {
                        pointFaceHit intersection1 = pointFaceHit(hit, startPoint);
                        intersection1.setBlockFaceFlag();
                        intersection1.setCutEdgeFlag();
                        intersection1.face1_ = this->faces_.found(tgtFace);
                        pts.append(hit.point());
                        hitMap.insert(hit.point(), intersection1);
                        countHits++;
                    }
                }
            }
        }
    }

    void backgroundBlock::shootBlockRays(
        const face &srcFace,
        const UList<point> &srcPoints,
        const face &tgtFace,
        const UList<point> &tgtPoints,
        List<point> &pts,
        HashTable<pointFaceHit, point> &hitMap)
    {
        int countHits = 0;
        for (label i = 0; i < srcFace.size(); ++i)
        {
            // Start point of the ray (vertex of source face)
            point startPoint = srcPoints[srcFace[i]];

            // Compute the edge vector (to the next vertex)
            vector edgeVec = srcFace.edge(i, srcPoints);

            // Calculate the edge distance as the magnitude of the edge vector
            scalar edgeDistance = mag(edgeVec);

            // Perform ray intersection with the target face
            pointHit hit = tgtFace.ray(
                startPoint,            // Start point of ray
                edgeVec,               // Direction vector (edge)
                tgtPoints,             // Points of the target face
                intersection::FULL_RAY // Use full-ray intersection
            );

            // If the ray hits the face, check the distance
            if (hit.hit())
            {
                // Calculate distance from the startPoint to the hit point
                scalar hitDistance = hit.distance();

                // Only include the hit if the distance is less than the edge distance
                if (hitDistance >= 0 && hitDistance <= edgeDistance)
                {
                    if (!hitMap.found(hit.point()))
                    {
                        pointFaceHit intersection1 = pointFaceHit(hit, startPoint);
                        intersection1.setEdgeFlag();
                        intersection1.face1_ = this->faces_.found(srcFace);
                        pts.append(hit.point());
                        hitMap.insert(hit.point(), intersection1);
                        countHits++;
                    }
                    else
                    {
                        pointFaceHit intersection1 = hitMap[hit.point()];
                        intersection1.face2_ = this->faces_.found(srcFace);
                    }
                }
            }
        }
    }

    void backgroundBlock::generateFace(List<point> pts, HashTable<pointFaceHit, point> &hitMap, face &outputFace, pointField &outputPoints, point &targetPoint)
    {
        Info << "Here" << endl;
        if (pts.size() < 3)
        {
            Info << "Insufficient points to form a face." << endl;
            return;
        }
        DynamicList<int> ls;
        ls.setSize(pts.size());

        for (int i = 0; i < pts.size(); i++)
        {
            ls[i] = i;
        }

        outputFace.setSize(pts.size());
        int count = 0;
        int current = 0;

        outputFace[count] = 0;
        count++;
        ls.remove(current);

        while (count < pts.size())
        {
            pointFaceHit hitI = hitMap[pts[current]];
            if (hitI.isInside_ == true)
            {
                for (int j = 0; j < ls.size(); j++)
                {
                    pointFaceHit hitJ = hitMap[pts[ls[j]]];
                    if (pts[current] == hitJ.previousPoint_)
                    {
                        current = ls[j];
                        outputFace[count] = current;
                        ls.remove(j);
                        count++;
                    }
                }
            }
            else if (hitI.isFace_ == true)
            {
                for (int j = 0; j < ls.size(); j++)
                {
                    pointFaceHit hitJ = hitMap[pts[ls[j]]];
                    if (hitI.face1_ == hitJ.face1_ || hitI.face1_ == hitJ.face2_)
                    {
                        current = ls[j];
                        outputFace[count] = current;
                        ls.remove(j);
                        count++;
                    }
                }
            }
            else if (hitI.isEdge_ == true)
            {
                for (int j = 0; j < ls.size(); j++)
                {
                    pointFaceHit hitJ = hitMap[pts[ls[j]]];
                    if (hitI.face1_ == hitJ.face1_ || hitI.face2_ == hitJ.face1_ || hitI.face1_ == hitJ.face2_ || hitI.face2_ == hitJ.face2_)
                    {
                        current = ls[j];
                        outputFace[count] = current;
                        ls.remove(j);
                        count++;
                    }
                }
            }
        }
        outputPoints = pts;
        vector normalVector = outputFace.areaNormal(pts);
        vector directionToTarget = targetPoint - outputPoints[outputFace[0]];

        if (Foam::dot(normalVector, directionToTarget) < 0)
        {
            outputFace = outputFace.reverseFace();
            // Info << "Face order reversed to match the target point direction." << endl;
        }
    }

    int backgroundBlock::countVerticesInsideBlock(const face &face1, const pointField &points) // Attention:: Need to work on it on edited blocks
    {
        int insideCount = 0;

        // Loop through each face in the provided face list
        for (label vertexIndex : face1)
        {
            const point &vertex = points[vertexIndex];

            // Check if the vertex lies within the bounding box of the block
            if (this->bounds_.contains(vertex))
            {
                ++insideCount;
            }
        }

        return insideCount;
    }

    int backgroundBlock::countVerticesInsideSurface(const faceList &faces, const pointField &points) // Attention:: to check
    {
        int insideCount = 0;

        // Loop through each face in the provided face list
        for (const face &face1 : faces)
        {
            // Loop through each vertex index in the current face
            for (label vertexIndex : face1)
            {
                const point &vertex = points[vertexIndex];

                // Ray-casting check: Cast a ray from the vertex in a direction (e.g., along the x-axis)
                int intersections = 0;
                vector rayDir(1, 0, 0); // Choose an arbitrary direction for ray-casting

                for (const face &blockFace : this->faces_)
                {
                    pointHit hit = blockFace.ray(vertex, rayDir, this->points_, intersection::FULL_RAY);

                    // If the ray intersects the face, increase the intersection count
                    if (hit.hit())
                    {
                        ++intersections;
                    }
                }

                // If the intersection count is odd, the point is inside the closed surface
                if (intersections % 2 == 1)
                {
                    ++insideCount;
                }
            }
        }

        return insideCount;
    }

    bool backgroundBlock::isPointInsideSurface(const point &p)
    {
        int intersections = 0;
        vector rayDir(1, 0, 0); // Arbitrary ray direction (e.g., x-axis)

        // Iterate over each face in the surface
        for (const face &blockFace : faces_)
        {
            // Cast a ray from the point p in the rayDir direction
            pointHit hit = blockFace.ray(p, rayDir, points_, intersection::FULL_RAY);

            // If the ray intersects the face, increment intersection count
            if (hit.hit() && hit.distance() >= 0)
            {
                ++intersections;
            }
        }

        // If the number of intersections is odd, the point is inside the surface
        return (intersections % 2 == 1);
    }

    void backgroundBlock::addPoints(const pointField &points, Foam::HashTable<label, point> pointMap, const pointField &mergedPoints)
    {
        for (const auto &pt : points)
        {
            if (!pointMap.found(pt))
            {
                pointMap.insert(pt, mergedPoints.size());
                mergedPoints.append(pt);
            }
        }
    }

}
