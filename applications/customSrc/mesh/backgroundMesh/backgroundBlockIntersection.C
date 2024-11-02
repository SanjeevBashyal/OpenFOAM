#include "backgroundMesh.H"

namespace Bashyal
{
    void backgroundBlock::intersectClosedSurface(const faceList &faces, const pointField &points, point insidePoint)
    {
        pointField mergedPoints;
        faceList mergedFaces;
        Foam::HashTable<label, point> pointMap;
        Foam::HashTable<label, face> faceMap;

        

        List<List<point>> blockFaceHitPoints;
        blockFaceHitPoints.setSize(faces_.size());

        HashTable<pointFaceHit, point> hitMap;

        // Loop through each face in cubeAggregate's face list
        for (int i = 0; i < faces.size(); i++)
        {
            pointField outputPoints;
            outputPoints.clear();
            face outputFace;
            outputFace.clear();

            this->intersectSurfaceFace(faces[i], i, points, outputFace, outputPoints, insidePoint, hitMap, blockFaceHitPoints);
            if (outputPoints.size()>0)
            {
                this->addPoints(outputPoints, pointMap, mergedPoints);
                this->addFace(outputPoints, outputFace, pointMap, mergedFaces, mergedPoints);
            }
            
        }

        for (int i = 0; i < faces_.size(); i++)
        {
            pointField outputPoints;
            outputPoints.clear();
            face outputFace;
            outputFace.clear();

            List<point> blockFaceHitPointsI = blockFaceHitPoints[i];
            this->generateBlockFace(blockFaceHitPointsI, hitMap, i, faces, points, outputFace, outputPoints);
            this->addPoints(outputPoints, pointMap, mergedPoints);
            this->addFace(outputPoints, outputFace, pointMap, mergedFaces, mergedPoints);
        }
        this->points_ = mergedPoints;
        this->faces_ = mergedFaces;
    }

    void backgroundBlock::generateBlockFace(List<point> &pts, HashTable<pointFaceHit, point> &hitMap, int faceIndex, const faceList &faces, const Foam::pointField &points, face &outputFace, pointField &outputPoints)
    {
        face faceI = faces_[faceIndex];

        DynamicList<int> ls;
        ls.setSize(pts.size());

        for (int i = 0; i < pts.size(); i++)
        {
            ls[i] = i;
        }

        int countVerticesInside = this->countVerticesInsideSurface(faceI.points(points_), faces, points);
        if (countVerticesInside == faceI.size())
        {
            outputFace.setSize(0);
            outputPoints.setSize(0);
            return;
        }

        int totalOutputPoints = pts.size() + faceI.size() - countVerticesInside;
        outputFace.setSize(totalOutputPoints);
        outputPoints.setSize(totalOutputPoints);

        int hitIndexer = 0;     // hold how many hits have been incorporated
        int vertexIndexer = 0;  // hold latest index of vertices on check
        int currentIndexer = 0; // hold latest index of outputPoints and outputFace

        bool initFlag = true;
        bool vertexFlag = false;
        bool edgeThroughVertexFlag = false;
        bool edgeThroughFaceFlag = false;
        bool faceFlag = false;

        while (currentIndexer < totalOutputPoints)
        {
            if (initFlag)
            {
                initFlag = false;

                point pt = points_[faceI[vertexIndexer]];
                if (!this->isPointInsideSurface(pt, faces, points))
                {
                    outputPoints[currentIndexer] = pt;
                    outputFace[currentIndexer] = currentIndexer;
                    currentIndexer++;

                    vertexFlag = true;
                    // searchFlag = true;
                }
                else
                {
                    vertexFlag = true;
                }
                // continue;
            }

            if (vertexFlag)
            {
                vertexFlag = false;
                point pt = outputPoints[outputPoints.size() - 1];
                for (int i = 0; i < ls.size(); i++)
                {
                    pointFaceHit hitI = hitMap[pts[ls[i]]];
                    if (hitI.isBlockEdge_ == true)
                    {
                        if (hitI.previousPoint_ == pt)
                        {
                            outputPoints[currentIndexer] = hitI.pt_;
                            outputFace[currentIndexer] = currentIndexer;
                            currentIndexer++;
                            edgeThroughVertexFlag = true;

                            ls.remove(i);
                            hitIndexer++;
                            break;
                        }
                    }
                }
            }

            if (edgeThroughVertexFlag)
            {
                edgeThroughVertexFlag = false;
                pointFaceHit hitPrevious = hitMap[outputPoints[outputPoints.size() - 1]];
                for (int i = 0; i < ls.size(); i++)
                {
                    pointFaceHit hitI = hitMap[pts[ls[i]]];
                    if (hitI.isBlockFace_ == true)
                    {
                        if (hitPrevious.cutFace1_ == hitI.cutFace1_ || hitPrevious.cutFace1_ == hitI.cutFace2_)
                        {
                            outputPoints[currentIndexer] = hitI.pt_;
                            outputFace[currentIndexer] = currentIndexer;
                            currentIndexer++;
                            faceFlag = true;

                            ls.remove(i);
                            hitIndexer++;
                            break;
                        }
                    }
                    else if (hitI.isBlockEdge_ == true)
                    {
                        if (hitPrevious.cutFace1_ == hitI.cutFace1_ || hitPrevious.cutFace1_ == hitI.cutFace2_)
                        {
                            outputPoints[currentIndexer] = hitI.pt_;
                            outputFace[currentIndexer] = currentIndexer;
                            currentIndexer++;
                            edgeThroughFaceFlag = true;

                            ls.remove(i);
                            hitIndexer++;
                        }
                    }
                }
            }

            if (faceFlag)
            {
                faceFlag = false;
                pointFaceHit hitPrevious = hitMap[outputPoints[outputPoints.size() - 1]];
                for (int i = 0; i < ls.size(); i++)
                {
                    pointFaceHit hitI = hitMap[pts[ls[i]]];
                    if (hitI.isBlockFace_ == true)
                    {
                        if (hitPrevious.cutFace1_ == hitI.cutFace1_ || hitPrevious.cutFace1_ == hitI.cutFace2_ || hitPrevious.cutFace2_ == hitI.cutFace1_ || hitPrevious.cutFace2_ == hitI.cutFace2_)
                        {
                            outputPoints[currentIndexer] = hitI.pt_;
                            outputFace[currentIndexer] = currentIndexer;
                            currentIndexer++;
                            faceFlag = true;

                            ls.remove(i);
                            hitIndexer++;
                            break;
                        }
                    }
                    else if (hitI.isBlockEdge_ == true)
                    {
                        if (hitPrevious.cutFace1_ == hitI.cutFace1_ || hitPrevious.cutFace1_ == hitI.cutFace2_ || hitPrevious.cutFace2_ == hitI.cutFace1_ || hitPrevious.cutFace2_ == hitI.cutFace2_)
                        {
                            outputPoints[currentIndexer] = hitI.pt_;
                            outputFace[currentIndexer] = currentIndexer;
                            currentIndexer++;
                            edgeThroughFaceFlag = true;

                            ls.remove(i);
                            hitIndexer++;
                        }
                    }
                }
            }

            if (edgeThroughFaceFlag)
            {
                edgeThroughFaceFlag = false;
                pointFaceHit hitPrevious = hitMap[outputPoints[outputPoints.size() - 1]];
                point prevPoint = hitPrevious.previousPoint_;

                label whichVertex = faces_[faceIndex].found(points_.found(prevPoint));
                vertexIndexer = whichVertex + 1;
                initFlag = true;
            }

            if (!(initFlag || vertexFlag || edgeThroughVertexFlag || faceFlag || edgeThroughFaceFlag))
            {
                vertexIndexer++;
                initFlag = true;
            }
        }
    }

    void backgroundBlock::intersectSurfaceFace(const face &faceI, const int faceIndex, const pointField &points, face &outputFace, pointField &outputPoints, point targetPoint, HashTable<pointFaceHit, point> &hitMap, List<List<point>> &blockFaceHitPoints)
    {
        int countVerticesInside = this->countVerticesInsideSurface(faceI.points(points), faces_, points_);
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

        List<point> pts;
        pts.clear();

        // Loop through each face in the backgroundBlock
        for (const face &blockFace : this->getFaces())
        {

            // Call findIntersections between faceI and the current blockFace
            this->findIntersections(
                faceI, faceIndex, blockFace, // Faces to intersect
                points, this->getPoints(),   // Points for each face
                pts,                         // Store intersections
                hitMap,
                blockFaceHitPoints);
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
        const int faceIndex,
        const face &face2,
        const UList<point> &points1,
        const UList<point> &points2,
        List<point> &pts,
        HashTable<pointFaceHit, point> &hitMap,
        List<List<point>> &blockFaceHitPoints)
    {
        // Use shootSurfaceRays method to shoot rays from face1 to face2 and vice versa
        shootSurfaceRays(face1, faceIndex, points1, face2, points2, pts, hitMap, blockFaceHitPoints);
        shootBlockRays(face2, points2, face1, faceIndex, points1, pts, hitMap, blockFaceHitPoints);
    }

    void backgroundBlock::shootSurfaceRays(
        const face &srcFace,
        const int faceIndex,
        const UList<point> &srcPoints,
        const face &tgtFace,
        const UList<point> &tgtPoints,
        List<point> &pts,
        HashTable<pointFaceHit, point> &hitMap,
        List<List<point>> &blockFaceHitPoints)
    {
        int countHits = 0;
        for (label i = 0; i < srcFace.size(); ++i)
        {
            // Start point of the ray (vertex of source face)
            point startPoint = srcPoints[srcFace[i]];

            if (this->contains(startPoint))
            {
                if (this->isPointInsideSurface(startPoint, this->faces_, this->points_))
                {
                    if (!hitMap.found(startPoint))
                    {
                        pointFaceHit intersection1 = pointFaceHit(startPoint);
                        intersection1.setBlockInsideFlag();
                        pts.append(startPoint);
                        hitMap.insert(startPoint, intersection1);
                    }
                    else
                    {
                        // Alert!! Need to update hitMap data       Hint: Use two start points and the problem will be solved
                        // Alert!! Need to connect inside Vertices  Hint: If blockFace, search for inside point
                        if (!pts.found(startPoint))
                        {
                            pts.append(startPoint);
                        }
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

                        label blockFaceIndex = this->faces_.found(tgtFace);
                        intersection1.blockFace1_ = blockFaceIndex;
                        intersection1.cutFace1_ = faceIndex;

                        pts.append(hit.point());
                        blockFaceHitPoints[blockFaceIndex].append(hit.point());
                        hitMap.insert(hit.point(), intersection1);
                        countHits++;
                    }
                    else
                    {
                        pointFaceHit intersection1 = hitMap[hit.point()];
                        intersection1.cutFace2_ = faceIndex;

                        label blockFaceIndex = this->faces_.found(tgtFace);
                        if (!blockFaceHitPoints[blockFaceIndex].found(hit.point()))
                        {
                            blockFaceHitPoints[blockFaceIndex].append(hit.point());
                        }
                        

                        if (!pts.found(hit.point()))
                        {
                            pts.append(hit.point());
                        }
                    }
                }
            }
        }
    }

    void backgroundBlock::shootBlockRays(
        const face &srcFace,
        const UList<point> &srcPoints,
        const face &tgtFace,
        const int faceIndex,
        const UList<point> &tgtPoints,
        List<point> &pts,
        HashTable<pointFaceHit, point> &hitMap,
        List<List<point>> &blockFaceHitPoints)
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
                        intersection1.setBlockEdgeFlag();
                        intersection1.setCutFaceFlag();

                        label blockFaceIndex = this->faces_.found(srcFace);
                        intersection1.blockFace1_ = blockFaceIndex;
                        intersection1.cutFace1_ = faceIndex;

                        pts.append(hit.point());
                        hitMap.insert(hit.point(), intersection1);
                        blockFaceHitPoints[blockFaceIndex].append(hit.point());
                        countHits++;
                    }
                    else
                    {
                        pointFaceHit intersection1 = hitMap[hit.point()];
                        label blockFaceIndex = this->faces_.found(srcFace);
                        intersection1.blockFace2_ = blockFaceIndex;

                        if (!blockFaceHitPoints[blockFaceIndex].found(hit.point()))
                        {
                            blockFaceHitPoints[blockFaceIndex].append(hit.point());
                        }

                        if (!pts.found(hit.point()))
                        {
                            pts.append(hit.point());
                        }
                    }
                }
            }
        }
    }

    void backgroundBlock::generateFace(List<point> pts, HashTable<pointFaceHit, point> &hitMap, face &outputFace, pointField &outputPoints, point &targetPoint)
    {
        // Info << "Here" << endl;
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
            if (hitI.isBlockInside_ == true)
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
                        break;
                    }
                }
            }
            else if (hitI.isBlockFace_ == true)
            {
                for (int j = 0; j < ls.size(); j++)
                {
                    pointFaceHit hitJ = hitMap[pts[ls[j]]];
                    if (hitI.blockFace1_ == hitJ.blockFace1_ || hitI.blockFace1_ == hitJ.blockFace2_)
                    {
                        current = ls[j];
                        outputFace[count] = current;
                        ls.remove(j);
                        count++;
                        break;
                    }
                }
            }
            else if (hitI.isBlockEdge_ == true)
            {
                for (int j = 0; j < ls.size(); j++)
                {
                    pointFaceHit hitJ = hitMap[pts[ls[j]]];
                    if (hitI.blockFace1_ == hitJ.blockFace1_ || hitI.blockFace2_ == hitJ.blockFace1_ || hitI.blockFace1_ == hitJ.blockFace2_ || hitI.blockFace2_ == hitJ.blockFace2_)
                    {
                        current = ls[j];
                        outputFace[count] = current;
                        ls.remove(j);
                        count++;
                        break;
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

    int backgroundBlock::countVerticesInsideBlockBounds(const face &face1, const pointField &points) // Attention:: Need to work on it on edited blocks
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

    int backgroundBlock::countVerticesInsideSurface(const pointField &pts, const faceList &faces, const pointField &points) // Attention:: to check
    {
        int insideCount = 0;

        // Loop through each face in the provided face list
        for (const point &pt : pts)
        {
            if (this->isPointInsideSurface(pt, faces, points))
            {
                ++insideCount;
            }
        }

        return insideCount;
    }

    bool backgroundBlock::isPointInsideSurface(const point &p, const faceList &faces, const pointField &points)
    {
        int intersections = 0;
        vector rayDir(1, 0, 0); // Arbitrary ray direction (e.g., x-axis)

        // Iterate over each face in the surface
        for (const face &blockFace : faces)
        {
            // Cast a ray from the point p in the rayDir direction
            pointHit hit = blockFace.ray(p, rayDir, points, intersection::FULL_RAY);

            // If the ray intersects the face, increment intersection count
            if (hit.hit() && hit.distance() >= 0)
            {
                ++intersections;
            }
        }

        // If the number of intersections is odd, the point is inside the surface
        return (intersections % 2 == 1);
    }

    void backgroundBlock::addPoints(const pointField &points, Foam::HashTable<label, point> &pointMap, pointField &mergedPoints)
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

    void backgroundBlock::addFace(const pointField &points, const face &faceI, Foam::HashTable<label, point> &pointMap, faceList &mergedFaces, pointField &mergedPoints)
    {

        // Adjust face points according to the mergedPoints
        Foam::face globalFace;
        for (const auto &pt : faceI)
        {
            // Retrieve the point coordinates using pt as an index into blockPoints
            const Foam::point &pointCoords = points[pt]; // Using blockPoints to get the actual coordinates

            // Add to pointMap_ if not already present, and retrieve the global index
            label globalPointIdx;
            if (pointMap.found(pointCoords))
            {
                globalPointIdx = pointMap[pointCoords];
            }
            else
            {
                // If not found, add the point and assign it a new global index
                globalPointIdx = mergedPoints.size();
                mergedPoints.append(pointCoords);
                pointMap.insert(pointCoords, globalPointIdx);
            }

            // Add global point index to the globalFace
            globalFace.append(globalPointIdx);
        }
        mergedFaces.append(globalFace);
    }

}
