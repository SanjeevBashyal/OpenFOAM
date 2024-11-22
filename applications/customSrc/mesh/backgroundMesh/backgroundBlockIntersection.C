#include "backgroundMesh.H"

namespace Bashyal
{
    void backgroundBlock::intersectClosedSurface(const faceList &faces, const pointField &points, point insidePoint, word stringPtr)
    {
        pointField mergedPoints;
        faceList mergedFaces;
        wordList mergePatches;
        wordList mergeStringPtrs;

        Foam::HashTable<label, point> pointMap;
        Foam::HashTable<label, face> faceMap;

        List<List<point>> cutFaceHitPoints;
        cutFaceHitPoints.setSize(faces.size());

        List<List<point>> blockFaceHitPoints;
        blockFaceHitPoints.setSize(faces_.size());

        HashTable<pointFaceHit, point> hitMap;

        // Loop through each face in cubeAggregate's face list
        for (int i = 0; i < faces.size(); i++)
        {
            this->intersectSurfaceFace(faces[i], i, points, insidePoint, hitMap, blockFaceHitPoints, cutFaceHitPoints[i]);
        }

        // Generation of Intersected Cut Face
        for (int i = 0; i < faces.size(); i++)
        {
            pointField outputPoints;
            outputPoints.clear();
            face outputFace;
            outputFace.clear();

            List<point> cutFaceHitPointsI = cutFaceHitPoints[i];

            if (cutFaceHitPointsI.size() > 0)
            {
                this->generateCutFace(cutFaceHitPointsI, hitMap, i, outputFace, outputPoints, insidePoint);
                this->addPoints(outputPoints, pointMap, mergedPoints);
                this->addFace(outputPoints, outputFace, pointMap, mergedFaces, mergedPoints);
                mergePatches.append("aggregate");
                mergeStringPtrs.append(stringPtr);
            }
        }

        // Generation of Intersected Block Face
        for (int i = 0; i < faces_.size(); i++)
        {
            pointField outputPoints;
            outputPoints.clear();
            face outputFace;
            outputFace.clear();

            List<point> blockFaceHitPointsI = blockFaceHitPoints[i];
            this->generateBlockFace(blockFaceHitPointsI, hitMap, i, faces, points, outputFace, outputPoints);
            if (outputPoints.size() > 2)
            {
                this->addPoints(outputPoints, pointMap, mergedPoints);
                this->addFace(outputPoints, outputFace, pointMap, mergedFaces, mergedPoints);
                mergePatches.append(patches_[i]);
                mergeStringPtrs.append(word::null);
            }
        }
        this->points_ = mergedPoints;
        this->faces_ = mergedFaces;
        this->patches_ = mergePatches;
        this->stringPtrs_ = mergeStringPtrs;

        if (this->points_.size() == 0)
        {
            this->ncells_ = 0;
        }
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
                }
                // continue;
            }

            if (vertexFlag)
            {
                vertexFlag = false;
                point pt = points_[faceI[vertexIndexer]];
                for (int i = 0; i < ls.size(); i++)
                {
                    pointFaceHit hitI = hitMap[pts[ls[i]]];
                    if (hitI.isBlockEdge_ == true)
                    {
                        point rayPoint;
                        if (hitI.blockFace1_ == faceIndex)
                        {
                            rayPoint = hitI.rayVertex1_;
                        }
                        else if (hitI.blockFace2_ == faceIndex)
                        {
                            rayPoint = hitI.rayVertex2_;
                        }
                        if (this->arePointsSame(rayPoint, pt, 1e-6))
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
                pointFaceHit hitPrevious = hitMap[outputPoints[currentIndexer - 1]];
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
                pointFaceHit hitPrevious = hitMap[outputPoints[currentIndexer - 1]];
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
                point prevPoint;
                edgeThroughFaceFlag = false;
                pointFaceHit hitPrevious = hitMap[outputPoints[currentIndexer - 1]];
                if (hitPrevious.blockFace1_ == faceIndex)
                {
                    prevPoint = hitPrevious.rayVertex1_;
                }
                else if ((hitPrevious.blockFace2_ == faceIndex))
                {
                    prevPoint = hitPrevious.rayVertex2_;
                }

                label whichVertex = faces_[faceIndex].find(points_.find(prevPoint));
                vertexIndexer = whichVertex + 1;
                initFlag = true;
            }

            if (!(initFlag || vertexFlag || edgeThroughVertexFlag || faceFlag || edgeThroughFaceFlag))
            {
                vertexIndexer++;
                initFlag = true;
                if (vertexIndexer >= faceI.size() && currentIndexer < totalOutputPoints)
                {
                    // pointFaceHit hitI = hitMap[pts[ls[0]]];
                    // outputPoints[currentIndexer] = hitI.pt_;
                    // outputFace[currentIndexer] = currentIndexer;
                    // currentIndexer++;
                    // faceFlag = true;

                    // ls.remove(0);
                    // hitIndexer++;
                    break;
                }
            }
        }
    }

    void backgroundBlock::intersectSurfaceFace(const face &faceI, const int faceIndex, const pointField &points, point targetPoint, HashTable<pointFaceHit, point> &hitMap, List<List<point>> &blockFaceHitPoints, List<point> &cutFaceHitPointsI)
    {
        // int countVerticesInside = this->countVerticesInsideSurface(faceI.points(points), faces_, points_);
        // if (countVerticesInside == faceI.size())
        // {
        //     face reversedFace = faceI.reverseFace();
        //     outputPoints = reversedFace.points(points);
        //     outputFace.setSize(reversedFace.size());
        //     for (label i = 0; i < outputFace.size(); ++i)
        //     {
        //         // Set the vertex index directly, assuming ordered sequentially
        //         outputFace[i] = i;
        //     }
        //     return;
        // }

        cutFaceHitPointsI.clear();

        // Loop through each face in the backgroundBlock
        for (const face &blockFace : this->getFaces())
        {

            // Call findIntersections between faceI and the current blockFace
            this->findIntersections(
                faceI, faceIndex, blockFace, // Faces to intersect
                points, this->getPoints(),   // Points for each face
                cutFaceHitPointsI,           // Store intersections
                hitMap,
                blockFaceHitPoints);
        }

        if (cutFaceHitPointsI.size() == 0)
        {
            return;
        }

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

                        // Alert!! Need to check for OR || can do problem
                        // Alert!! Need mesh and face refiner for smaller intersections
                        if (!this->foundOnList(pts, startPoint, 1e-6))
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
                    // if (!hitMap.found(hit.point()))
                    point hitPt;
                    if (!this->pointOnList(hitMap.toc(), hit.point(), hitPt, 1e-6))
                    {
                        pointFaceHit intersection1 = pointFaceHit(hit, startPoint);
                        intersection1.setBlockFaceFlag();
                        intersection1.setCutEdgeFlag();

                        label blockFaceIndex = this->faces_.find(tgtFace);
                        intersection1.blockFace1_ = blockFaceIndex;
                        intersection1.cutFace1_ = faceIndex;

                        pts.append(hit.point());
                        blockFaceHitPoints[blockFaceIndex].append(hit.point());
                        hitMap.insert(hit.point(), intersection1);
                        countHits++;
                    }
                    else
                    {
                        pointFaceHit &intersection1 = hitMap[hitPt];
                        intersection1.rayVertex2_ = startPoint;
                        intersection1.cutFace2_ = faceIndex;

                        label blockFaceIndex = this->faces_.find(tgtFace);
                        if (!this->foundOnList(blockFaceHitPoints[blockFaceIndex], hitPt, 1e-6))
                        {
                            blockFaceHitPoints[blockFaceIndex].append(hitPt);
                        }

                        if (!this->foundOnList(pts, hitPt, 1e-6))
                        {
                            pts.append(hitPt);
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
                    // if (!hitMap.found(hit.point()))
                    point hitPt;
                    if (!this->pointOnList(hitMap.toc(), hit.point(), hitPt, 1e-6))
                    {
                        pointFaceHit intersection1 = pointFaceHit(hit, startPoint);
                        intersection1.setBlockEdgeFlag();
                        intersection1.setCutFaceFlag();

                        label blockFaceIndex = this->faces_.find(srcFace);
                        intersection1.blockFace1_ = blockFaceIndex;
                        intersection1.cutFace1_ = faceIndex;

                        pts.append(hit.point());
                        hitMap.insert(hit.point(), intersection1);
                        blockFaceHitPoints[blockFaceIndex].append(hit.point());
                        countHits++;
                    }
                    else
                    {
                        pointFaceHit &intersection1 = hitMap[hitPt];
                        label blockFaceIndex = this->faces_.find(srcFace);
                        intersection1.blockFace2_ = blockFaceIndex;
                        intersection1.rayVertex2_ = startPoint;

                        if (!this->foundOnList(blockFaceHitPoints[blockFaceIndex], hitPt, 1e-6))
                        {
                            blockFaceHitPoints[blockFaceIndex].append(hitPt);
                        }

                        if (!this->foundOnList(pts, hitPt, 1e-6))
                        {
                            pts.append(hitPt);
                        }
                    }
                }
            }
        }
    }

    void backgroundBlock::generateCutFace(List<point> pts, HashTable<pointFaceHit, point> &hitMap, int faceIndex, face &outputFace, pointField &outputPoints, point &targetPoint)
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
                    if (this->arePointsSame(pts[current], hitJ.rayVertex1_, 1e-6) || this->arePointsSame(pts[current], hitJ.rayVertex2_, 1e-6))
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
                    if (hitJ.isBlockInside_)
                    {
                        if (this->arePointsSame(hitI.rayVertex1_, hitJ.pt_, 1e-6) || this->arePointsSame(hitI.rayVertex2_, hitJ.pt_, 1e-6))
                        {
                            current = ls[j];
                            outputFace[count] = current;
                            ls.remove(j);
                            count++;
                            break;
                        }
                    }
                    else if (hitI.blockFace1_ == hitJ.blockFace1_ || hitI.blockFace1_ == hitJ.blockFace2_)
                    {
                        current = ls[j];
                        outputFace[count] = current;
                        ls.remove(j);
                        count++;
                        break;
                    }
                    else if (this->arePointsSame(hitI.rayVertex1_, hitJ.rayVertex1_, 1e-6) && this->arePointsSame(hitI.rayVertex2_, hitJ.rayVertex2_, 1e-6))
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

    bool backgroundBlock::pointOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, Foam::point &outputPoint, const double tolerance)
    {
        for (const Foam::point &pt : pointList)
        {
            // Calculate the distance between the points
            double distance = Foam::mag(pt - checkPoint);

            // Check if the distance is within the specified tolerance
            if (distance <= tolerance)
            {
                outputPoint = pt;
                return true;
            }
        }
        return false;
    }

    bool backgroundBlock::foundOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, const double tolerance)
    {
        for (const Foam::point &pt : pointList)
        {
            // Calculate the distance between the points
            double distance = Foam::mag(pt - checkPoint);

            // Check if the distance is within the specified tolerance
            if (distance <= tolerance)
            {
                return true;
            }
        }
        return false;
    }

    bool backgroundBlock::arePointsSame(const Foam::point &point1, const Foam::point &point2, const double tolerance)
    {
        double distance = Foam::mag(point1 - point2);

        // Check if the distance is within the specified tolerance
        if (distance <= tolerance)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

}
