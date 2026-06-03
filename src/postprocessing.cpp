// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "postprocessing.h"
#include "bvhpoint.h"
#include "marchingcube.h"
#include "parallel.h"
#include "marchingtriangle.h"

void postProcess(std::vector<Vec3> &curve, const SDF &sdf, float spacing, size_t iterationsCount) {
    if (!curve.size()) return;

    float growSpeed = 1.25f;
    float smoothness = 1.0f;
    float strength = 0.5f;

    auto triangles = marchingCube([&](size_t x, size_t y, size_t z) { return sdf.getField()[sdf.idx3To1(x, y, z)]; },
                                  {sdf.getSize(0), sdf.getSize(1), sdf.getSize(2)}, sdf.getCellSize(), true);
                                  
    size_t vertsCount = 0;
    for (const auto &t : triangles) {
        vertsCount = std::max(vertsCount, t[0].idx + 1);
        vertsCount = std::max(vertsCount, t[1].idx + 1);
        vertsCount = std::max(vertsCount, t[2].idx + 1);
    }
    std::vector<Vec3> borderVertex(vertsCount);
    for (const auto &t : triangles) {
        borderVertex[t[0].idx] = t[0].v;
        borderVertex[t[1].idx] = t[1].v;
        borderVertex[t[2].idx] = t[2].v;
    }

    for (size_t iter = 0; iter < iterationsCount; ++iter) {
        printProgress(float(iter+1)/iterationsCount);

        // Create a reference frame on the curve.
        std::vector<Vec3> points(borderVertex);
        points.insert(points.end(), curve.begin(), curve.end());
        BVHPoint bvhPoint(points);
        
        // Grow the curve.
        std::vector<Vec3> oldCurve(curve);
        std::vector<Vec3> curveBack(curve.size());
        Parallel::For(0, curve.size(), [&](size_t i) {
            // Create a reference frame on the curve.
            Vec3 p = curve[i];
            size_t prevInd = i ? i - 1 : curve.size() - 1;
            size_t nextInd = i < curve.size()-1 ? i + 1 : 0;
            Vec3 T = (curve[nextInd]-curve[prevInd]).normalize();
            Vec3 N = T.normalToTangent().normalize();
            Vec3 B = N.cross(T).normalize();

            // Get all nearby points.
            float maxRadius = growSpeed*spacing;
            auto neighbors = bvhPoint.pointsInSphereIndex(p, maxRadius);

            // Iterate over the directions orthogonal to the curve.
            size_t directionsCount = 16;
            Vec3 sum(0.0f);
            float weightSum = 0.0f;
            for (size_t dInd = 0; dInd < directionsCount/2; ++dInd) {
                float theta = float(dInd)/(directionsCount/2)*M_PI;
                Vec3 direction = cosf(theta)*B + sinf(theta)*N;

                // Find the closest points in the both the positive and negative direction.
                float distanceNegative = -maxRadius, distancePositive = maxRadius;
                float weightNegative = 1.0f, weightPositive = 1.0f;
                float deltaNegative = distanceNegative + 0.5f*spacing, deltaPositive = distancePositive - 0.5f*spacing;
                for (size_t j : neighbors) {
                    if (i == j) continue;
                    Vec3 q = points[j];
                    float distance = (q-p).length();

                    // Use the tangent distance to define closest as this filters points along the tangent directions.
                    float tangentDistance = distance*distance/direction.dot(q-p);

                    if (tangentDistance > 0 && tangentDistance < distancePositive) {
                        // The spacing with the border has to be halved.
                        bool isClosestABorder = j < borderVertex.size();
                        // Only move half the distance as the other point will likely move too, unless is a border.
                        deltaPositive = isClosestABorder ? distance - 0.5f*spacing : 0.5f*(distance - spacing);
                        // Increase the weight for the boder to create a stronger energy barrier.
                        weightPositive = isClosestABorder ? 4.0f : 1.0f;
                        distancePositive = tangentDistance;
                    } else if (tangentDistance < 0 && tangentDistance > distanceNegative) {
                        bool isClosestABorder = j < borderVertex.size();
                        deltaNegative = isClosestABorder ? -distance + 0.5f*spacing : -0.5f*(distance - spacing);
                        weightNegative = isClosestABorder ? 4.0f : 1.0f;
                        distanceNegative = tangentDistance;
                    }
                }

                // Find the correct position to have the right spacing.
                Vec3 correctNegative = (p + deltaNegative * direction);
                Vec3 correctPositive = (p + deltaPositive * direction);

                // Average the contribution of the varius direction.
                sum += weightNegative*correctNegative;
                sum += weightPositive*correctPositive;
                weightSum += weightNegative;
                weightSum += weightPositive;
            }

            // Update the curve with the average.
            curveBack[i] = sum/weightSum;
        });

        // Smooth the curve.
        float lambda = 1.0f/(smoothness + 1e-6);
        Parallel::For(0, curve.size(), [&](size_t i) {
            Vec3 a = lambda*curveBack[i];
            Vec3 b = curveBack[(i+1)%curve.size()];
            Vec3 c = curveBack[i ? i-1 : curve.size()-1];
            curve[i] = (a+b+c)/(lambda+2.0f);
        });

        // Reduce impact
        for (size_t i = 0; i < curve.size(); ++i) {
            curve[i] = strength*curve[i] + (1-strength)*oldCurve[i];
        }

        // Resample the curve.
        curve = resampleCycle(curve, 0.25f*spacing);
    }

    for (auto &point : curve) {
        point += sdf.getOrigin();
    }
}