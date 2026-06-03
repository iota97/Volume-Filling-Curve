// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "stitcher.h"
#include "parallel.h"
#include "bvhpointtagged.h"
#include "marchingtriangle.h"
#include <iostream>
#include <fstream>

void saveCurveToPLY(const char *path, const std::vector<Vec3> &curve) {
    size_t vertsCount = curve.size();
    std::vector<std::array<uint8_t, 3>> colors(vertsCount);
    double totalLength = 0;
    for (size_t i = 0; i < vertsCount; ++i) {
        totalLength += (curve[i] - curve[(i+1) % vertsCount]).length();
    }
    double length = 0;
    for (size_t i = 0; i < vertsCount; ++i) {
        length += (curve[i] - curve[(i+1) % vertsCount]).length();
        colors[i] = turbo(length/totalLength);
    }

    std::ofstream plyFile(path);
    plyFile << "ply\nformat ascii 1.0\nelement vertex " << vertsCount << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            << "element edge " << vertsCount << "\nproperty int vertex1\nproperty int vertex2\n"
            << "end_header\n";
    
    for (size_t i = 0; i < vertsCount; ++i) {
        auto v = curve[i];
        plyFile << v.x << " " << v.y << " " << v.z << " " 
            << uint32_t(colors[i][0]) << " "  << uint32_t(colors[i][1]) << " " << uint32_t(colors[i][2]) << std::endl;
    }

    if (curve.size()) {
        for (size_t i = 0; i < curve.size()-1; ++i) {
            plyFile << i << " " << i+1 << std::endl;
        }
        plyFile << curve.size()-1 << " " << 0 << std::endl;
    }
}

std::vector<Vec3> stitch(std::vector<std::vector<Vec3>> cycles, float spacing) {
    if (!cycles.size()) return std::vector<Vec3>();
    float searchRadius = 2.5f;
    float startNumber = cycles.size();
    while (cycles.size() > 1) {
        printProgress((startNumber-cycles.size()+2)/startNumber);
        std::vector<Vec3Tagged> cyclePoints;
        for (size_t i = 0; i < cycles.size(); ++i) {
            for (size_t j = 0; j < cycles[i].size(); ++j) {
                cyclePoints.push_back({cycles[i][j], i, j});
            }
        }
        BVHPointTagged bvh(cyclePoints);

        size_t shortestCycle = 0;
        size_t shortestSize = cycles[0].size();
        for (size_t i = 1; i < cycles.size(); ++i) {
            if (cycles[i].size() < shortestSize) {
                shortestSize = cycles[i].size();
                shortestCycle = i;
            }
        }

        struct Stitch {
            size_t a0[2];
            size_t a1[2];
            size_t b0[2];
            size_t b1[2];
        };
        Stitch stitch{};

        float shortestDistance = std::numeric_limits<float>::infinity();
        for (size_t i = 0; i < shortestSize; ++i) {
            size_t idxA = i;
            size_t idxB = (i+4) % shortestSize;
            Vec3 posA = cycles[shortestCycle][idxA];
            Vec3 posB = cycles[shortestCycle][idxB];
            Vec3Tagged nearestA = bvh.nearestPointDifferent(posA, shortestCycle, searchRadius*spacing);
            if (nearestA.cycleID == size_t(-1)) continue;

            size_t nearestBIdx0 = (nearestA.vertsID + 4) % cycles[nearestA.cycleID].size();
            size_t nearestBIdx1 = nearestA.vertsID < 4 ? cycles[nearestA.cycleID].size()-4 : nearestA.vertsID-4;
            size_t nearestBIdx = (posB-cycles[nearestA.cycleID][nearestBIdx0]).length2() < (posB-cycles[nearestA.cycleID][nearestBIdx1]).length2() ? nearestBIdx0 : nearestBIdx1;
            Vec3Tagged nearestB = { cycles[nearestA.cycleID][nearestBIdx], nearestA.cycleID, nearestBIdx };

            float distance = (posA-nearestA.v).length() + (posB-nearestB.v).length();
            if (distance < shortestDistance) {
                shortestDistance = distance;
                stitch = {{shortestCycle, idxA}, {shortestCycle, idxB}, {nearestA.cycleID, nearestA.vertsID}, {nearestB.cycleID, nearestB.vertsID}};
            }
        }

        if (std::isinf(shortestDistance)) {
            //std::cout << "Warning: could not stitch, increasing search radius." << std::endl;
            searchRadius *= 2.0f;
            continue;
        }

        std::vector<Vec3> newCycle;
        if (stitch.a0[1] < stitch.a1[1]) {
            //a1 -> max
            for (size_t i = stitch.a1[1]; i < cycles[stitch.a1[0]].size(); ++i) {
                newCycle.push_back(cycles[stitch.a1[0]][i]);
            }
            //0 -> a0
            for (size_t i = 0; i <= stitch.a0[1]; ++i) {
                newCycle.push_back(cycles[stitch.a0[0]][i]);
            }
        } else {
            // a1 -> a0
            for (size_t i = stitch.a1[1]; i <= stitch.a0[1]; ++i) {
                newCycle.push_back(cycles[stitch.a1[0]][i]);
            }
        }

        size_t distInc = stitch.b0[1] < stitch.b1[1] ? stitch.b1[1] - stitch.b0[1] : cycles[stitch.b0[0]].size() - stitch.b0[1] + stitch.b1[1];
        size_t distDec = stitch.b0[1] > stitch.b1[1] ? stitch.b0[1] - stitch.b1[1] : cycles[stitch.b0[0]].size() - stitch.b1[1] + stitch.b0[1];
        if (distInc > distDec) {
            if (stitch.b1[1] < stitch.b0[1]) {
                // b0 -> max
                for (size_t i = stitch.b0[1]; i < cycles[stitch.b0[0]].size(); ++i) {
                    newCycle.push_back(cycles[stitch.b0[0]][i]);
                }
                //0 -> b1
                for (size_t i = 0; i <= stitch.b1[1]; ++i) {
                    newCycle.push_back(cycles[stitch.b1[0]][i]);
                }
            } else {
                // b0 -> b1
                for (size_t i = stitch.b0[1]; i <= stitch.b1[1]; ++i) {
                    newCycle.push_back(cycles[stitch.b1[0]][i]);
                }
            }
        } else {
            if (stitch.b1[1] > stitch.b0[1]) {
                // b0 -> 0
                for (size_t i = stitch.b0[1]; i >= 0; --i) {
                    newCycle.push_back(cycles[stitch.b0[0]][i]);
                    if (!i) break;
                }
                // max -> b1
                for (size_t i = cycles[stitch.b1[0]].size()-1; i >= stitch.b1[1]; --i) {
                    newCycle.push_back(cycles[stitch.b1[0]][i]);
                    if (!i) break;
                }
            } else {
                // b0 -> b1
                for (size_t i = stitch.b0[1]; i >= stitch.b1[1]; --i) {
                    newCycle.push_back(cycles[stitch.b1[0]][i]);
                    if (!i) break;
                }
            }

        }

        cycles.push_back(newCycle);
        cycles.erase(std::next(cycles.begin(), stitch.a0[0]));
        cycles.erase(std::next(cycles.begin(), stitch.b0[0] < stitch.a0[0] ? stitch.b0[0] : stitch.b0[0]-1));
    }

    return resampleCycle(cycles[0], 0.25f*spacing);
}