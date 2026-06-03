// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "biphasor.h"
#include "marchingcube.h"
#include "marchingtriangle.h"
#include "phasor.h"
#include <chrono>
#include <iostream>

std::vector<std::vector<Vec3>> getIsolines(const SDF &sdf, float spacing, const std::vector<Vec3> &normals, const std::vector<Vec3> &bitangent) {
    auto startTime = std::chrono::steady_clock::now();
    std::cout << "\nNormal Phasor..." << std::endl;
    std::vector<float> phasesNormal = computePhases(sdf, normals, spacing, true);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nBitangent Phasor..." << std::endl;
    std::vector<float> phasesBitangent = computePhases(sdf, bitangent, spacing, false);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nMarching Cubes..." << std::endl;
    float cellSize = sdf.getCellSize();
    std::array<size_t, 3> domainSize{sdf.getSize(0), sdf.getSize(1), sdf.getSize(2)};
    //std::swap(phasesNormal, phasesBitangent);
    auto valueNormal = [&](size_t x, size_t y, size_t z) { 
        return cosf(phasesNormal[sdf.idx3To1(x, y, z)]);
    };
    auto triangles = marchingCube(valueNormal, domainSize, cellSize);
    
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nMarching Triangles..." << std::endl;

    auto weight = [&] (Vec3 pos, size_t cellIndex) {
        return (std::max(0.0f, cellSize - (pos-sdf.getPosition(cellIndex)).length()))/cellSize;
    };

    auto valueBitangent = [&](const Vec3 &pos) {
        float val = 0, wSum = 0;
        auto xyz = sdf.idx1To3(sdf.getIndex(pos + sdf.getOrigin()));
        uint32_t minX = xyz[0] > 0 ? xyz[0]-1 : xyz[0];
        uint32_t maxX = xyz[0] < sdf.getSize(0)-1 ? xyz[0]+1 : xyz[0];
        uint32_t minY = xyz[1] > 0 ? xyz[1]-1 : xyz[1];
        uint32_t maxY = xyz[1] < sdf.getSize(1)-1 ? xyz[1]+1 : xyz[1];
        uint32_t minZ = xyz[2] > 0 ? xyz[2]-1 : xyz[2];
        uint32_t maxZ = xyz[2] < sdf.getSize(2)-1 ? xyz[2]+1 : xyz[2];

        for (uint32_t x = minX; x <= maxX; ++x) {
            for (uint32_t y = minY; y <= maxY; ++y) {
                for (uint32_t z = minZ; z <= maxZ; ++z) {
                    float w = weight(pos + sdf.getOrigin(), sdf.idx3To1({x, y, z}));
                    val += w*cosf(phasesBitangent[sdf.idx3To1({x, y, z})]);
                    wSum += w;
                }
            }
        }
        return val/wSum;
    };
    auto segments = marchingTriangle(valueBitangent, triangles);
    auto cycles = connectCycles(segments, 0.25f*spacing);

    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    //saveTriangleToViridisPLY("../data/surfaces.ply", triangles, sdf.getOrigin(), valueBitangent);
    //saveCyclesToPLY("../data/lines.ply", cycles, sdf.getOrigin());

    return cycles;
}