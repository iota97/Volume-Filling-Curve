// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include "stlloader.h"
#include "sdf.h"
#include "biphasor.h"
#include "direction.h"
#include "stitcher.h"
#include "parallel.h"
#include "postprocessing.h"
#include "tangent.h"

int main (int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " [mesh.stl] [spacing] [curve.ply] <tangent_field> <normal_field>" << std::endl;
        std::cout << "Tangent fields:" << std::endl;
        std::cout << "    0 -> Constant x-direction [default]" << std::endl;
        std::cout << "    1 -> Constant y-direction" << std::endl;
        std::cout << "    2 -> Constant z-direction" << std::endl;
        std::cout << "    3 -> Two constant piecewise directions" << std::endl;
        std::cout << "    4 -> Radial 1" << std::endl;
        std::cout << "    5 -> Radial 2" << std::endl;
        std::cout << "    6 -> Othogonal to the mesh boundary" << std::endl;
        std::cout << "    7 -> Parallel to the mesh boundary" << std::endl;
        std::cout << "    8 -> Random" << std::endl;
        std::cout << "    9 -> Fancy" << std::endl;
        std::cout << "Normal fields (Figure 17):" << std::endl;
        std::cout << "    0 -> Smoothest [default]" << std::endl;
        std::cout << "    1 -> Parallel to boundary" << std::endl;
        exit(1);
    }

    const char *meshPath = argv[1];
    float spacing = std::atof(argv[2]);
    size_t mode = argc > 4 ? std::atoi(argv[4]) : 0;
    bool parallel = argc > 5 && std::atoi(argv[5]);

    auto startTime = std::chrono::steady_clock::now();
    std::cout << "STL Loading..." << std::endl;
    StlLoader mesh(meshPath);
    printProgress(1.0f);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;
    
    startTime = std::chrono::steady_clock::now();
    std::cout << "\nSDF Generation..." << std::endl;
    SDF sdf(mesh.getTriangleSoup(), 0.25f*spacing);
    std::cout << "Grid size: [" << sdf.getSize(0) << ", " << sdf.getSize(1) << ", " << sdf.getSize(2) << "]" << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;
    auto startTimeTotal = std::chrono::steady_clock::now();

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nTangents Generation..." << std::endl;
    std::vector<Vec3> tangents = generateTangents(sdf, mode);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nFrame Generation..." << std::endl;
    std::vector<Vec3> normals;
    std::vector<Vec3> bitangents;
    generateDirections(sdf, tangents, bitangents, normals, parallel);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    auto cycles = getIsolines(sdf, spacing, normals, bitangents);

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nStitching..." << std::endl;
    auto curve = stitch(cycles, spacing);
    //saveCurveToPLY("../data/stitched.ply", curve);

    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nPost Processing..." << std::endl;
    postProcess(curve, sdf, spacing, 8);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;

    startTime = std::chrono::steady_clock::now();
    std::cout << "\nSaving Curve..." << std::endl;
    saveCurveToPLY(argv[3], curve);
    printProgress(1.0f);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()/1000.0 << " s" << std::endl;
    std::cout << "\nTotal time (excluding SDF generation): " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTimeTotal).count()/1000.0 << " s" << std::endl;

    return 0;
}