// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "stlloader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <set>

StlLoader::StlLoader(const char *path) {
    std::ifstream fileBinary(path, std::ifstream::ate | std::ios::binary);
    size_t fileSize = fileBinary.tellg(); 
    uint32_t triangleCount;
    fileBinary.seekg(0, std::ios::beg);
    const char *headerASCII = "solid";
    char header[5];
    fileBinary.read(header, 5);
    if (memcmp(header, headerASCII, 5) == 0) {
        fileBinary.seekg(80);
        fileBinary.read((char*)&triangleCount, 4);
        if (84 + (triangleCount * 50) != fileSize) {
            triangleCount = uint32_t(-1);
        }
    } else {
        fileBinary.seekg(80);
        fileBinary.read((char*)&triangleCount, 4);
        if (84 + (triangleCount * 50) != fileSize) {
            std::cerr << "Invalid STL" << std::endl;
            return;
        }
    }
    std::unordered_map<Vec3, size_t> verticesMap;
    size_t usedVerts = 0;
    float x, y, z;
    Vec3 vert;
    if (triangleCount == uint32_t(-1)) {
        fileBinary.close();
        std::ifstream fileASCII(path);
        std::string line;
        std::istringstream lineStream;
        while (std::getline(fileASCII, line)) {
            size_t idx = line.find("outer loop");
            if (idx != std::string::npos) {
                std::array<size_t, 3> triangle;
                std::array<Vec3, 3> triVerts;
                for (uint8_t i = 0; i < 3; ++i) {
                    while (std::getline(fileASCII, line) && line.find("vertex") == std::string::npos);
                    lineStream = std::istringstream(line.substr(line.find("vertex")+6));
                    lineStream >> x >> y >> z;
                    vert = Vec3(x, y, z);
                    triVerts[i] = vert;
                    if (verticesMap.count(vert) == 0) {
                        verticesMap.insert({vert, triangle[i] = usedVerts++});
                    } else {
                        triangle[i] = verticesMap[vert];
                    }
                }
                if ((triVerts[0]-triVerts[2]).cross(triVerts[1]-triVerts[2]).length() > std::numeric_limits<float>::min())
                    triangles.push_back(triangle);
            }
        }
        fileASCII.close();
    } else {
        for (size_t i = 0; i < triangleCount; ++i) {
            fileBinary.seekg(84+i*50+12, std::ios::beg);
            std::array<size_t, 3> triangle;
            std::array<Vec3, 3> triVerts;
            for (uint8_t j = 0; j < 3; ++j) {
                fileBinary.read((char*)&x, 4);
                fileBinary.read((char*)&y, 4);
                fileBinary.read((char*)&z, 4);
                vert = Vec3(x, y, z);
                triVerts[j] = vert;
                if (verticesMap.count(vert) == 0) {
                    verticesMap.insert({vert, triangle[j] = usedVerts++});
                } else {
                    triangle[j] = verticesMap[vert];
                }
            }
            if ((triVerts[0]-triVerts[2]).cross(triVerts[1]-triVerts[2]).length() > std::numeric_limits<float>::min())
                triangles.push_back(triangle);
        }
        fileBinary.close();
    }
    vertices.resize(verticesMap.size());
    for (const auto &[key, value] : verticesMap) {
        vertices[value] = key;
    }

    std::set<std::array<Vec3, 3>> triangleSoupMap;
    for (const auto &tri : triangles) {
        triangleSoupMap.insert({vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]});
    }

    for (const auto &tri : triangleSoupMap) {
        triangleSoup.push_back(tri);
    }
}
