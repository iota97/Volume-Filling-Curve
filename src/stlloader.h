// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "mathutils.h"
#include <vector>
#include <array>

class StlLoader {
public:
    StlLoader(const char *path);

    const std::vector<Vec3> &getVertices() const { return vertices; }
    const std::vector<std::array<size_t, 3>> &getTriangles() const { return triangles; }
    const std::vector<std::array<Vec3, 3>> &getTriangleSoup() const { return triangleSoup; }
    
private:
    std::vector<Vec3> vertices;
    std::vector<std::array<size_t, 3>> triangles;
    std::vector<std::array<Vec3, 3>> triangleSoup;
};