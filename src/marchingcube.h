// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include <vector>
#include <functional>
#include "mathutils.h"

std::vector<std::array<Vec3Idx, 3>> marchingCube(const std::function<float(size_t, size_t, size_t)> &value, std::array<size_t, 3> domainSize, float cellSize, float quiet = false);

void saveTriangleToPLY(const char *path, const  std::vector<std::array<Vec3Idx, 3>> &triangles, const Vec3 &offset);

void saveTriangleToViridisPLY(const char *path, const  std::vector<std::array<Vec3Idx, 3>> &triangles, const Vec3 &offset, const std::function<float(Vec3)> &value);