// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include <vector>
#include <functional>
#include "mathutils.h"

std::vector<std::array<Vec3Idx, 2>> marchingTriangle(const std::function<float(Vec3)> &value, std::vector<std::array<Vec3Idx, 3>> triangles, float quiet = false);
std::vector<std::vector<Vec3>> connectCycles(const std::vector<std::array<Vec3Idx, 2>> &segments, float spacing);
std::vector<Vec3> resampleCycle(const std::vector<Vec3> &cycle, float spacing);

void saveCyclesToPLY(const char *path, const std::vector<std::vector<Vec3>> &cycles, const Vec3 &offset);
void saveSegmentToPLY(const char *path, const std::vector<std::array<Vec3Idx, 2>> &segments);