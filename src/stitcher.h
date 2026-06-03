// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "mathutils.h"
#include <vector>

std::vector<Vec3> stitch(std::vector<std::vector<Vec3>> cycles, float spacing);

void saveCurveToPLY(const char *path, const std::vector<Vec3> &curve);