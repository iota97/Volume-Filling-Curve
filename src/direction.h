// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "sdf.h"

void generateDirections(const SDF &sdf, const std::vector<Vec3> &tangents, std::vector<Vec3> &bitangents, std::vector<Vec3> &normals, bool parallel = true);