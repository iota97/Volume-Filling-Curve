// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "sdf.h"

std::vector<std::vector<Vec3>> getIsolines(const SDF &sdf, float spacing, const std::vector<Vec3> &normals, const std::vector<Vec3> &bitangent);