// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "sdf.h"

void postProcess(std::vector<Vec3> &curve, const SDF &sdf, float spacing, size_t iterationsCount = 8);