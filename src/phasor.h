// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "sdf.h"

std::vector<float> computePhases(const SDF &sdf, const std::vector<Vec3> &directions, float spacing, bool border = true);