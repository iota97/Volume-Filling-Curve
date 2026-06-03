// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "mathutils.h"
#include <vector>
#include <array>

class SDF {
public:
    SDF(const std::vector<std::array<Vec3, 3>> &triangles, float cellSize, uint32_t border = 2, uint32_t raysCount = 5, bool quiet = false);
    SDF(const char* path);

    void saveToBin(const char *path) const;
    void saveSliceToPLY(const char *path) const;
    std::array<uint32_t, 3> getSize() const { return size; }
    uint32_t getSize(uint32_t i) const { return size[i]; }
    float getCellSize() const { return cellSize; }
    const std::vector<float> &getField() const { return field; }
    std::vector<float> &getField() { return field; }
    Vec3 getPosition(size_t i) const;
    size_t getIndex(const Vec3 &position) const;
    float getVal(const Vec3 &position) const;
    Vec3 getGrad(const Vec3 &position) const;
    Vec3 getGrad(size_t i) const;
    float getCurvature(const Vec3 &position) const;
    Vec3 getGrad(const std::array<uint32_t, 3> &xyz) const;
    float getVal(size_t x, size_t y, size_t z) const { return field[idx3To1(x, y, z)]; }
    size_t idx3To1(size_t x, size_t y, size_t z) const;
    size_t idx3To1(const std::array<uint32_t, 3> &c) const;
    std::array<uint32_t, 3> idx1To3(size_t i) const;
    Vec3 getOrigin() const { return origin; }

private:
    float weight(Vec3 pos, size_t cellIndex) const;
    float sdfVal(Vec3 pos, size_t cellIndex) const;

    std::vector<float> field;
    std::array<uint32_t, 3> size;
    float cellSize;
    Vec3 origin;
};