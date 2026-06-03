// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "direction.h"
#include "grid.h"

void generateDirections(const SDF &sdf, const std::vector<Vec3> &tangents, std::vector<Vec3> &bitangents, std::vector<Vec3> &normals, bool parallel) {
    normals.resize(tangents.size());
    bitangents.resize(tangents.size());

    struct Data {
        Vec3 normal{ 0, 0, 1 };
        Vec3 tangent{ 1, 0, 0 };
    };

    GridHierarchy<Data> grid({sdf.getSize(0), sdf.getSize(1), sdf.getSize(2)}, sdf.getCellSize(), Vec3(sdf.getOrigin()));
    Parallel::For(0, tangents.size(), [&](size_t i) {
       grid.getData(i)->normal = sdf.getGrad(i);
       grid.getData(i)->tangent = tangents[i];
    });

    auto smooth = [&](size_t i, GridHierarchy<Data> *grid) {
        Vec3 avg(0);
        for (size_t n : grid->getNeighborhood(i)) {
            Vec3 d = grid->getData(n)->normal;
            avg += (d.dot(avg) < 0 ? -1 : 1)*d;
        }
        grid->getData(i)->normal = avg.normalize();

        if (fabsf(grid->getData(i)->tangent.dot(grid->getData(i)->normal)) > 0.99f) {
            grid->getData(i)->normal = grid->getData(i)->tangent.normalToTangent().normalize();
        }
        Vec3 bitangent = grid->getData(i)->normal.cross(grid->getData(i)->tangent).normalize();
        grid->getData(i)->normal = grid->getData(i)->tangent.cross(bitangent).normalize();
    };

    auto prolong = [&](size_t i, GridHierarchy<Data> *gridDense, GridHierarchy<Data> *gridCoarse) {
        size_t n = gridDense->denseToCoarse(i);
        gridDense->getData(i)->normal = gridCoarse->getData(n)->normal;
    };

    auto restrictFunc = [&](size_t i, GridHierarchy<Data> *gridDense, GridHierarchy<Data> *gridCoarse) {
        Vec3 avg(0);
        for (size_t n : gridCoarse->coarseToDense(i)) {
            Vec3 d = gridDense->getData(n)->tangent;
            avg += (d.dot(avg) < 0 ? -1 : 1)*d;
        }
        gridCoarse->getData(i)->tangent = avg.normalize();
        gridCoarse->getData(i)->normal = gridCoarse->getData(i)->tangent.normalToTangent().normalize();
    };

    if (parallel) {
        grid.smooth(smooth, 8, false);
    } else {
        grid.optimize(smooth, restrictFunc, prolong);
    }

    Parallel::For(0, tangents.size(), [&](size_t i) {
       normals[i] = grid.getData(i)->normal.normalize();
       bitangents[i] = tangents[i].cross(normals[i]).normalize();
    });
}