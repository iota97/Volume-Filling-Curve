// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "phasor.h"
#include "grid.h"

std::vector<float> computePhases(const SDF &sdf, const std::vector<Vec3> &directions, float spacing, bool border) {
    struct Data {
        Vec3 normal{ 0,0,1 };
        float phase{};
        bool frozen{};
        bool outside{};
    };

    GridHierarchy<Data> grid({sdf.getSize(0), sdf.getSize(1), sdf.getSize(2)}, sdf.getCellSize(), Vec3(sdf.getOrigin()));

    Parallel::For(0, grid.getTotalSize(), [&](size_t i) {
        float sdfVal = sdf.getField()[i];
        Vec3 sdfGrad = sdf.getGrad(grid.idx1To3(i));
        float cellDiag = 0.375f*spacing*sqrt(3);
        bool isBoundary = sdfVal < 0.5f*cellDiag && sdfVal > -0.5f*cellDiag;
        grid.getData(i)->normal = directions[i].normalize();
        grid.getData(i)->frozen = border && isBoundary;
        grid.getData(i)->outside = sdfVal >= 0.0f;
        grid.getData(i)->phase = M_PI*sdfVal/spacing * (directions[i].dot(sdfGrad) < 0 ? 1 : -1);
        if (sdfVal >= 0.0f) grid.getData(i)->phase = 0;
    });

    auto restrictFunc = [&](size_t i, GridHierarchy<Data> *gridDense, GridHierarchy<Data> *gridCoarse) {
        float x = 0, y = 0;
        bool frozen = false;
        bool outside = true;
        Vec3 norm(0);

        for (size_t n : gridCoarse->coarseToDense(i)) {
            if (gridDense->getData(n)->outside) continue;
            norm += (norm.dot(gridDense->getData(n)->normal) < 0 ? -1 : 1) * gridDense->getData(n)->normal;
        }
        gridCoarse->getData(i)->normal = norm.normalize();
            
        for (size_t n : gridCoarse->coarseToDense(i)) {
            if (gridDense->getData(n)->outside) continue;
            outside = false;
            if (gridDense->getData(n)->frozen) {
                frozen = true;
                Vec3 d_i = gridCoarse->getData(i)->normal;
                Vec3 d_j = gridDense->getData(n)->normal;
                float f_i = 0.5f/spacing;
                float f_j = 0.5f/spacing;
                alignPhases(gridCoarse->getPosition(i), gridDense->getPosition(n), d_i, d_j, f_i, f_j, gridDense->getData(n)->phase, x, y);
            }
        }

        gridCoarse->getData(i)->phase = atan2f(y, x);
        gridCoarse->getData(i)->frozen = frozen;
        gridCoarse->getData(i)->outside = outside;
    };

    auto prolong = [&](size_t i, GridHierarchy<Data> *gridDense, GridHierarchy<Data> *gridCoarse) {
        if (gridDense->getData(i)->outside || gridDense->getData(i)->frozen) return;
        size_t n = gridDense->denseToCoarse(i);
        Vec3 d_i = gridDense->getData(i)->normal;
        Vec3 d_j = gridCoarse->getData(n)->normal;
        float f_i = 0.5f/spacing;
        float f_j = 0.5f/spacing;
        float x = 0, y = 0;
        alignPhases(gridDense->getPosition(i), gridCoarse->getPosition(n), d_i, d_j, f_i, f_j, gridCoarse->getData(n)->phase, x, y);
        gridDense->getData(i)->phase = atan2f(y, x);
    };

    auto smooth = [&](size_t i, GridHierarchy<Data> *grid) {
        if (grid->getData(i)->frozen || grid->getData(i)->outside) return;
        float x = 0, y = 0;
        for (size_t n : grid->getNeighborhood(i)) {
            if (grid->getData(n)->outside) continue;
            Vec3 d_i = grid->getData(i)->normal;
            Vec3 d_j = grid->getData(n)->normal;
            float f_i = 0.5f/spacing;
            float f_j = 0.5f/spacing;
            alignPhases(grid->getPosition(i), grid->getPosition(n), d_i, d_j, f_i, f_j, grid->getData(n)->phase, x, y);
        }
        grid->getData(i)->phase = atan2f(y, x);
    };

    grid.optimize(smooth, restrictFunc, prolong);

    std::vector<float> phases(grid.getTotalSize());
    Parallel::For(0, grid.getTotalSize(), [&](size_t i) {
        phases[i] = grid.getData(i)->phase;
    });
    return phases;
}

