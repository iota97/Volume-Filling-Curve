// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "tangent.h"
#include "grid.h"

std::vector<Vec3> generateTangents(const SDF &sdf, size_t mode) {
    std::vector<Vec3> tangents(sdf.getSize(0)*sdf.getSize(1)*sdf.getSize(2));
    GridHierarchy<Vec3> grid({sdf.getSize(0), sdf.getSize(1), sdf.getSize(2)}, sdf.getCellSize(), Vec3(sdf.getOrigin()));
    Parallel::For(0, tangents.size(), [&](size_t i) {
        Vec3 p = sdf.getPosition(i);
        Vec3 e = Vec3(sdf.getSize(0)-1, sdf.getSize(1)-1, sdf.getSize(2)-1)*sdf.getCellSize();
        switch (mode) {
            default:
            case 0:
                *grid.getData(i) = Vec3(1, 0, 0);
                break;
            case 1:
                *grid.getData(i) = Vec3(0, 1, 0);
                break;
            case 2:
                *grid.getData(i) = Vec3(0, 0, 1);
                break;
            case 3:
                *grid.getData(i) = p.x < 0.5f*e.x ? Vec3(1, 0, 0) : Vec3(0, 0, 1);
                break;
            case 4:
                *grid.getData(i) = Vec3(p.x-0.5f*e.x, 0, p.z-0.5f*e.z).normalize();
                break;
            case 5:
                *grid.getData(i) = Vec3(0.5f*e.x-p.x, 0, p.z-0.5f*e.z).normalize();
                break;
            case 6:
                *grid.getData(i) = sdf.getGrad(i);
                break;
            case 7:
                *grid.getData(i) = sdf.getGrad(i).normalToTangent();
                break;
            case 8:
                *grid.getData(i) = randomDirectionFromUint(i);
                break;
            case 9:
                *grid.getData(i) = Vec3(sinf(20.0f*p.x/e.x), 0, cosf(20.0f*p.x/e.x));
                break;
        }
    });

    size_t smoothIter = mode >= 6 && mode <= 7 ? 32 : 1;
    if (mode == 8) smoothIter = 256;
    grid.smooth([&](size_t i, GridHierarchy<Vec3> *grid) {
        Vec3 avg(0);
        for (size_t n : grid->getNeighborhood(i)) {
            Vec3 d = *grid->getData(n);
            avg += (d.dot(avg) < 0 ? -1 : 1)*d;
        }
        *grid->getData(i) = avg.normalize();
    }, smoothIter);

    Parallel::ForProgress(0, tangents.size(), [&](size_t i) {
        tangents[i] = *grid.getData(i);
    });

    // grid.saveVec3ToAttributePLY("../data/tangents.ply", [&](size_t i, const GridHierarchy<Vec3> *grid, bool &save) {
    //     save = grid->idx1To3(i)[0] % 10 == 0 && grid->idx1To3(i)[2] % 10 == 0;
    //     save &= sdf.getField()[i] < 0 && grid->idx1To3(i)[1] == size_t(0.35f*grid->getSize()[1]);
    //     return *grid->getData(i);
    // });

    return tangents;
}