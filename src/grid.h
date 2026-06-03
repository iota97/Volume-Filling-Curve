// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once

#include "mathutils.h"
#include "parallel.h"
#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>

template <class T>
class GridHierarchy {
public:
    GridHierarchy(const std::array<uint32_t, 3> &size, float cellSize, Vec3 origin = Vec3(0)) : size(size), cellSize(cellSize), origin(origin) {
        data.resize(size_t(size[0])*size[1]*size[2]);
        Downsample();
    }

    ~GridHierarchy() { delete coarse; }
    GridHierarchy(const GridHierarchy &) = delete;
    GridHierarchy &operator=(const GridHierarchy &) = delete;

    GridHierarchy(const GridHierarchy &&other) {
        data = other.data;
        size = other.size;
        cellSize = other.cellSize;
        coarse = other.coarse;
        dense = other.dense;
        other.coarse = nullptr;
        other.dense = nullptr;
    }

    GridHierarchy &operator=(GridHierarchy &&other) {
        if (this != &other) {
            data = other.data;
            size = other.size;
            cellSize = other.cellSize;
            coarse = other.coarse;
            dense = other.dense;
            other.coarse = nullptr;
            other.dense = nullptr;
        }
        return *this;
    }

    size_t getLodsCount() const {
        size_t count = 1;
        GridHierarchy* curr = coarse;
        while (curr) {
            curr = curr->coarse;
            count++;
        }
        return count;
    }

    const T &getData(size_t x, size_t y, size_t z) const {
        return data[idx3To1(x, y, z)];
    }

    size_t getTotalSize() const {
        return size[0]*size[1]*size[2];
    }

    T &getData(size_t x, size_t y, size_t z) {
        return data[idx3To1(x, y, z)];
    }

    Vec3 getPosition(size_t i) const {
        auto coord = idx1To3(i);
        return cellSize*Vec3(coord[0]+0.5, coord[1]+0.5, coord[2]+0.5) + origin;
    }

    struct NeighborhoodIter {
        const GridHierarchy<T> *grid;
        size_t i;
        uint8_t j;
        size_t p;
        NeighborhoodIter(const GridHierarchy<T> *grid, size_t i) : grid(grid), i(i), j(0) {
            ++(*this);
        }

        NeighborhoodIter(const GridHierarchy<T> *grid, size_t i, size_t j) : grid(grid), i(i), j(j) {}

        NeighborhoodIter &operator++() {
            do {
                auto coord = grid->idx1To3(i);
                switch (j++) {
                case 0:
                    p = coord[0] < grid->size[0]-1 ? grid->idx3To1(coord[0]+1, coord[1], coord[2]) : size_t(-1);
                    break;
                case 1:
                    p = coord[0] > 0 ? grid->idx3To1(coord[0]-1, coord[1], coord[2]) : size_t(-1);
                    break;
                case 2:
                    p = coord[1] < grid->size[1]-1 ? grid->idx3To1(coord[0], coord[1]+1, coord[2]) : size_t(-1);
                    break;
                case 3:
                    p = coord[1] > 0 ? grid->idx3To1(coord[0], coord[1]-1, coord[2]) : size_t(-1);
                    break;
                case 4:
                    p = coord[2] < grid->size[2]-1 ? grid->idx3To1(coord[0], coord[1], coord[2]+1) : size_t(-1);
                    break;
                case 5:
                    p = coord[2] > 0 ? grid->idx3To1(coord[0], coord[1], coord[2]-1) : size_t(-1);
                    break;
                default:
                    break;
                }
            } while (j < 7 && p == size_t(-1));
            return *this;
        }

        size_t operator*() const {
            return p;
        }

        bool operator!=(const NeighborhoodIter &o) const {
            return i != o.i || j != o.j;
        }
    };

    struct Neighborhood {
        Neighborhood(size_t i, const GridHierarchy<T> *grid) : grid(grid), i(i) {}
        const GridHierarchy<T> *grid;
        size_t i;
        NeighborhoodIter begin() const { return NeighborhoodIter(grid, i); }
        NeighborhoodIter end() const { return NeighborhoodIter(grid, i, 7); }
    };

    Neighborhood getNeighborhood(size_t i) const {
        return Neighborhood(i, this);
    };

    typedef std::function<void(size_t, GridHierarchy<T> *)> InitFunction;
    typedef std::function<void(size_t, GridHierarchy<T> *)> SmoothFunction;
    typedef std::function<bool(size_t, GridHierarchy<T> *)> SmoothConvergenceFunction;
    typedef std::function<void(size_t, GridHierarchy<T> *, GridHierarchy<T> *)> RestrictFunction;
    typedef std::function<void(size_t, GridHierarchy<T> *, GridHierarchy<T> *)> ProlongFunction;

    void init(const InitFunction &initFunc) {
        Parallel::For(0, data.size(), [this, &initFunc](size_t i) {
            initFunc(i, this);
        });
    }

    void optimize(const SmoothFunction &smoothFunc, const RestrictFunction &restrictFunc, const ProlongFunction &prolongFunc, uint32_t iterationCount = 32, bool quiet = false) {
        if (coarse) {
            restrict(restrictFunc);
            coarse->optimize(smoothFunc, restrictFunc, prolongFunc, iterationCount, true);
            prolong(prolongFunc);
            smooth(smoothFunc, iterationCount, quiet);
        }
    }

    T *getData(size_t i) { return &data[i]; }
    std::vector<T> *getData() { return &data; }
    const T *getData(size_t i) const { return &data[i]; }
    const std::vector<T> *getData() const { return &data; }

    typedef std::function<float(size_t, const GridHierarchy<T> *)> FloatFunction;
    typedef std::function<std::array<uint8_t, 3>(size_t, const GridHierarchy<T>*, bool& save)> ColorFunction;
    typedef std::function<Vec3(size_t, const GridHierarchy<T> *, bool &save)> Vec3Function;

    void saveFloatToBin(const char *path, const FloatFunction& floatFunc) const {
        std::vector<float> floatDat(data.size());
        Parallel::For(0, data.size(), [this, &floatDat, &floatFunc](size_t i) {
            floatDat[i] = floatFunc(i, this);
        });
        std::ofstream sdfFile(path, std::ios::out | std::ios::binary);
        sdfFile.write((char*)&size[0], sizeof(uint32_t));
        sdfFile.write((char*)&size[1], sizeof(uint32_t));
        sdfFile.write((char*)&size[2], sizeof(uint32_t));
        sdfFile.write((char*)&cellSize, sizeof(float));
        sdfFile.write((char*)floatDat.data(), floatDat.size()*sizeof(float));
    }

    void saveVec3ToPLY(const char *path, const Vec3Function& vec3Func) const {
        std::vector<Vec3> pos;
        std::vector<Vec3> norm;

        for (size_t i = 0; i < data.size(); ++i) {
            bool add = true;
            Vec3 normal = vec3Func(i, this, add);
            if (add) {
                pos.push_back(getPosition(i));
                norm.push_back(0.5f*cellSize * normal);
            }
        }
        
        std::ofstream plyFile(path);
        plyFile << "ply\nformat ascii 1.0\nelement vertex " << 2*pos.size() << "\n"
                << "property float x\nproperty float y\nproperty float z\n"
                << "element edge " << pos.size() << "\nproperty int vertex1\nproperty int vertex2\n"
                << "end_header\n";
        
        for (size_t i = 0; i < pos.size() ; ++i) {
            plyFile << pos[i].x << " " << pos[i].y << " " << pos[i].z << std::endl;
            plyFile << pos[i].x + norm[i].x << " " << pos[i].y + norm[i].y << " " << pos[i].z + norm[i].z << std::endl;
        }

        for (size_t i = 0; i < pos.size(); ++i) {
            plyFile << 2*i << " " << 2*i+1 << std::endl;
        }
    }

    void saveVec3ToAttributePLY(const char *path, const Vec3Function& vec3Func) const {
        std::vector<Vec3> pos;
        std::vector<Vec3> norm;

        for (size_t i = 0; i < data.size(); ++i) {
            bool add = true;
            Vec3 normal = vec3Func(i, this, add);
            if (add) {
                pos.push_back(getPosition(i));
                norm.push_back(normal);
            }
        }
        
        std::ofstream plyFile(path);
        plyFile << "ply\nformat ascii 1.0\nelement vertex " << pos.size() << "\n"
                << "property float x\nproperty float y\nproperty float z\n"
                << "property float nx\nproperty float ny\nproperty float nz\n"
                << "end_header\n";
        
        for (size_t i = 0; i < pos.size() ; ++i) {
            plyFile << pos[i].x << " " << pos[i].y << " " << pos[i].z << " " << norm[i].x << " " << norm[i].y << " " << norm[i].z << std::endl;
        }
    }

    size_t getIndex(const Vec3 &position) const {
        Vec3 pos = position - origin;
        return idx3To1({uint32_t(floorf(pos.x/cellSize)), uint32_t(floorf(pos.y/cellSize)), uint32_t(floorf(pos.z/cellSize))});
    }

    void saveColorToPLY(const char* path, const ColorFunction& colorFunc, const ProlongFunction& prolongFunc, size_t samples) const {
        GridHierarchy<T> *superGrid = new GridHierarchy<T>();
        size_t scale = 1 << samples;
        superGrid->size = { uint32_t(size[0] * scale), uint32_t(size[1] * scale), uint32_t(size[2] * scale) };
        superGrid->data.resize(superGrid->size[0] * superGrid->size[1] * superGrid->size[2]);
        superGrid->cellSize = cellSize / scale;
        superGrid->dense = nullptr;
        superGrid->coarse = nullptr;
        superGrid->origin = origin;
        GridHierarchy<T>* curr = superGrid;
        for (size_t i = 0; i < samples; ++i) {
            curr->coarse = new GridHierarchy<T>();
            curr->coarse->size = { uint32_t(ceil(curr->size[0] * 0.5)), uint32_t(ceil(curr->size[1] * 0.5)), uint32_t(ceil(curr->size[2] * 0.5)) };
            curr->coarse->data.resize(curr->coarse->size[0] * curr->coarse->size[1] * curr->coarse->size[2]);
            curr->coarse->cellSize = curr->cellSize * 2;
            curr->coarse->dense = const_cast<GridHierarchy<T> *>(curr);
            curr->coarse->origin = origin;
            curr = curr->coarse;
        }

        Parallel::For(0, curr->data.size(), [curr, this](size_t i) {
            curr->data[i] = data[i];
        });

        for (size_t i = 0; i < samples; ++i) {
            Parallel::For(0, curr->dense->data.size(), [curr, &prolongFunc](size_t i) {
                prolongFunc(i, curr->dense, curr);
            });
            curr = curr->dense;
        }

        auto &grid = *curr;
        std::vector<Vec3> pos;
        std::vector<std::array<uint32_t, 4>> faces;
        std::vector<std::array<uint8_t, 3>> colors;
        for (size_t x = 0; x < grid.getSize()[0] - 1; ++x) {
            for (size_t y = 0; y < grid.getSize()[1] - 1; ++y) {
                for (size_t z = 0; z < grid.getSize()[2] - 1; ++z) {
                    auto _000 = grid.idx3To1({ uint32_t(x), uint32_t(y) , uint32_t(z) });
                    auto _001 = grid.idx3To1({ uint32_t(x), uint32_t(y) , uint32_t(z) + 1 });
                    auto _010 = grid.idx3To1({ uint32_t(x), uint32_t(y) + 1 , uint32_t(z) });
                    auto _011 = grid.idx3To1({ uint32_t(x), uint32_t(y) + 1 , uint32_t(z) + 1 });
                    auto _100 = grid.idx3To1({ uint32_t(x) + 1, uint32_t(y) , uint32_t(z) });
                    auto _101 = grid.idx3To1({ uint32_t(x) + 1, uint32_t(y) , uint32_t(z) + 1 });
                    auto _110 = grid.idx3To1({ uint32_t(x) + 1, uint32_t(y) + 1 , uint32_t(z) });
                    auto _111 = grid.idx3To1({ uint32_t(x) + 1, uint32_t(y) + 1 , uint32_t(z) + 1 });
                    bool save;
                    colorFunc(_000, &grid, save); if (!save) continue;
                    colorFunc(_001, &grid, save); if (!save) continue;
                    colorFunc(_010, &grid, save); if (!save) continue;
                    colorFunc(_011, &grid, save); if (!save) continue;
                    colorFunc(_100, &grid, save); if (!save) continue;
                    colorFunc(_101, &grid, save); if (!save) continue;
                    colorFunc(_110, &grid, save); if (!save) continue;
                    colorFunc(_111, &grid, save); if (!save) continue;

                    for (uint8_t f = 0; f < 6; ++f) {
                        size_t a, b, c, d;
                        switch (f) {
                        case 0:
                            a = _000, b = _001, c = _011, d = _010;
                            break;
                        case 1:
                            a = _100, b = _101, c = _111, d = _110;
                            break;
                        case 2:
                            a = _001, b = _101, c = _100, d = _000;
                            break;
                        case 3:
                            a = _011, b = _111, c = _110, d = _010;
                            break;
                        case 4:
                            a = _000, b = _010, c = _110, d = _100;
                            break;
                        default:
                            a = _001, b = _011, c = _111, d = _101;
                            break;
                        }

                        size_t curr = pos.size();
                        pos.push_back(grid.getPosition(a));
                        pos.push_back(grid.getPosition(b));
                        pos.push_back(grid.getPosition(c));
                        pos.push_back(grid.getPosition(d));
                        colors.push_back(colorFunc(a, &grid, save));
                        colors.push_back(colorFunc(b, &grid, save));
                        colors.push_back(colorFunc(c, &grid, save));
                        colors.push_back(colorFunc(d, &grid, save));

                        faces.push_back({ uint32_t(curr + 0) , uint32_t(curr + 1) , uint32_t(curr + 2) , uint32_t(curr + 3) });
                    }
                }
            }
        }

        delete superGrid;

        uint32_t usedVerts = 0;
        std::unordered_map<Vec3, uint32_t> vertsID;
        std::unordered_set<uint32_t> usedVert;
        for (size_t i = 0; i < pos.size(); ++i) {
            if (!vertsID.count(pos[i])) {
                vertsID[pos[i]] = usedVerts++;
                usedVert.insert(i);
            }
        }

        std::set<std::array<uint32_t, 4>> usedFaces;
        for (size_t i = 0; i < faces.size(); ++i) {
            uint32_t a = vertsID[pos[faces[i][0]]];
            uint32_t b = vertsID[pos[faces[i][1]]];
            uint32_t c = vertsID[pos[faces[i][2]]];
            uint32_t d = vertsID[pos[faces[i][3]]];
            std::array<uint32_t, 4> face{ a, b, c, d };
            usedFaces.insert(face);
        }

        std::ofstream plyFile(path);
        plyFile << "ply\nformat binary_little_endian 1.0\nelement vertex " << usedVerts << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            << "element face " << usedFaces.size() << "\nproperty list uchar uint vertex_indices\n"
            << "end_header\n";
        plyFile.close();
        plyFile = std::ofstream(path, std::ios::out | std::ios::binary | std::ios::app);

        for (size_t i = 0; i < pos.size(); ++i) {
            if (usedVert.count(i)) {
                plyFile.write(reinterpret_cast<const char*>(&pos[i].x), 4);
                plyFile.write(reinterpret_cast<const char*>(&pos[i].y), 4);
                plyFile.write(reinterpret_cast<const char*>(&pos[i].z), 4);
                plyFile.write(reinterpret_cast<const char*>(&colors[i][0]), 1);
                plyFile.write(reinterpret_cast<const char*>(&colors[i][1]), 1);
                plyFile.write(reinterpret_cast<const char*>(&colors[i][2]), 1);
            }
        }

        for (const auto &face : usedFaces) {
            uint8_t four = 4;
            plyFile.write(reinterpret_cast<const char*>(&four), 1);
            plyFile.write(reinterpret_cast<const char*>(&face[0]), 4);
            plyFile.write(reinterpret_cast<const char*>(&face[1]), 4);
            plyFile.write(reinterpret_cast<const char*>(&face[2]), 4);
            plyFile.write(reinterpret_cast<const char*>(&face[3]), 4);
        }
    }

    std::array<uint32_t, 3> idx1To3(size_t i) const {
        return {uint32_t(i/size[1]/size[2]), uint32_t((i/size[2])%size[1]), uint32_t(i % size[2])};
    }

    size_t idx3To1(size_t x, size_t y, size_t z) const {
        return x*size[1]*size[2] + y*size[2] + z;
    }

    size_t idx3To1(const std::array<uint32_t, 3> &c) const {
        return size_t(c[0])*size[1]*size[2] + size_t(c[1])*size[2] + size_t(c[2]);
    }

    std::vector<size_t> coarseToDense(size_t i) const {
        std::vector<size_t> neigh;
        auto coarseCoord = idx1To3(i);
        neigh.push_back(dense->idx3To1(2*coarseCoord[0], 2*coarseCoord[1], 2*coarseCoord[2]));
        if (2*coarseCoord[0]+1 < dense->size[0]) {
            neigh.push_back(dense->idx3To1(2*coarseCoord[0]+1, 2*coarseCoord[1], 2*coarseCoord[2]));
            if (2*coarseCoord[1]+1 < dense->size[1]) {
                neigh.push_back(dense->idx3To1(2*coarseCoord[0]+1, 2*coarseCoord[1]+1, 2*coarseCoord[2]));
                if (2*coarseCoord[2]+1 < dense->size[2]) {
                    neigh.push_back(dense->idx3To1(2*coarseCoord[0]+1, 2*coarseCoord[1]+1, 2*coarseCoord[2]+1));
                }
            }
            if (2*coarseCoord[2]+1 < dense->size[2]) {
                neigh.push_back(dense->idx3To1(2*coarseCoord[0]+1, 2*coarseCoord[1], 2*coarseCoord[2]+1));
            }
        }
        if (2*coarseCoord[1]+1 < dense->size[1]) {
            neigh.push_back(dense->idx3To1(2*coarseCoord[0], 2*coarseCoord[1]+1, 2*coarseCoord[2]));
            if (2*coarseCoord[2]+1 < dense->size[2]) {
                neigh.push_back(dense->idx3To1(2*coarseCoord[0], 2*coarseCoord[1]+1, 2*coarseCoord[2]+1));
            }
        }
        if (2*coarseCoord[2]+1 < dense->size[2]) {
            neigh.push_back(dense->idx3To1(2*coarseCoord[0], 2*coarseCoord[1], 2*coarseCoord[2]+1));
        }
        return neigh;
    }

    size_t denseToCoarse(size_t i) const {
        auto coarseDense = idx1To3(i);
        return coarse->idx3To1(coarseDense[0]/2, coarseDense[1]/2, coarseDense[2]/2);
    }

    const std::array<uint32_t, 3> & getSize() const { return size; }
    float getCellSize() const { return cellSize; }

    void smooth(const SmoothFunction &smoothFunc, uint32_t iterationCount = 32, bool quiet = true) {
        for (uint32_t iter = 0; iter < 2*iterationCount; ++iter) {
            if (!quiet) {
                printProgress(float(iter)/(2*iterationCount-1));
            }
            Parallel::For(0, size[0]*size[1]*size[2], [this, &smoothFunc, iter](size_t i) {
                auto coord = idx1To3(i);
                if ((coord[0] + coord[1] + coord[2]) % 2 == iter % 2) {
                    smoothFunc(i, this);
                }
            });
        }
    }

    void smoothConvergence(const SmoothConvergenceFunction &smoothFunc) {
        for (uint32_t iter = 0;; ++iter) {
            bool changed = Parallel::ForAny(0, size[0]*size[1]*size[2], [this, &smoothFunc, iter](size_t i) {
                auto coord = idx1To3(i);
                if ((coord[0] + coord[1] + coord[2]) % 2 == iter % 2) {
                    return smoothFunc(i, this);
                }
                return false;
            });
            if (!changed) break;
        }
    }

private:
    GridHierarchy() {}

    void Downsample() {
        if (size[0]+size[1]+size[2] > 3) {
            coarse = new GridHierarchy<T>();
            coarse->size = {uint32_t(ceil(size[0]*0.5)), uint32_t(ceil(size[1]*0.5)), uint32_t(ceil(size[2]*0.5))};
            coarse->data.resize(coarse->size[0]*coarse->size[1]*coarse->size[2]);
            coarse->cellSize = cellSize*2;
            coarse->dense = this;
            coarse->origin = origin;
            coarse->Downsample();
        }
    }

    void restrict(const RestrictFunction &restrictFunc) {
        Parallel::For(0, coarse->data.size(), [this, &restrictFunc](size_t i) {
            restrictFunc(i, this, coarse);
        });
    }

    void prolong(const ProlongFunction &prolongFunc) {
        Parallel::For(0, data.size(), [this, &prolongFunc](size_t i) {
            prolongFunc(i, this, coarse);
        });
    }

    std::vector<T> data;
    std::array<uint32_t, 3> size;
    float cellSize;
    Vec3 origin;
    GridHierarchy<T> *coarse = nullptr;
    GridHierarchy<T> *dense = nullptr;
};