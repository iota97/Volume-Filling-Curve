// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "sdf.h"
#include "bvhtriangle.h"
#include "parallel.h"

SDF::SDF(const std::vector<std::array<Vec3, 3>> &triangles, float cellSize, uint32_t border, uint32_t raysCount, bool quiet) : cellSize(cellSize), origin(-float(border*cellSize)) {
    BVHTriangle bvh(triangles);
    Vec3 minAABB, maxAABB;
    bvh.getAABB(minAABB, maxAABB);
    Vec3 deltaAABB = maxAABB-minAABB;
    Vec3 offset = minAABB-border*cellSize;
    size[0] = uint32_t(ceilf(deltaAABB.x/cellSize)+2*border);
    size[1] = uint32_t(ceilf(deltaAABB.y/cellSize)+2*border);
    size[2] = uint32_t(ceilf(deltaAABB.z/cellSize)+2*border);
    field.resize(size_t(size[0])*size[1]*size[2]);

    Parallel::ForGrainProgress(0, field.size(), [&](size_t i) {
        size_t x = i/size[1]/size[2], y = (i/size[2])%size[1], z = i % size[2];
        Vec3 point = cellSize*Vec3(x, y, z)+0.5*cellSize + offset;
        Vec3 closest = bvh.closestPoint(point);
        Vec3 delta = point-closest;
        uint32_t outside = 0;
        uint32_t seed = i;
        for (uint32_t j = 0; j < raysCount; ++j) {
            outside += bvh.intersectionCount(point, randomDirectionFromUint(seed)) % 2;
            seed = randomUintFromUint(seed);
        }
        float sign = outside > raysCount/2 ? -1 : 1;
        field[i] = sign * delta.length();
    }, 4096, quiet);
}

SDF::SDF(const char* path) {
    std::ifstream sdfFile(path, std::ios::binary);
    sdfFile.read((char*)&size[0], sizeof(uint32_t));
    sdfFile.read((char*)&size[1], sizeof(uint32_t));
    sdfFile.read((char*)&size[2], sizeof(uint32_t));
    sdfFile.read((char*)&cellSize, sizeof(float));
    field.resize(size[0]*size[1]*size[2]);
    sdfFile.read((char*)field.data(), field.size()*sizeof(float));
    origin = Vec3(-2.0f*cellSize);
}

void SDF::saveToBin(const char *path) const {
    std::ofstream sdfFile(path, std::ios::out | std::ios::binary);
    sdfFile.write((char*)&size[0], sizeof(uint32_t));
    sdfFile.write((char*)&size[1], sizeof(uint32_t));
    sdfFile.write((char*)&size[2], sizeof(uint32_t));
    sdfFile.write((char*)&cellSize, sizeof(float));
    sdfFile.write((char*)field.data(), field.size()*sizeof(float));
}

Vec3 SDF::getPosition(size_t i) const {
    std::array<uint32_t, 3> coord{uint32_t(i/size[1]/size[2]), uint32_t((i/size[2])%size[1]), uint32_t(i % size[2])};
    return cellSize*Vec3(coord[0]+0.5, coord[1]+0.5, coord[2]+0.5) + origin;
}

size_t SDF::getIndex(const Vec3 &position) const {
    Vec3 pos = position - origin;
    return idx3To1({uint32_t(floorf(pos.x/cellSize)), uint32_t(floorf(pos.y/cellSize)), uint32_t(floorf(pos.z/cellSize))});
}

float SDF::weight(Vec3 pos, size_t cellIndex) const {
    return (std::max(0.0f, cellSize - (pos-getPosition(cellIndex)).length()))/cellSize;
}

float SDF::sdfVal(Vec3 pos, size_t cellIndex) const {
    float sdfVal = 0, wSum = 0;
    auto xyz = idx1To3(cellIndex);
    uint32_t minX = xyz[0] > 0 ? xyz[0]-1 : xyz[0];
    uint32_t maxX = xyz[0] < size[0]-1 ? xyz[0]+1 : xyz[0];
    uint32_t minY = xyz[1] > 0 ? xyz[1]-1 : xyz[1];
    uint32_t maxY = xyz[1] < size[1]-1 ? xyz[1]+1 : xyz[1];
    uint32_t minZ = xyz[2] > 0 ? xyz[2]-1 : xyz[2];
    uint32_t maxZ = xyz[2] < size[2]-1 ? xyz[2]+1 : xyz[2];

    for (uint32_t x = minX; x <= maxX; ++x) {
        for (uint32_t y = minY; y <= maxY; ++y) {
            for (uint32_t z = minZ; z <= maxZ; ++z) {
                float w = weight(pos, idx3To1({x, y, z}));
                sdfVal += w*field[idx3To1({x, y, z})];
                wSum += w;
            }
        }
    }
    return sdfVal/wSum;
}

float SDF::getVal(const Vec3 &position) const {
    return sdfVal(position, getIndex(position));
}

std::array<uint32_t, 3> SDF::idx1To3(size_t i) const {
    return {uint32_t(i/size[1]/size[2]), uint32_t((i/size[2])%size[1]), uint32_t(i % size[2])};
}

size_t SDF::idx3To1(const std::array<uint32_t, 3> &c) const {
    return size_t(c[0])*size[1]*size[2] + size_t(c[1])*size[2] + size_t(c[2]);
}

size_t SDF::idx3To1(size_t x, size_t y, size_t z) const {
    return x*size[1]*size[2] + y*size[2] + z;
}

Vec3 SDF::getGrad(const Vec3 &position) const {
    auto xyz = idx1To3(getIndex(position));
    return getGrad(xyz);
}

Vec3 SDF::getGrad(size_t i) const {
    return getGrad(idx1To3(i));
}

Vec3 SDF::getGrad(const std::array<uint32_t, 3>& xyz) const {
    std::array<uint32_t, 3> xm1yz = xyz[0] > 0 ? std::array<uint32_t, 3>{xyz[0] - 1, xyz[1], xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xp1yz = xyz[0] < getSize()[0] - 1 ? std::array<uint32_t, 3>{xyz[0] + 1, xyz[1], xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xym1z = xyz[1] > 0 ? std::array<uint32_t, 3>{xyz[0], xyz[1] - 1, xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xyp1z = xyz[1] < getSize()[1] - 1 ? std::array<uint32_t, 3>{xyz[0], xyz[1] + 1, xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xyzm1 = xyz[2] > 0 ? std::array<uint32_t, 3>{xyz[0], xyz[1], xyz[2] - 1} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xyzp1 = xyz[2] < getSize()[2] - 1 ? std::array<uint32_t, 3>{xyz[0], xyz[1], xyz[2] + 1} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    float dx = xm1yz[0] == uint32_t(-1) || xp1yz[0] == uint32_t(-1) ? getCellSize() : 2.0f * getCellSize();
    float dy = xym1z[0] == uint32_t(-1) || xyp1z[0] == uint32_t(-1) ? getCellSize() : 2.0f * getCellSize();
    float dz = xyzm1[0] == uint32_t(-1) || xyzp1[0] == uint32_t(-1) ? getCellSize() : 2.0f * getCellSize();
    if (xm1yz[0] == uint32_t(-1)) xm1yz = xyz;
    if (xp1yz[0] == uint32_t(-1)) xp1yz = xyz;
    if (xym1z[0] == uint32_t(-1)) xym1z = xyz;
    if (xyp1z[0] == uint32_t(-1)) xyp1z = xyz;
    if (xyzm1[0] == uint32_t(-1)) xyzm1 = xyz;
    if (xyzp1[0] == uint32_t(-1)) xyzp1 = xyz;

    float dfdx = (field[idx3To1(xp1yz)] - field[idx3To1(xm1yz)]) / dx;
    float dfdy = (field[idx3To1(xyp1z)] - field[idx3To1(xym1z)]) / dy;
    float dfdz = (field[idx3To1(xyzp1)] - field[idx3To1(xyzm1)]) / dz;

    return Vec3(dfdx, dfdy, dfdz).normalize();
}

float SDF::getCurvature(const Vec3& position) const {
    auto xyz = idx1To3(getIndex(position));
    std::array<uint32_t, 3> xm1yz = xyz[0] > 0 ? std::array<uint32_t, 3>{xyz[0] - 1, xyz[1], xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xp1yz = xyz[0] < getSize()[0] - 1 ? std::array<uint32_t, 3>{xyz[0] + 1, xyz[1], xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xym1z = xyz[1] > 0 ? std::array<uint32_t, 3>{xyz[0], xyz[1] - 1, xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xyp1z = xyz[1] < getSize()[1] - 1 ? std::array<uint32_t, 3>{xyz[0], xyz[1] + 1, xyz[2]} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xyzm1 = xyz[2] > 0 ? std::array<uint32_t, 3>{xyz[0], xyz[1], xyz[2] - 1} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    std::array<uint32_t, 3> xyzp1 = xyz[2] < getSize()[2] - 1 ? std::array<uint32_t, 3>{xyz[0], xyz[1], xyz[2] + 1} : std::array<uint32_t, 3>{ uint32_t(-1), uint32_t(-1), uint32_t(-1) };
    if (xm1yz[0] == uint32_t(-1)) xm1yz = xyz;
    if (xp1yz[0] == uint32_t(-1)) xp1yz = xyz;
    if (xym1z[0] == uint32_t(-1)) xym1z = xyz;
    if (xyp1z[0] == uint32_t(-1)) xyp1z = xyz;
    if (xyzm1[0] == uint32_t(-1)) xyzm1 = xyz;
    if (xyzp1[0] == uint32_t(-1)) xyzp1 = xyz;

    float x = (1.0f - getGrad(xm1yz).dot(getGrad(xp1yz))) * 0.5f;
    float y = (1.0f - getGrad(xym1z).dot(getGrad(xyp1z))) * 0.5f;
    float z = (1.0f - getGrad(xyzm1).dot(getGrad(xyzp1))) * 0.5f;

    return (x + y + z) / 3.0f;
}

void SDF::saveSliceToPLY(const char *path) const {
    std::vector<Vec3> pos;
    std::vector<std::array<uint32_t, 4>> faces;
    std::vector<std::array<uint32_t, 2>> edges;
    std::vector<float> values;

    size_t y = size[1]/2;
    for (size_t x = 0; x < size[0]-1; ++x) {
        for (size_t z = 0; z < size[2]-1; ++z) {
            uint32_t curr = pos.size();
            pos.push_back(getPosition(idx3To1(x, y, z)));
            pos.push_back(getPosition(idx3To1(x, y, z+1)));
            pos.push_back(getPosition(idx3To1(x+1, y, z+1)));
            pos.push_back(getPosition(idx3To1(x+1, y, z)));
            faces.push_back({curr+0, curr+1, curr+2, curr+3});

            for (uint32_t i = 0; i < 4; ++i) {
                values.push_back(field[getIndex(pos[curr+i])]);
            }
        }
    }

    float maxVal = -std::numeric_limits<float>::infinity(), minVal = std::numeric_limits<float>::infinity();
    for (float val : values) {
        maxVal = std::max(maxVal, val);
        minVal = std::min(minVal, val);
    }

    for (auto &face : faces) {
        for (float bias = -maxVal; bias <= -minVal; bias += 4.0f*cellSize) {
            uint8_t a = 0, b = 1, c = 2;
            for (uint8_t i = 0; i < 2; ++i) {
                uint32_t i0 = face[a], i1 = face[b], i2 = face[c];
                float v0 = values[face[a]] + bias, v1 = values[face[b]] + bias, v2 = values[face[c]] + bias;
                uint32_t a0, a1, b0, b1;
                uint8_t mask = ((v0 >= 0) << 2) | ((v1 >= 0) << 1) | ((v2 >= 0) << 0);
                switch (mask) {
                    case 0b001:
                    case 0b110:
                        a0 = i0; a1 = i2; b0 = i1, b1 = i2;
                        break;
                    case 0b010:
                    case 0b101:
                        a0 = i0; a1 = i1; b0 = i1, b1 = i2;
                        break;
                    case 0b011:
                    case 0b100:
                        a0 = i0; a1 = i1; b0 = i0, b1 = i2;
                        break;
                    default:
                        a = 0, b = 2, c = 3;
                        continue;
                }

                float va0 = values[a0] + bias, va1 = values[a1] + bias, vb0 = values[b0] + bias, vb1 = values[b1] + bias;
                Vec3 pa = pos[a0]-va0/(va1-va0)*(pos[a1]-pos[a0]);
                Vec3 pb = pos[b0]-vb0/(vb1-vb0)*(pos[b1]-pos[b0]);
                uint32_t curr = pos.size();
                pos.push_back(pa);
                pos.push_back(pb);
                values.push_back(0);
                values.push_back(0);
                edges.push_back({curr, curr+1});
                a = 0, b = 2, c = 3;
            }
        }
    }

    for (float &val : values) {
        val = (val-minVal)/(maxVal-minVal);
    }

    uint32_t usedVerts = 0;
    std::unordered_map<Vec3, uint32_t> vertsID;
    std::unordered_set<uint32_t> usedVert;
    for (size_t i = 0; i < pos.size(); ++i) {
        if (!vertsID.count(pos[i])) {
            vertsID[pos[i]] = usedVerts++;
            usedVert.insert(i);
        }
    }

    for (size_t i = 0; i < faces.size(); ++i) {
        faces[i][0] = vertsID[pos[faces[i][0]]];
        faces[i][1] = vertsID[pos[faces[i][1]]];
        faces[i][2] = vertsID[pos[faces[i][2]]];
        faces[i][3] = vertsID[pos[faces[i][3]]];
    }

    for (size_t i = 0; i < edges.size(); ++i) {
        edges[i][0] = vertsID[pos[edges[i][0]]];
        edges[i][1] = vertsID[pos[edges[i][1]]];
    }

    std::ofstream plyFile(path);
    plyFile << "ply\nformat binary_little_endian 1.0\nelement vertex " << usedVerts << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            << "element face " << faces.size() << "\nproperty list uchar uint vertex_indices\n"
            << "element edge " << edges.size() << "\nproperty int vertex1\nproperty int vertex2\n"
            << "end_header\n";
    plyFile.close();
    plyFile = std::ofstream(path, std::ios::out | std::ios::binary | std::ios::app);

    for (size_t i = 0; i < pos.size(); ++i) {
        if (!usedVert.count(i)) continue;
        plyFile.write(reinterpret_cast<const char*>(&pos[i].x), 4);
        plyFile.write(reinterpret_cast<const char*>(&pos[i].y), 4);
        plyFile.write(reinterpret_cast<const char*>(&pos[i].z), 4);
        auto color = turbo(values[i]);
        plyFile.write(reinterpret_cast<const char*>(&color[0]), 1);
        plyFile.write(reinterpret_cast<const char*>(&color[1]), 1);
        plyFile.write(reinterpret_cast<const char*>(&color[2]), 1);
    }

    for (const auto &face : faces) {
        uint8_t four = 4;
        plyFile.write(reinterpret_cast<const char*>(&four), 1);
        plyFile.write(reinterpret_cast<const char*>(&face[0]), 4);
        plyFile.write(reinterpret_cast<const char*>(&face[1]), 4);
        plyFile.write(reinterpret_cast<const char*>(&face[2]), 4);
        plyFile.write(reinterpret_cast<const char*>(&face[3]), 4);
    }

    for (const auto &edge : edges) {
        plyFile.write(reinterpret_cast<const char*>(&edge[0]), 4);
        plyFile.write(reinterpret_cast<const char*>(&edge[1]), 4);
    }
}