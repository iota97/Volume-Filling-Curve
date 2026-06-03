// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#include "marchingtriangle.h"
#include <algorithm>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <iostream>
#include <fstream>

static std::vector<std::vector<std::array<uint8_t, 2>>> cases {
    {},
    {{0, 1}, {0, 2}},
    {{0, 1}, {1, 2}},
    {{1, 2}, {2, 0}},
    {{1, 2}, {2, 0}},
    {{0, 1}, {1, 2}},
    {{0, 1}, {0, 2}},
    {},
};

std::vector<std::array<Vec3Idx, 2>> marchingTriangle(const std::function<float(Vec3)> &value, std::vector<std::array<Vec3Idx, 3>> triangles, float quiet) {
    std::vector<std::array<Vec3Idx, 2>> segments;
    std::map<std::array<size_t, 2>, size_t> idxMap;
    for (size_t i = 0; i < triangles.size(); ++i) {
        if (!quiet) printProgress(double(i+1)/triangles.size());
        const auto &tris = triangles[i];
        size_t idx = 0;
        for (uint8_t i = 0; i < 3; ++i) {
            float val = value(tris[i].v);
            idx |= (val >= 0) << i;
        }

        std::vector<std::array<uint8_t, 2>> &edges = cases[idx];
        if (edges.size() == 2) {
            Vec3 pa0 = tris[edges[0][0]].v, pa1 = tris[edges[0][1]].v;
            if (pa1 < pa0) std::swap(pa0, pa1);
            Vec3 pb0 = tris[edges[1][0]].v, pb1 = tris[edges[1][1]].v;
            if (pb1 < pb0) std::swap(pb0, pb1);

            float va0 = value(pa0), va1 = value(pa1);
            float vb0 = value(pb0), vb1 = value(pb1);
            if (fabsf(va0) < 1e-3) va0 = va0 >= 0 ? 1e-3 : -1e-3; 
            if (fabsf(va1) < 1e-3) va1 = va1 >= 0 ? 1e-3 : -1e-3; 
            if (fabsf(vb0) < 1e-3) vb0 = vb0 >= 0 ? 1e-3 : -1e-3; 
            if (fabsf(vb1) < 1e-3) vb1 = vb1 >= 0 ? 1e-3 : -1e-3;

            std::array<size_t, 2> ca{tris[edges[0][0]].idx, tris[edges[0][1]].idx};
            std::array<size_t, 2> cb{tris[edges[1][0]].idx, tris[edges[1][1]].idx};
            if (ca[1] < ca[0]) std::swap(ca[0], ca[1]);
            if (cb[1] < cb[0]) std::swap(cb[0], cb[1]);
            size_t idxA, idxB;
            if (idxMap.count(ca)) {
                idxA = idxMap[ca];
            } else {
                idxA = idxMap.size();
                idxMap[ca] = idxA; 
            }

            if (idxMap.count(cb)) {
                idxB = idxMap[cb];
            } else {
                idxB = idxMap.size();
                idxMap[cb] = idxB; 
            }

            Vec3Idx va{ pa0-va0/(va1-va0)*(pa1-pa0), idxA};
            Vec3Idx vb{ pb0-vb0/(vb1-vb0)*(pb1-pb0), idxB};

            segments.push_back({va, vb});
        }
    }
    return segments;
}

std::vector<std::vector<Vec3>> connectCycles(const std::vector<std::array<Vec3Idx, 2>> &segments, float spacing) {
    size_t vertsCount = 0;
    std::vector<std::vector<Vec3>> cycles;

    for (const auto &s : segments) {
        vertsCount = std::max(vertsCount, s[0].idx + 1);
        vertsCount = std::max(vertsCount, s[1].idx + 1);
    }

    std::vector<Vec3> verts(vertsCount);
    for (const auto &s : segments) {
        verts[s[0].idx] = s[0].v;
        verts[s[1].idx] = s[1].v;
    }

    std::unordered_map<size_t, std::array<size_t, 2>> connectivity;
    std::unordered_set<size_t> usedVerts;
    for (const std::array<Vec3Idx, 2> &segment : segments) {
        if (segment[0].idx == segment[1].idx) continue;
        usedVerts.insert(segment[0].idx);
        usedVerts.insert(segment[1].idx);
        if (connectivity.count(segment[0].idx)) {
            connectivity[segment[0].idx][1] = segment[1].idx;
        } else {
            connectivity[segment[0].idx] = {segment[1].idx, size_t(-1)};
        }

        if (connectivity.count(segment[1].idx)) {
            connectivity[segment[1].idx][1] = segment[0].idx;
        } else {
            connectivity[segment[1].idx] = {segment[0].idx, size_t(-1)};
        }
    }

    for (const auto &adj : connectivity) {
        if (adj.second[1] ==  size_t(-1)) {
            std::cerr << "Fatal Error: open cycle detected!" << std::endl;
            return cycles;
        }
    }

    size_t totalVerts = usedVerts.size();
    while (!usedVerts.empty()) {
        std::vector<size_t> cycleIdx;
        size_t current = *usedVerts.begin();
        usedVerts.erase(current);
        cycleIdx.push_back(current);
        while (true) {
            const std::array<size_t, 2> &neighbors = connectivity[current];
            size_t next = cycleIdx.size() < 2 || cycleIdx[cycleIdx.size()-2] == neighbors[1] ? neighbors[0] : neighbors[1];
            if (next == cycleIdx[0]) {
                break;
            }
            cycleIdx.push_back(next);
            usedVerts.erase(next);
            current = next;
            if (cycleIdx.size() > totalVerts) {
                std::cerr << "Warning: broken cycle detected!" << std::endl;
                cycleIdx.clear();
                break;
            }
        }

        std::vector<Vec3> cycle;
        for (size_t idx : cycleIdx) {
            cycle.push_back(verts[idx]);
        }
        cycle = resampleCycle(cycle, spacing);
        if (cycle.size() > 16) {
            cycles.push_back(cycle);
        }
    }

    return cycles;
}


std::vector<Vec3> resampleCycle(const std::vector<Vec3> &cycle, float spacing) {
    std::vector<Vec3> newCycle;
    newCycle.push_back(cycle[0]);
    for (size_t i = 1; i < cycle.size()+1; ++i) {
        float segDistance = (cycle[i % cycle.size()] - newCycle.back()).length();
        if (segDistance <= 0.5f*spacing) continue;
        if (segDistance >= spacing) {
            float t = spacing/segDistance;
            Vec3 pos = t*cycle[i % cycle.size()] + (1-t) * newCycle.back();
            newCycle.push_back(pos);
            --i;
        }
    }
    return newCycle;
}

void saveSegmentToPLY(const char *path, const std::vector<std::array<Vec3Idx, 2>> &segments) {
    size_t vertsCount = 0;
    for (const auto &s : segments) {
        vertsCount = std::max(vertsCount, s[0].idx + 1);
        vertsCount = std::max(vertsCount, s[1].idx + 1);
    }

    std::vector<Vec3> verts(vertsCount);
    for (const auto &s : segments) {
        verts[s[0].idx] = s[0].v;
        verts[s[1].idx] = s[1].v;
    }

    std::ofstream plyFile(path);
    plyFile << "ply\nformat ascii 1.0\nelement vertex " << vertsCount << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "element edge " << segments.size() << "\nproperty int vertex1\nproperty int vertex2\n"
            << "end_header\n";
        
    for (const auto &v : verts) {
        plyFile << v.x << " " << v.y << " " << v.z << std::endl;
    }

    for (const auto &s : segments) {
        plyFile << s[0].idx << " " << s[1].idx << std::endl;
    }
}

void saveCyclesToPLY(const char *path, const std::vector<std::vector<Vec3>> &cycles, const Vec3 &offset) {
    size_t vertsCount = 0;
    for (const auto &cycle : cycles) {
        vertsCount += cycle.size();
    }

    std::ofstream plyFile(path);
    plyFile << "ply\nformat ascii 1.0\nelement vertex " << vertsCount << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "element edge " << vertsCount << "\nproperty int vertex1\nproperty int vertex2\n"
            << "end_header\n";
    
    for (const auto &cycle : cycles) {
        for (const auto v : cycle) {
            plyFile << v.x + offset.x << " " << v.y + offset.y << " " << v.z + offset.z << std::endl;
        }
    }

    size_t off = 0;
    for (const auto &cycle : cycles) {
        for (size_t i = 0; i < cycle.size()-1; ++i) {
            plyFile << off+i << " " << off+i+1 << std::endl;
        }
        plyFile << off + cycle.size()-1 << " " << off << std::endl;
        off += cycle.size();
    }
}