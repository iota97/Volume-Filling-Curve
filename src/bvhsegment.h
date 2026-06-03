// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#include "mathutils.h"
#include <vector>
#include <array>
#include <unordered_set>
#include <iostream>
#include <fstream>

class BVHSegment {
public:
    BVHSegment(const std::vector<std::array<Vec3, 2>> &edges) {
        edgeID.resize(edges.size());
        for (size_t i = 0; i < edgeID.size(); ++i) {
            edgeID[i] = i;
        }
        edge.resize(edges.size());
        for (size_t i = 0; i < edge.size(); ++i) {
            Vec3 centroid = (edges[i][0]+edges[i][1])*0.5f;
            edge[i] = {edges[i][0], edges[i][1], centroid};
        }
        build();
    }

    Vec3 closestPoint(const Vec3 &point) const {
        Vec3 closest = edge[0].a;
        float minDistance2 = std::numeric_limits<float>::infinity();

        std::vector<size_t> stack{0};
        while (stack.size()) {
            size_t idx = stack.back();
            stack.pop_back();
            while (true) {
                if (!sphereAABBIntersection(point, minDistance2 + 1e-4, nodes[idx].minAABB, nodes[idx].maxAABB)) break;
                if (nodes[idx].isLeaf()) {
                   for (size_t offset = nodes[idx].firstLeft, last = offset + nodes[idx].edgeCount; offset < last; ++offset) {
                        size_t edgeIdx = edgeID[offset];
                        Vec3 edgeClosest = closestPointOnSegment(point, edge[edgeIdx].a, edge[edgeIdx].b);
                        float distance2 = (edgeClosest-point).length2();
                        if (distance2 < minDistance2) {
                            minDistance2 = distance2;
                            closest = edgeClosest;
                        }
                    }
                    break;
                }
                size_t leftIdx = nodes[idx].firstLeft;
                size_t rightIdx = leftIdx + 1;

                if (pointAABBIntersection(point, nodes[leftIdx].minAABB, nodes[leftIdx].maxAABB)) {
                    idx = leftIdx;
                    stack.push_back(rightIdx);
                } else if (pointAABBIntersection(point, nodes[rightIdx].minAABB, nodes[rightIdx].maxAABB)) {
                    idx = rightIdx;
                    stack.push_back(leftIdx);
                } else {
                    float leftDist2 = pointAABBDistance2(point, nodes[leftIdx].minAABB, nodes[leftIdx].maxAABB);
                    float rightDist2 = pointAABBDistance2(point, nodes[rightIdx].minAABB, nodes[rightIdx].maxAABB);
                    idx = leftDist2 < rightDist2 ? leftIdx : rightIdx;
                    stack.push_back(leftDist2 < rightDist2 ? rightIdx : leftIdx);
                }
            }
        }
        return closest;
    }

private:
    struct Node {
        Vec3 minAABB;
        Vec3 maxAABB;
        size_t edgeCount;
        size_t firstLeft;
        bool isLeaf() const { return edgeCount > 0; }
    };

    struct edge {
        Vec3 a, b, centroid;
    };

    std::vector<Node> nodes;
    std::vector<edge> edge;
    std::vector<size_t> edgeID;

    void build() {
        nodes.push_back(Node());
        nodes[0].firstLeft = 0;
        nodes[0].edgeCount = edgeID.size();
        updateLastBounds();
        subdivideNodes();
    }

    void updateLastBounds() {
        nodes.back().minAABB = Vec3(std::numeric_limits<float>::infinity());
        nodes.back().maxAABB = Vec3(-std::numeric_limits<float>::infinity());

        for (size_t offset = nodes.back().firstLeft, last = offset + nodes.back().edgeCount; offset < last; ++offset) {
            size_t edgeIdx = edgeID[offset];
            for (uint8_t axis = 0; axis < 3; ++axis) {
                nodes.back().minAABB[axis] = std::min(nodes.back().minAABB[axis], edge[edgeIdx].a[axis]);
                nodes.back().minAABB[axis] = std::min(nodes.back().minAABB[axis], edge[edgeIdx].b[axis]);
                nodes.back().maxAABB[axis] = std::max(nodes.back().maxAABB[axis], edge[edgeIdx].a[axis]);
                nodes.back().maxAABB[axis] = std::max(nodes.back().maxAABB[axis], edge[edgeIdx].b[axis]);
            }
        }
    }

    void subdivideNodes() {
        std::vector<size_t> stack{0};
        while (stack.size()) {
            size_t idx = stack.back();
            stack.pop_back();

            if (nodes[idx].edgeCount > 2) {
                Vec3 extent = nodes[idx].maxAABB - nodes[idx].minAABB;
                uint8_t axis = 0;
                for (uint8_t i = 1; i < 3; ++i) {
                    if (extent[i] > extent[axis])
                        axis = i;
                }

                float split = 0;
                for (size_t i = 0; i < nodes[idx].edgeCount; ++i) {
                    split += edge[edgeID[nodes[idx].firstLeft+i]].centroid[axis];
                }
                split /= nodes[idx].edgeCount;

                size_t i = nodes[idx].firstLeft;
                size_t j = i + nodes[idx].edgeCount - 1;
                while (i <= j) {
                    if (edge[edgeID[i]].centroid[axis] < split) {
                        i += 1;
                    } else {
                        std::swap(edgeID[i], edgeID[j]);
                        if (!j) break;
                        j -= 1;
                    }
                }

                size_t leftCount = i - nodes[idx].firstLeft;
                if (leftCount > 0 && leftCount < nodes[idx].edgeCount) {
                    Node leftNode;
                    leftNode.firstLeft = nodes[idx].firstLeft;
                    leftNode.edgeCount = leftCount;
                    nodes[idx].firstLeft = nodes.size();
                    stack.push_back(nodes.size());

                    Node rightNode;
                    rightNode.firstLeft = i;
                    rightNode.edgeCount = nodes[idx].edgeCount - leftCount;
                    nodes[idx].edgeCount = 0;
                    stack.push_back(nodes.size()+1);

                    nodes.push_back(rightNode);
                    updateLastBounds();
                    nodes.push_back(leftNode);
                    updateLastBounds();
                }
            }
        }
    }

    static Vec3 closestPointOnSegment(const Vec3 &p, const Vec3 &a, const Vec3 &b) {
        Vec3 edge = b-a;
        Vec3 deltaStart = p-a;
        float l = deltaStart.dot(edge);
        float s = edge.length2();
        if (l > 0 && l < s) {
            return a + l/s * edge;
        } else {
            if (l <= 0) {
                return a;
            } else {
                return b;
            }
        }
    }

    static bool pointAABBIntersection(const Vec3 &point, const Vec3 &minAABB, const Vec3 &maxAABB) {
        for (uint8_t axis = 0; axis < 3; ++axis) {
            if (point[axis] < minAABB[axis] || point[axis] > maxAABB[axis]) 
                return false;
        }
        return true;
    }

    static bool sphereAABBIntersection(const Vec3 &center, float radius2, const Vec3 &minAABB, const Vec3 &maxAABB) {
        float dist2 = 0.0f;
        for (uint8_t axis = 0; axis < 3; ++axis) {
            float m = minAABB[axis]-center[axis], M = center[axis]-maxAABB[axis];
            dist2 += (m > 0)*m*m + (M > 0)*M*M;
        }
        return dist2 <= radius2;
    }

    static float pointAABBDistance2(const Vec3 &point, const Vec3 &minAABB, const Vec3 &maxAABB) {
        Vec3 closest(point);
        for (uint8_t axis = 0; axis < 3; ++axis) {
            if (point[axis] < minAABB[axis]) {
                closest[axis] = minAABB[axis];
            } else if (point[axis] > maxAABB[axis]) {
                closest[axis] = maxAABB[axis];
            }
        }
        return (point-closest).length2();
    }
};