// MIT License
// Copyright 2026 Giovanni Cocco and Inria

#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <limits>
#include <unordered_map>
#include <ostream>
#include <array>

struct Vec3 {
    float x, y, z;
    Vec3() = default;
    Vec3(float v) : x(v), y(v), z(v) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    Vec3(const float *data, uint32_t idx) : x(data[3*idx+0]), y(data[3*idx+1]), z(data[3*idx+2]) {}
    float &operator[](uint8_t axis) { return axis ? (axis == 1 ? y : z) : x; }
    float operator[](uint8_t axis) const { return axis ? (axis == 1 ? y : z) : x; }
    Vec3 operator+(const Vec3& o) const { Vec3 res; res.x = x + o.x; res.y = y + o.y; res.z = z + o.z; return res; }
    void operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; }
    void operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; }
    Vec3 operator-(const Vec3& o) const { Vec3 res; res.x = x - o.x; res.y = y - o.y; res.z = z - o.z; return res; }
    Vec3 operator-() const { Vec3 res; res.x = -x; res.y = -y; res.z = -z; return res; }
    Vec3 operator*(float s) const { Vec3 res; res.x = x*s; res.y = y*s; res.z = z*s; return res; }
    Vec3 operator/(float s) const { Vec3 res; res.x = x/s; res.y = y/s; res.z = z/s; return res; }
    float dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
    float length() const { return sqrt(x*x + y*y + z*z); }
    float length2() const { return x*x + y*y + z*z; }
    Vec3 normalize() const { Vec3 res; float l = length(); res.x = x/l; res.y = y/l; res.z = z/l; return l > 0 ? res : Vec3(1, 0, 0); }
    Vec3 project(const Vec3& n) const { return n * dot(n); }
    Vec3 orthonormalize(const Vec3& n) const { return (*this - project(n)).normalize(); }
    Vec3 normalToTangent() const { float s = z > 0.0f ? 1.0f : -1.0f, a = -1.0f/(s+z), b = x*y*a; return Vec3(1.0f+s*x*x*a, s*b, -s*x); }
    Vec3 cross(const Vec3 &rhs) const { return Vec3{y * rhs.z - z * rhs.y, -(x * rhs.z - z * rhs.x), x * rhs.y - y * rhs.x}; }
    float normSqr() const { return dot(*this); }
    bool operator==(const Vec3& o) const { return x == o.x && y == o.y && z == o.z; }
    bool operator!=(const Vec3& o) const { return x != o.x || y != o.y || z != o.z; }
    bool operator<(const Vec3& o) const { if (x < o.x) return true; if (x > o.x) return false; if (y < o.y) return true; if (y > o.y) return false; if (z < o.z) return true; return false; }
    Vec3 clamp(const Vec3& min, const Vec3& max) { return Vec3{ std::min(max.x, std::max(min.x, x)), std::min(max.y, std::max(min.y, y)), std::min(max.z, std::max(min.z, z)) }; }
};
inline Vec3 operator*(float s, const Vec3& v) { return v*s; }

struct Mat3 {
    Vec3 rows[3];
    Mat3(float a0, float a1, float a2, float b0, float b1, float b2, float c0, float c1, float c2) { rows[0] = Vec3(a0, a1, a2); rows[1] = Vec3(b0, b1, b2); rows[2] = Vec3(c0, c1, c2); }
};

inline Vec3 operator*(const Mat3 &m, const Vec3& v) { 
    return Vec3(m.rows[0].dot(v), m.rows[1].dot(v), m.rows[2].dot(v)); 
}

inline Vec3 operator*(const Vec3& v, const Mat3 &m) { 
    return m.rows[0]*v.x + m.rows[1]*v.y + m.rows[2]*v.z;
}

template <> 
struct std::hash<Vec3> {
    std::size_t operator()(const Vec3& v) const {
    return ((std::hash<float>()(v.x)
             ^ (std::hash<float>()(v.y) << 1)) >> 1)
             ^ (std::hash<float>()(v.z) << 1);
    }
};

inline std::ostream& operator<<(std::ostream& stream, const Vec3& v) {
    return stream << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

inline uint32_t randomUintFromUint(uint32_t v) {
    uint32_t state = 747796405*v + 2891336453;
    uint32_t word = ((state >> ((state >> 28) + 4)) ^ state) * 277803737;
    return (word >> 22) ^ word;
}

inline uint32_t randomUintFromUint(uint32_t v, uint32_t c) {
    for (uint32_t i = 0; i < c; ++i) {
        v = randomUintFromUint(v);
    }
    return v;
}

inline float randomFloatFromUint(uint32_t v) {
    return float(randomUintFromUint(v))/std::numeric_limits<uint32_t>::max();
}

inline Vec3 randomDirectionFromUint(uint32_t v) {
    return (2*Vec3(randomFloatFromUint(v), randomFloatFromUint(13*v+5), randomFloatFromUint(17*v+7))-1).normalize();
}

inline void printProgress(double percentage) {
    static uint8_t last = 0xff;
    if (uint8_t(percentage*100.0f) != last) {
        last = uint8_t(percentage*100.0f);
        const char* bar = "||||||||||||||||||||||||||||||||||||||||";
        const uint32_t barWidth = std::strlen(bar);
        uint32_t val = (uint32_t) (percentage * 100);
        uint32_t lpad = (uint32_t) (percentage * barWidth);
        uint32_t rpad = barWidth - lpad;
        printf("\r%3d%% [%.*s%*s]", val, lpad, bar, rpad, "");
        fflush(stdout);
    }
    if (percentage == 1.0) {
        printf("\n");
        last = 0xff;
    }
}

inline Vec3 barycentricCoord(const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c) {
    Vec3 v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = 1.0f/(d00 * d11 - d01 * d01);
    float y = (d11 * d20 - d01 * d21) * denom;
    float z = (d00 * d21 - d01 * d20) * denom;
    float x = 1.0f - y - z;
    return {x, y, z};
}

inline Vec3 triangleNormal(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
    return (a-c).cross(b-c).normalize();
}

inline float cellSizeFromLayerHeight(float h) { 
    return h*0.375f;
}

inline void alignPhases(const Vec3 &p_i, const Vec3 &p_j, const Vec3 &d_i, const Vec3 &d_j, float f_i, float f_j, float phi_j, float &x, float &y, float w = 1.0f) {
    Vec3 p_ij = (p_i+p_j)*0.5f;
    float a = 2.0f*float(M_PI)*f_i * (p_ij-p_i).dot(d_i);
    float b = 2.0f*float(M_PI)*f_j * (p_ij-p_j).dot(d_j) + phi_j;
    float angle = d_i.dot(d_j) > 0.0f ? b - a : -b - a;
    x += w*fabsf(d_i.dot(d_j))*cosf(angle);
    y += w*fabsf(d_i.dot(d_j))*sinf(angle);
}

inline std::array<uint8_t, 3> viridis(float t) {
    static const float c0[3]{ 0.2777273272234177, 0.005407344544966578, 0.3340998053353061 };
    static const float c1[3]{ 0.1050930431085774, 1.404613529898575, 1.384590162594685 };
    static const float c2[3]{ -0.3308618287255563, 0.214847559468213, 0.09509516302823659 };
    static const float c3[3]{ -4.634230498983486, -5.799100973351585, -19.33244095627987 };
    static const float c4[3]{ 6.228269936347081, 14.17993336680509, 56.69055260068105 };
    static const float c5[3]{ 4.776384997670288, -13.74514537774601, -65.35303263337234 };
    static const float c6[3]{ -5.435455855934631, 4.645852612178535, 26.3124352495832 };

    uint8_t r = 255 * (c0[0] + t * (c1[0] + t * (c2[0] + t * (c3[0] + t * (c4[0] + t * (c5[0] + t * c6[0]))))));
    uint8_t g = 255 * (c0[1] + t * (c1[1] + t * (c2[1] + t * (c3[1] + t * (c4[1] + t * (c5[1] + t * c6[1]))))));
    uint8_t b = 255 * (c0[2] + t * (c1[2] + t * (c2[2] + t * (c3[2] + t * (c4[2] + t * (c5[2] + t * c6[2]))))));
    return std::array<uint8_t, 3>{r, g, b};
}

inline std::array<uint8_t, 3> turbo(float t) {
    static const float c0[3]{ 0.1140890109226559, 0.06288340699912215, 0.2248337216805064 };
    static const float c1[3]{ 6.716419496985708, 3.182286745507602, 7.571581586103393 };
    static const float c2[3]{ -66.09402360453038, -4.9279827041226, -10.09439367561635 };
    static const float c3[3]{ 228.7660791526501, 25.04986699771073, -91.54105330182436 };
    static const float c4[3]{ -334.8351565777451, -69.31749712757485, 288.5858850615712 };
    static const float c5[3]{ 218.7637218434795, 67.52150567819112, -305.2045772184957 };
    static const float c6[3]{ -52.88903478218835, -21.54527364654712, 110.5174647748972 };

    uint8_t r = 255 * std::max(0.0f, std::min(1.0f, c0[0] + t * (c1[0] + t * (c2[0] + t * (c3[0] + t * (c4[0] + t * (c5[0] + t * c6[0])))))));
    uint8_t g = 255 * std::max(0.0f, std::min(1.0f, c0[1] + t * (c1[1] + t * (c2[1] + t * (c3[1] + t * (c4[1] + t * (c5[1] + t * c6[1])))))));
    uint8_t b = 255 * std::max(0.0f, std::min(1.0f, c0[2] + t * (c1[2] + t * (c2[2] + t * (c3[2] + t * (c4[2] + t * (c5[2] + t * c6[2])))))));
    return std::array<uint8_t, 3>{r, g, b};
}

inline std::array<uint8_t, 3> twilight(float t) {
    static const float c0[3] { 0.996106, 0.851653, 0.940566 };
    static const float c1[3] { -6.529620, -0.183448, -3.940750 };
    static const float c2[3] { 40.899579, -7.894242, 38.569228 };
    static const float c3[3] { -155.212979, 4.404793, -167.925730 };
    static const float c4[3] { 296.687222, 24.084913, 315.087856 };
    static const float c5[3] { -261.270519, -29.995422, -266.972991 };
    static const float c6[3] { 85.335349, 9.602600, 85.227117 };

    uint8_t r = 255 * (c0[0] + t * (c1[0] + t * (c2[0] + t * (c3[0] + t * (c4[0] + t * (c5[0] + t * c6[0]))))));
    uint8_t g = 255 * (c0[1] + t * (c1[1] + t * (c2[1] + t * (c3[1] + t * (c4[1] + t * (c5[1] + t * c6[1]))))));
    uint8_t b = 255 * (c0[2] + t * (c1[2] + t * (c2[2] + t * (c3[2] + t * (c4[2] + t * (c5[2] + t * c6[2]))))));
    return std::array<uint8_t, 3>{r, g, b};
}

struct Vec3Idx {
    Vec3 v;
    size_t idx;
};