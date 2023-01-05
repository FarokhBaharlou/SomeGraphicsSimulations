#pragma once

#include <cstdint>

template<typename T>
class Vec3
{
public:
    Vec3() : x(0), y(0), z(0) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3 operator*(const T& r) const { return Vec3(x * r, y * r, z * r); }
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    template<typename U>
    Vec3 operator/(const Vec3<U>& v) const { return Vec3(x / v.x, y / v.y, z / v.z); }
    friend Vec3 operator/(const T r, const Vec3& v)
    {
        return Vec3(r / v.x, r / v.y, r / v.z);
    }
    const T& operator[](size_t i) const { return (&x)[i]; }
    T& operator[](size_t i) { return (&x)[i]; }
    T length2() const { return x * x + y * y + z * z; }
    friend Vec3 operator*(const float& r, const Vec3& v)
    {
        return Vec3(v.x * r, v.y * r, v.z * r);
    }
    friend std::ostream& operator<<(std::ostream& os, const Vec3<T>& v)
    {
        os << v.x << " " << v.y << " " << v.z << std::endl; return os;
    }
    T x, y, z;
};

template<typename T>
Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
    return Vec3<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

template<typename T>
T dot(const Vec3<T>& va, const Vec3<T>& vb)
{
    return va.x * vb.x + va.y * vb.y + va.z * vb.z;
}

template<typename T>
void normalize(Vec3<T>& vec)
{
    T len2 = vec.length2();
    if (len2 > 0)
    {
        T invLen = 1 / sqrt(len2);
        vec.x *= invLen, vec.y *= invLen, vec.z *= invLen;
    }
}

template<typename T>
class Matrix44
{
public:
    Matrix44() {}
    Matrix44(T m00, T m01, T m02, T m03,
             T m10, T m11, T m12, T m13,
             T m20, T m21, T m22, T m23,
             T m30, T m31, T m32, T m33)
    {
        m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
        m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
        m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
        m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
    }
    // ?????
    Matrix44 inverse() const { Matrix44 matInv = *this; return matInv; }
    T* operator[](size_t i) { return &m[i][0]; }
    const T* operator[](size_t i) const { return &m[i][0]; }
    T m[4][4] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
};

template<typename T>
void matVecMult(const Matrix44<T>& m, Vec3<T>& v)
{
    Vec3<T> vt;
    vt.x = v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
    vt.y = v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
    vt.z = v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2];

    v = vt;
}

template<typename T>
void matPointMult(const Matrix44<T>& m, Vec3<T>& p)
{
    Vec3<T> pt;
    pt.x = p.x * m[0][0] + p.y * m[1][0] + p.z * m[2][0] + m[3][0];
    pt.y = p.x * m[0][1] + p.y * m[1][1] + p.z * m[2][1] + m[3][1];
    pt.z = p.x * m[0][2] + p.y * m[1][2] + p.z * m[2][2] + m[3][2];
    T w = p.x * m[0][3] + p.y * m[1][3] + p.z * m[2][3] + m[3][3];
    if (w != 1)
    {
        pt.x /= w, pt.y /= w, pt.z /= w;
    }

    p = pt;
}

using Vec3f = Vec3<float>;
using Vec3b = Vec3<bool>;
using Vec3i = Vec3<int32_t>;
using Vec3ui = Vec3<uint32_t>;
using Matrix44f = Matrix44<float>;