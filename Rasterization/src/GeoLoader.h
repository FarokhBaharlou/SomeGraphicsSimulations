#pragma once

#include <fstream>
#include <sstream>
#include "Utility.h"

void loadGeoFile(const char* file, uint32_t& numFaces, std::unique_ptr<Vec3f[]>& verts, std::unique_ptr<Vec2f[]>& st, std::unique_ptr<uint32_t[]>& vertsIndex)
{
    std::ifstream ifs;
    try
    {
        ifs.open(file);
        if (ifs.fail()) throw;
        std::stringstream ss;
        ss << ifs.rdbuf();
        ss >> numFaces;
        uint32_t vertsIndexArraySize = 0;
        // reading face index array
        for (uint32_t i = 0; i < numFaces; ++i)
        {
            uint32_t tmp;
            ss >> tmp; //faceIndex[i];
            vertsIndexArraySize += tmp; //faceIndex[i];
        }
        vertsIndex = std::unique_ptr<uint32_t[]>(new uint32_t[vertsIndexArraySize]);
        uint32_t vertsArraySize = 0;
        // reading vertex index array
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i)
        {
            ss >> vertsIndex[i];
            if (vertsIndex[i] > vertsArraySize) vertsArraySize = vertsIndex[i];
        }
        vertsArraySize += 1;
        // reading vertices
        verts = std::unique_ptr<Vec3f[]>(new Vec3f[vertsArraySize]);
        for (uint32_t i = 0; i < vertsArraySize; ++i)
        {
            ss >> verts[i].x >> verts[i].y >> verts[i].z;
        }
        // reading normals
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i)
        {
            Vec3f normal;
            ss >> normal.x >> normal.y >> normal.z;
        }
        // reading st coordinates
        st = std::unique_ptr<Vec2f[]>(new Vec2f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i)
        {
            ss >> st[i].x >> st[i].y;
        }
    }
    catch (...)
    {
        ifs.close();
    }
    ifs.close();
}