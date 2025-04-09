#pragma once

#include "src/CGALTypedefs.h"

#include <unordered_map>

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "./FaceFiber.h"
#include "./DisjointSet.h"

#include "./CGALTypedefs.h"

class Data
{
  public:
    Data() {}

    // Min/max range coordinates
    float minF, maxF;
    float minG, maxG;

    // Min/max domain coordinates
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;

    std::vector<float> currentFiberPoint;
    void printSheetHistogram();

    // Compute the min/max F, G and X, Y, Z coordinates
    void computeMinMaxRangeDomainCoordinates();

    // Which of the tets in the tetrahedra array contain a fiber
    std::vector<bool> tetsWithFibers;

    std::string longnameF, longnameG, units;

    std::vector<std::vector<size_t>> tetrahedra;
    std::vector<std::vector<float>> vertexDomainCoordinates;

    std::vector<float> vertexCoordinatesF;
    std::vector<float> vertexCoordinatesG;

    std::string fibersFile = "./fibers.vtp";
    std::vector<FaceFiberPoint> faceFibers;
    void saveFibers();

    void generatefFaceFibersForSheet(const int sheetId, const int numberOfFiberPoints);
    void generatefFaceFibersForSheets(const int sheetOutputCount, const int numberOfFiberPoints, const std::string);

    void computeTetExitPoints(const float, const float, const std::vector<float> = {1,1,1});
    void computeTetExitPointsNew(const float, const float, const std::vector<float> = {1,1,1});
void computeTetExitPointsNewNew(const float, const float, const bool, const int reebSheetIdOnly = -1, const std::vector<float> = {1,1,1});



    void readData(const std::string&, const float&);
    void readDataVTU(const std::string&, const float&);



    // Reeb stuff relted content
    Arrangement_2 arr;

    // Ideally we do want this with a proper data structure with a STAR, then there's no write conflits
    // Make sure to always keep the edge (u, v) such that u < v in index value,
    std::map<std::pair<int, int>, std::set<int>> upperLink;
    std::map<std::pair<int, int>, std::set<int>> lowerLink;



    //
    // <Reeb space related stuff>
    //

    // The points of the line segments that define the arrangement, does not include new intersection points
    std::vector<Point_2> arrangementPoints;

    // The inverse map of arrangementPoints, returns the index of a point
    std::map<Point_2, int> arrangementPointsIdices;

    // The integer ID of every face in the arrangement
    std::map<Arrangement_2::Face_const_handle, int> arrangementFacesIdices;
    std::vector<Arrangement_2::Face_const_handle> arrangementIndexToFace;

    // Tell us which triangles (as sets of IDs) are connected (part of a tetrahedron)
    std::set<std::pair<std::set<int>, std::set<int>>> connectedTriangles;



    //std::set<std::pair<std::set<int>, std::vector<int>>> adjacentTriangles;
    //std::map<std::set<int>, std::vector<std::set<int>>> adjacentTriangles;

    // Triangle to index
    std::vector<std::set<int>> indexToTriangle;

    // Index to triangle
    //std::map<std::set<int>, int> triangleToIndex;
    std::unordered_map<std::set<int>, int, MyHash<std::set<int>>> triangleToIndex;

    std::vector<std::vector<int>> adjacentTrianglesIndex;


    // The jacobi type of each edge 0 for definite, 1 for regular and n > 1 for indefinite of type n
    std::unordered_map<std::pair<int, int>, int, MyHash<std::pair<int, int>>> jacobiType;

    // The number of fiber components for each preimage graph, more of a utility thing
    //std::vector<int> arrangementFiberComponents;



    // Big memory usage is here


    // The connected components of the preimage graph for each face in the arrangement
    std::vector<DisjointSet<int>> preimageGraphs;

    // The actual Reeb space, map from a connected component of a preimage graph to a sheet.
    DisjointSet<std::pair<int, int>> reebSpace;

    // For each faceID we ahve a number of triangle seeds, which are given in fiberComponentId
    std::vector<std::vector<std::pair<int, int>>> fiberSeeds;

    std::unordered_map<int, CartesianPolygon_2> sheetPolygon;

    std::map<int, double> sheetArea;

    std::unordered_set<int> incompleteSheets;

    // There's also this vector in plotwidget
    //std::vector<QPolygon> arrangementPolygons;
    

    // The vertices of the correspondence graph H are pairs of int <faceID, fiber component ID>
    // I'd like to map between them and int indices for the disjoint set
    //std::vector<std::pair<int, int>> indexToVertexH;
    //std::unordered_map<std::pair<int, int>, int, MyHash<std::pair<int, int>>> vertexHtoIndex;

    // Maps the IDs of the Reeb space sheets to consequitive integers, useful for colouring things
    std::unordered_map<int, int> sheetToColour;

    // Data structure to query the arrangement
    // Given a point, which face is it in?
    std::unique_ptr<Point_location> pl;  // nullptr by default

    void sortVertices();
    void printMesh();


    //
    // </Reeb space related stuff>
    //

    // Colour map for Reeb space sheets and fibers
    const std::vector<std::vector<float>> fiberColours = {
        {1.0f, 0.0f, 0.0f, 0.392f},    // Vivid Red
        {0.0f, 1.0f, 0.0f, 0.392f},    // Bright Green
        {0.0f, 0.0f, 1.0f, 0.392f},    // Pure Blue
        {1.0f, 1.0f, 0.0f, 0.392f},    // Bright Yellow
        {0.0f, 1.0f, 1.0f, 0.392f},    // Cyan
        {1.0f, 0.0f, 1.0f, 0.392f},    // Magenta
        {0.58f, 0.0f, 1.0f, 0.392f},   // Deep Purple
        {0.0f, 0.45f, 0.7f, 0.392f},   // Ocean Blue
        {1.0f, 0.5f, 0.0f, 0.392f},    // Orange
        {0.0f, 0.6f, 0.5f, 0.392f},    // Teal
        {1.0f, 0.84f, 0.0f, 0.392f},   // Gold
        {0.85f, 0.4f, 0.55f, 0.392f},  // Mauve
        {0.4f, 0.8f, 1.0f, 0.392f},    // Sky Blue
        {0.2f, 0.8f, 0.2f, 0.392f},    // Leaf Green
        {0.9f, 0.3f, 0.3f, 0.392f},    // Coral Red
        {0.6f, 0.6f, 0.0f, 0.392f}     // Olive
    };
};
