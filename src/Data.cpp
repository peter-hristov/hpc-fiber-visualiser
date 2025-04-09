#include "Data.h"
#include "./Timer.h"
#include "./DisjointSet.h"
#include "src/CGALTypedefs.h"
#include "src/FaceFiber.h"

#include <CGAL/enum.h>
#include <filesystem>
#include <cassert>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <queue>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <random>
#include <ranges>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>

#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

using namespace std;

double randomPerturbation(double epsilon) {
    static std::mt19937 gen(std::random_device{}()); // Seeded generator
    std::uniform_real_distribution<double> dist(-epsilon, epsilon);
    return dist(gen);
}

void Data::saveFibers()
{
    std::cout << "Saving fibers in " << this->fibersFile << std::endl;
    //std::cout << "The fiber has size " << this->faceFibers.size() << std::endl;  

    // 1. Create the points
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto idArray = vtkSmartPointer<vtkIntArray>::New();
    auto colourArray = vtkSmartPointer<vtkFloatArray>::New();

    idArray->SetName("SheetId");
    idArray->SetNumberOfComponents(1);

    colourArray->SetName("Colour");
    colourArray->SetNumberOfComponents(3);

    for (const FaceFiberPoint &p : this->faceFibers)
    {
        // Insert points and corresponding IDs
        points->InsertNextPoint(p.point[0], p.point[1], p.point[2]);
        idArray->InsertNextValue(p.sheetId);

        const vector<float> sheetColour = this->fiberColours[this->sheetToColour[p.sheetId] % this->fiberColours.size()];
        float color[3] = {sheetColour[0], sheetColour[1], sheetColour[2]};
        colourArray->InsertNextTuple(color);
    }

    // 3. Create the cells (wrap polyline in cell array)
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 1 ; i < this->faceFibers.size() ; i+=2)
    {
        //cout << i << " " << this->faceFibers[i-1].sheetId << " " << this->faceFibers[i].sheetId << endl;
        //cout << i << " " << this->faceFibers[i-1].triangleId << " " << this->faceFibers[i].triangleId << endl;
        //printf("(%f, %f, %f) -> (%f, %f, %f)\n", this->faceFibers[i-1].point[0], this->faceFibers[i-1].point[1], this->faceFibers[i-1].point[2], this->faceFibers[i].point[0], this->faceFibers[i].point[1], this->faceFibers[i].point[2]);

        if (this->faceFibers[i-1].sheetId == this->faceFibers[i].sheetId)
        {
            // One edge segment
            auto polyLine = vtkSmartPointer<vtkPolyLine>::New();
            polyLine->GetPointIds()->SetNumberOfIds(2);
            polyLine->GetPointIds()->SetId(0, i-1);
            polyLine->GetPointIds()->SetId(1, i);

            cells->InsertNextCell(polyLine);
        }
    }

    // 4. Create the polydata object
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(cells);

    // 5. Attach the VertexID array to the point data
    polyData->GetPointData()->AddArray(idArray);
    polyData->GetPointData()->AddArray(colourArray);
    polyData->GetPointData()->SetScalars(colourArray);  // optional: for coloring

    // 6. Write to .vtp file (XML format)
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(this->fibersFile.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}

void Data::generatefFaceFibersForSheets(const int sheetOutputCount, const int numberOfFiberPoints, const std::string folderPath)
{
    namespace fs = std::filesystem;

    //string folderPath = "./sheetFibers";
    fs::path folderPathFs(folderPath);
    if (!fs::exists(folderPathFs)) 
    {
        fs::create_directory(folderPathFs);
    }

    for (const auto &[sheetId, colourId] : this->sheetToColour)
    {
        if (this->incompleteSheets.contains(sheetId))
        {
            printf("Skipping fiber %d, it's incomplete.",  sheetId);
        }

        if (colourId > sheetOutputCount || this->incompleteSheets.contains(sheetId))
        {
            continue;
        }

        std::cout << "-------------------------------------------------------------------------------------------- Generating fibers for sheet " << sheetId << "..." << std::endl;
        this->generatefFaceFibersForSheet(sheetId, numberOfFiberPoints);

        //std::cout << "Saving fibers..." << std::endl;
        this->fibersFile = folderPathFs.string() + "/fibers_" + std::to_string(sheetId) + ".vtp";
        this->saveFibers();
        this->faceFibers.clear();
    }

    this->fibersFile = "./fibers.vtp";

}

void Data::generatefFaceFibersForSheet(const int sheetId, const int numberOfFiberPoints)
{
    CartesianPolygon_2 &polygon = this->sheetPolygon[sheetId];

    if (polygon.size() == 0)
    {
        return;
    }


    // Compute the controid so that we can pull all verties towards it
    CartesianPoint centroid = CGAL::centroid(polygon.vertices_begin(), polygon.vertices_end());

    // To make sure we don't write a comma at the end of the array

    vector<vector<float>> fiberPoints;

    // If need only one, get it at the center
    if (numberOfFiberPoints == 1)
    {
        fiberPoints.push_back({(float)centroid.x(), (float)centroid.y()});
        this->computeTetExitPointsNewNew(fiberPoints[0][0], fiberPoints[0][1], false, sheetId);

        return;
    }


    // If we need more, sample along the boundary
    for (const CartesianPoint &point : polygon) 
    {
        // Get point from CGAL (and convert to double )
        float u = point.x();
        float v = point.y();

        // Interpolate closer to the centroid to make sure we are in the sheet ( if the sheet is "convex enough")
        float alpha = 0.2;
        u = (1 - alpha) * u + alpha * centroid.x();
        v = (1 - alpha) * v + alpha * centroid.y();

        fiberPoints.push_back({u, v});
        //fiberPoints.push_back({centroid.x(), centroid.y()});
    }


    // Calculate step size we only want some of the fiber points, not all
    double step = static_cast<double>(fiberPoints.size() - 1) / (numberOfFiberPoints - 1);

    for (int i = 0; i < numberOfFiberPoints; ++i) 
    {
        int index = static_cast<int>(i * step);
        this->computeTetExitPointsNewNew(fiberPoints[index][0], fiberPoints[index][1], false, sheetId);
    }
}


void Data::printSheetHistogram()
{
    Timer::start();
    if (this->currentFiberPoint.size() == 0)
    {
        return;
    }

    assert (this->currentFiberPoint.size() == 2);

    set<int> intersectedSheets;
    CartesianLine line(CartesianPoint(0, this->currentFiberPoint[0]), CartesianPoint(1, this->currentFiberPoint[1])); // Line through (0, 1) and (1, 0)
    for (const auto &[sheetId, polygon] : this->sheetPolygon)
    {
        for (const CartesianSegment& segment : polygon.edges()) 
        {
            if (CGAL::do_intersect(segment, line)) 
            {
                intersectedSheets.insert(sheetId);
            }
        }
    }

    std::cout << "Intersected sheets: ";
    for (const auto &sheetId : intersectedSheets)
    {
        std::cout << sheetId << " (a = )" << this->sheetArea[sheetId] << std::endl;

    }
    std::cout << "\n";
    Timer::stop("Computed face intersected by the line  :");
}


void Data::computeTetExitPointsNewNew(const float u, const float v, const bool clearFibers, const int reebSheetIdOnly, const std::vector<float> color)
{
    if (true == clearFibers)
    {
        this->faceFibers.clear();
    }

    //
    // Get the ID of the face we are intersecting
    //

    //Timer::start();

    // Store the current fiber point
    this->currentFiberPoint = {u, v};

    // The query point (u, v)
    Point_2 query_point(u, v);


    // Locate the point in the arrangement
    CGAL::Object result = this->pl->locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    int currentFaceID = 0;

    if (CGAL::assign(face, result)) 
    {
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } else 
    {
        assert(false);
    }
    //Timer::stop("Computed active arrangement face       :");



    //Timer::start();

    // The sizes of these data structures are linear in the size of the fiber, not an issue
    std::queue<int> bfsQueue;
    // This also acts as the visited array
    std::unordered_map<int, int> triangleSheetId;
    // Needed to close loops closed fibers, hack because we are not visiting tets, but triangles
    std::unordered_set<pair<int, int>, MyHash<pair<int, int>>> activeAdjacentTrianglesConnected;
    // Cache barycentric coordintes, they are expensive to compute
    std::unordered_map<int, std::array<double, 3>> triangleBarycentricCoordinates;

    if (reebSheetIdOnly == -1)
    {
        std::cout << "There are " << this->fiberSeeds[currentFaceID].size() << " fiber components with (sheet IDs, sorted IDs): ";
    }

    //vector<int> sheetIds;
    for (const auto &[triangleId, fiberComponentId] : this->fiberSeeds[currentFaceID])
    {
        const int sheetId = this->reebSpace.findTriangle({currentFaceID, fiberComponentId});

        if (reebSheetIdOnly == -1 || sheetId == reebSheetIdOnly)
        {
            bfsQueue.push(triangleId);
            //cout << "Adding seed \n";
        }

        //const int sheetId = this->reebSpace.findTriangle({currentFaceID, fiberComponentId});
        //triangleColour[triangleId] = this->sheetToColour[sheetId] % this->fiberColours.size();
        triangleSheetId[triangleId] = sheetId;

        //sheetIds.push_back(sheetId);

        if (reebSheetIdOnly == -1)
        {
            printf("(%d, %d) ", sheetId, this->sheetToColour[sheetId]);
        }
    }
    //std::cout << std::endl;


    // Define query point
    CartesianPoint P(u, v);

    while (false == bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        const int currentSheeId = triangleSheetId[currentTriangleId];
        bfsQueue.pop();

        const vector<float> sheetColour = this->fiberColours[this->sheetToColour[currentSheeId] % this->fiberColours.size()];

        const set<int> triangleUnpacked = this->indexToTriangle[currentTriangleId];
        const vector<int> triangleIndices = std::vector<int>(triangleUnpacked.begin(), triangleUnpacked.end());

        std::array<double, 3> barycentricCoordinatesCurrent;

        if (triangleBarycentricCoordinates.contains(currentTriangleId))
        {
            barycentricCoordinatesCurrent = triangleBarycentricCoordinates[currentTriangleId];
        }
        else
        {
            const CartesianPoint A(this->vertexCoordinatesF[triangleIndices[0]], this->vertexCoordinatesG[triangleIndices[0]]);
            const CartesianPoint B(this->vertexCoordinatesF[triangleIndices[1]], this->vertexCoordinatesG[triangleIndices[1]]);
            const CartesianPoint C(this->vertexCoordinatesF[triangleIndices[2]], this->vertexCoordinatesG[triangleIndices[2]]);
            CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesCurrent.begin());
            triangleBarycentricCoordinates[currentTriangleId] = barycentricCoordinatesCurrent;
        }

        // Sanity check
        assert(barycentricCoordinatesCurrent[0] > 0 && barycentricCoordinatesCurrent[1] > 0 && barycentricCoordinatesCurrent[2] > 0);

        // Look at the neighbours
        for (const int &neighbourTriagleId : this->adjacentTrianglesIndex[currentTriangleId])
        {
            if (neighbourTriagleId == currentTriangleId)
            {
                continue;
            }

            // We can't skip visited neighbours, because there may be a fiber between us (completing a circle)
            //if (triangleColour.contains(neighbourTriagle)) { continue; }

            const set<int> triangle2Unpacked = this->indexToTriangle[neighbourTriagleId];
            const vector<int> triangle2Indices = std::vector<int>(triangle2Unpacked.begin(), triangle2Unpacked.end());



            // The neighbour is active if we've already seen it
            bool isActive = triangleSheetId.contains(neighbourTriagleId);

            // Or if the image of the triangle contains the fiber points
            // We use a fast test to avoid having to use barycentric coordinates all the time
            if (false == isActive)
            {
                CartesianPoint A(this->vertexCoordinatesF[triangle2Indices[0]], this->vertexCoordinatesG[triangle2Indices[0]]);
                CartesianPoint B(this->vertexCoordinatesF[triangle2Indices[1]], this->vertexCoordinatesG[triangle2Indices[1]]);
                CartesianPoint C(this->vertexCoordinatesF[triangle2Indices[2]], this->vertexCoordinatesG[triangle2Indices[2]]);

                std::vector<CartesianPoint> triangle = {A, B, C};
                const auto result = CGAL::bounded_side_2(triangle.begin(), triangle.end(), P);

                isActive = (result == CGAL::ON_BOUNDED_SIDE);
            }

            // Determine if the triangle is active
            if (isActive)
            {
                // Only add the neighbour if we have not already visited it
                if (false == triangleSheetId.contains(neighbourTriagleId))
                {
                    // BFS things
                    bfsQueue.push(neighbourTriagleId);
                    triangleSheetId[neighbourTriagleId] = currentSheeId;
                }

                // Even if we have aleady added a neighbour, maybe there still isn't a fiber between us (for finishing loops)

                // At this point, we know that both us and we neighbour are active, is there already a fiber between us? Then skip
                if (activeAdjacentTrianglesConnected.contains({currentTriangleId, neighbourTriagleId}))
                {
                    continue;
                }
                else
                {
                    activeAdjacentTrianglesConnected.insert({currentTriangleId, neighbourTriagleId});
                    activeAdjacentTrianglesConnected.insert({neighbourTriagleId, currentTriangleId});
                }




                // Compute barycentric coordinates for drawing
                std::array<double, 3> barycentricCoordinatesNeighbour;
                if (triangleBarycentricCoordinates.contains(neighbourTriagleId))
                {
                    barycentricCoordinatesNeighbour = triangleBarycentricCoordinates[neighbourTriagleId];
                }
                else
                {
                    CartesianPoint A(this->vertexCoordinatesF[triangle2Indices[0]], this->vertexCoordinatesG[triangle2Indices[0]]);
                    CartesianPoint B(this->vertexCoordinatesF[triangle2Indices[1]], this->vertexCoordinatesG[triangle2Indices[1]]);
                    CartesianPoint C(this->vertexCoordinatesF[triangle2Indices[2]], this->vertexCoordinatesG[triangle2Indices[2]]);

                    CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesNeighbour.begin());
                    triangleBarycentricCoordinates[neighbourTriagleId] = barycentricCoordinatesNeighbour;
                }

                


                //
                // Add a fiber segment
                //
                FaceFiberPoint fb(barycentricCoordinatesCurrent[0], barycentricCoordinatesCurrent[1], {
                        this->vertexDomainCoordinates[triangleIndices[0]],
                        this->vertexDomainCoordinates[triangleIndices[1]],
                        this->vertexDomainCoordinates[triangleIndices[2]],
                        },
                        sheetColour);
                fb.sheetId = currentSheeId;
                fb.triangleId = currentTriangleId;
                this->faceFibers.push_back(fb);

                FaceFiberPoint fb2(barycentricCoordinatesNeighbour[0], barycentricCoordinatesNeighbour[1], {
                        this->vertexDomainCoordinates[triangle2Indices[0]],
                        this->vertexDomainCoordinates[triangle2Indices[1]],
                        this->vertexDomainCoordinates[triangle2Indices[2]],
                        },
                        sheetColour);
                fb2.sheetId = currentSheeId;
                fb2.triangleId = neighbourTriagleId;
                this->faceFibers.push_back(fb2);

                //printf("Adding fiber between %d -> %d\n", currentTriangleId, neighbourTriagleId);
            }
        }
    }

    //Timer::stop("Computed fiber in                      :");
}

void Data::computeTetExitPointsNew(const float u, const float v, const std::vector<float> color)
{
    //this->faceFibers.clear();

    //
    // Get the ID of the face we are intersecting
    //

    //Timer::start();

    // The query point (u, v)
    Point_2 query_point(u, v);


    // Locate the point in the arrangement
    CGAL::Object result = this->pl->locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    int currentFaceID = 0;

    if (CGAL::assign(face, result)) 
    {
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } else 
    {
        assert(false);
    }
    //Timer::stop("Computed active arrangement face       :");



    int i = -1;
    for (const auto &[triangle, triangleId] : this->preimageGraphs[currentFaceID].data)
    {
        i++;
        int j = -1;
        for (const auto&[triangle2, triangleId2] : this->preimageGraphs[currentFaceID].data)
        {
            j++;
            if (j <= i) { continue; }

            const set<int> triangleUnpacked = this->indexToTriangle[triangle];
            const set<int> triangle2Unpacked = this->indexToTriangle[triangle2];

            if (this->connectedTriangles.contains({triangleUnpacked, triangle2Unpacked}))
            {
                const int componentID = this->preimageGraphs[currentFaceID].find(triangleId);
                //const int pairToHIndex = this->vertexHtoIndex[{currentFaceID, componentID}];
                const int sheetID = this->reebSpace.findTriangle({currentFaceID, componentID});
                const int sheetColourID = this->sheetToColour[sheetID] % this->fiberColours.size();
                const vector<float> sheetColour = this->fiberColours[sheetColourID];

                //
                // Get the IDs and barycentri coordinates for the first point
                //
                vector<int> vertexIds;

                for(const int &vertexId : triangleUnpacked)
                {
                    vertexIds.push_back(vertexId);
                }

                CartesianPoint A(this->vertexCoordinatesF[vertexIds[0]], this->vertexCoordinatesG[vertexIds[0]]);
                CartesianPoint B(this->vertexCoordinatesF[vertexIds[1]], this->vertexCoordinatesG[vertexIds[1]]);
                CartesianPoint C(this->vertexCoordinatesF[vertexIds[2]], this->vertexCoordinatesG[vertexIds[2]]);

                // Define query point
                CartesianPoint P(u, v);

                std::array<double, 3> coordinates;
                CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, coordinates.begin());

                assert(coordinates[0] >= 0 && coordinates[1] >= 0 && coordinates[2] >= 0);

                FaceFiberPoint fb(coordinates[0], coordinates[1], {
                        this->vertexDomainCoordinates[vertexIds[0]],
                        this->vertexDomainCoordinates[vertexIds[1]],
                        this->vertexDomainCoordinates[vertexIds[2]],
                        },
                        sheetColour);
                this->faceFibers.push_back(fb);




                // Get the IDs and barycentri coordinates for the second point
                vector<int> vertexIds2;

                for(const int &vertexId : triangle2Unpacked)
                {
                    vertexIds2.push_back(vertexId);
                }


                // Define triangle vertices
                CartesianPoint A2(this->vertexCoordinatesF[vertexIds2[0]], this->vertexCoordinatesG[vertexIds2[0]]);
                CartesianPoint B2(this->vertexCoordinatesF[vertexIds2[1]], this->vertexCoordinatesG[vertexIds2[1]]);
                CartesianPoint C2(this->vertexCoordinatesF[vertexIds2[2]], this->vertexCoordinatesG[vertexIds2[2]]);

                std::array<double, 3> coordinates2;
                CGAL::Barycentric_coordinates::triangle_coordinates_2(A2, B2, C2, P, coordinates2.begin());
                assert(coordinates2[0] >= 0 && coordinates2[1] >= 0 && coordinates2[2] >= 0);

                FaceFiberPoint fb2(coordinates2[0], coordinates2[1], {
                        this->vertexDomainCoordinates[vertexIds2[0]],
                        this->vertexDomainCoordinates[vertexIds2[1]],
                        this->vertexDomainCoordinates[vertexIds2[2]],
                        },
                        sheetColour);

                this->faceFibers.push_back(fb2);
            }
        }
    }
}







void Data::computeTetExitPoints(const float u, const float v, const std::vector<float> color)
{
    this->faceFibers.clear();
    this->tetsWithFibers = vector<bool>(this->tetrahedra.size(), false);

    //
    // Get the ID of the face we are intersecting
    //

    // The query point (u, v)
    Point_2 query_point(u, v);

    // Locate the point in the arrangement
    CGAL::Object result = this->pl->locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    int currentFaceID = 0;

    if (CGAL::assign(face, result)) 
    {
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } else 
    {
        assert(false);
    }

    // For every tet, compute the two exit points
    for(size_t tetId = 0 ; tetId < this->tetrahedra.size(); tetId++)
    {
        const auto tet = this->tetrahedra[tetId];

        // For every triangle in every tet, get the fiber point in it
        for(int i = 0 ; i < 4 ; i++)
        {
            for(int j = i + 1 ; j < 4 ; j++)
            {
                for(int k = j + 1 ; k < 4 ; k++)
                {
                    float x1 = this->vertexCoordinatesF[tet[i]];
                    float y1 = this->vertexCoordinatesG[tet[i]];

                    float x2 = this->vertexCoordinatesF[tet[j]];
                    float y2 = this->vertexCoordinatesG[tet[j]];

                    float x3 = this->vertexCoordinatesF[tet[k]];
                    float y3 = this->vertexCoordinatesG[tet[k]];


                    const float xmin = std::min({x1, x2, x3});
                    const float xmax = std::max({x1, x2, x3});
                    const float ymin = std::min({y1, y2, y3});
                    const float ymax = std::max({y1, y2, y3});

                    // This triangle is not relevant, point is outside the bounding box
                    if (u < xmin || u > xmax || v < ymin || v > ymax) 
                    {
                        continue;
                    }


                    float det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

                    float alpha = ((y2 - y3) * (u - x3) + (x3 - x2) * (v - y3)) / det;
                    float beta = ((y3 - y1) * (u - x3) + (x1 - x3) * (v - y3)) / det;
                    float gamma = 1 - alpha - beta;

                    // Are we inside the triangle. We exclude the 0 and 1 because weird things happen there
                    if (alpha > 0 && beta > 0 && gamma > 0 && alpha < 1 && beta < 1 && gamma < 1)
                    {
                        //printf("In triangle %ld, %ld, %ld in tet %ld.\n", tet[i], tet[j], tet[k], tetId);
                        //printf("In triangle (%f, %f) | (%f, %f) | (%f, %f) comparing with point (%f, %f) and alpha = %f, betta = %f, gamma = %f.\n", x1, y1, x2, y2, x3, y3, x, y, alpha, betta, gamma);

                        //const std::set<int> triangle({tet[i], tet[j], tet[k]});

                        const int triangleVertexA = tet[i];
                        const int triangleVertexB = tet[j];
                        const int triangleVertexC = tet[k];

                        //printf("Current triangle (%d, %d, %d)\n", triangleVertexA, triangleVertexB, triangleVertexC);


                        //for (const auto [key, value] : this->faceDisjointSets[currentFaceID].data)
                        //{
                        //cout << "Triangle ";
                        //for (const auto v : key)
                        //{
                        //cout << v << " ";
                        //}

                        //cout << " with root " << this->faceDisjointSets[currentFaceID].find(value) << " ( " << this->faceDisjointSets[currentFaceID].findTriangle(key) << ") " << endl;

                        //}


                        //for (const auto &[key, value] : this->reebSpace.data)
                        //{
                        //printf("Face ID = %d, fiber component root = %d, SheetID = %d\n", key.first, key.second, this->reebSpace.findTriangle(key));

                        //}


                        // Which sheets does this fiber belong to?
                        // 1. Triangle -> Face ComponentID
                        const int triangleID = this->triangleToIndex[std::set<int>({triangleVertexA, triangleVertexB, triangleVertexC})];
                        const int componentID = this->preimageGraphs[currentFaceID].findTriangle(triangleID);
                        //printf("The face ID is %d and the component ID is = %d\n", currentFaceID, componentID);

                        // 2. Fac ComponentID -> Reeb Space Sheet
                        //const int pairToHIndex = this->vertexHtoIndex[{currentFaceID, componentID}];
                        //const int sheetID = this->reebSpace.findTriangle(pairToHIndex);
                        const int sheetID = this->reebSpace.findTriangle({currentFaceID, componentID});

                        //printf("The Sheet ID is = %d\n", sheetID);

                        const int sheetColourID = this->sheetToColour[sheetID];

                        // 3. Get the colou of the sheet
                        const vector<float> sheetColour = this->fiberColours[sheetColourID];


                        FaceFiberPoint fb(alpha, beta, {
                                this->vertexDomainCoordinates[tet[i]],
                                this->vertexDomainCoordinates[tet[j]],
                                this->vertexDomainCoordinates[tet[k]],
                                },
                                sheetColour);

                        this->faceFibers.push_back(fb);
                        this->tetsWithFibers[tetId] = true;
                    }
                }
            }
        }
    }
}

void
Data::computeMinMaxRangeDomainCoordinates()
{
    // Compute the min/max domain coordinates
    this->minX = this->vertexDomainCoordinates[0][0];
    this->maxX = this->vertexDomainCoordinates[0][0];

    this->minY = this->vertexDomainCoordinates[0][1];
    this->maxY = this->vertexDomainCoordinates[0][1];

    this->minZ = this->vertexDomainCoordinates[0][2];
    this->maxZ = this->vertexDomainCoordinates[0][2];

    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->minX = std::min(this->minX, this->vertexDomainCoordinates[i][0]);
        this->maxX = std::max(this->maxX, this->vertexDomainCoordinates[i][0]);

        this->minY = std::min(this->minY, this->vertexDomainCoordinates[i][1]);
        this->maxY = std::max(this->maxY, this->vertexDomainCoordinates[i][1]);

        this->minZ = std::min(this->minZ, this->vertexDomainCoordinates[i][2]);
        this->maxZ = std::max(this->maxZ, this->vertexDomainCoordinates[i][2]);
    }


    // Compute the min/max range coordinates
    this->minF = this->vertexCoordinatesF[0];
    this->maxF = this->vertexCoordinatesF[0];

    this->minG = this->vertexCoordinatesG[0];
    this->maxG = this->vertexCoordinatesG[0];

    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->minF = std::min(this->minF, this->vertexCoordinatesF[i]);
        this->maxF = std::max(this->maxF, this->vertexCoordinatesF[i]);

        this->minG = std::min(this->minG, this->vertexCoordinatesG[i]);
        this->maxG = std::max(this->maxG, this->vertexCoordinatesG[i]);
    }

    // Add some padding to the range coordinates for better visibility
    //this->minF -= .2;
    //this->maxF += .2;
    //this->minG -= .2;
    //this->maxG += .2;

    this->minF -= 0.1 * (this->maxF - this->minF);
    this->maxF += 0.1 * (this->maxF - this->minF);
    this->minG -= 0.1 * (this->maxG - this->minG);
    this->maxG += 0.1 * (this->maxG - this->minG);
}

void Data::sortVertices()
{
    vector<pair<CartesianPoint, int>> pointWithIndices(this->vertexCoordinatesG.size());

    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        pointWithIndices[i] = {CartesianPoint(this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]), i};
    }

    //cout << "Before sorting" << endl;
    //this->printMesh();

    std::sort(pointWithIndices.begin(), pointWithIndices.end(),
            [](const pair<CartesianPoint, int> &a, const pair<CartesianPoint, int>& b) {
            return CGAL::compare_xy(a.first, b.first) == CGAL::SMALLER;
            });

    //cout << "Printing points... " << endl;
    //for (int i = 0 ; i < pointWithIndices.size() ; i++)
    //{
        //printf("%d - (%f, %f) - %d\n", i, pointWithIndices[i].first.x(), pointWithIndices[i].first.y(), pointWithIndices[i].second);
    //}


    // Set up the inverse index search
    vector<int> meshIDtoSortIndex(pointWithIndices.size());
    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        meshIDtoSortIndex[pointWithIndices[i].second] = i;
    }

    //cout << "Swapping..." << endl;
    //for (int i = 0 ; i < pointWithIndices.size() ; i++)
    //{
        //printf("%d -> %d\n", i, meshIDtoSortIndex[i]);

    //}

    //
    // Now we can swap things around
    //

    // Set up copies of the originals for the swap, otherwise editin in place causes errors
    std::vector<std::vector<float>> vertexDomainCoordinatesOriginal = this->vertexDomainCoordinates;
    std::vector<float> vertexCoordinatesFOriginal = this->vertexCoordinatesF;
    std::vector<float> vertexCoordinatesGOriginal = this->vertexCoordinatesG;

    // Swap tet indices
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        for (int j = 0 ; j < this->tetrahedra[i].size() ; j++)
        {
            this->tetrahedra[i][j] = meshIDtoSortIndex[this->tetrahedra[i][j]];
        }
    }

    // Swap domain positions
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->vertexDomainCoordinates[i] = vertexDomainCoordinatesOriginal[pointWithIndices[i].second];
    }
    
    // Swap range positions
    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->vertexCoordinatesF[i] = vertexCoordinatesFOriginal[pointWithIndices[i].second];
        this->vertexCoordinatesG[i] = vertexCoordinatesGOriginal[pointWithIndices[i].second];
    }
    //cout << "After sorting..." << endl;
    //this->printMesh();

}

void Data::printMesh()
{
    // Print vertex domain coordinates
    cout << "Vertex domain coordinates: " << endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f, %f)\n", i, this->vertexDomainCoordinates[i][0], this->vertexDomainCoordinates[i][1], this->vertexDomainCoordinates[i][2]);
    }
    // Print vertex range coordinates
    cout << "Vertex range coordinates: " << endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f)\n", i, this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]);
    }

    // Print tetrahedra
    cout << "Tetrahedra: " << endl;
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        printf("%d - (%ld, %ld, %ld, %ld)\n", i, this->tetrahedra[i][0], this->tetrahedra[i][1], this->tetrahedra[i][2], this->tetrahedra[i][3]);
    }
}
void Data::readDataVTU(const string &filename, const float &perturbationEpsilon)
{
    // Read the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();


    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    int numVertices = mesh->GetPoints()->GetNumberOfPoints(); 
    int numTets = mesh->GetNumberOfCells();

    // Initialize all the data arrays
    this->vertexCoordinatesF = std::vector<float>(numVertices, 0);
    this->vertexCoordinatesG = std::vector<float>(numVertices, 0);
    this->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    this->vertexDomainCoordinates = std::vector<std::vector<float>>(numVertices, {0, 0, 0});



    // Print vertices
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    //std::cout << "Vertices:\n";
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        //std::cout << "Vertex " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";

        this->vertexDomainCoordinates[i][0] = p[0];
        this->vertexDomainCoordinates[i][1] = p[1];
        this->vertexDomainCoordinates[i][2] = p[2];
    }

    // Print tetrahedra
    //std::cout << "\nTetrahedra:\n";
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell* cell = mesh->GetCell(i);
        if (cell->GetNumberOfPoints() == 4) { // Tetrahedron check
            //std::cout << "Tetrahedron " << i << ": ";
            for (vtkIdType j = 0; j < 4; j++) {
                //std::cout << cell->GetPointId(j) << " ";
                this->tetrahedra[i][j] = cell->GetPointId(j);
            }
            //std::cout << "\n";
        }
    }

    // Print vertex data arrays
    //std::cout << "\nVertex Data Arrays:\n";
    vtkPointData* pointData = mesh->GetPointData();

    assert(pointData->GetNumberOfArrays() >= 2);

    vtkDataArray* fDataArray = pointData->GetArray(1);
    vtkDataArray* gDataArray = pointData->GetArray(0);

    assert(fDataArray->GetNumberOfTuples() == numVertices);
    assert(gDataArray->GetNumberOfTuples() == numVertices);

    for (vtkIdType i = 0; i < fDataArray->GetNumberOfTuples(); i++) 
    {
        this->vertexCoordinatesF[i] = fDataArray->GetTuple1(i);
    }

    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        this->vertexCoordinatesG[i] = gDataArray->GetTuple1(i);
    }

    // Add numerical perturbation, taken from Tierny 2027 Jacobi Fiber Surfaces
    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        //const float e = randomPerturbation(1e-8);
        //const float iFloat = static_cast<float>(i);

        // Tierny 2017 - Jacobi Fiber Surfaces Version - not that useful
        //this->vertexCoordinatesF[i] += iFloat * e;
        //this->vertexCoordinatesG[i] += iFloat * e * iFloat * e;

        this->vertexCoordinatesF[i] += randomPerturbation(perturbationEpsilon);
        this->vertexCoordinatesG[i] += randomPerturbation(perturbationEpsilon);
    }

    // Now we can sort
    this->sortVertices();

    // Compute the ranges for a bounding box in the range.
    this->computeMinMaxRangeDomainCoordinates();
}

void
Data::readData(const string &filename, const float &perturbationEpsilon)
{
    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
    if (false == dataFile.is_open()) { throw "Could not open data file."; }

    // Read in data in a string and skip the comments
    string rawStringData;
    string myline;
    while (dataFile) {
        std::getline (dataFile, myline);
        if (myline[0] == '#')
        {
            //std::cout << myline << '\n';
        }
        else
        {
            rawStringData += " " + myline;
        }
    }

    // Set up the inputstream from the string
    std::istringstream dataStream(rawStringData);

    // Read in the number of vertices and tets
    int numVertices, numTets;
    dataStream >> numVertices >> numTets;

    // Initialize all the data arrays
    this->vertexCoordinatesF = std::vector<float>(numVertices, 0);
    this->vertexCoordinatesG = std::vector<float>(numVertices, 0);
    this->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    this->vertexDomainCoordinates = std::vector<std::vector<float>>(numVertices, {0, 0, 0});

    // Read in the domain coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexDomainCoordinates[i][0];
        dataStream >> this->vertexDomainCoordinates[i][1];
        dataStream >> this->vertexDomainCoordinates[i][2];
    }

    // Read in the range coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexCoordinatesF[i];
        dataStream >> this->vertexCoordinatesG[i];
    }
    
    // Read in the tetrahedron configuration
    for  (int i = 0 ; i < numTets ; i++)
    {
        dataStream >> this->tetrahedra[i][0];
        dataStream >> this->tetrahedra[i][1];
        dataStream >> this->tetrahedra[i][2];
        dataStream >> this->tetrahedra[i][3];
    }

    this->computeMinMaxRangeDomainCoordinates();

    this->sortVertices();
}
