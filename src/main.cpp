#include <set>
#include <filesystem>
#include <string>
#include <unordered_map>

#include "./Timer.h"
#include "./ReebSpace.h"
#include "./CGALTypedefs.h"
#include "./DisjointSet.h"
#include "./utility/CLI11.hpp"

#include <CGAL/Union_find.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellData.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>



using namespace std;
namespace fs = std::filesystem;


void writePolygons(Data *data, string filename)
{

    // Create the multiblock dataset
    auto multiBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    int counter = 0;

    for (const auto &[sheetId, polygon] : data->sheetPolygon)
    {
        const vector<float> sheetColour = data->fiberColours[data->sheetToColour[sheetId] % data->fiberColours.size()];


        //if (true)
        if (data->incompleteSheets.contains(sheetId))
        {
            cout << "Writing sheet " << sheetId << endl;

            auto points = vtkSmartPointer<vtkPoints>::New();
            auto polygonsCells = vtkSmartPointer<vtkCellArray>::New();

            auto cellColors = vtkSmartPointer<vtkFloatArray>::New();
            cellColors->SetNumberOfComponents(3);  // RGB
            cellColors->SetName("Colors");

            auto cellSheetIds = vtkSmartPointer<vtkIntArray>::New();
            cellSheetIds->SetNumberOfComponents(1);
            cellSheetIds->SetName("SheetId");

            // Loop through all faces to see which ones are in the sheet
            for (auto f = data->arr.faces_begin(); f != data->arr.faces_end(); ++f) 
            {
                const int currentFaceID = data->arrangementFacesIdices[f];

                // This adds one polygon
                for (const auto &[triangleId, fiberComponentId] : data->fiberSeeds[currentFaceID])
                {

                    const int componentSheetId = data->reebSpace.findTriangle({currentFaceID, fiberComponentId});

                    // Now we can add the polygon
                    if (componentSheetId == sheetId)
                    {
                        // Get the starting index for these new points
                        vtkIdType startIndex = points->GetNumberOfPoints();

                        typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
                        typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
                        do {
                            typename Arrangement_2::Halfedge_const_handle e = curr;

                            // Get point from CGAL (and convert to double )
                            const float u = CGAL::to_double(e->source()->point().x());
                            const float v = CGAL::to_double(e->source()->point().y());

                            points->InsertNextPoint(u, v, 0.0);
                        } while (++curr != circ);


                        // Write the face as a polygon
                        auto poly = vtkSmartPointer<vtkPolygon>::New();
                        poly->GetPointIds()->SetNumberOfIds(polygon.size());
                        for (vtkIdType j = 0; j < static_cast<vtkIdType>(polygon.size()); ++j)
                        {
                            poly->GetPointIds()->SetId(j, startIndex + j);
                        }

                        polygonsCells->InsertNextCell(poly);

                        // Assign per-cell data
                        cellColors->InsertNextTypedTuple(sheetColour.data());
                        cellSheetIds->InsertNextValue(sheetId);
                    }
                }
            }

            // Create one final polydata object
            auto polyData = vtkSmartPointer<vtkPolyData>::New();
            polyData->SetPoints(points);
            polyData->SetPolys(polygonsCells);
            polyData->GetCellData()->AddArray(cellColors);
            polyData->GetCellData()->AddArray(cellSheetIds);
            polyData->GetCellData()->SetScalars(cellColors);

            multiBlock->SetBlock(counter++, polyData);
        }
        else
        {
            //continue;
            // Create points and cell array
            auto points = vtkSmartPointer<vtkPoints>::New();
            auto polygonsCells = vtkSmartPointer<vtkCellArray>::New();

            auto cellColors = vtkSmartPointer<vtkFloatArray>::New();
            cellColors->SetNumberOfComponents(3);  // RGB
            cellColors->SetName("Colors");

            auto cellSheetIds = vtkSmartPointer<vtkIntArray>::New();
            cellSheetIds->SetNumberOfComponents(1);  // RGB
            cellSheetIds->SetName("SheetId");

            // Add points to the points object
            vtkIdType startIndex = points->GetNumberOfPoints();
            for (const CartesianPoint &point : polygon) 
            {
                points->InsertNextPoint(point.x(), point.y(), 0.0);
            }

            // Create the polygon and add it to the cell array
            auto poly = vtkSmartPointer<vtkPolygon>::New();
            poly->GetPointIds()->SetNumberOfIds(polygon.size());
            for (vtkIdType j = 0; j < static_cast<vtkIdType>(polygon.size()); ++j)
            {
                poly->GetPointIds()->SetId(j, startIndex + j);
            }
            polygonsCells->InsertNextCell(poly);

            // Add the color for this polygon
            cellColors->InsertNextTypedTuple(sheetColour.data());
            cellSheetIds->InsertNextValue(sheetId);


            // Create the polydata
            auto polyData = vtkSmartPointer<vtkPolyData>::New();
            polyData->SetPoints(points);
            polyData->SetPolys(polygonsCells);
            polyData->GetCellData()->AddArray(cellSheetIds);  // Add the integer array
            polyData->GetCellData()->AddArray(cellColors);
            polyData->GetCellData()->SetScalars(cellColors);


            // Add the polydata as a block in the multiblock dataset
            multiBlock->SetBlock(counter++, polyData);
        }




        // Write to VTP file
        //auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        //writer->SetFileName(filename.c_str());
        //writer->SetInputData(polyData);
        //writer->Write();

    }


    // Write the multiblock dataset to a VTK file
    auto writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(multiBlock);
    writer->Write();


}


int main(int argc, char* argv[])
{

    // Parse the command line arguments
    CLI::App cliApp("Fiber Visualiser");

    string filename;
    cliApp.add_option("--file, -f", filename, "Input data filename. Has to be either .txt of .vti.")->required();

    bool performanceRun = false;
    cliApp.add_flag("--performanceRun, -p", performanceRun, "Only compute the Reeb space, no graphics..");

    bool discardFiberSeeds = false;
    cliApp.add_flag("--discardPreimageGraphs, -d", discardFiberSeeds, "Discard the seeds for generating fibers based on sheets, discard to save a bit of memory (not too much).");

    float perturbationEpsilon = 1e-2;
    cliApp.add_option("--epsilon, -e", perturbationEpsilon, "Strength of the numerial perturbation in the range [-e, e].");

    string outputSheetPolygonsFilename;
    cliApp.add_option("--outputSheetPolygons, -o", outputSheetPolygonsFilename, "Filename where to output the coordinates of the polygons that represent each sheet.");

    int fiberSampling = 1;
    cliApp.add_option("--fiberSampling, -s", fiberSampling, "When saving fibers per component, how many do we save. Default is to save the centroid, otherwise sample along the boundary.");

    int sheetOutputCount = 10;
    cliApp.add_option("--sheetOutputCount", sheetOutputCount, "How many sheets to sample for automatic feature extraction.");

    string outputSheetFibersFolder;
    cliApp.add_option("--outputSheetFibersFolder", outputSheetFibersFolder, "Folder in which to ouput fiber for each sheet.");

    //string outputFibersFilename = "./fibers.vtp";
    //cliApp.add_option("--outputFibers", outputSheetPolygonsFilename, "Filename where to save the visible fiber components. Must be .vtp");

    CLI11_PARSE(cliApp, argc, argv);

    // For convenience
    if (performanceRun == true)
    {
        discardFiberSeeds = true;
    }

    fs::path filePath(filename);
    
    if (!fs::exists(filePath)) 
    {
        std::cerr << "Error: File does not exist: " << filename << std::endl;
        return 0;
    }

    // Set up the data class
    Data* data = new Data();

    std::string extension = filePath.extension().string();
    if (extension == ".vtu") 
    {
        data->readDataVTU(filename, perturbationEpsilon);
    } 
    else if (extension == ".txt") 
    {
        data->readData(filename, perturbationEpsilon);
    } 
    else 
    {
        std::cerr << "Error: Unsupported file type: " << extension << std::endl;
    }

    // Compute the 2D arrangement
    ReebSpace::computeArrangement(data);

    //std::cout << "Press Enter to continue...";
    //std::cin.get();  // waits for Enter key

    Timer::start();
    ReebSpace::computeUpperLowerLink(data);
    Timer::stop("Computed upper and lower link          :");

    Timer::start();
    ReebSpace::computeTriangleAdjacency(data);
    Timer::stop("Computed triangle adjacency            :");


    //cout << "Triangles to Index " << endl;
    //for (const auto &[triangle, triangleId] : data->triangleToIndex)
    //{
        //cout << "Triangle = ";
        //for (const auto v : triangle)
        //{
            //cout << v << " ";
        //}
        //cout << "  ID = " << triangleId << endl;
    //}

    //cout << "\n\nIndex to triangle" << endl;
    //for (int i = 0 ; i < data->indexToTriangle.size() ; i++)
    //{
        //cout << "  ID = " << i << " ";

        //cout << "Triangle = ";
        //for (const auto v : data->indexToTriangle[i])
        //{
            //cout << v << " ";
        //}

        //cout << endl;
    //}


    Timer::start();
    ReebSpace::testTraverseArrangement(data);
    Timer::stop("Computed empty traversal               :");

    Timer::start();
    ReebSpace::computePreimageGraphs(data, discardFiberSeeds);
    Timer::stop("Computed {G_F} and H                   :");

    //Timer::start();
    //ReebSpace::computeCorrespondenceGraph(data);
    //Timer::stop("Computed H                             :");

    std::cout << "Postprocessing..." << std::endl;
    //Timer::start();
    ReebSpace::computeReebSpacePostprocess(data);
    //Timer::stop("Computed RS(f) Postprocess             :");

    //std::cout << "Press Enter to continue...";
    //std::cin.get();  // waits for Enter key

    // Save all the polygons
    if (false == outputSheetPolygonsFilename.empty())
    {
        fs::path filePathOutput(outputSheetPolygonsFilename);

        // Write to the file
        std::ofstream outFile(filePathOutput);
        if (!outFile) 
        {
            std::cerr << "Error: Could not open file for writing: " << filePathOutput << std::endl;
            return 1;
        }

        outFile << data->sheetPolygon.size() << std::endl;

        for (const auto &[sheetId, polygon] : data->sheetPolygon)
        {
            outFile << "SheetId = " << sheetId << std::endl;


            // Compute the controid so that we can pull all verties towards it
            CartesianPoint centroid = CGAL::centroid(polygon.vertices_begin(), polygon.vertices_end());

            // To make sure we don't write a comma at the end of the array
            int pointsWritten = 0;

            outFile << "[";
            for (const CartesianPoint &point : polygon) 
            {
                // Get point from CGAL (and convert to double )
                float u = point.x();
                float v = point.y();

                // Interpolate closer to the centroid
                float alpha = 0.5;
                u = (1 - alpha) * u + alpha * centroid.x();
                v = (1 - alpha) * v + alpha * centroid.x();

                outFile << u << ", " << v << ", " << 0;
                if (pointsWritten < polygon.size() - 1)
                {
                    outFile << ", ";

                }

                pointsWritten++;
            }
            outFile << "]" << std::endl;

        }

        outFile.close();
    }

    if (performanceRun == true)
    {
        return 0;
    }


    Timer::start();
    data->pl = std::make_unique<Point_location>(data->arr);
    Timer::stop("Arrangement search structure           :");

    if (false == outputSheetFibersFolder.empty())
    {
        data->generatefFaceFibersForSheets(sheetOutputCount, fiberSampling, outputSheetFibersFolder);
    }


    writePolygons(data, "./polygons.vtm");

    delete data;

    // return to caller
    return 0;
} // main()
