#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;

int main() {
    // 1. Create a CGAL segment
    Segment_2 seg(Point_2(0, 0), Point_2(1, 1));

    // 2. Convert to VTK
    auto points = vtkSmartPointer<vtkPoints>::New();
    vtkIdType id1 = points->InsertNextPoint(CGAL::to_double(seg.source().x()), CGAL::to_double(seg.source().y()), 0.0);
    vtkIdType id2 = points->InsertNextPoint(CGAL::to_double(seg.target().x()), CGAL::to_double(seg.target().y()), 0.0);

    auto line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, id1);
    line->GetPointIds()->SetId(1, id2);

    auto lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell(line);

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("segment.vtp");
    writer->SetInputData(polyData);
    writer->Write();

    return 0;
}
