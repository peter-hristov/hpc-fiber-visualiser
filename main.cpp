
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPolyData.h>
#include <vtkVertex.h>
#include <vtkXMLPolyDataWriter.h>

#include <vector>
#include <set>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Segment_2                                       Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;

void write_vtk_vtp(const std::vector<Segment_2>& segments,
                   const std::vector<Point_2>& intersections,
                   const std::string& filename) {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    auto vertices = vtkSmartPointer<vtkCellArray>::New();

    std::map<Point_2, vtkIdType> point_id_map;
    vtkIdType next_id = 0;

    // Add endpoints of segments
    for (const auto& seg : segments) {
        for (const auto& pt : { seg.source(), seg.target() }) {
            if (point_id_map.count(pt) == 0) {
                vtkIdType id = points->InsertNextPoint(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), 0.0);
                point_id_map[pt] = id;
            }
        }
        vtkIdType id1 = point_id_map[seg.source()];
        vtkIdType id2 = point_id_map[seg.target()];

        auto line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, id1);
        line->GetPointIds()->SetId(1, id2);
        lines->InsertNextCell(line);
    }

    // Add intersection points
    for (const auto& pt : intersections) {
        if (point_id_map.count(pt) == 0) {
            vtkIdType id = points->InsertNextPoint(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), 0.0);
            point_id_map[pt] = id;
        }
        vtkIdType id = point_id_map[pt];

        auto vertex = vtkSmartPointer<vtkVertex>::New();
        vertex->GetPointIds()->SetId(0, id);
        vertices->InsertNextCell(vertex);
    }

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    polyData->SetVerts(vertices);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    writer->Write();

    std::cout << "Wrote " << filename << "\n";
}

int main() {
    Segment_2 seg1(Point_2(0, 0), Point_2(1, 1));
    Segment_2 seg2(Point_2(0, 1), Point_2(1, 0));

    Arrangement_2 arr;
    insert(arr, seg1);
    insert(arr, seg2);

    std::vector<Segment_2> segments = { seg1, seg2 };
    std::set<Point_2> intersections;

    // CGAL computes all vertices and we collect ones that are intersection points
    for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        if (!vit->is_isolated() && vit->degree() > 2) {
            intersections.insert(vit->point());
        }
    }

    write_vtk_vtp(segments, {intersections.begin(), intersections.end()}, "curves.vtp");
    return 0;
}
