#include <CGAL/trace.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

int main(int argc, char* argv[])	{
	// Poisson options
	FT sm_angle = 20.0; // Min triangle angle in degrees.
	FT sm_radius; // Max triangle size w.r.t. point set average spacing.
	FT sm_distance = 0.375; // Surface Approximation error w.r.t. point set average spacing.
	int to_average_spacing = 6;
	std::string input_path = "input.xyz";
	
	std::istringstream ss;
	
	ss.clear();
	ss.str(argv[1]);
	ss >> sm_angle;
	
	ss.clear();
	ss.str(argv[2]);
	ss >> sm_radius;
	
	ss.clear();
	ss.str(argv[3]);
	ss >> sm_distance;
	
	ss.clear();
	ss.str(argv[4]);
	ss >> to_average_spacing;
	
	input_path = argv[5];
	
	std::cout<<"sm_angle: "<<sm_angle<<std::endl;
	std::cout<<"sm_radius: "<<sm_radius<<std::endl;
	std::cout<<"sm_distance: "<<sm_distance<<std::endl;
	std::cout<<"to_average_spacing: "<<to_average_spacing<<std::endl;
	std::cout<<"input_path: "<<input_path<<std::endl;
	
	


	// Reads the point set file in points[].
	// Note: read_xyz_points_and_normals() requires an iterator over points
	// + property maps to access each point's position and normal.
	// The position property map can be omitted here as we use iterators over Point_3 elements.
	PointList points;
	std::ifstream stream(input_path.c_str());
	if (!stream ||
		!CGAL::read_xyz_points_and_normals(
							stream,
							std::back_inserter(points),
							CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())))	{
		std::cerr << "Error: cannot read file " + input_path << std::endl;
		return EXIT_FAILURE;
	}

	// Creates implicit function from the read points using the default solver.
	// Note: this method requires an iterator over points
	// + property maps to access each point's position and normal.
	// The position property map can be omitted here as we use iterators over Point_3 elements.
	Poisson_reconstruction_function function(
							points.begin(), points.end(),
							CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()));
	std::cout << "Info: Poisson_reconstruction_function done" << std::endl;

	// Computes the Poisson indicator function f()
	// at each vertex of the triangulation.
	if ( ! function.compute_implicit_function() )	{
		std::cout << "Error: ! function.compute_implicit_function()" << std::endl;
		return EXIT_FAILURE;
	}

	// Computes average spacing
	FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), to_average_spacing /* knn = 1 ring */);
	std::cout << "Info: CGAL::compute_average_spacing() done" << std::endl;
	std::cout<<"average spacing: "<<average_spacing<<std::endl;

	// Gets one point inside the implicit surface
	// and computes implicit function bounding sphere radius.
	Point inner_point = function.get_inner_point();
	Sphere bsphere = function.bounding_sphere();
	FT radius = std::sqrt(bsphere.squared_radius());

	// Defines the implicit surface: requires defining a
	// conservative bounding sphere centered at inner point.
	FT sm_sphere_radius = 5.0 * radius;
	FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
	Surface_3 surface(function,
					Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
					sm_dichotomy_error/sm_sphere_radius);

	std::cout << "Info: Surface_3 surface() done" << std::endl;

	// Defines surface mesh generation criteria
	CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
														sm_radius * average_spacing,  // Max triangle size
														sm_distance * average_spacing); // Approximation error

	std::cout << "Info: CGAL::Surface_mesh_default_criteria_3<STr> criteria() done" << std::endl;

	// Generates surface mesh with manifold option
	STr tr; // 3D Delaunay triangulation for surface mesh generation
	C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
	CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
							surface,                              // implicit surface
							criteria,                             // meshing criteria
							CGAL::Manifold_with_boundary_tag());  // require manifold mesh
	std::cout << "Info: CGAL::make_surface_mesh() done" << std::endl;

	if(tr.number_of_vertices() == 0)	{
		//return EXIT_FAILURE;
		std::cout << "Error: tr.number_of_vertices() == 0" << std::endl;
		return false;
	}

	std::ofstream out("output.off");
	CGAL::output_surface_facets_to_off(out, c2t3);

    //return EXIT_SUCCESS;
    std::cout << "Info: CGAL triangulation done" << std::endl;
    return true;
}
