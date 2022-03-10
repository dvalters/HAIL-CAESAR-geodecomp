#include <libgeodecomp/misc/apitraits.h>
#include <libgeodecomp/geometry/coord.h>

#include <typemaps.h>
#include <apitraits.h>

#ifndef CELL_H
#define CELL_H


class Cell
{
  friend class Typemaps;
  friend int main(int argc, char **argv);
  
  // refactor - Should replace these defines with type alias declarations (= C++11 template typedef)
  // refactor - check that grid orientation makes sense (write test)
#define thisCell_old neighborhood[LibGeoDecomp::Coord<2>( 0, 0)]
#define west_old neighborhood[LibGeoDecomp::Coord<2>(-1, 0)]
#define east_old neighborhood[LibGeoDecomp::Coord<2>( 1, 0)]
#define north_old neighborhood[LibGeoDecomp::Coord<2>( 0, -1)]
#define south_old neighborhood[LibGeoDecomp::Coord<2>( 0,  1)]
  
  class API :
    public LibGeoDecomp::APITraits::HasStencil<LibGeoDecomp::Stencils::VonNeumann<2,1> >,
    public LibGeoDecomp::APITraits::HasCubeTopology<2>,
    public LibGeoDecomp::APITraits::HasCustomMPIDataType<Cell>,
    public LibGeoDecomp::APITraits::HasOpaqueMPIDataType<Cell>
  {};
  

public:
  enum CellType : int {INTERNAL=0,	                                     \
		       EDGE_WEST=1, EDGE_NORTH=2, EDGE_EAST=3, EDGE_SOUTH=4, \
		       CORNER_NW=5, CORNER_NE=6, CORNER_SE=7, CORNER_SW=8,   \
		       NODATA=9};

  // Grid variables are instances variables of Cell class (each grid cell has its own value)
  CellType celltype;
  double elevation;
  double water_depth;
  double qx;
  double qy;

  // Overall simulation parameters are 
  static double water_depth_erosion_threshold;
  static double no_data_value;
  static double edgeslope;
  static double hflow_threshold;
  static double mannings;
  static double froude_limit;
  static double time_factor;
  static double courant_number;
  static double maxdepth;
  constexpr static double gravity = 9.81;
  
  
  
  explicit Cell(CellType celltype_in, double elevation_in, double water_depth_in, \
	      double qx_in, double qy_in);
  
  void set_global_timefactor();
  double set_local_timefactor();
  void froude_check(double &q, double hflow);
  void update_q(const double &q_old, double &q_new, double hflow, double tempslope, double local_time_factor);
  
  template<typename COORD_MAP> void initialise_grid_value_updates(const COORD_MAP& neighborhood);
  template<typename COORD_MAP> void update(const COORD_MAP& neighborhood, unsigned nanoStep);
  template<typename COORD_MAP> void catchment_waterinputs(const COORD_MAP& neighborhood);
  template<typename COORD_MAP> void catchment_water_input_and_hydrology(const COORD_MAP& neighborhood, double local_time_factor);
  template<typename COORD_MAP> void water_flux_out(const COORD_MAP& neighborhood);
  template<typename COORD_MAP> void flow_route_x(const COORD_MAP& neighborhood);
  template<typename COORD_MAP> void flow_route_y(const COORD_MAP& neighborhood);
  template<typename COORD_MAP> void depth_update(const COORD_MAP& neighborhood);
  template<typename COORD_MAP> void update_qx(const COORD_MAP& neighborhood, double hflow, double tempslope, double local_time_factor);
  template<typename COORD_MAP> void update_qy(const COORD_MAP& neighborhood, double hflow, double tempslope, double local_time_factor);
  template<typename COORD_MAP> void update_water_depth(const COORD_MAP& neighborhood, double east_qx, double south_qy, double local_time_factor);
  template<typename COORD_MAP> void discharge_check(const COORD_MAP& neighborhood, double &q, double neighbour_water_depth, double Delta);
  
};








#endif


