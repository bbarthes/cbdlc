
#include "default_defines.h"
#include "kernel_list.h"

static char *kernels_name[KERNEL_TAB_SIZE] = {
  "gravity", // Gravity
  "eating", // Eating
  "growing", // Growing
  "null_kernel", // reset_int (NOT USED)
  "null_kernel", // box count_all (NOT USED)
  "null_kernel", // box_count (NOT_USED)
  "null_kernel", // scan (NOT USED)
  "null_kernel", // scan_down (NOT USED)
  "null_kernel", // copy (NOT USED)
  "null_kernel", // box_sort_all (NOT_USED)
  "null_kernel", // box_sort (NOT_USED)
  "null_kernel", // box_force (NOT_USED)
  "lennard_jones", // force
  "border_collision",  // bounce
  "update_position", // update_position
#ifdef _SPHERE_MODE_
  "update_vertices", // update_vertices
#else
  "update_vertice", // update_vertices
#endif

  "zero_speed", // zero_speed
  "atom_collision", // collision
  "null_kernel", // NULL 
};



char *kernel_name(unsigned kernel_num)
{
  return kernels_name[kernel_num];
}
