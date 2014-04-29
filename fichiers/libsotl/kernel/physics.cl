#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable

#define ALRND(a,n) (((size_t)(n)+(a)-1)&(~(size_t)((a)-1)))

#define ALIGN_MASK (~(15UL))

#if USE_DOUBLE == 0
#define coord_t float3
#define calc_t float
#else
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define coord_t double3 
#define calc_t double
#endif

static inline coord_t load3coord (__global calc_t * pos, int offset)
{
    coord_t f;

    f.x = *pos; pos += offset;
    f.y = *pos; pos += offset;
    f.z = *pos;

    return f;
}

static inline void store3coord (__global calc_t * pos, coord_t f, int offset)
{
    *pos = f.x; pos += offset;
    *pos = f.y; pos += offset;
    *pos = f.z;
}

static inline void inc3coord (__global calc_t * pos, coord_t f, int offset)
{
    *pos += f.x; pos += offset;
    *pos += f.y; pos += offset;
    *pos += f.z;
}

static inline calc_t squared_dist(coord_t a, coord_t b)
{
  coord_t f = a - b;
  return dot(f, f);
}

static inline calc_t lj_squared_v(calc_t r2)
{
    calc_t rr2 = 1.0 / r2;
    calc_t r6;

    r6 = LENNARD_SIGMA * LENNARD_SIGMA * rr2;
    r6 = r6 * r6 * r6;

    calc_t de = 24 * LENNARD_EPSILON * rr2 * (2.0f * r6 * r6 - r6);

    return de;
}

__kernel
void update_position(__global calc_t *pos, __global calc_t *spd,
                     __constant calc_t *min, __constant calc_t *max,
                     unsigned offset)
{
    const unsigned gid = get_global_id(0);
    coord_t my_spd = load3coord(spd + gid, offset);

    inc3coord (pos + gid, my_spd, offset);
}

// This kernel is executed with offset x 3 threads
//
__kernel
void zero_speed(__global calc_t *speed)
{
  int index = get_global_id(0);

  speed[index] = 0.0;
}

// This kernel updates the Vertex Buffer Object according to positions of atoms
// It is executed by natoms * 3 threads
//
__kernel
void update_vertice(__global float* vbo, __global calc_t * pos, unsigned offset)
    // global = 3 * nb_vertices, local = 1
{
    unsigned index = get_global_id (0);
    unsigned axe = index % 3;
    unsigned no_atom = index / 3;

    vbo[index] = pos[no_atom + axe * offset];
}


__kernel
void border_collision (__global calc_t * pos, __global calc_t * speed,
		  __constant calc_t * min, __constant calc_t * max,
		  calc_t radius, unsigned natoms, unsigned offset)
{
  // TODO
}


static inline void freeze_atom(__global calc_t * speed, int offset)
{
    *speed = 0.0; speed += offset;
    *speed = 0.0; speed += offset;
    *speed = 0.0;
}

static void check_collision (__global calc_t * pos, __global calc_t * speed,
			     calc_t radius, unsigned index,
			     unsigned natoms, unsigned offset)
{
    coord_t mypos;

    mypos = load3coord (pos + index, offset);

    for (unsigned n = index + 1; n < natoms; n++)
    {
        coord_t other = load3coord (pos + n, offset);
        if (distance (mypos, other) <= 2 * radius)
        {
	  // Frozen Bubble mode
	  freeze_atom (speed + index, offset);
	  freeze_atom (speed + n, offset);
        }
    }
}

// This kernel is executed with one thread per atom
//
__kernel
void atom_collision (__global calc_t * pos, __global calc_t * speed,
		     calc_t radius, unsigned natoms, unsigned offset)
{
    unsigned index = get_global_id (0);

    if(index < natoms)
      check_collision (pos, speed, radius, index, natoms, offset);
}

// This kernel is executed with one thread per atom
//
__kernel
void gravity (__global calc_t * pos, __global calc_t * speed, calc_t g,
	      unsigned natoms, unsigned offset)
{
  // TODO
}

// This kernel is executed with one thread per atom
//
__kernel
void lennard_jones (__global calc_t * pos,
		       __global calc_t * speed,
		       unsigned natoms, unsigned offset)
{
  // TODO
}


__kernel void null_kernel (void)
{}

// The following kernels only work under SPHERE_MODE

// This kernel is executed with 3 * vertices_per_atom threads
__kernel
void eating(__global float *model, float dy)
{
    int index = get_global_id(0);
    int N = get_global_size(0);

    // y coordinate?
    dy *= ((index % 3) == 1) ? 1.0f : 0.0f;
    // vertex belongs to upper part of sphere?
    dy *= (index >= N/2) ? 1.0f : -1.0f;

    model[index] += dy;
}

// This kernel is executed with 3 * vertices_per_atom threads
__kernel
void growing(__global float *model, float factor)
{
    unsigned index = get_global_id(0);
    unsigned N = get_global_size(0);

    // we only modify ghosts' shape, so we add N
    model[index + N] *= factor;
}

// This kernel updates vertices in whole vbo using vertices for 2 atom "models" : pacman & ghost
// It is executed by total number of vertices * 3
//
__kernel
void update_vertices(__global float *vbo, __global float *pos,
		     unsigned offset,
		     __global float *model, unsigned vertices_per_atom)
{
    unsigned index = get_global_id(0);
    unsigned tpa = 3 * vertices_per_atom;
    unsigned slice = (index % 3) * offset;
    __global float *p = pos + slice;
    unsigned idx = index / tpa;

    vbo[index] = model[index % (2 * tpa)] + p[idx];
}
