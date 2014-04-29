#include "default_defines.h"
#include "global_definitions.h"
#include "device.h"
#include "openmp.h"
#include "sotl.h"

#ifdef HAVE_LIBGL
#include "vbo.h"
#endif

#include <stdio.h>
#include <omp.h>

static int *atom_state = NULL;

#ifdef HAVE_LIBGL

#define SHOCK_PERIOD  50
float couleur[24][3] = {
	{1.0,1.0,1.0},
	{1.0,1.0,0.0},
	{1.0,0.0,0.0},
	{1.0,0.0,1.0},
	{0.0,0.0,1.0},
	{0.0,1.0,1.0},
	{0.0,1.0,0.0},
	
	{0.5,0.5,0.5},
	{0.5,0.5,0.0},
	{0.5,0.0,0.0},
	{0.5,0.0,0.5},
	{0.0,0.0,0.5},
	{0.0,0.5,0.5},
	{0.0,0.5,0.0},
	
	{0.25,0.25,0.25},
	{0.25,0.25,0.0},
	{0.25,0.0,0.0},
	{0.25,0.0,0.25},
	{0.0,0.0,0.25},
	{0.0,0.25,0.25},
	{0.0,0.5,0.0},
	
	{1.0,0.0,0.25},
	{1.0,0.25,0.0},
	{0.0,0.25,1.0}
	
	};


// Update OpenGL Vertex Buffer Object
//
static void omp_update_vbo (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  //sotl_domain_t *domain = &dev->domain;
#pragma omp for
  for (unsigned n = 0; n < set->natoms; n++) {
    vbo_vertex[n*3 + 0] = set->pos.x[n];
    vbo_vertex[n*3 + 1] = set->pos.y[n];
    vbo_vertex[n*3 + 2] = set->pos.z[n];

	
    // Atom color depends on z coordinate
    {
      //float ratio = (float)atom_state[n];
      //(set->pos.z[n] - domain->min_ext[2]) / (domain->max_ext[2] - domain->min_ext[2]);

      vbo_color[n*3 + 0] = couleur[atom_state[n]][0];//(1.0 - ratio) * atom_color[0].R + ratio * 1.0;
      vbo_color[n*3 + 1] = couleur[atom_state[n]][1];//(1.0 - ratio) * atom_color[0].G + ratio * 0.0;
      vbo_color[n*3 + 2] = couleur[atom_state[n]][2];//(1.0 - ratio) * atom_color[0].B + ratio * 0.0;
      //atom_state[n]--;
    }
  }
}
#endif

// Update positions of atoms by adding (dx, dy, dz)
//
static void omp_move (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  #pragma omp for
  for (unsigned n = 0; n < set->natoms; n++) {
    set->pos.x[n] += set->speed.dx[n];
    set->pos.y[n] += set->speed.dy[n];
    set->pos.z[n] += set->speed.dz[n];
  }
}

// Apply gravity force
//
static void omp_gravity (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  const calc_t g = 0.005;

#pragma omp for
  for (unsigned n = 0; n < set->natoms; n++) {
  
  	set->speed.dy[n] -= g;
  }
}

static void omp_bounce (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  sotl_domain_t *domain = &dev->domain;

  //TODO
  
  #pragma omp for
  for (unsigned n = 0; n < set->natoms; n++) {
  
  	for(int i = 0; i <3 ; i++)
  	{
  		if(set->pos.x[n+set->offset*i]< domain->min_ext[i])
  		{
  			set->speed.dx[n+set->offset*i]*= -1;
  			//atom_state[n] = SHOCK_PERIOD;
  		}
  		if(set->pos.x[n+set->offset*i] > domain->max_ext[i])
  		{
  			set->speed.dx[n+set->offset*i]*= -1;
  			//atom_state[n] = SHOCK_PERIOD;
		}
  	
  	}
  }
}

static calc_t squared_distance (sotl_atom_set_t *set, unsigned p1, unsigned p2)
{
  calc_t *pos1 = set->pos.x + p1,
    *pos2 = set->pos.x + p2;

  calc_t dx = pos2[0] - pos1[0],
         dy = pos2[set->offset] - pos1[set->offset],
         dz = pos2[set->offset*2] - pos1[set->offset*2];

  return dx * dx + dy * dy + dz * dz;
}

static calc_t lennard_jones (calc_t r2)
{
  calc_t rr2 = 1.0 / r2;
  calc_t r6;

  r6 = LENNARD_SIGMA * LENNARD_SIGMA * rr2;
  r6 = r6 * r6 * r6;

  return 24 * LENNARD_EPSILON * rr2 * (2.0f * r6 * r6 - r6);
}

static void swap(sotl_atom_set_t *set, int i, int j){
	float pt;
	float spt;
	int st;
	for(int k=0;k<3;k++){
		//#pragma omp critical (pos)
		{
		pt = set->pos.x[set->offset * k + j];
		set->pos.x[set->offset * k + j] = set->pos.x[set->offset * k + i];
		set->pos.x[set->offset * k + i] = pt;
		}
		
		//#pragma omp critical (speed)
		{
		spt = set->speed.dx[set->offset * k + j];
		set->speed.dx[set->offset * k + j] = set->speed.dx[set->offset * k + i];
		set->speed.dx[set->offset * k + i] = spt;
		}
	}
	//#pragma omp critical (state)
	{
	st = atom_state[j];
	atom_state[j] = atom_state[i];
	atom_state[i] = st;
	}
}

static void sortBubble(sotl_atom_set_t *set){ // from wikipedia
    int n = set->natoms;
    int swaped = 1;
    while(swaped == 1){
		swaped = 0;
		/*float lol;
	    float kil;
	    int oldI;*/
		for(int i = 1 ; i < n; i++){
		   float old = set->pos.x[set->offset * 2 + (i-1)];
		   float cur = set->pos.x[set->offset * 2 + i];
		   
		   if(old > cur){
			   swap(set, i, i-1);
			   /*for(int k=0;k<3;k++){
					lol = set->pos.x[set->offset * k + (i-1)];
					set->pos.x[set->offset * k + (i-1)] = set->pos.x[set->offset * k + i];
					set->pos.x[set->offset * k + i] = lol;
					
					kil = set->speed.dx[set->offset * k + (i-1)];
					set->speed.dx[set->offset * k + (i-1)] = set->speed.dx[set->offset * k + i];
					set->speed.dx[set->offset * k + i] = kil;
				}
				oldI = atom_state[i-1];
				atom_state[i-1] = atom_state[i];
				atom_state[i] = oldI;*/
				swaped = 1;
			}
		}
       
    }
}



static void oddEvenSort(sotl_atom_set_t *set){ // from wikipedia
	
	//#pragma omp single
	//{
	int sorted = 0;
	int n = set->natoms;

	
	//#pragma omp master
	//{
	while(!sorted)
	{
		sorted=1;
		#pragma omp for
		for(int i = 1; i < n-1; i += 2)
		{
			float next = set->pos.x[set->offset * 2 + (i+1)];
			float cur = set->pos.x[set->offset * 2 + i];
			if(cur > next)
			{
				swap(set, i, i+1);
				sorted = 0;
			}
		}
		#pragma omp barrier
	 
		#pragma omp for
		for(int i = 0; i < n-1; i += 2)
		{
			float next = set->pos.x[set->offset * 2 + (i+1)];
			float cur = set->pos.x[set->offset * 2 + i];
			if(cur > next)
			{
				swap(set, i, i+1);
				sorted = 0;
			}
		}
		#pragma omp barrier
	}
}


static void sortAtomBox(unsigned *atom_box_pos, unsigned *atom_box_size, unsigned offset_box, sotl_atom_set_t *set){
	int n = set->natoms;
	//#pragma omp for 
	for(int i = 0 ; i < n; i++){
		int posTmp[3];
		/*for(int j=0;j<3;j++){
			posTmp[i] = set->pos.x[set->offset * j + i] / LENNARD_SQUARED_CUTOFF; //position in box axis X or Y or Z
		}*/
		atom_box_pos[posTmp[0] + offset_box * 1 * posTmp[1] + offset_box * 2 * posTmp[2]][atom_box_size] = i;
		atom_box_size[posTmp[0] + offset_box * 1 * posTmp[1] + offset_box * 2 * posTmp[2]]++;
	}
}

static void boxForce(sotl_device_t *dev){
	//#pragma omp single
	//{
		
		unsigned *atom_box_pos = NULL;
		unsigned *atom_box_size = NULL;
		sotl_atom_set_t *set = &dev->atom_set;
		sotl_domain_t *domain = &dev->domain;
		int n = set->natoms;
		
		unsigned offset_box = domain->max_ext[0]/LENNARD_SQUARED_CUTOFF; // work only on cube !!!
		unsigned nb_box = offset_box * 3;
		atom_box_pos = malloc(set->natoms * nb_box * sizeof(int)); // Each box can contains all the atoms
		atom_box_size = malloc(nb_box * sizeof(int));
	//}
	
	//#pragma omp for 
	for(int i = 0 ; i < n+1; i++){ // init size to 0 for each box
			atom_box_size[i] = 0;
	}
	
	sortAtomBox(atom_box_pos, atom_box_size, offset_box, set);
	
	//#pragma omp for 
	/*for(unsigned boxNb = 0; boxNb < nb_box; boxNb++) {
		//int *box_pos = atom_box_pos[boxNb];
		unsigned box_size = atom_box_size[boxNb];
		
		for (unsigned current = 0; current < box_size; current++) {
			calc_t force[3] = { 0.0, 0.0, 0.0 };
			
			for(int oX = -1; oX < 2; oX++){ // compare with other box (27)
				for(int oY = -1; oY < 2; oY++){
					for(int oZ = -1; oZ < 2; oZ++){
						unsigned idBox = (boxNb + oX + oY * offset_box + oZ * offset_box);
						if(idBox >= 0 && idBox < nb_box){ // if no array overflow 
							//int *box_pos_other = atom_box_pos[idBox];
							for (unsigned other = 0; other < atom_box_size[idBox]; other++){
								if (current != other) {
									calc_t sq_dist = squared_distance (set, current, other);
									if (sq_dist < LENNARD_SQUARED_CUTOFF) {
										calc_t intensity = lennard_jones (sq_dist);

										force[0] += intensity * (set->pos.x[current] - set->pos.x[other]);
										force[1] += intensity * (set->pos.x[set->offset + current] -
													 set->pos.x[set->offset + other]);
										force[2] += intensity * (set->pos.x[set->offset * 2 + current] -
													 set->pos.x[set->offset * 2 + other]);
									}
								}
							}
						}
					}
				}
			}
			atom_state[current] = omp_get_thread_num();
			set->speed.dx[current] += force[0];
			set->speed.dx[set->offset + current] += force[1];
			set->speed.dx[set->offset * 2 + current] += force[2];
		}
	}*/
	free(atom_box_pos);
	free(atom_box_size);
//}
}

static void omp_forceOrg(sotl_device_t *dev){
  sotl_atom_set_t *set = &dev->atom_set;
  #pragma omp for 
  for (unsigned current = 0; current < set->natoms; current++) {
    calc_t force[3] = { 0.0, 0.0, 0.0 };
	
    for (unsigned other = 0; other < set->natoms; other++)
      if (current != other) {
	calc_t sq_dist = squared_distance (set, current, other);

	if (sq_dist < LENNARD_SQUARED_CUTOFF) {
	  calc_t intensity = lennard_jones (sq_dist);

	  force[0] += intensity * (set->pos.x[current] - set->pos.x[other]);
	  force[1] += intensity * (set->pos.x[set->offset + current] -
				   set->pos.x[set->offset + other]);
	  force[2] += intensity * (set->pos.x[set->offset * 2 + current] -
				   set->pos.x[set->offset * 2 + other]);
	}

      }

    atom_state[current] = omp_get_thread_num();
    set->speed.dx[current] += force[0];
    set->speed.dx[set->offset + current] += force[1];
    set->speed.dx[set->offset * 2 + current] += force[2];
  }
	
}


static void omp_force (sotl_device_t *dev)
{
  //sotl_atom_set_t *set = &dev->atom_set;
  
  #pragma omp single
  boxForce(dev);
  
    //*
  //sort table
  //#pragma omp single
  //sortBubble(set);
  
  //#pragma omp single
  //oddEvenSort(set);
  
  //compute
  /*
  #pragma omp for
  for (unsigned current = 0; current < set->natoms; current++) {
    calc_t force[3] = { 0.0, 0.0, 0.0 };
    for (unsigned other = 0; other < set->natoms; other++){
		//if(set->pos.x[set->offset * 2 + other] > (set->pos.x[set->offset * 2 + current] + LENNARD_SQUARED_CUTOFF)) break;
		
		if (current != other) {
			calc_t sq_dist = squared_distance (set, current, other);
			if (sq_dist < LENNARD_SQUARED_CUTOFF) {
				calc_t intensity = lennard_jones (sq_dist);

				force[0] += intensity * (set->pos.x[current] - set->pos.x[other]);
				force[1] += intensity * (set->pos.x[set->offset + current] -
							 set->pos.x[set->offset + other]);
				force[2] += intensity * (set->pos.x[set->offset * 2 + current] -
							 set->pos.x[set->offset * 2 + other]);
			}
		}
	}
	atom_state[current] = omp_get_thread_num();
	set->speed.dx[current] += force[0];
	set->speed.dx[set->offset + current] += force[1];
	set->speed.dx[set->offset * 2 + current] += force[2];
  
	}//*/
  
//old stuff
//*
  //omp_forceOrg(dev);

  //*/
}


// Main simulation function
//
void omp_one_step_move (sotl_device_t *dev)
{

#pragma omp parallel default(shared)
{
  // Apply gravity force
  //
  if (gravity_enabled)
    omp_gravity (dev);

  // Compute interactions between atoms
  //
  if (force_enabled)
    omp_force (dev);

  // Bounce on borders
  //
  if(borders_enabled)
    omp_bounce (dev);

  // Update positions
  //
  omp_move (dev);

#ifdef HAVE_LIBGL
  // Update OpenGL position
  //
  if (dev->display)
    omp_update_vbo (dev);
#endif
}
}

void omp_init (sotl_device_t *dev)
{
#ifdef _SPHERE_MODE_
  sotl_log(ERROR, "Sequential implementation does currently not support SPHERE_MODE\n");
  exit (1);
#endif

  borders_enabled = 1;

  dev->compute = SOTL_COMPUTE_OMP; // dummy op to avoid warning
}

void omp_alloc_buffers (sotl_device_t *dev)
{
  atom_state = calloc(dev->atom_set.natoms, sizeof(int));
  printf("natoms: %d\n", dev->atom_set.natoms);
}

void omp_finalize (sotl_device_t *dev)
{
  free(atom_state);

  dev->compute = SOTL_COMPUTE_SEQ; // dummy op to avoid warning
}
