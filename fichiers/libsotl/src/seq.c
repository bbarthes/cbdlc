
#include "default_defines.h"
#include "global_definitions.h"
#include "device.h"
#include "seq.h"
#include "sotl.h"

#ifdef HAVE_LIBGL
#include "vbo.h"
#endif

#include <stdio.h>

static int *atom_state = NULL;

#ifdef HAVE_LIBGL

#define SHOCK_PERIOD  50

// Update OpenGL Vertex Buffer Object
//

static void seq_update_vbo (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  sotl_domain_t *dom = &dev->domain;

  for (unsigned n = 0; n < set->natoms; n++) {
    vbo_vertex[n*3 + 0] = set->pos.x[n];
    vbo_vertex[n*3 + 1] = set->pos.y[n];
    vbo_vertex[n*3 + 2] = set->pos.z[n];

    // Atom color depends on z coordinate
   // if(atom_state[n]  > 0)
    {
            float ratio =  (float)atom_state[n]/dom->total_boxes;
      //float ratio = (set->pos.z[n] - domain->min_ext[2]) / (domain->max_ext[2] - domain->min_ext[2]);

      vbo_color[n*3 + 0] = (1.0 - ratio) * atom_color[0].R + ratio * 1.0;
      vbo_color[n*3 + 1] = (1.0 - ratio) * atom_color[0].G + ratio * 0.5;
      vbo_color[n*3 + 2] = (1.0 - ratio) * atom_color[0].B + ratio * 0.0;
     // printf("%d\n",atom_state[n]);
      //atom_state[n]--;
    }
  }
}
#endif

typedef struct box {
	
	int * nbAtomToBox; // nombre d'atome par box
	int * preNbAtomToBox; // nombre prefix d'atome par box
	calc_t * swapPosx;
	calc_t * swapSpeedx;
	int * swapState;
	//unsigned * numboxAtom; //numero de la box pour un atom.
} box_t;

static box_t * boxSort = NULL;

// Update positions of atoms by adding (dx, dy, dz)
//
static void seq_move (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;

  for (unsigned n = 0; n < set->natoms; n++) {
    set->pos.x[n] += set->speed.dx[n];
    set->pos.y[n] += set->speed.dy[n];
    set->pos.z[n] += set->speed.dz[n];
  }
}

// Apply gravity force
//
static void seq_gravity (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  const calc_t g = 0.005;

  for (unsigned n = 0; n < set->natoms; n++) {
  
  	set->speed.dy[n] -= g;
  }
}

static void seq_bounce (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  sotl_domain_t *domain = &dev->domain;

  for (unsigned n = 0; n < set->natoms; n++) {
  	for(int i = 0; i <3 ; i++)
  	{
  		if(set->pos.x[n+set->offset*i]<= domain->min_ext[i])
  		{
  			//TODO pour éviter d'avoir des atome qui ne sont pas dans des boite on force les cordonné au plus/moins = max_ext/min_ext
  			set->pos.x[n+set->offset*i] = domain->min_ext[i];
  			
  			if(set->speed.dx[n+set->offset*i] < 0)
  			{set->speed.dx[n+set->offset*i]*= -1;
            //atom_state[n] = SHOCK_PERIOD;
  			}
  		}
  		if(set->pos.x[n+set->offset*i] >= domain->max_ext[i])
  		{
  			set->pos.x[n+set->offset*i] = domain->max_ext[i];
  		
  			if(set->speed.dx[n+set->offset*i] > 0){
				set->speed.dx[n+set->offset*i]*= -1;
                //atom_state[n] = SHOCK_PERIOD;
  			}
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

static void seq_force_old (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;

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

    set->speed.dx[current] += force[0];
    set->speed.dx[set->offset + current] += force[1];
    set->speed.dx[set->offset * 2 + current] += force[2];
  }
}


static void switchPtr ( sotl_atom_set_t *set)
{
	int * spState =  atom_state;
	calc_t * spPosx = set->pos.x;
	calc_t * spSpeedx = set->speed.dx;
	
	atom_state = boxSort->swapState;
	
	set->pos.x = boxSort->swapPosx;
	set->pos.y = set->pos.x + set->offset;
	set->pos.z = set->pos.y + set->offset;

	set->speed.dx = boxSort->swapSpeedx;
	set->speed.dy = set->speed.dx + set->offset;
	set->speed.dz = set->speed.dy + set->offset;
	
	boxSort->swapState = spState;
	boxSort->swapPosx = spPosx;
	boxSort->swapSpeedx = spSpeedx;


}

static void resTab(int * tab , unsigned sizeofTab)
{
    for(unsigned i = 0 ; i < sizeofTab; i++)
    {
        tab[i]= 0;
    }
}

static void prefixTab(int * source , int * dest, unsigned size)
{
    dest[0] = 0;
    for(unsigned i =0 ; i < size; i++)
    {
        dest[i+1] = dest[i] + source[i];
    }
}

static void moveAtomBox(int nAtom, int numbox, sotl_atom_set_t *set)
{

    int newPos = boxSort->preNbAtomToBox[numbox] + boxSort->nbAtomToBox[numbox];

    for(unsigned j = 0 ; j < 3 ; j++)
    {
        int offset = set->offset * j;
        boxSort->swapPosx[newPos + offset] = set->pos.x[nAtom + offset];
        boxSort->swapSpeedx[newPos + offset] = set->speed.dx[nAtom +offset];
    }
    boxSort->swapState[newPos] = atom_state[nAtom];
    boxSort->nbAtomToBox[numbox]++;
}

static void sortAtomBox( sotl_domain_t *dom , sotl_atom_set_t *set)
{
	int n = set->natoms;
	
	
    free(boxSort->nbAtomToBox);
    boxSort->nbAtomToBox = NULL;
    boxSort->nbAtomToBox = atom_set_box_count(dom,set);

    prefixTab(boxSort->nbAtomToBox, boxSort->preNbAtomToBox, dom->total_boxes);
	resTab(boxSort->nbAtomToBox, dom->total_boxes);

	for(int i = 0 ; i < n; i++)
	{
        int pos =atom_get_num_box(dom,set->pos.x[i], set->pos.y[i], set->pos.z[i],BOX_SIZE_INV);
        moveAtomBox(i,pos,set);
        boxSort->swapState[i] = pos;

	}

	switchPtr(set);

}

int get_num_box(const sotl_domain_t *dom, const int x, const int y,
                     const int z)
{
    int box_x, box_y, box_z;
    int box_id;
    int rrc  = BOX_SIZE_INV;

    box_x = x * rrc;
    box_y = y * rrc;
    box_z = z * rrc;

    box_id =  box_z * dom->boxes[0] * dom->boxes[1] +
              box_y * dom->boxes[0] +
              box_x;

    if(box_id < 0 && (unsigned)box_id < dom->total_boxes)
        return -1;
		
    return box_id;
}


static void seq_force_cube (sotl_device_t *dev)
{
	sotl_atom_set_t *set = &dev->atom_set;
	sotl_domain_t *dom = &dev->domain;
	
	sortAtomBox(dom,set);
	


	
    //*
    unsigned currentBox;
    int otherBox;
    calc_t sq_dist;

	for(int x = 0; x < dom->boxes[0]; x++)
        for(int y = 0; y < dom->boxes[1]; y++)
			for(int z = 0; z < dom->boxes[2]; z++)
			{
                currentBox = get_num_box(dom,x,y,z);
				
				for(int xadj = -1; xadj < 2; xadj++)
					for(int yadj = -1; yadj < 2; yadj++)
						for(int zadj = -1; zadj < 2; zadj++)
						{

							int xother = x +xadj;
							int yother = y +yadj;
							int zother = z +zadj;
                            otherBox = get_num_box(dom,xother,yother,zother);
                            if(otherBox > -1)
							{

                              //  printf("current %d , other %d  \n", currentBox,otherBox);

								for (unsigned current = boxSort->preNbAtomToBox[currentBox]; current < boxSort->preNbAtomToBox[currentBox+1]; current++) 
								{
                                  //  printf("bite %d\n",current);
									calc_t force[3] = { 0.0, 0.0, 0.0 };

									for (unsigned other = boxSort->preNbAtomToBox[otherBox]; other < boxSort->preNbAtomToBox[otherBox+1]; other++)
									{
                                      //  printf("coucou %d\n",other);
										if (current != other) 
										{
                                            sq_dist = squared_distance (set, current, other);

											if (sq_dist < LENNARD_SQUARED_CUTOFF) 
											{
											  calc_t intensity = lennard_jones (sq_dist);

											  force[0] += intensity * (set->pos.x[current] - set->pos.x[other]);
											  force[1] += intensity * (set->pos.x[set->offset + current] -
														   set->pos.x[set->offset + other]);
											  force[2] += intensity * (set->pos.x[set->offset * 2 + current] -
														   set->pos.x[set->offset * 2 + other]);
											}

										}
									}
									set->speed.dx[current] += force[0];
									set->speed.dx[set->offset + current] += force[1];
									set->speed.dx[set->offset * 2 + current] += force[2];
								}
							}
						}
			}
            //*/ // seq_force_old(dev);
	
}

static void seq_force (sotl_device_t *dev)
{
    //seq_force_old(dev);
    seq_force_cube(dev);
}


// Main simulation function
//
void seq_one_step_move (sotl_device_t *dev)
{
  // Apply gravity force
  //
  if (gravity_enabled)
    seq_gravity (dev);

  // Compute interactions between atoms
  //
  if (force_enabled)
    seq_force (dev);

  // Bounce on borders
  //
  if(borders_enabled)
    seq_bounce (dev);

  // Update positions
  //
  seq_move (dev);

#ifdef HAVE_LIBGL
  // Update OpenGL position
  //
  if (dev->display)
    seq_update_vbo (dev);
#endif
}

void seq_init (sotl_device_t *dev)
{
#ifdef _SPHERE_MODE_
  sotl_log(ERROR, "Sequential implementation does currently not support SPHERE_MODE\n");
  exit (1);
#endif

  borders_enabled = 1;

  dev->compute = SOTL_COMPUTE_SEQ; // dummy op to avoid warning
}

void seq_alloc_buffers (sotl_device_t *dev)
{
    sotl_domain_t *dom = &dev->domain;
    sotl_atom_set_t *set = &dev->atom_set;

	boxSort = malloc (sizeof (box_t));
		
    boxSort->nbAtomToBox = malloc(dom->total_boxes * sizeof(int));
    boxSort->preNbAtomToBox = malloc((dom->total_boxes +1) * sizeof(int));
	
	boxSort->swapPosx = malloc(sizeof(calc_t) * set->offset * 3);
	boxSort->swapSpeedx = malloc(sizeof(calc_t) * set->offset * 3);
	boxSort->swapState = malloc(set->natoms * sizeof(int));
	//boxSort->numboxAtom = malloc(set->natoms * sizeof(unsigned));
	
  atom_state = calloc(dev->atom_set.natoms, sizeof(int));
  printf("natoms: %d\n", dev->atom_set.natoms);
}

void seq_finalize (sotl_device_t *dev)
{
	free(boxSort->nbAtomToBox);
	free(boxSort->preNbAtomToBox);
	free(boxSort->swapPosx);
	free(boxSort->swapSpeedx);
	free(boxSort->swapState);
	//free(boxSort->numboxAtom);
	free(boxSort);
	
	
  free(atom_state);

  dev->compute = SOTL_COMPUTE_SEQ; // dummy op to avoid warning
}
