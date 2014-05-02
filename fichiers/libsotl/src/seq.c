
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


typedef struct box {
	unsigned nbBox; //nombre de box
	unsigned nbBoxLigne; //nombre de box sur un axe
	int * nbAtomToBox; // nombre d'atome par box
	int * preNbAtomToBox; // nombre prefix d'atome par box
	calc_t * swapPosx;
	calc_t * swapSpeedx;
	int * swapState;
	unsigned * numboxAtom; //numero de la box pour un atom.
} box_t;

static box_t * boxSort = NULL;

// Update OpenGL Vertex Buffer Object
//
static void seq_update_vbo (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  //sotl_domain_t *domain = &dev->domain;

  for (unsigned n = 0; n < set->natoms; n++) {
    vbo_vertex[n*3 + 0] = set->pos.x[n];
    vbo_vertex[n*3 + 1] = set->pos.y[n];
    vbo_vertex[n*3 + 2] = set->pos.z[n];

    // Atom color depends on z coordinate
    if(atom_state[n]  > 0)
    {
      float ratio =  (float)atom_state[n]/ SHOCK_PERIOD; //(set->pos.z[n] - domain->min_ext[2]) / (domain->max_ext[2] - domain->min_ext[2]);
//	float ratior = (float)atom_state[n]/boxSort->nbBox;
	
      vbo_color[n*3 + 0] = (1.0 - ratio) * atom_color[0].R + ratio * 1.0;
      vbo_color[n*3 + 1] = (1.0 - ratio) * atom_color[0].G + ratio * 0.0;
      vbo_color[n*3 + 2] = (1.0 - ratio) * atom_color[0].B + ratio * 0.5;
      atom_state[n]--;
    }
  }
}
#endif

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

  //TODO
  for (unsigned n = 0; n < set->natoms; n++) {
  
  
  	for(int i = 0; i <3 ; i++)
  	{
  		if(set->pos.x[n+set->offset*i]<= domain->min_ext[i])
  		{
  		
  			//TODO pour éviter d'avoir des atome qui ne sont pas dans des boite on force les cordonné au plus/moins = max_ext/min_ext
  			set->pos.x[n+set->offset*i] = domain->min_ext[i];
  			
  			if(set->speed.dx[n+set->offset*i] < 0)
  			{set->speed.dx[n+set->offset*i]*= -1;
  			atom_state[n] = SHOCK_PERIOD;
  			}
  		}
  		if(set->pos.x[n+set->offset*i] >= domain->max_ext[i])
  		{
  		
  			set->pos.x[n+set->offset*i] = domain->max_ext[i];
  		
  			if(set->speed.dx[n+set->offset*i] > 0){
  			set->speed.dx[n+set->offset*i]*= -1;
  			atom_state[n] = SHOCK_PERIOD;
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

static void sortBubble(sotl_atom_set_t *set)
{ // from wikipedia
	unsigned n = set->natoms;
	int swap = 1;

	float swapPos;
	float swapSpeed;
	int swapState;
	float old ;
	float cur;
	int offset;
	
	while(swap == 1)
	{
		swap = 0;
		
		for(unsigned i = 1 ; i < n; i++)
		{
			old = set->pos.x[set->offset * 2 + (i-1)];
			cur = set->pos.x[set->offset * 2 + i];
			if(old > cur)
			{
				for(int k=0;k<3;k++)
				{
					offset = set->offset * k;
					
					swapPos = set->pos.x[offset + (i-1)];
					set->pos.x[offset + (i-1)] = set->pos.x[offset + i];
					set->pos.x[offset + i] = swapPos;
					
					swapSpeed = set->speed.dx[offset + (i-1)];
					set->speed.dx[offset + (i-1)] = set->speed.dx[offset + i];
					set->speed.dx[offset + i] = swapSpeed;
				}
				
				//*
				swapState = atom_state[i-1];
				atom_state[i-1] = atom_state[i];
				atom_state[i] = swapState;
				
				//*/
				
				swap = 1;
			}
		}
	}
}

static unsigned xyzToNumBox(unsigned x, unsigned y, unsigned z)
{
	//return x + (LENNARD_SQUARED_CUTOFF * y) + (2 * LENNARD_SQUARED_CUTOFF * z) ;
	//return x + (LENNARD_SQUARED_CUTOFF * y) + LENNARD_SQUARED_CUTOFF * z ;
	return x + y * boxSort->nbBoxLigne + z * ( boxSort->nbBoxLigne * boxSort->nbBoxLigne);
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


//TODO fonction de dev inutile remplacer par switchPtr
static void copiePtr( sotl_atom_set_t *set)
{

	for (unsigned n = 0; n < set->natoms; n++ )
	{
		atom_state[n] = boxSort->swapState[n];
		for(int j = 0 ; j < 3 ; j++)
		{
			set->pos.x[n + j * set->offset] = boxSort->swapPosx[n + j * set->offset];
			set->speed.dx[n + j * set->offset] = boxSort->swapSpeedx[n + j * set->offset];
		}
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


static void sortAtomBox(sotl_atom_set_t *set){
	int n = set->natoms;
	
	resTab(boxSort->nbAtomToBox, boxSort->nbBox);
	resTab(boxSort->preNbAtomToBox, boxSort->nbBox+1);
	
	unsigned * numbox = boxSort->numboxAtom;
	unsigned posTmp[3];
	
	for(int i = 0 ; i < n; i++)
	{
	
		for(int j=0;j<3;j++)
		{
			posTmp[j] = (unsigned)(set->pos.x[set->offset * j + i] / LENNARD_SQUARED_CUTOFF); //position in box axis X or Y or Z
		
		}
		numbox[i] = xyzToNumBox(posTmp[0],posTmp[1],posTmp[2]);
		if(numbox[i] < boxSort-> nbBox)
		
		boxSort->nbAtomToBox[numbox[i]]++;
		
	}
	
	prefixTab(boxSort->nbAtomToBox, boxSort->preNbAtomToBox, boxSort->nbBox);
	resTab(boxSort->nbAtomToBox, boxSort->nbBox);
	
	for(int i = 0 ; i < n; i++)
	{
		
		//TODO si l'atome sort du cub mon ne peut plus le calculer
		if(numbox[i] < boxSort-> nbBox)
			moveAtomBox(i,numbox[i],set);
		
	}
	
	switchPtr(set);

}

static void boxForceTest(sotl_device_t *dev)
{

	sotl_atom_set_t *set = &dev->atom_set;
	sotl_domain_t *domain = &dev->domain;

	sortAtomBox(set);
	int cubeSize = 3;
	unsigned numBox,numBoxAdj;
	//float maxSize = LENNARD_SQUARED_CUTOFF * boxSort->nbBoxLigne;
	printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
	int patacoin = 0;
	for(int x = 0; x < cubeSize; x++)
		for(int y = 0; y < cubeSize; y++)
			for(int z = 0; z < cubeSize; z++)
			{
				//numBox = xyzToNumBox(x,y,z);
				for(int xadj = -1; xadj < 2; xadj++)
					for(int yadj = -1; yadj < 2; yadj++)
						for(int zadj = -1; zadj < 2; zadj++)
						{
							int xboxadj = x +xadj;
							int yboxadj = y +yadj;
							int zboxadj = z +zadj;
							if(xboxadj > -1 && xboxadj < cubeSize  && yboxadj > -1 && yboxadj < cubeSize  && zboxadj > -1 && zboxadj < cubeSize)
							{
								patacoin++;
								printf("main : x %d, y %d, z %d\n", x, y , z);
								printf("other : x %d, y %d, z %d\n", xboxadj, yboxadj , zboxadj);
								//numBoxAdj = xyzToNumBox(xboxadj,yboxadj,zboxadj);
								/*
								for (unsigned current = boxSort->preNbAtomToBox[numBox]; current < boxSort->preNbAtomToBox[numBox+1]; current++) 
								{
									calc_t force[3] = { 0.0, 0.0, 0.0 };

									for (unsigned other = boxSort->preNbAtomToBox[numBoxAdj]; other < boxSort->preNbAtomToBox[numBoxAdj+1]; other++)
									{
										if (current != other) 
										{
											calc_t sq_dist = squared_distance (set, current, other);

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
								}*/
							}
						}
			}
			printf("pa %d \n", patacoin); 
			exit(0);
}

static void boxForce(sotl_device_t *dev)
{boxForceTest(dev);
/*
	sotl_atom_set_t *set = &dev->atom_set;
	sotl_domain_t *domain = &dev->domain;

	sortAtomBox(set);
	
	unsigned numBox,numBoxAdj;
	//float maxSize = LENNARD_SQUARED_CUTOFF * boxSort->nbBoxLigne;
	for(int x = 0; x < boxSort->nbBoxLigne; x++)
		for(int y = 0; y < boxSort->nbBoxLigne; y++)
			for(int z = 0; z < boxSort->nbBoxLigne; z++)
			{
				numBox = xyzToNumBox(x,y,z);
				for(int xadj = -1; xadj < 2; xadj++)
					for(int yadj = -1; yadj < 2; yadj++)
						for(int zadj = -1; zadj < 2; zadj++)
						{
							int xboxadj = x +xadj;
							int yboxadj = y +yadj;
							int zboxadj = z +zadj;
							if(xboxadj > -1 && xboxadj < boxSort->nbBoxLigne  && yboxadj > -1 && yboxadj < boxSort->nbBoxLigne  && zboxadj > -1 && zboxadj < boxSort->nbBoxLigne)
							{
								numBoxAdj = xyzToNumBox(xboxadj,yboxadj,zboxadj);
								
								for (unsigned current = boxSort->preNbAtomToBox[numBox]; current < boxSort->preNbAtomToBox[numBox+1]; current++) 
								{
									calc_t force[3] = { 0.0, 0.0, 0.0 };

									for (unsigned other = boxSort->preNbAtomToBox[numBoxAdj]; other < boxSort->preNbAtomToBox[numBoxAdj+1]; other++)
									{
										if (current != other) 
										{
											calc_t sq_dist = squared_distance (set, current, other);

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
	*/
	/*
	for(unsigned boxNb = 0; boxNb < boxSort->nbBox; boxNb++) {
		//int *box_pos = atom_box_pos[boxNb];
		unsigned box_size = boxSort->nbAtomToBox[boxNb];
		
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
			//atom_state[current] = omp_get_thread_num();
			set->speed.dx[current] += force[0];
			set->speed.dx[set->offset + current] += force[1];
			set->speed.dx[set->offset * 2 + current] += force[2];
		}
	}
	//*/
}


static void seq_force (sotl_device_t *dev)
{
/*	sotl_atom_set_t *set = &dev->atom_set;*/
	boxForce(dev);

	//sortBubble(set);

/*
	for (unsigned current = 0; current < set->natoms; current++) 
	{
		calc_t force[3] = { 0.0, 0.0, 0.0 };

		for (unsigned other = 0; other < set->natoms; other++)
		{
			if (current != other) 
			{
				calc_t sq_dist = squared_distance (set, current, other);

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
	//*/
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
	sotl_atom_set_t *set = &dev->atom_set;
	sotl_domain_t *domain = &dev->domain;
	
	
	unsigned nbBoxLigne = (domain->max_ext[0]/LENNARD_SQUARED_CUTOFF) +1; // work only on cube !!!
	
	boxSort = malloc (sizeof (box_t));
	boxSort->nbBoxLigne = nbBoxLigne;
	boxSort->nbBox = nbBoxLigne *nbBoxLigne * nbBoxLigne ;
	
	
	boxSort->nbAtomToBox = malloc(boxSort->nbBox * sizeof(int));
	boxSort->preNbAtomToBox = malloc((boxSort->nbBox +1) * sizeof(int));
	
	boxSort->swapPosx = malloc(sizeof(calc_t) * set->offset * 3);
	boxSort->swapSpeedx = malloc(sizeof(calc_t) * set->offset * 3);
	boxSort->swapState = malloc(set->natoms * sizeof(int));
	boxSort->numboxAtom = malloc(set->natoms * sizeof(unsigned));
	
  atom_state = calloc(set->natoms, sizeof(int));
  printf("natoms: %d\n", set->natoms);
}

void seq_finalize (sotl_device_t *dev)
{
  free(atom_state);
  free(boxSort->nbAtomToBox);
  free(boxSort->preNbAtomToBox);
  free(boxSort->swapPosx);
  free(boxSort->swapSpeedx);
  free(boxSort->swapState);
  free(boxSort);
  free(boxSort->numboxAtom);

  dev->compute = SOTL_COMPUTE_SEQ; // dummy op to avoid warning
}
