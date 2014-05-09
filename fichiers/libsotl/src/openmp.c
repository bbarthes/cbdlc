#include "default_defines.h"
#include "global_definitions.h"
#include "device.h"
#include "openmp.h"
#include "sotl.h"
#include <assert.h>

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
  sotl_domain_t *dom = &dev->domain;
  #pragma omp for
  for (unsigned n = 0; n < set->natoms; n++) {
    vbo_vertex[n*3 + 0] = set->pos.x[n];
    vbo_vertex[n*3 + 1] = set->pos.y[n];
    vbo_vertex[n*3 + 2] = set->pos.z[n];

	
    // Atom color depends on z coordinate
    {
      //float ratio = (float)atom_state[n];
      //(set->pos.z[n] - domain->min_ext[2]) / (domain->max_ext[2] - domain->min_ext[2]);
/*
      vbo_color[n*3 + 0] = couleur[atom_state[n]][0];//(1.0 - ratio) * atom_color[0].R + ratio * 1.0;
      vbo_color[n*3 + 1] = couleur[atom_state[n]][1];//(1.0 - ratio) * atom_color[0].G + ratio * 0.0;
      vbo_color[n*3 + 2] = couleur[atom_state[n]][2];//(1.0 - ratio) * atom_color[0].B + ratio * 0.0;
      //atom_state[n]--;*/
      
      float ratio =  (float)atom_state[n]/dom->total_boxes;
      vbo_color[n*3 + 0] = (1.0 - ratio) * atom_color[0].R + ratio * 1.0;
      vbo_color[n*3 + 1] = (1.0 - ratio) * atom_color[0].G + ratio * 0.5;
      vbo_color[n*3 + 2] = (1.0 - ratio) * atom_color[0].B + ratio * 0.0;
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

#define CUBE_SENTCIL_SIZE  81

static int cube_stencil[CUBE_SENTCIL_SIZE] = { // 27 * 3
                            -1,-1,-1,
                            0,-1,-1,
                            1,-1,-1,
                            -1,0,-1,
                            0,0,-1,
                            1,0,-1,
                            -1,1,-1,
                            0,1,-1,
                            1,1,-1,

                            -1,-1,0,
                            0,-1,0,
                            1,-1,0,
                            -1,0,0,
                            0,0,0,
                            1,0,0,
                            -1,1,0,
                            0,1,0,
                            1,1,0,

                            -1,-1,1,
                            0,-1,1,
                            1,-1,1,
                            -1,0,1,
                            0,0,1,
                            1,0,1,
                            -1,1,1,
                            0,1,1,
                            1,1,1
                             };

static box_t * boxSort = NULL;

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
		//printf("omp_gravity thread %d\n",omp_get_thread_num());
  	set->speed.dy[n] -= g;
  }
}

static void omp_bounce (sotl_device_t *dev)
{
    sotl_atom_set_t *set = &dev->atom_set;
    sotl_domain_t *domain = &dev->domain;


    #pragma omp for
    for (unsigned n = 0; n < set->natoms; n++)
    {
        for(int i = 0; i <3 ; i++)
        {
            if(set->pos.x[n+set->offset*i]<= domain->min_ext[i])
            {
                //TODO pour éviter d'avoir des atome qui ne sont pas dans des boite on force les cordonné au plus/moins = max_ext/min_ext
                //set->pos.x[n+set->offset*i] = domain->min_ext[i];

                if(set->speed.dx[n+set->offset*i] < 0)
                {
                    set->speed.dx[n+set->offset*i]*= -1;
                }
            }
            if(set->pos.x[n+set->offset*i] >= domain->max_ext[i])
            {
                //set->pos.x[n+set->offset*i] = domain->max_ext[i];

                if(set->speed.dx[n+set->offset*i] > 0)
                {
                    set->speed.dx[n+set->offset*i]*= -1;
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


/* début trie en z*/

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

static void sortZAtomBox(unsigned *atom_box_pos, unsigned *atom_box_size, unsigned offset_box, sotl_atom_set_t *set){
    int n = set->natoms;
    //#pragma omp for
    for(int i = 0 ; i < n; i++){
        int posTmp[3];
        //*
          for(int j=0;j<3;j++){
            posTmp[i] = set->pos.x[set->offset * j + i] / LENNARD_SQUARED_CUTOFF; //position in box axis X or Y or Z
        }
        //*/
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

    sortZAtomBox(atom_box_pos, atom_box_size, offset_box, set);

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

/* fin trie en z*/


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


/*début trie en boite*/

static void omp_switchPtr ( sotl_atom_set_t *set)
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

static void omp_resTab(int * tab , unsigned sizeofTab)
{
		//#pragma omp for
    for(unsigned i = 0 ; i < sizeofTab; i++)
    {
        tab[i]= 0;
    }
}

static void omp_prefixTab(int * source , int * dest, unsigned size)
{
    dest[0] = 0;
    //#pragma omp for
    for(unsigned i =0 ; i < size; i++)
    {
        dest[i+1] = dest[i] + source[i];
    }
}

static void omp_moveAtomBox(int nAtom, int numbox, sotl_atom_set_t *set)
{

    int newPos = boxSort->preNbAtomToBox[numbox] + boxSort->nbAtomToBox[numbox];

    assert(numbox >= 0);

		//#pragma omp for
    for(unsigned j = 0 ; j < 3 ; j++)
    {
        int offset = set->offset * j;
        boxSort->swapPosx[newPos + offset] = set->pos.x[nAtom + offset];
        boxSort->swapSpeedx[newPos + offset] = set->speed.dx[nAtom +offset];
    }
    
    boxSort->nbAtomToBox[numbox]++;
}

static void omp_sortAtomBox( sotl_domain_t *dom , sotl_atom_set_t *set)
{
	int n = set->natoms;
	int box_id;
	
	#pragma omp single
	{
		//printf("omp_sortAtomBox sing thread %d\n",omp_get_thread_num());
		free(boxSort->nbAtomToBox);
		boxSort->nbAtomToBox = NULL;
		boxSort->nbAtomToBox = atom_set_box_count(dom,set);

	omp_prefixTab(boxSort->nbAtomToBox, boxSort->preNbAtomToBox, dom->total_boxes);
	omp_resTab(boxSort->nbAtomToBox, dom->total_boxes);
	}
	
	#pragma omp barrier
	
	#pragma omp for private(box_id)
	for(int i = 0 ; i < n; i++)
	{
		//printf("omp_sortAtomBox thread %d\n",omp_get_thread_num());
		box_id =atom_get_num_box(dom,set->pos.x[i], set->pos.y[i], set->pos.z[i],BOX_SIZE_INV);
		omp_moveAtomBox(i,box_id,set);
		boxSort->swapState[i] = box_id;

	}
	
	#pragma omp single
	omp_switchPtr(set);
	
	#pragma omp barrier
}

int omp_get_num_box(const sotl_domain_t *dom, const int x, const int y, const int z)
{
    if(x  <0 || y <0 || z <0 )
        return -1;

    int box_x, box_y, box_z;
    int box_id;
    int rrc  = BOX_SIZE_INV;

    box_x = x * rrc;
    box_y = y * rrc;
    box_z = z * rrc;

    box_id =  box_z * dom->boxes[0] * dom->boxes[1] +
              box_y * dom->boxes[0] +
              box_x;

    if(box_id < 0)
        return -1;
    if((unsigned) box_id >= dom->total_boxes)
		return -1;
		
    return box_id;
}

void omp_get_pos_box(const sotl_domain_t *dom, int nbBox, unsigned* pos)
{
    if(nbBox < 0 || pos == NULL){
        return;
    }

    int oY = dom->boxes[0];
    int oZ = oY * dom->boxes[1];

    int box_x, box_y, box_z;
    int rrc  = BOX_SIZE_INV;

    box_z =  nbBox/oZ;
    box_x =  nbBox%oZ;
    box_y =  box_x/oY;
    box_x =  box_x%oY;

    pos[0] = box_x/rrc;
    pos[1] = box_y/rrc;
    pos[2] = box_z/rrc;
}


static void omp_computeForce(sotl_atom_set_t *set, int currentBox, int otherBox){

    calc_t sq_dist;
    calc_t force[3];

    for (int current = boxSort->preNbAtomToBox[currentBox];
    current < boxSort->preNbAtomToBox[currentBox+1];
    current++)
    {

        force[0] = 0.0;
        force[1] = 0.0;
        force[2] = 0.0;

        for (int other = boxSort->preNbAtomToBox[otherBox];
        other < boxSort->preNbAtomToBox[otherBox+1];
        other++)
        {

            if (current != other ) //&& current != NULL && other != NULL)//current other pas négatif car unsigned
            {
                assert(current >= 0);
                assert(other >= 0);
                sq_dist = squared_distance (set, current, other);

                if (sq_dist < LENNARD_SQUARED_CUTOFF)
                {
                    calc_t intensity = lennard_jones (sq_dist);

                    force[0] += intensity * (set->pos.x[current] - set->pos.x[other]);
                    force[1] += intensity * (set->pos.x[set->offset + current] - set->pos.x[set->offset + other]);
                    force[2] += intensity * (set->pos.x[set->offset * 2 + current] - set->pos.x[set->offset * 2 + other]);
                }
            }
        }
        set->speed.dx[current] += force[0];
        set->speed.dx[set->offset + current] += force[1];
        set->speed.dx[set->offset * 2 + current] += force[2];
  }

}

/*
void checkGetNum(sotl_domain_t *dom, int currentBox){
    // CHECK
    // printf("CHECK BEGIN\n");
    unsigned pos[3] ={0,0,0};
    get_pos_box(dom, currentBox, pos);
    if(x!=pos[0])   printf("fuck X %d,%d\n",x,pos[1]);
    if(y!=pos[1])   printf("fuck Y\n");
    if(z!=pos[2])   printf("fuck Z\n");
    //printf("CHECK END\n");
    // END CHECK

}*/

static void omp_seq_force_cube (sotl_device_t *dev)
{
	sotl_atom_set_t *set = &dev->atom_set;
	sotl_domain_t *dom = &dev->domain;
	//printf("omp_seq_force_cube thread %d\n",omp_get_thread_num());
	
  omp_sortAtomBox(dom,set);
	

    int otherBox;

		#pragma omp for private(otherBox)
    for(unsigned currentBox = 0; currentBox < dom->total_boxes; currentBox++)
    {
				unsigned pos[3] = {0,0,0};
        omp_get_pos_box(dom, currentBox, pos);


        if(boxSort->nbAtomToBox[currentBox] >0)
        {
            for(int pos_s = 0; pos_s < CUBE_SENTCIL_SIZE; pos_s+=3)
            {
                int xother = pos[0] +cube_stencil[pos_s];
                int yother = pos[1] +cube_stencil[pos_s+1];;
                int zother = pos[2] +cube_stencil[pos_s+2];

                otherBox = omp_get_num_box(dom,xother,yother,zother);

                if(otherBox > -1 && boxSort->nbAtomToBox[otherBox] > 0)
                {
                    omp_computeForce(set, currentBox, otherBox);
                }
            }
        }
    }

	
}



static void omp_force (sotl_device_t *dev)
{
  //sotl_atom_set_t *set = &dev->atom_set;
  
  //#pragma omp single
  //boxForce(dev);
  
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
/*
  omp_forceOrg(dev);

  //*/
  
  
  omp_seq_force_cube (dev);
}


// Main simulation function
//
void omp_one_step_move (sotl_device_t *dev)
{
	//omp_set_num_threads(16);
	
	#pragma omp parallel //default(shared)
	{
		//printf("thread %d\n",omp_get_thread_num());
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

void omp_finalize (sotl_device_t *dev)
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
