/*
  CSC D18 - Path Tracer code.

  Derived from the ray tracer starter code. Most function 
  names are identical, though in practice the implementation
  should be much simpler!

  You only need to modify or add code in sections
  clearly marked "TO DO" - remember to check what
  functionality is actually needed for the corresponding
  assignment!

  Last updated: Aug. 2017   - F.J.E.
*/

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:		
* 2) Student Name:		
*
* 1) Student number:
* 2) Student number:
* 
* 1) UtorID
* 2) UtorID
* 
* We hereby certify that the work contained here is our own
*
* ____________________             _____________________
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils_path.h"			// <-- This includes PathTracer.h
#ifdef _OPENMP
#include <omp.h>
#endif
//#define __USE_IS			// Use importance sampling for diffuse materials
//#define __USE_ES			// Use explicit light sampling
//#define __DEBUG			// <-- Use this to turn on/off debugging output

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct textureNode *texture_list;
unsigned long int NUM_RAYS;
int MAX_DEPTH;

#include "buildScene.c"			// Import scene definition

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // Inputs:
 //   *ray    -  A pointer to the ray being traced
 //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
 //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
 //              projection
 // Outputs:
 //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
 //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
 //              this ray (this is required so you can do the shading)
 //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
 //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
 //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

 // Os is intentionally unused: self-intersections are rejected by the lambda > EPSILON check instead.
 // Skipping Os entirely would break cases where a ray must hit its source object again (e.g. refraction
 // exiting a sphere, self-shadowing, or non-convex objects).
 struct object3D *o;
 struct point3D tp, tn;
 double tl, ta, tb;

 *lambda=-1;
 *obj=NULL;
 for (o=object_list; o!=NULL; o=o->next)
 {
  o->intersect(o,ray,&tl,&tp,&tn,&ta,&tb);
  if (tl>EPSILON && (*lambda<0 || tl<*lambda))
  {
   *lambda=tl;
   *obj=o;
   *p=tp;
   *n=tn;
   *a=ta;
   *b=tb;
  }
 }
}

void PathTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os, int CEL)
{
 // Trace one light path through the scene.
 //
 // Parameters:
 //   *ray   -  A pointer to the ray being traced
 //   depth  -  Current recursion depth for recursive raytracing
 //   *col   - Pointer to an RGB colour structure so you can return the object colour
 //            at the intersection point of this ray with the closest scene object.
 //   *Os    - 'Object source' is a pointer to the object from which the ray 
 //            originates so you can discard self-intersections due to numerical
 //            errors. NULL for rays originating from the center of projection. 
 
 double lambda;			// Lambda at intersection
 double a,b;			// Texture coordinates
 struct object3D *obj;		// Pointer to object at intersection
 struct point3D p;		// Intersection point
 struct point3D n;		// Normal at intersection
 double R,G,B;			// Handy in case you need to keep track of some RGB colour value
 double dice;			// Handy to keep a random value
 struct ray3D *next_ray;	// For the new ray to be used in recursive calls
 
 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return black (no light coming into pixel from this path).
 {
  col->R=ray->Ir;	// These are accumulators, initialized at 0. Whenever we find a source of light these
  col->G=ray->Ig;	// get incremented accordingly. At the end of the recursion, we return whatever light
  col->B=ray->Ib;	// we accumulated into these three values.
  return;
 }

 findFirstHit(ray,&lambda,Os,&obj,&p,&n,&a,&b);
 if (obj==NULL)
 {
  col->R=ray->Ir;
  col->G=ray->Ig;
  col->B=ray->Ib;
  return;
 }

 if (obj->isLightSource)
 {
  if (CEL)
  {
   ray->Ir+=ray->R*obj->col.R;
   ray->Ig+=ray->G*obj->col.G;
   ray->Ib+=ray->B*obj->col.B;
  }
  col->R=ray->Ir;
  col->G=ray->Ig;
  col->B=ray->Ib;
  return;
 }

 // Make the normal face the incoming ray
 if (dot(&n,&ray->d)>0) scaleVector(-1,&n);

 // Orthonormal basis u, v perpendicular to n (axis must not be parallel to n)
 struct point3D axis={1,0,0,1};
 if (fabs(n.px)>.9)
 {
  axis.px=0;
  axis.py=1;
 }
 struct point3D *u=cross(&n,&axis);
 normalize(u);
 struct point3D *v=cross(&n,u);

 // Uniform hemisphere sample around n
 double costheta=rand01();
 double sintheta=sqrt(1-costheta*costheta);
 double phi=2*PI*rand01();
 struct point3D d;
 d.px=costheta*n.px+sintheta*(cos(phi)*u->px+sin(phi)*v->px);
 d.py=costheta*n.py+sintheta*(cos(phi)*u->py+sin(phi)*v->py);
 d.pz=costheta*n.pz+sintheta*(cos(phi)*u->pz+sin(phi)*v->pz);
 d.pw=1;
 free(u);
 free(v);

 next_ray=newRay(&p,&d);
 next_ray->R=ray->R*obj->col.R*2*costheta;
 next_ray->G=ray->G*obj->col.G*2*costheta;
 next_ray->B=ray->B*obj->col.B*2*costheta;
 next_ray->Ir=ray->Ir;
 next_ray->Ig=ray->Ig;
 next_ray->Ib=ray->Ib;
 PathTrace(next_ray,depth+1,col,obj,1);
 free(next_ray);
}

int main(int argc, char *argv[])
{
 // Main function for the path tracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;		// Will hold the final image
 struct view *cam;		// Camera and view for this scene
 int sx;			// Size of the  image
 int num_samples;		// Number of samples to use per pixel
 char output_name[1024];	// Name of the output file for the .ppm image file
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for pixels
 int i,j,k;			// Counters for pixel coordinates and samples
 double *rgbIm;			// Image is now double precision floating point since we
				// will be accumulating brightness differences with a 
				// wide dynamic range
 struct object3D *obj;		// Will need this to process lightsource weights
 double *wght;			// Holds weights for each pixel - to provide log response
 double pct,wt;
 
 time_t t1,t2;
 FILE *f;
				
 if (argc<5)
 {
  fprintf(stderr,"PathTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: PathTracer size rec_depth num_samples output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   num_samples = Number of samples per pixel\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 num_samples=atoi(argv[3]);
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 fprintf(stderr,"Number of samples = %d\n",num_samples);
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 texture_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 wght=(double *)calloc(sx*sx,sizeof(double));
 if (!im||!wght)
 {
  fprintf(stderr,"Unable to allocate memory for image\n");
  exit(0);
 }
 else rgbIm=(double *)im->rgbdata;
 for (i=0;i<sx*sx;i++) *(wght+i)=1.0;
 
 buildScene();		// Create a scene. 
 
 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center
 e.px=0;
 e.py=0;
 e.pz=-15;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin.
 g.px=0-e.px;
 g.py=0-e.py;
 g.pz=0-e.pz;
 g.pw=1;
 // In this case, the camera is looking along the world Z axis, so
 // vector w should end up being [0, 0, -1]

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=1;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list, texture_list);
  deleteImage(im);
  exit(0);
 }

 du=cam->wsize/sx;		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/sx;		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix:\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 // Update light source weights - will give you weights for each light source that add up to 1
 obj=object_list;
 pct=0;
 while (obj!=NULL)
 {
  if (obj->isLightSource)
   pct+=obj->LSweight;
  obj=obj->next;
 }
 obj=object_list;
 while (obj!=NULL)
 {
  if (obj->isLightSource)
  {
   obj->LSweight/=pct;
  }
  obj=obj->next;
 }
 fprintf(stderr,"\n");

 NUM_RAYS=0;

 t1=time(NULL);

#ifdef _OPENMP
#pragma omp parallel
 seedThreadRNG(omp_get_thread_num());
#else
 seedThreadRNG(0);
#endif

 fprintf(stderr,"Rendering pass... ");
 for (k=0; k<num_samples; k++)
 {
  fprintf(stderr,"%d/%d, ",k,num_samples);
#pragma omp parallel for schedule(dynamic,1) private(i,j,pc,wt,ray,col,d)
  for (j=0;j<sx;j++)		// For each of the pixels in the image
  {
   for (i=0;i<sx;i++)
   {
    // Random sample within the pixel's area
    pc.px=(cam->wl+((i+rand01())*du));
    pc.py=(cam->wt+((j+rand01())*dv));
    pc.pz=cam->f;
    pc.pw=1;

    // Convert image plane sample coordinates to world coordinates
    matVecMult(cam->C2W,&pc);

    // Now compute the ray direction
    memcpy(&d,&pc,sizeof(struct point3D));
    subVectors(&cam->e,&d);		// Direction is d=pc-e
    normalize(&d);

    // Create a ray and do the raytracing for this pixel.
    ray=newRay(&pc,&d);

    if (ray!=NULL)
    {
     wt=*(wght+i+(j*sx));
     PathTrace(ray,1, &col,NULL,1);
     (*(rgbIm+((i+(j*sx))*3)+0))+=col.R*pow(2,-log(wt));
     (*(rgbIm+((i+(j*sx))*3)+1))+=col.G*pow(2,-log(wt));
     (*(rgbIm+((i+(j*sx))*3)+2))+=col.B*pow(2,-log(wt));
     wt+=col.R;
     wt+=col.G;
     wt+=col.B;
     *(wght+i+(j*sx))=wt;
     free(ray);   
    }
   } // end for i
  } // end for j  
  if (k%100==0)  dataOutput(rgbIm,sx,&output_name[0]);  		// Update output image every 100 passes
 } // End for k 
 t2=time(NULL);

#pragma omp parallel reduction(+:NUM_RAYS)
 NUM_RAYS+=thread_rays;

 fprintf(stderr,"\nDone!\n");

 dataOutput(rgbIm,sx,&output_name[0]);
 
 fprintf(stderr,"Total number of rays created: %ld\n",NUM_RAYS);
 fprintf(stderr,"Rays per second: %f\n",(double)NUM_RAYS/(double)difftime(t2,t1));

 // Exit section. Clean up and return.
 cleanup(object_list,texture_list);			// Object and texture lists
 deleteImage(im);					// Rendered image
 free(cam);						// camera view
 free(wght);
 exit(0);
}

