
#include "pmd.h"
#include "quaternion.h"

#ifdef USE_VTK
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointSource.h>
#include <vtkProperty.h>
#include <vtkCommand.h>
#endif

#include <iostream>
#include <pthread.h>
#include <float.h>

#include <stdio.h>

#define PI 3.14159

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

#ifdef USE_VTK
vtkPoints* AtomVisualize::pts = 0;
vtkCellArray* AtomVisualize::conn = 0;
vtkPolyData* AtomVisualize::poly = 0;
vtkUnsignedCharArray* AtomVisualize::colors = 0;
vtkPointSource* AtomVisualize::pointSource = 0;
vtkPolyDataMapper* AtomVisualize::inputMapper = 0;
vtkActor* AtomVisualize::inputActor = 0;
vtkRenderer* AtomVisualize::ren1 = 0;
vtkRenderWindow* AtomVisualize::renWin = 0;
vtkRenderWindowInteractor* AtomVisualize::renWinInteractor = 0;
KeyPressInteractorStyle* AtomVisualize::style = 0;
vtkTimerCallback2* AtomVisualize::cb = 0;
#else
int AtomVisualize::mouse_old_x = 0; 
int AtomVisualize::mouse_old_y = 0;
int AtomVisualize::mouse_buttons = 0;
float AtomVisualize::zoomFactor = 1.0f;
int AtomVisualize::displayList = -1;
#endif

simflat_t* AtomVisualize::sim = 0;
bool AtomVisualize::newDataAvailable = false;
float AtomVisualize::minCentro = 0.0f;
float AtomVisualize::maxCentro = 0.0f;
float AtomVisualize::minX = 0.0f; 
float AtomVisualize::maxX = 0.0f; 
float AtomVisualize::minY = 0.0f; 
float AtomVisualize::maxY = 0.0f; 
float AtomVisualize::minZ = 0.0f; 
float AtomVisualize::maxZ = 0.0f;
Quaternion* AtomVisualize::q = 0;


#ifdef USE_VTK

class vtkTimerCallback2 : public vtkCommand
{
  public:
    static vtkTimerCallback2 *New()
    {
      vtkTimerCallback2 *cb = new vtkTimerCallback2;
      cb->TimerCount = 0;
      return cb;
    }
 
    virtual void Execute(vtkObject *caller, unsigned long eventId, void * vtkNotUsed(callData))
    {
      if (vtkCommand::TimerEvent == eventId)
      {
        ++this->TimerCount;
      }
      //std::cout << this->TimerCount << std::endl;
      vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(caller);
      viz->renderData();
      iren->Render(); 
    }
 
  private:
    int TimerCount;
  public:
    AtomVisualize* viz;
};

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static KeyPressInteractorStyle* New();
    vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);
    AtomVisualize* atomVisualizer;
 
    virtual void OnKeyPress() 
    {
      // Get the keypress
      vtkRenderWindowInteractor *rwi = this->Interactor;
      std::string key = rwi->GetKeySym();
 
      // Output the key that was pressed
      //std::cout << "Pressed " << key << std::endl;
 
      // Handle a "normal" key
      if (key == "u")
      {
        atomVisualizer->renderData(); 
        atomVisualizer->renWinInteractor->Render();
      }
      
      // Forward events
      vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
};
vtkStandardNewMacro(KeyPressInteractorStyle);

#else

void AtomVisualize::keyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'u' : case 'U' : 
      displayList = -1; 
      glutPostRedisplay();
      break;
  };
}

void AtomVisualize::idle()
{
}

void AtomVisualize::mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) mouse_buttons |= 1<<button;
  else if (state == GLUT_UP) mouse_buttons = 0;

  mouse_old_x = x;
  mouse_old_y = y;
  glutPostRedisplay();
}

void AtomVisualize::motion(int x, int y)
{
  float dx = x - mouse_old_x;
  float dy = y - mouse_old_y;

  if (mouse_buttons == 1)
  {
    Quaternion newRotX;
    newRotX.setEulerAngles(-0.2*dx*PI/180.0, 0.0, 0.0);
    q->mul(newRotX);

    Quaternion newRotY;
    newRotY.setEulerAngles(0.0, 0.0, -0.2*dy*PI/180.0);
    q->mul(newRotY);
  }
  else if (mouse_buttons == 4)
  {
    zoomFactor += dy/100.0f;
  }
 
  mouse_old_x = x;
  mouse_old_y = y;
  glutPostRedisplay();
}

#endif


AtomVisualize::AtomVisualize()
{
#ifdef USE_VTK
  pts = vtkPoints::New();
  conn = vtkCellArray::New();
  poly = vtkPolyData::New();
  colors = vtkUnsignedCharArray::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  inputMapper = vtkPolyDataMapper::New();
  inputActor = vtkActor::New();
  inputActor->GetProperty()->SetPointSize(3.0f);
 
  ren1 = vtkRenderer::New();
  ren1->SetBackground(0.1, 0.2, 0.4);
  ren1->AddActor(inputActor);

  renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren1);
  renWin->SetSize(1024, 1024);

  renWinInteractor = vtkRenderWindowInteractor::New();
  renWinInteractor->SetRenderWindow(renWin);

  style = KeyPressInteractorStyle::New();
  renWinInteractor->SetInteractorStyle(style);
  style->atomVisualizer = this;
  style->SetCurrentRenderer(ren1);

#ifdef VIZ_ANIMATION
  cb = vtkTimerCallback2::New();
  cb->viz = this;
  renWinInteractor->AddObserver(vtkCommand::TimerEvent, cb);
#endif

#endif
  q = new Quaternion();
  newDataAvailable = false;
}


void AtomVisualize::interact()
{
#ifdef USE_VTK
  renWinInteractor->Initialize(); 
  int timerId = renWinInteractor->CreateRepeatingTimer(100);
  renWinInteractor->Start();
#else
  glutDisplayFunc(renderData);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutIdleFunc(idle);
  glutMainLoop();
#endif
}


void AtomVisualize::initialize(simflat_t *atoms)
{
  sim = NULL;
  sim = blankSimulation(atoms->pot);
  for (int i=0; i<3; i++) sim->bounds[i] = atoms->bounds[i];
  for (int i=0; i<3; i++) sim->boxsize[i] = atoms->boxsize[i];
  allocDomains(sim);

#ifndef USE_VTK
  int argc = 0; char **argv = 0; glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(1024, 1024);
  glutCreateWindow("cruft"); 
#endif 
}


void AtomVisualize::updateData(simflat_t *atoms) 
{
  int result = pthread_mutex_trylock( &mutex1 );
  if (result != 0) return;

  copyDomains(sim, atoms);

  pthread_mutex_unlock( &mutex1 );
  newDataAvailable = true;

#ifdef VIZ_ANIMATION
#ifndef USE_VTK
  displayList = -1;
  glutPostRedisplay();
#endif
#endif
}


void AtomVisualize::renderData()
{
    pthread_mutex_lock( &mutex1 );
    //if (n > 50000) n = 50000;

#ifdef USE_VTK 
    int prevPoints = pts->GetNumberOfPoints();
    vtkIdType acnt = 0;
#else
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    float cameraFOV = 45.0f;
    float zNear = 1.0f;
    float zFar = 1000.0f;
    float aspectRatio = 1.0f;
    float z1 = 1.5*((maxX-minX)/2.0)/tan(PI*cameraFOV/(180.0f*2.0f));
    float z2 = 1.5*((maxY-minY)/2.0)/tan(PI*cameraFOV/(180.0f*2.0f));
    float z3 = 1.5*maxZ;
    float cameraZ = zoomFactor*std::max(z1, std::max(z2, z3));

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(cameraFOV, aspectRatio, zNear, zFar);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0f, 0.0f, cameraZ, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    glPushMatrix();

    glTranslatef(-(maxX-minX)/2.0f, -(maxY-minY)/2.0f, -(maxZ-minZ)/2.0f);

    float rotationMatrix[16];
    q->getRotMat(rotationMatrix);
    glMultMatrixf(rotationMatrix);

    float centerX = (maxX-minX)/2.0f;  float centerY = (maxY-minY)/2.0f;  float centerZ = (maxZ-minZ)/2.0f;
    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    float offsetX = matrix[0]*centerX + matrix[1]*centerY + matrix[2]*centerZ; 
    float offsetY = matrix[4]*centerX + matrix[5]*centerY + matrix[6]*centerZ;
    float offsetZ = matrix[8]*centerX + matrix[9]*centerY + matrix[10]*centerZ;
    offsetX = centerX - offsetX; offsetY = centerY - offsetY; offsetZ = centerZ - offsetZ;
    glTranslatef(-offsetX, -offsetY, -offsetZ);

    bool creatingDisplayList = false;
    if (displayList == -1)
    {
      displayList = glGenLists(1);      
      glNewList(displayList, GL_COMPILE_AND_EXECUTE);
      creatingDisplayList = true;
      minX = minY = minZ = FLT_MAX;
      maxX = maxY = maxZ = FLT_MIN;
      glPointSize(3.0);
      glBegin(GL_POINTS);
    }
    else
    {
      glCallList(displayList);
      glPopMatrix();
      glutSwapBuffers();
      pthread_mutex_unlock( &mutex1 );
      return;
    }
#endif

    if (newDataAvailable)
    {
      computeCentrosymmetry(sim);
      newDataAvailable = false;
    }

    for(int iBox=0; iBox<sim->nboxes; iBox++) 
    { 
      int ii;
      int ioff;
      if ( ! sim->natoms[iBox]) continue;
      for(ioff=MAXATOMS*iBox,ii=0; ii<sim->natoms[iBox]; ii++,ioff++) 
      {
        unsigned char color[3];  colorMap(minCentro, maxCentro, sim->centro[ioff], color); 

#ifdef USE_VTK
        if (acnt < prevPoints)
        {
          pts->SetPoint(acnt, sim->dcenter[iBox][0] + sim->r[ioff][0], sim->dcenter[iBox][1] + sim->r[ioff][1], sim->dcenter[iBox][2] + sim->r[ioff][2]);         
          colors->SetTupleValue(acnt, color);
        }
        else
        {
          pts->InsertPoint(acnt, sim->dcenter[iBox][0] + sim->r[ioff][0], sim->dcenter[iBox][1] + sim->r[ioff][1], sim->dcenter[iBox][2] + sim->r[ioff][2]);   
          colors->InsertNextTupleValue(color);
          conn->InsertNextCell(1, &acnt);   
        }   
        acnt++;
#else
        glColor3ub(color[0], color[1], color[2]);
        float x = sim->dcenter[iBox][0] + sim->r[ioff][0];  if (x > maxX) maxX = x;  if (x < minX) minX = x;
        float y = sim->dcenter[iBox][1] + sim->r[ioff][1];  if (y > maxY) maxY = y;  if (y < minY) minY = y;
        float z = sim->dcenter[iBox][2] + sim->r[ioff][2];  if (z > maxZ) maxZ = z;  if (z < minZ) minZ = z;
        glVertex3f(x, y, z);
#endif
      }
    }

#ifdef USE_VTK
    poly->SetPoints(pts);  poly->SetVerts(conn);  poly->GetPointData()->SetScalars(colors);
    poly->Modified();

    inputMapper->SetInput(poly);
    inputActor->SetMapper(inputMapper);

    poly->Update();
    inputMapper->Update();
#else
    glEnd();
    if (creatingDisplayList) glEndList();
    glPopMatrix();
    glutSwapBuffers();
#endif

    pthread_mutex_unlock( &mutex1 );
}


void AtomVisualize::computeCentrosymmetry(simflat_t *atoms)
{
  int             i, j, i1, i2, k, kmax;
  real_t          dx, dy, dz, rsq, rmax;
 
  std::cout << "Computing centrosymmetry" << std::endl;
  std::cout << "Stats: " << atoms->nboxes << " " << atoms->nbx[0] << " " << atoms->nbx[1] << " " << atoms->nbx[2] << std::endl;
  int cntr = 0;

  // Within-box nearest neighbors
  for(int iBox=0; iBox<atoms->nboxes; iBox++) 
  { 
    if (iBox % 100000 == 0) printf("Computing centrosymmetry for box %d of %d, %d elements\n", iBox, atoms->nboxes, atoms->natoms[iBox]);
    int ii, jj;
    int ioff, joff;
    if ( ! atoms->natoms[iBox]) continue;

    for(ioff=MAXATOMS*iBox,ii=0; ii<atoms->natoms[iBox]; ii++,ioff++) 
    {
      Neighbor12 nlist;
      int Nbors = 0;
      real_t centro = 0.0;

      for(joff=MAXATOMS*(iBox),jj=0; jj<atoms->natoms[iBox]; jj++,joff++)
      {
        if (ioff == joff) continue;

        dx = atoms->r[joff][0] - atoms->r[ioff][0];
        dy = atoms->r[joff][1] - atoms->r[ioff][1];
        dz = atoms->r[joff][2] - atoms->r[ioff][2];
        rsq = dx * dx + dy * dy + dz * dz;

        // Determine where to store pt2 in pt1's neighbor list 
        if (Nbors < 12) i1 = Nbors++;
        else {  // already have 12 neighbors, need to check distances 
	  kmax = 0;
	  rmax = nlist[0].rsq;
	  for (k = 1; k < 12; k++)
	    if (nlist[k].rsq > rmax) {
	      kmax = k;
	      rmax = nlist[k].rsq;
	    }
	  if (rsq > rmax) i1 = 100;   // ie. do not store 
	  else i1 = kmax;
        }
        if (i1 < 12) {
	  nlist[i1].rsq = rsq;
	  nlist[i1].r[0] = dx;
	  nlist[i1].r[1] = dy;
	  nlist[i1].r[2] = dz;
          nlist[i1].index = joff;
        }
      }


      // Nearest neighbors in adjacent cells
      int bx, by, bz, nx, ny, nz, nxc, nyc, nzc, b1, b2;
  
      bz = iBox / (atoms->nbx[0]*atoms->nbx[1]);
      by = (iBox - bz*atoms->nbx[0]*atoms->nbx[1]) / (atoms->nbx[0]);
      bx = iBox - bz*atoms->nbx[0]*atoms->nbx[1] - by*atoms->nbx[0];

      for (nxc=bx-1; nxc<=bx+1; nxc++)
      {
        for (nyc=by-1; nyc<=by+1; nyc++)
        {
          for (nzc=bz-1; nzc<=bz+1; nzc++)
          {
            if ((nxc == bx) && (nyc == by) && (nzc == bz)) continue;
            //if ((nx < 0) || (ny < 0) || (nz < 0)) continue;
            //if ((nx >= atoms->nbx[0]) || (ny >= atoms->nbx[1]) || (nz >= atoms->nbx[2])) continue;

            float xoffset = 0.0;  float yoffset = 0.0;  float zoffset = 0.0;
            nx = nxc; ny = nyc; nz = nzc;
            if (nx < 0) { nx = atoms->nbx[0]-1; xoffset = -atoms->bounds[0]; }  
            if (nx >= atoms->nbx[0]) { nx = 0; xoffset = atoms->bounds[0]; }
            if (ny < 0) { ny = atoms->nbx[1]-1; yoffset = -atoms->bounds[1]; }
            if (ny >= atoms->nbx[1]) { ny = 0; yoffset = atoms->bounds[1]; }
            if (nz < 0) { nz = atoms->nbx[2]-1; zoffset = -atoms->bounds[2]; }  
            if (nz >= atoms->nbx[2]) { nz = 0; zoffset = atoms->bounds[2]; }
              
            int b1 = bz*atoms->nbx[0]*atoms->nbx[1] + by*atoms->nbx[0] + bx;
            int b2 = nz*atoms->nbx[0]*atoms->nbx[1] + ny*atoms->nbx[0] + nx;
        
            for (joff=MAXATOMS*b2,jj=0; jj<atoms->natoms[b2]; jj++,joff++)
            {
              dx = (xoffset + atoms->dcenter[b2][0] + atoms->r[joff][0]) - (atoms->dcenter[b1][0] + atoms->r[ioff][0]);
              dy = (yoffset + atoms->dcenter[b2][1] + atoms->r[joff][1]) - (atoms->dcenter[b1][1] + atoms->r[ioff][1]);
              dz = (zoffset + atoms->dcenter[b2][2] + atoms->r[joff][2]) - (atoms->dcenter[b1][2] + atoms->r[ioff][2]);
              rsq = dx * dx + dy * dy + dz * dz;

              // Determine where to store pt2 in pt1's neighbor list 
              if (Nbors < 12) i1 = Nbors++;
              else   // already have 12 neighbors, need to check distances 
              {
	        kmax = 0;
	        rmax = nlist[0].rsq;
	        for (k = 1; k < 12; k++)
                {
	          if (nlist[k].rsq > rmax) 
                  {
	            kmax = k;
	            rmax = nlist[k].rsq;
	          }
                }

	        if (rsq > rmax) i1 = 100;   // ie. do not store 
	        else i1 = kmax;
              }
              if (i1 < 12) 
              {
	        nlist[i1].rsq = rsq;
	        nlist[i1].r[0] = dx;
	        nlist[i1].r[1] = dy;
	        nlist[i1].r[2] = dz;
                nlist[i1].index = joff;
              }           
            }
          }
        }
      }

      
      // Compute centrosymmetry from neighbors
      real_t sum, summin, cmin = 1000.0, cmax = -1.0;
      int kmin;
      if (Nbors < 12) atoms->centro[ioff] = 1000.0;  // surface atoms, etc.
      else
      {
#ifdef APPROX_CENTROSYM
        // "Lazy" but faster centrosymmetry computation
        atoms->centro[ioff] = 0.0;
        for (j = 0; j < 11; j++)
        {
	  if (nlist[j].rsq > 0.0)   // j not paired up yet 
          {
	    // look for neighbor k most nearly opposite from j 
	    summin = 1000.0;
	    for (k = j+1; k < 12; k++)
            {
	      if (nlist[k].rsq > 0.0)  // k not paired up yet 
              {
	        sum = (nlist[j].r[0] + nlist[k].r[0])*
	    	      (nlist[j].r[0] + nlist[k].r[0]) +
		      (nlist[j].r[1] + nlist[k].r[1]) *
		      (nlist[j].r[1] + nlist[k].r[1]) +
		      (nlist[j].r[2] + nlist[k].r[2]) *
		      (nlist[j].r[2] + nlist[k].r[2]);
	        if (sum < summin) {kmin = k; summin = sum;}
	      }
            }

	    // Paired up atoms j and kmin, add to centrosym parameter and mark 
	    atoms->centro[ioff] += summin;
	    nlist[j].rsq = nlist[kmin].rsq = -1.0;
	  }
        }
#else
        // Slow but more accurate centrosymmetry computation
        atoms->centro[ioff] = 1000.0;     
        int plist[12];
        for (int f=0; f<12; f++) plist[f] = -1;
        for (int p0=0; p0<11; p0++)
        {
          int n = 0; int p = p0 + 1;
          plist[n] = p;  plist[p] = n;
          for (int p1=0; p1<9; p1++)
          {
            int n = 0; while (plist[n] >= 0) n++;
            int p = n+1;  int c = 0;
            while (c <= p1)
            {
              if (plist[p] == -1) c++;
              if (c <= p1) p++;
            }
            plist[n] = p;  plist[p] = n;
            for (int p2=0; p2<7; p2++)
            {
              int n = 0; while (plist[n] >= 0) n++;
              int p = n+1;  int c = 0;
              while (c <= p2)
              {
                if (plist[p] == -1) c++;
                if (c <= p2) p++;
              }
              plist[n] = p;  plist[p] = n;
              for (int p3=0; p3<5; p3++)
              {
                int n = 0; while (plist[n] >= 0) n++;
                int p = n+1;  int c = 0;
                while (c <= p3)
                {
                  if (plist[p] == -1) c++;
                  if (c <= p3) p++;
                }
                plist[n] = p;  plist[p] = n;
                for (int p4=0; p4<3; p4++)
                {
                  int n = 0; while (plist[n] >= 0) n++;
                  int p = n+1;  int c = 0;
                  while (c <= p4)
                  {
                    if (plist[p] == -1) c++;
                    if (c <= p4) p++;
                  }
                  plist[n] = p;  plist[p] = n;
                  for (int p5=0; p5<1; p5++)
                  {
                    int n = 0; while (plist[n] >= 0) n++;
                    int p = n+1;  int c = 0;
                    while (c <= p5)
                    {
                      if (plist[p] == -1) c++;
                      if (c <= p5) p++;
                    }
                    plist[n] = p;  plist[p] = n;

                    float testCentro = 0.0;
                    int tlist[12];
                    for (int t=0; t<12; t++) tlist[t] = plist[t]; 
                    for (int pi=0; pi<12; pi++)
                    {
                      if (tlist[pi] >= 0)
                      {
                        testCentro += (nlist[pi].r[0] + nlist[plist[pi]].r[0])*
	      	                      (nlist[pi].r[0] + nlist[plist[pi]].r[0]) +
	   	                      (nlist[pi].r[1] + nlist[plist[pi]].r[1]) *
		                      (nlist[pi].r[1] + nlist[plist[pi]].r[1]) +
		                      (nlist[pi].r[2] + nlist[plist[pi]].r[2]) *
		                      (nlist[pi].r[2] + nlist[plist[pi]].r[2]);
                        tlist[plist[pi]] = -1;  tlist[pi] = -1; 
                      } 
                    }
                    if (testCentro < atoms->centro[ioff]) atoms->centro[ioff] = testCentro;
                    plist[n] = -1;  plist[p] = -1;
                  }
                  plist[n] = -1;  plist[p] = -1;
                }
                plist[n] = -1;  plist[p] = -1;
              }
              plist[n] = -1;  plist[p] = -1;
            }
            plist[n] = -1;  plist[p] = -1;
          }
          plist[n] = -1;  plist[p] = -1;
        }
#endif
        //if (cntr < 5)
          //printf("Value %d: %f %f %f %f\n", ioff, atoms->dcenter[iBox][0] + atoms->r[ioff][0], atoms->dcenter[iBox][1] + atoms->r[ioff][1], 
            //                                      atoms->dcenter[iBox][2] + atoms->r[ioff][2], atoms->centro[ioff]); 
        //cntr++;
      }  
    }
  }

  // Find range
  minCentro = 1000.0;  maxCentro = 0.0;
  for(int iBox=0; iBox<atoms->nboxes; iBox++) 
  { 
    int ii, ioff;
    if ( ! atoms->natoms[iBox]) continue;

    for(ioff=MAXATOMS*iBox,ii=0; ii<atoms->natoms[iBox]; ii++,ioff++)
    {
      if (atoms->centro[ioff] < minCentro) minCentro = atoms->centro[ioff];
      if (atoms->centro[ioff] > maxCentro) maxCentro = atoms->centro[ioff];
    }
  }
  std::cout << "Range: " << minCentro << " " << maxCentro << std::endl;
  if (maxCentro > 100.0) maxCentro = 100.0;

  printf("Finished computing centrosymmetry\n");
}


void AtomVisualize::colorMap(float min, float max, float value, unsigned char* color)
{
  // HSV rainbow for height field, stolen form Manta
  const float V = 0.7f, S = 1.0f;
  float H = (1.0f - (value - min) / (max - min));

  if (H < 0.0f) H = 0.0f;
  else if (H > 1.0f) H = 1.0f;
  H *= 4.0f;

  float i = floor(H);
  float f = H - i;

  float p = V * (1.0 - S);
  float q = V * (1.0 - S * f);
  float t = V * (1.0 - S * (1 - f));

  float R, G, B;
  if (i == 0.0) { R = V; G = t; B = p; } 
  else if (i == 1.0) { R = q; G = V; B = p; } 
  else if (i == 2.0) { R = p; G = V; B = t; } 
  else if (i == 3.0) { R = p; G = q; B = V; } 
  else if (i == 4.0) { R = t; G = p; B = V; } 
  else { R = V; G = p; B = q; }
 
  color[0] = 255*R; color[1] = 255*G; color[2] = 255*B;
  //std::cout << min << " " << max << " " << value << " " << R << " " << G << " " << B << std::endl;
};

