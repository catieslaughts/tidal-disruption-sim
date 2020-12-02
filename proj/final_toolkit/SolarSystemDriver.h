#ifndef __SolarSystemDriver_h__
#define __SolarSystemDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "Particles.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLParticles.h"
#include "math.h"
#include "InClassDemoDriver.h"
#include <fstream>
#include <sstream>

#define PI 3.14159265
#define meters_per_AU 1.496e11 //meters per AU
#define km_per_AU 1.496e8 //km per AU
#define Ro_per_AU 215.032 //solar radii per AU
#define kg_per_Mo 1.9891e30 //kg per solar mass
#define s_per_year 3.154e7

template<int d> class SolarSystemDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.005;

	// set draw_axes = true in ../viewer/src/OpenGLViewer.h, on line 30

	Particles<d> particles;
	Array<Point> points;
	OpenGLSegmentMesh* opengl_trace[8]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

	OpenGLSphere* sun;
	real sun_mass = (real)1;//in solar masses
	real g = (real)1.3271244e20*pow(s_per_year,2)/pow(meters_per_AU,3);//in AU^3/year^2/solar mass
	int n = 8;
	real sun_size_multiplyer =100.; //changes radii to make viewing easier
	real planet_size_multiplyer=100.;



	//in solar masses:
	real masses[8] = {3.3e23/kg_per_Mo,4.8e24/kg_per_Mo,5.9e24/kg_per_Mo,6.4e23/kg_per_Mo,1.8e27/kg_per_Mo,5.6e26/kg_per_Mo,8.6e25/kg_per_Mo,10.2e25/kg_per_Mo};
	//in AU:
	real distances[8] = {0.39, 0.723, 1, 1.524, 5.203, 9.539, 19.18, 30.06};
	//in AU/year:
	real velocities[8] = {47.3*s_per_year/km_per_AU, 35.0*s_per_year/km_per_AU, 29.8*s_per_year/km_per_AU, 24.1*s_per_year/km_per_AU, 13.1*s_per_year/km_per_AU, 9.7*s_per_year/km_per_AU, 6.8*s_per_year/km_per_AU, 5.4*s_per_year/km_per_AU};
	//in AU:
	real radii[8] = {4878/km_per_AU, 12104/km_per_AU, 12756/km_per_AU, 6794/km_per_AU, 142984/km_per_AU, 120536/km_per_AU, 51118/km_per_AU, 49532/km_per_AU};
	Array<OpenGLSphere*> suns;								////spheres

public:
	virtual void Initialize()
	{
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		////initialize sphere
		sun=Add_Interactive_Object<OpenGLSphere>();
		sun->pos=Vector3(0,0,0);
		sun->radius=(real)sun_size_multiplyer/Ro_per_AU; //in AU
		Set_Color(sun,OpenGLColor(1.,.75,.0,.5));
		sun->Set_Data_Refreshed();
		sun->Initialize();

		////initialize planet points
		real rad = 0.1;
		int ind = 0;
		particles.Resize(n);
		points.resize(n);

		for(int i=0;i<n;i++){
			particles.V(i) = Vector3(.0,.0,velocities[i]);
			particles.M(i) = masses[i];
			particles.X(i) = Vector3(distances[i],0,0);
			points[i].Set_Radius(radii[i]*planet_size_multiplyer);
			points[i].Initialize(this);
			points[i].Sync_Data(particles.X(i));
			points[i].Set_Color(.0,.0,.75);
		}

		////set lights
		auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
		OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
		OpenGLUbos::Update_Lights_Ubo();

		////initialize a segment mesh to visualize the trace
		for(int i=0; i<n; i++){
			opengl_trace[i]=Add_Interactive_Object<OpenGLSegmentMesh>();
			opengl_trace[i]->mesh.Elements().resize(1);
			opengl_trace[i]->mesh.Vertices().resize(1);
			opengl_trace[i]->mesh.Vertices()[0]=particles.X(i);
			opengl_trace[i]->mesh.elements[0]=Vector2i(0,0);
			opengl_trace[i]->Set_Data_Refreshed();
			opengl_trace[i]->Initialize();
		}
	}

	void Particle_Force_Accumulation(){
		for(int i=0;i<n;i++){
			particles.F(i) = Vector3::Zero();
		}

		for(int i=0; i<n;i++){
			real r_sun = (sun->pos - particles.X(i)).norm();
			Vector3 unit_sun = (sun->pos - particles.X(i))/r_sun;

			particles.F(i) += (g*sun_mass*particles.M(i)/pow(r_sun,2))*unit_sun;

			for(int j=i+1; j<n;j++){
				real r = (particles.X(j) - particles.X(i)).norm();
				Vector3 norm = (particles.X(j) - particles.X(i)).normalized();

				Vector3 grav_f = norm * g * particles.M(i) * particles.M(j) / pow(r,2);
				particles.F(i) += grav_f;
				particles.F(j) += -grav_f;
			}
		}
	}

	virtual void Advance(const real dt){
		Particle_Force_Accumulation();
		int tot = 0;
		for(int i = 0; i < particles.Size(); i++){
			particles.V(i) += particles.F(i) * dt / particles.M(i);
			particles.X(i) += particles.V(i) * dt;

			if (particles.X(i)[1]<=-5){
				particles.V(i)[1]=0;
			}
		}
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		////update and sync data for sun
		sun->Set_Data_Refreshed();

		//std::cout << opengl_trace[2]->mesh.Vertices();

		for(int i=0;i<n;i++){
			points[i].Sync_Data(particles.X(i));
			//points[i].Set_Color(1,.5,0);

			opengl_trace[i]->mesh.Vertices().push_back(particles.X(i));
			int m=(int)opengl_trace[i]->mesh.Elements().size();
			opengl_trace[i]->mesh.Elements().push_back(Vector2i(m-1,m-2));
			opengl_trace[i]->Set_Data_Refreshed();
		}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		// fluid.Advance(dt);
		Advance(dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	////Keyboard interaction
	virtual void Initialize_Common_Callback_Keys()
	{
		OpenGLViewer::Initialize_Common_Callback_Keys();
		Bind_Callback_Key('a',&Keyboard_Event_A_Func,"press A");
	}

	virtual void Keyboard_Event_A()
	{
		std::cout<<"A pressed"<<std::endl;
	}
	Define_Function_Object(SolarSystemDriver,Keyboard_Event_A);

protected:
	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif
