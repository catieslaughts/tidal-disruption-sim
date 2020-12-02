#ifndef __BlackholeCenteredDriver_h__
#define __BlackholeCenteredDriver_h__
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

template<int d> class BlackholeCenteredDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.025;

	// set draw_axes = true in ../viewer/src/OpenGLViewer.h, on line 30

	int n;
	Particles<d> particles;
	Array<Point> points;
	OpenGLSegmentMesh* opengl_trace=nullptr;

	OpenGLSphere* blackhole;
	real blackhole_mass = (real)300;

	// real g = (real)1.3271244e20;
	// real g = (real) 1e10;
	real g = 0.4;

	/*******INITIAL CONDITIONS*********/

	real color_k = 40;
	real color_xo = 0.0136;
	Vector3 vo = Vector3(-3.,.0,7.);
	Vector3 xo = Vector3::Unit(0) * 6;

	/*********END INITIAL CONDITIONS**********/
	//TODO try color_k = 85;

public:
	virtual void Initialize()
	{
		// fluid.Initialize();
		OpenGLViewer::Initialize();
	}

	//synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		//initialize sphere
		blackhole=Add_Interactive_Object<OpenGLSphere>();
		blackhole->pos=VectorD::Zero();
		blackhole->radius=(real).1;
		Set_Color(blackhole,OpenGLColor(.2,.2,.2,.1));
		blackhole->Set_Data_Refreshed();
		blackhole->Initialize();

		// read in location data for star points
		std::ifstream pointfile ("starpoints.txt");
		std::string line;
		Array<Vector3> newpoints;
		if(pointfile.is_open()){
			while(getline(pointfile,line)){
				Array<real> newline;
				std::stringstream ss(line);
				real i;
				while (ss >> i){
					newline.push_back(i);
					if (ss.peek() == ','){
						ss.ignore();
					}
				}
				if(newline.size() == 3){
					newpoints.push_back(Vector3(newline[0], newline[1], newline[2]));
				}
			}
		}

		//initialize star points
		n = newpoints.size();
		real rad = 0.1;
		particles.Resize(n);
		points.resize(n);
		for(int i=0;i<n;i++){
			particles.M(i) = (real)1 / (real)n;
			particles.V(i) = vo;
			particles.X(i) = xo + newpoints[i] * rad;

			// sync particles with display points
			points[i].Set_Radius(.02);
			points[i].Initialize(this);
			points[i].Sync_Data(particles.X(i));
			points[i].Set_Color(1,.5,0);
		}

		//set lights
		auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
		OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
		OpenGLUbos::Update_Lights_Ubo();

		////initialize a segment mesh to visualize the trace
		opengl_trace=Add_Interactive_Object<OpenGLSegmentMesh>();
		opengl_trace->mesh.Elements().resize(1);
		opengl_trace->mesh.Vertices().resize(1);
		opengl_trace->mesh.Vertices()[0]=particles.X(n-1);
		opengl_trace->mesh.elements[0]=Vector2i(0,0);
		opengl_trace->Set_Data_Refreshed();
		opengl_trace->Initialize();
	}

	Vector3 Color_Map(real val){
		real f = 1/(1+exp(-color_k*(val - color_xo))); //logistic
		int r, g, b;
		r = g = b = 0;

		//rainbow
		real a = (1 - f) / 0.25;
		int x = (int)floor(a);
		int y = (int)floor(255*(a - x));
		switch(x){
			case 0:{
				r = 255;g = y;b = 0;
			}break;
			case 1:{
				r = 255 - y;g = 255;b = 0;
			}break;
			case 2:{
				r = 0;g = 255;b = y;
			}break;
			case 3:{
				r = 0;g = 255 - y;b = 255;
			}break;
			case 4:{
				r = 0;g = 0;b = 255;
			}break;
		}

		//yellow to red
		// real a = 1 - f;
		// int y = floor(255*a);
		// r = 255;g = y;b = 0;

		return Vector3((real)r,(real)g,(real)b) / 255;
	}

	void Particle_Force_Accumulation(){
		for(int i=0;i<n;i++){
			particles.F(i) = Vector3::Zero();
		}

		for(int i=0; i<n;i++){
			real r_bh = (blackhole->pos - particles.X(i)).norm();
			Vector3 unit_bh = (blackhole->pos - particles.X(i)).normalized();
			Vector3 bh_f = (blackhole_mass*particles.M(i)/pow(r_bh,2))*unit_bh;
			particles.F(i) += bh_f;
			particles.C(i) = bh_f.norm();

			for(int j=i+1; j<n;j++){
				real r = (particles.X(j) - particles.X(i)).norm();
				if(r >= 0.21 && r < 10){
					Vector3 norm = (particles.X(j) - particles.X(i)).normalized();
					Vector3 grav_f = norm * g * particles.M(i) * particles.M(j) / pow(r,2);
					particles.F(i) += grav_f;
					particles.F(j) += -grav_f;
				}
			}
		}
	}

	virtual void Advance(const real dt){
		Particle_Force_Accumulation();
		for(int i = 0; i < particles.Size(); i++){
			particles.V(i) += particles.F(i) * dt / particles.M(i);
			particles.X(i) += particles.V(i) * dt;
		}
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		//update and sync data for blackhole
		blackhole->Set_Data_Refreshed();

		//update and synch data for sun particles
		for(int i=0;i<n;i++){
			points[i].Sync_Data(particles.X(i));
			Vector3 c = Color_Map(particles.C(i));
			points[i].Set_Color(c[0], c[1], c[2]);
		}

		opengl_trace->mesh.Vertices().push_back(particles.X(n-1));
		int m=(int)opengl_trace->mesh.Elements().size();
		opengl_trace->mesh.Elements().push_back(Vector2i(m-1,m-2));
		opengl_trace->Set_Data_Refreshed();
	}

	//update simulation and visualization for each time step
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

	//Keyboard interaction
	virtual void Initialize_Common_Callback_Keys()
	{
		OpenGLViewer::Initialize_Common_Callback_Keys();
		Bind_Callback_Key('a',&Keyboard_Event_A_Func,"press A");
	}

	virtual void Keyboard_Event_A()
	{
		std::cout<<"A pressed"<<std::endl;
	}
	Define_Function_Object(BlackholeCenteredDriver,Keyboard_Event_A);

protected:
	//Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif
