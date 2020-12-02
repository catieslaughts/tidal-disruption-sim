//////////////////////////////////////////////////////////////////////////
////Dartmouth Physical Computing Starter Code
////http://www.dartmouth.edu/~boolzhu/cosc89.18.html
//////////////////////////////////////////////////////////////////////////

#ifndef __Particles_h__
#define __Particles_h__
#include "Common.h"

template<int d> class Particles
{using VectorD=Vector<real,d>;
public:
	////attributes
	ArrayPtr<VectorD> x;		////position
	ArrayPtr<VectorD> v;		////velocity
	ArrayPtr<VectorD> f;		////force
	ArrayPtr<real> m;			////mass
	ArrayPtr<real> c;			////color
	ArrayPtr<real> r;			////radius
	ArrayPtr<real> p;			////pressure
	ArrayPtr<real> den;			////density
	ArrayPtr<int> idx;			////index, for rigid body

	//////////////////////////////////////////////////////////////////////////
	////common functions
	Particles()
	{
		if(x==nullptr)x.reset(new Array<VectorD>());	
		if(v==nullptr)v.reset(new Array<VectorD>());	
		if(f==nullptr)f.reset(new Array<VectorD>());	
		if(m==nullptr)m.reset(new Array<real>());	
		if(c==nullptr)c.reset(new Array<real>());	
		if(r==nullptr)r.reset(new Array<real>());	
		if(p==nullptr)p.reset(new Array<real>());	
		if(den==nullptr)den.reset(new Array<real>());	
		if(idx==nullptr)idx.reset(new Array<int>());
	}

	void Resize(const int size)
	{
		x->resize((size_type)size,VectorD::Zero());
		v->resize((size_type)size,VectorD::Zero());
		f->resize((size_type)size,VectorD::Zero());
		m->resize((size_type)size,(real)0);
		c->resize((size_type)size,(real)0);
		r->resize((size_type)size,(real)0);
		p->resize((size_type)size,(real)0);
		den->resize((size_type)size,(real)0);
		idx->resize((size_type)size,0);
	}

	int Add_Element()
	{
		x->push_back(VectorD::Zero());
		v->push_back(VectorD::Zero());
		f->push_back(VectorD::Zero());
		m->push_back((real)0);
		c->push_back((real)0);
		r->push_back((real)0);
		p->push_back((real)0);
		den->push_back((real)0);
		idx->push_back(0);
		return (int)x->size()-1;
	}

	int Size() const {return (int)(*x).size();}

	//////////////////////////////////////////////////////////////////////////
	////functions for separate attributes
	////functions for x
	VectorD& X(const int i)
	{return (*x)[i];}
	
	const VectorD& X(const int i) const 
	{return (*x)[i];}

	Array<VectorD>* X()
	{return x.get();}

	const Array<VectorD>* X() const 
	{return x.get();}
	
	ArrayPtr<VectorD> XPtr()
	{return x;}
	
	const ArrayPtr<VectorD> XPtr() const
	{return x;}
	
	Array<VectorD>& XRef()
	{return *x;}

	const Array<VectorD>& XRef() const 
	{return *x;}

	//////////////////////////////////////////////////////////////////////////
	////functions for v
	VectorD& V(const int i)
	{return (*v)[i];}
	
	const VectorD& V(const int i) const 
	{return (*v)[i];}

	Array<VectorD>* V()
	{return v.get();}

	const Array<VectorD>* V() const 
	{return v.get();}
	
	ArrayPtr<VectorD> VPtr()
	{return v;}
	
	const ArrayPtr<VectorD> VPtr() const
	{return v;}
	
	Array<VectorD>& VRef()
	{return *v;}

	const Array<VectorD>& VRef() const 
	{return *v;}

	//////////////////////////////////////////////////////////////////////////
	////functions for f
	VectorD& F(const int i)
	{return (*f)[i];}
	
	const VectorD& F(const int i) const 
	{return (*f)[i];}

	Array<VectorD>* F()
	{return f.get();}

	const Array<VectorD>* F() const 
	{return f.get();}
	
	ArrayPtr<VectorD> FPtr()
	{return f;}
	
	const ArrayPtr<VectorD> FPtr() const
	{return f;}
	
	Array<VectorD>& FRef()
	{return *f;}

	const Array<VectorD>& FRef() const 
	{return *f;}

	//////////////////////////////////////////////////////////////////////////
	////functions for m
	real& M(const int i)
	{return (*m)[i];}
	
	const real& M(const int i) const 
	{return (*m)[i];}

	Array<real>* M()
	{return m.get();}

	const Array<real>* M() const 
	{return m.get();}
	
	ArrayPtr<real> MPtr()
	{return m;}
	
	const ArrayPtr<real> MPtr() const
	{return m;}
	
	Array<real>& MRef()
	{return *m;}

	const Array<real>& MRef() const 
	{return *m;}

	//////////////////////////////////////////////////////////////////////////
	////functions for c
	real& C(const int i)
	{return (*c)[i];}
	
	const real& C(const int i) const 
	{return (*c)[i];}

	Array<real>* C()
	{return c.get();}

	const Array<real>* C() const 
	{return c.get();}
	
	ArrayPtr<real> CPtr()
	{return c;}
	
	const ArrayPtr<real> CPtr() const
	{return c;}
	
	Array<real>& CRef()
	{return *c;}

	const Array<real>& CRef() const 
	{return *c;}

	//////////////////////////////////////////////////////////////////////////
	////functions for r
	real& R(const int i)
	{return (*r)[i];}
	
	const real& R(const int i) const 
	{return (*r)[i];}

	Array<real>* R()
	{return r.get();}

	const Array<real>* R() const 
	{return r.get();}
	
	ArrayPtr<real> RPtr()
	{return r;}
	
	const ArrayPtr<real> RPtr() const
	{return r;}
	
	Array<real>& RRef()
	{return *r;}

	const Array<real>& RRef() const 
	{return *r;}

	//////////////////////////////////////////////////////////////////////////
	////functions for p
	real& P(const int i)
	{return (*p)[i];}
	
	const real& P(const int i) const 
	{return (*p)[i];}

	Array<real>* P()
	{return p.get();}

	const Array<real>* P() const 
	{return p.get();}
	
	ArrayPtr<real> PPtr()
	{return p;}
	
	const ArrayPtr<real> PPtr() const
	{return p;}
	
	Array<real>& PRef()
	{return *p;}

	const Array<real>& PRef() const 
	{return *p;}

	//////////////////////////////////////////////////////////////////////////
	////functions for den
	real& D(const int i)
	{return (*den)[i];}
	
	const real& D(const int i) const 
	{return (*den)[i];}

	Array<real>* D()
	{return den.get();}

	const Array<real>* D() const 
	{return den.get();}
	
	ArrayPtr<real> DPtr()
	{return den;}
	
	const ArrayPtr<real> DPtr() const
	{return den;}
	
	Array<real>& DRef()
	{return *den;}

	const Array<real>& DRef() const 
	{return *den;}

	//////////////////////////////////////////////////////////////////////////
	////functions for idx
	int& I(const int i)
	{return (*idx)[i];}
	
	const int& I(const int i) const 
	{return (*idx)[i];}

	Array<int>* I()
	{return idx.get();}

	const Array<int>* I() const 
	{return idx.get();}
	
	ArrayPtr<int> IPtr()
	{return idx;}
	
	const ArrayPtr<int> IPtr() const
	{return idx;}
	
	Array<int>& IRef()
	{return *idx;}

	const Array<int>& IRef() const 
	{return *idx;}
};
#endif
