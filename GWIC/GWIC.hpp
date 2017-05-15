#include <gecode/int/branch.hh>
#include <gecode/int.hh>
#include <gecode/GWIC/head.hh>

namespace Gecode{  
 
	  
	//extern IntArgs Aa,Bb,Cc,Dd,reF;

	extern int nNum, nSym, _myc, _psize,nSymparcom;
	extern IntArgs PartA,re_index; 
	extern IntArgs MyB;
	extern int size1,size2;
	extern IntArgs ReA;
	extern int e_p_myb;
	//extern IntArgs WatNum,VarWNum,FailNum;
	extern int ChoH;
	extern int dosize;	
	extern int x_size;
	template<class View>
	 class GWIC : public Propagator {
	 protected: 
		ViewArray<View> x;
		Int::IntView y;

		 
		int index;//the index of the symmetry
		int var;//first pointer.pos
		int val;//first pointer.val
		int pos_a;//the position of first pointer
		int pos_b;//the position of rhs assignment
		int re;
		 
	 public: 
		//post
		GWIC(Space& home,  ViewArray<View>& x0, Int::IntView y,int r);
		//copy
		GWIC(Space& home, bool share, GWIC<View>& p);
		virtual Propagator* copy(Space& home, bool share);
		//cost
		virtual PropCost cost(const Space&, const ModEventDelta&) const;
		//propagation
		virtual ExecStatus propagate(Space& home, const ModEventDelta&);
		//post
		static ExecStatus post(Space& home, ViewArray<View>& x0,Int::IntView y, int r);
		//dispose
		virtual size_t dispose(Space& home);
		int symGoal(Space& home,int re);
	};

	template<class View>
	int
	GWIC<View>::symGoal(Space& home,int re) {
		
			 
			 
	 };
		
	// posting
	template<class View>
	inline
	GWIC<View>::GWIC(Space& home,   ViewArray<View>& x0,Int::IntView y0, int r) : Propagator(home), x(x0),y(y0){
		//initialize 
		index=r;
		var=val=pos_a=-1;
		pos_b=0;
		re=0;
		//x.subscribe(home,*this,Int::PC_INT_DOM); 
		y.subscribe(home,*this,Int::PC_INT_DOM);  
	} 
	template<class View>
	ExecStatus
	GWIC<View>::post(Space& home, ViewArray<View>& x0,Int::IntView y0,  int r) {
		(void) new (home) GWIC(home,x0,y0,r); 
		return ES_OK; 
	} 
	// disposal 
	template<class View>
	forceinline size_t
	GWIC<View>::dispose(Space& home) {
		 
		 
		ChoH=-1;
		(void) Propagator::dispose(home);
		return sizeof(*this); 
	} 
	// copying
	template<class View>
    forceinline
    GWIC<View>::GWIC(Space& home, bool share, GWIC<View>& p) : Propagator(home,share,p),
														index(p.index),var(p.var),val(p.val),pos_a(p.pos_a),pos_b(p.pos_b),re(p.re){
												 
		 
		x.update(home,share,p.x); 
		y.update(home,share,p.y); 
		
	 
	} 
	template<class View> Propagator*  GWIC<View>::copy(Space& home, bool share) {
		return new (home) GWIC(home,share,*this); 
	} 
	// cost computation 
	template<class View>
	PropCost
	GWIC<View>::cost(const Space&, const ModEventDelta&) const {
		return PropCost::linear(PropCost::LO,int(100000000000));
	} 
	
	 
	// propagation 
	template<class View>
	ExecStatus
	GWIC<View>::propagate(Space& home, const ModEventDelta&) { 
	 
		int mre=re;
		re=0;
		 
		if(var==-1 &&	val==-1 && pos_a+1!=_psize)//need to find the first pointer
		{	  
			int s;
			for(s=pos_a+1;s<_psize;s++)
			{
				 
				Index_class g=(*_symmetries)(index,PartA[s*2],PartA[1+s*2]);
				if(g.index==PartA[s*2]&&g.val==PartA[1+s*2]){
					int p=pos_b;
					for(;p<e_p_myb&& MyB[p+2]<s+1;p+=3)
					{	
						Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
						if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
						if(x[re_index[index*x_size+g.index]].in(g.val))
						{	
							//std::cout<<"prune2 "<<index<<":"<<MyB[p]<<","<<MyB[p+1]<<"\t"<<g.index<<","<<g.val<<"\n";
							GECODE_ME_CHECK(x[re_index[index*x_size+g.index]].nq(home,g.val));
							ReA[re++]=g.index;
							ReA[re++]=g.val;
						}
					}
					pos_b=p;
					continue;
				}
				if(!x[re_index[index*x_size+g.index]].in(g.val))//now the symmetry is broken
				{	
					 
					if(re>0) GECODE_ME_CHECK(y.nq(home,y.min()));  
					for(int i=0;i<re;i+=2)
					{	 

						MyB[e_p_myb++]=ReA[i]; MyB[e_p_myb++]=ReA[i+1]; MyB[e_p_myb++]=_psize-1;
					} 
					y.cancel(home,*this,Int::PC_INT_DOM);
					return home.ES_SUBSUMED(*this);
					 
				}
				if(x[re_index[index*x_size+g.index]].in(g.val)&&!x[re_index[index*x_size+g.index]].assigned())//now find the first pointer
				{
				
					var=g.index;
					val=g.val;
					pos_a=s; 
					 
					y.cancel(home,*this,Int::PC_INT_DOM);
					//add extra checks in order to prune the 
					int ns=0;
					int p=pos_b;
					int t;
					for(t=s+1;t<_psize;t++)
					{
						 
							 
						for(;p<e_p_myb&& MyB[p+2]<t;p+=3)
						{	 
							Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
							if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
							if(x[re_index[index*x_size+g.index]].in(g.val))
							{	
								ns=1;
								break;
							}
						}
						if(ns==1) break;
					 
						Index_class g=(*_symmetries)(index,PartA[2*t],PartA[1+2*t]);
						if(g.index==PartA[2*t]&&g.val==PartA[1+2*t]) continue;
						if(!x[re_index[index*x_size+g.index]].in(g.val))//now the symmetry is broken
						{	
							ns=-1; 
							break;
						}
						
		 
					}
					if(ns==-1) 
					{
						//std::cout<<"prune\t"; 
						if(re>0) GECODE_ME_CHECK(y.nq(home,y.min()));  
						for(int i=0;i<re;i+=2)
						{	 

							MyB[e_p_myb++]=ReA[i]; MyB[e_p_myb++]=ReA[i+1]; MyB[e_p_myb++]=_psize-1;
						} 
						 
						return home.ES_SUBSUMED(*this);
					}
					else
						x[re_index[index*x_size+var]].subscribe(home,*this,Int::PC_INT_VAL); 
					break;
				}
				else
				{
					int p=pos_b;
					for(;p<e_p_myb&& MyB[p+2]<s+1;p+=3)
					{	 
						Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
						if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
						if(x[re_index[index*x_size+g.index]].in(g.val))
						{	
							GECODE_ME_CHECK(x[re_index[index*x_size+g.index]].nq(home,g.val));
							ReA[re++]=g.index;
							ReA[re++]=g.val;
						}
					}
					pos_b=p;
				}
					
			}
			 
			if(s==_psize) pos_a=s-1;
			
		}
		 

		else 
		{	if( var!=-1 &&	val!=-1)//not subsumed yet and have constraints
			{	
				if(!x[re_index[index*x_size+var]].in(val))//the symmetry is broken 
				{	
					x[re_index[index*x_size+var]].cancel(home,*this,Int::PC_INT_VAL);  
					if(re>0) GECODE_ME_CHECK(y.nq(home,y.min()));  
					for(int i=0;i<re;i+=2)
					{	MyB[e_p_myb++]=ReA[i]; MyB[e_p_myb++]=ReA[i+1]; MyB[e_p_myb++]=_psize-1;
					}
					
					return home.ES_SUBSUMED(*this);
				}	
				if(x[re_index[index*x_size+var]].in(val)&&x[re_index[index*x_size+var]].assigned())//need to find a new first watched literal
				{	
					int p=pos_b;
					for(;p<e_p_myb&& MyB[p+2]<pos_a+1;p+=3)
					{	
						Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
						if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
						if(x[re_index[index*x_size+g.index]].in(g.val))
						{	
							
							GECODE_ME_CHECK(x[re_index[index*x_size+g.index]].nq(home,g.val));
							ReA[re++]=g.index;
							ReA[re++]=g.val;
						}
					}
					pos_b=p;
					x[re_index[index*x_size+var]].cancel(home,*this,Int::PC_INT_VAL);
					
					int s;
					
					for(s=pos_a+1;s<_psize;s++)
					{
						
						Index_class g=(*_symmetries)(index,PartA[2*s],PartA[1+2*s]);
						if(g.index==PartA[2*s]&&g.val==PartA[1+2*s]){
							int p=pos_b;
							for(;p<e_p_myb&& MyB[p+2]<s+1;p+=3)
							{	Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
								if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
								if(x[re_index[index*x_size+g.index]].in(g.val))
								{	
									GECODE_ME_CHECK(x[re_index[index*x_size+g.index]].nq(home,g.val));
									ReA[re++]=g.index;
									ReA[re++]=g.val;
								}
							}
							pos_b=p;
							continue;
						}
						if(!x[re_index[index*x_size+g.index]].in(g.val))//now the symmetry is broken
						{	
							for(int p=pos_b;p<e_p_myb&& MyB[p+2]<s;p+=3)
							{	Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
								if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
								if(x[re_index[index*x_size+g.index]].in(g.val))
								{	
									GECODE_ME_CHECK(x[re_index[index*x_size+g.index]].nq(home,g.val));
									ReA[re++]=g.index;
									ReA[re++]=g.val;
								}
							}
							if(re>0) GECODE_ME_CHECK(y.nq(home,y.min()));  
							for(int i=0;i<re;i+=2)
							{	MyB[e_p_myb++]=ReA[i]; MyB[e_p_myb++]=ReA[i+1]; MyB[e_p_myb++]=_psize-1;
							}
								 	 
							return home.ES_SUBSUMED(*this);
						}
						if(x[re_index[index*x_size+g.index]].in(g.val)&&!x[re_index[index*x_size+g.index]].assigned())//now find the first pointer
						{
							var=g.index;
							val=g.val;
							pos_a=s; 
							int ns=0;
							int t;
							int p=pos_b;
							for(t=s+1;t<_psize;t++)
							{
								for(;p<e_p_myb&& MyB[p+2]<t;p+=3)
								{	Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
									if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
									if(x[re_index[index*x_size+g.index]].in(g.val))
									{	
										ns=1;
										break;
									}
								}
								if(ns==1) break;
								
								Index_class g=(*_symmetries)(index,PartA[2*t],PartA[1+2*t]);
								if(g.index==PartA[2*t]&&g.val==PartA[1+2*t]) continue;
								if(!x[re_index[index*x_size+g.index]].in(g.val))//now the symmetry is broken
								{	
									ns=-1; 
									break;
								}

							}
							if(ns==-1) 
							{
								
								if(re>0) GECODE_ME_CHECK(y.nq(home,y.min()));  
								for(int i=0;i<re;i+=2)
								{	 
						
									MyB[e_p_myb++]=ReA[i]; MyB[e_p_myb++]=ReA[i+1]; MyB[e_p_myb++]=_psize-1;
								} 
							 
								return home.ES_SUBSUMED(*this);
							}
							else
								x[re_index[index*x_size+var]].subscribe(home,*this,Int::PC_INT_VAL);
							break;
						}
						else
						{	p=pos_b;
							for(;p<e_p_myb&& MyB[p+2]<s+1;p+=3)
							{	Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
								if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
								if(x[re_index[index*x_size+g.index]].in(g.val))
								{	
									GECODE_ME_CHECK(x[re_index[index*x_size+g.index]].nq(home,g.val));
									ReA[re++]=g.index;
									ReA[re++]=g.val;
								}
							}
							pos_b=p;
						}
					 
					}
					 
					//-----------------------------------------------------------------effect all nogoods whose LHSs are before s
					
					
					 
					 
					if(s==_psize)
					{
						//std::cout<<"set the first pointer\n";
						var=-1;
						val=-1;
						pos_a=s-1;
						y.subscribe(home,*this,Int::PC_INT_DOM);
					}
					
				}
				
			
				
			}
			else if(pos_b!=e_p_myb)
			{	int p=pos_b;
				for(;p<e_p_myb;p+=3)
				{	Index_class g=(*_symmetries)(index,MyB[p],MyB[p+1]);
					if(g.index==MyB[p]&&g.val==MyB[p+1]) continue;
					if(x[re_index[index*x_size+g.index]].in(g.val))
					{	
						GECODE_ME_CHECK(x[re_index[index*x_size+g.index]].nq(home,g.val));
						 
						ReA[re++]=g.index;
						ReA[re++]=g.val;
					}
				}
				pos_b=p;
			}
			
			 
			
			
		}
		//std::cout<<"\n"; 
		if(re>0) GECODE_ME_CHECK(y.nq(home,y.min()));  
		for(int i=0;i<re;i+=2)
		{	 
			
			MyB[e_p_myb++]=ReA[i]; MyB[e_p_myb++]=ReA[i+1]; MyB[e_p_myb++]=_psize-1;
		} 
		
	 	
	 
		return ES_FIX;

	
	} 


	void gwic(Space& home, const IntVarArgs& x,Int::IntView y, int r) { 
		// constraint post function 
		ViewArray<Int::IntView> xv(home,x); 
		if (GWIC<Int::IntView>::post(home,xv,y,r) != ES_OK)
			home.fail(); 
	} 
}