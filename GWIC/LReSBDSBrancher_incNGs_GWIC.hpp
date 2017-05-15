
//***********************************************************************//

#include <gecode/int/branch.hh>
#include <gecode/int.hh>
#include <gecode/GWIC/GWIC.hpp>

namespace Gecode {
 


	//the following are used for cs_sbds

	int s_index;
	IntArgs MVH;
	IntArgs Aa,Bb,Cc,Dd,E,F,G,H,reF;
	//IntArgs WatNum,VarWNum,FailNum;
	IntArgs re_index,subArray;
	int ChoH;
	int nNum, nSym,_myc,_psize,nSymparcom;
	int dosize,x_size;
	IntArgs PartA; 
	IntArgs MyB;
	int size1;
	int size2;
	IntArgs ReA;
	int e_f_myb;
	int e_p_myb;
	int ct;
	Index_class _symmetries (int id, int index, int val){
	   
	 if(reF[id]==0)
	 {	if(val==Aa[id]) val=Bb[id];
		else if(val==Bb[id]) val=Aa[id];}
	 else  
	 if(reF[id]==1)
	 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
		else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;}
	 else if(reF[id]==2)
	 {	if(index%nNum==Aa[id]-1) index=index+Bb[id]-Aa[id];
		else if(index%nNum==Bb[id]-1) index=index-Bb[id]+Aa[id];}
	 else if (reF[id]==3)
	 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
		else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
		if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
		else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
	 }
	 else if (reF[id]==4)
	 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
		else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
		if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
		else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
		if(index>=nNum*(E[id]-1)&&index<nNum*E[id]) index=index+(F[id]-E[id])*nNum;
		else if(index>=nNum*(F[id]-1)&&index<nNum*F[id]) index=index-(F[id]-E[id])*nNum;
	 }
	 else if (reF[id]==5)
	 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
		else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
		if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
		else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
		if(index%nNum==G[id]-1) index=index+H[id]-G[id];
		else if(index%nNum==H[id]-1) index=index-H[id]+G[id];
	 }
	 else if (reF[id]==6)
	 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
		else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
		if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
		else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
		if(index>=nNum*(E[id]-1)&&index<nNum*E[id]) index=index+(F[id]-E[id])*nNum;
		else if(index>=nNum*(F[id]-1)&&index<nNum*F[id]) index=index-(F[id]-E[id])*nNum;
		if(index%nNum==G[id]-1) index=index+H[id]-G[id];
		else if(index%nNum==H[id]-1) index=index-H[id]+G[id];
	 }
	 return Index_class(index,val); 
	}  
	  


//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//*************   branch ***************************//



 template<class View, int n, class Val, unsigned int a>
class LReSBDSBrancher : public ViewValBrancher<View,n,Val,a> {
    typedef typename ViewBrancher<View,n>::BranchFilter BranchFilter;
  public:
     
	int start;
 
	int p_size;
	
	int p_myb;
	IntVar y;
  protected:
    /// Constructor for cloning \a b
    LReSBDSBrancher(Space& home, bool share, LReSBDSBrancher& b);
    /// Constructor for creation
	
	
    LReSBDSBrancher(Home home, 
                 ViewArray<View>& x, IntVar y,
                 ViewSel<View>* vs[n], 
                 ValSelCommitBase<View,Val>* vsc,
                 BranchFilter bf,IntVarValPrint vvp);
  public:
    /// Return choice
    virtual const Choice* choice(Space& home);
    /// Return choice
    virtual const Choice* choice(const Space& home, Archive& e);
    /// Perform commit for choice \a c and alternative \a b
    virtual ExecStatus commit(Space& home, const Choice& c, unsigned int b);
    /// Perform cloning
    virtual Actor* copy(Space& home, bool share);
    /// Perform dispose
    virtual size_t dispose(Space& home);
    /// Delete brancher and return its size
    static BrancherHandle post(Home home,
                               ViewArray<View>& x,IntVar y, 
                               ViewSel<View>* vs[n],
                               ValSelCommitBase<View,Val>* vsc,
                               BranchFilter bf,IntVarValPrint vvp);
  };

  template<class View, int n, class Val, unsigned int a> 
  LReSBDSBrancher<View,n,Val,a>
  ::LReSBDSBrancher(Home home, ViewArray<View>& x, IntVar yv,
                 ViewSel<View>* vs[n],
                 ValSelCommitBase<View,Val>* vsc,
                 BranchFilter bf,IntVarValPrint vvp)
    : ViewValBrancher<View,n,Val,a>(home, x, vs, vsc, bf,vvp),y(yv) 
  {		
		
		PartA=IntArgs(2*x.size());
		p_size=0;
		start=0;
		p_myb=0;
		home.notice(*this, AP_DISPOSE);
		 
	  
  }

  template<class View, int n, class Val, unsigned int a>
  forceinline BrancherHandle
  LReSBDSBrancher<View,n,Val,a>::
  post(Home home, ViewArray<View>& x, IntVar y,
       ViewSel<View>* vs[n], ValSelCommitBase<View,Val>* vsc,
       BranchFilter bf,IntVarValPrint vvp) {
		return *new (home) LReSBDSBrancher<View,n,Val,a>(home,x,y,vs,vsc,bf,vvp);
  }

  template<class View, int n, class Val, unsigned int a>
  forceinline
  LReSBDSBrancher<View,n,Val,a>::
  LReSBDSBrancher(Space& home, bool shared, LReSBDSBrancher<View,n,Val,a>& b)
    : ViewValBrancher<View,n,Val,a>(home,shared,b) {
	  
		size=b.p_size;
		start=b.start;
		p_myb=b.p_myb;
		p_myb=e_p_myb;
		y.update(home,shared,b.y); 
		 
  }
  
  template<class View, int n, class Val, unsigned int a>
  Actor*
  LReSBDSBrancher<View,n,Val,a>::copy(Space& home, bool shared) {
     
    return new (home) LReSBDSBrancher<View,n,Val,a>(home,shared,*this);
  }


  // Compute choice
  template<class View, int n, class Val, unsigned int a>
  const Choice*
  LReSBDSBrancher<View,n,Val,a>::choice(Space& home) { //std::cout<<"choice \n";
	return ViewValBrancher<View,n,Val,a>::choice(home);
 
  }

 template<class View, int n, class Val, unsigned int a>
  const Choice*
  LReSBDSBrancher<View,n,Val,a>::choice(const Space& home, Archive& e) {

  
    return ViewValBrancher<View,n,Val,a>::choice(home,e);
  } 
template<class View, int n, class Val, unsigned int a>
  size_t
  LReSBDSBrancher<View,n,Val,a>::dispose(Space& home) {
    home.ignore(*this,AP_DISPOSE);
    (void) ViewValBrancher<View,n,Val,a>::dispose(home);
    return sizeof(LReSBDSBrancher<View,n,Val,a>);
  }
  template<class View, int n, class Val, unsigned int a>
  ExecStatus
  LReSBDSBrancher<View,n,Val,a>
  ::commit(Space& home, const Choice& c, unsigned int b) {
  
		const PosValChoice<Val>& pvc
		  = static_cast<const PosValChoice<Val>&>(c);
		int pos = pvc.pos().pos;
		int val = pvc.val();
		_myc=1-b;

		//std::cout<<bb<<"\n";
		
		if(b==0)
		{      	
			p_myb=e_p_myb;
			 
			 //ajust symmetries
			 PartA[p_size*2]=pos;
			 PartA[p_size*2+1]=val;
			 p_size++;
			 _psize=p_size;
			 e_p_myb=p_myb;
			 
			 GECODE_ME_CHECK(this->x[pos].eq(home,val));
			 StatusStatistics c;
			 home.status(c); 
			 
				   
        }
		else
		{		
			//symmetry breaking constraint
			rel(home,y,IRT_GR,y.min());
			
			MyB[p_myb++]=pos; MyB[p_myb++]=val; MyB[p_myb++]=p_size-1;
			 
			e_p_myb=p_myb;
			 
			_psize=p_size;
			GECODE_ME_CHECK(this->x[pos].nq(home,val));
			StatusStatistics c;
			home.status(c);  
				 
        }
		
  

		return ES_OK;
	}
  
//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//**************************************************//
//*************   post ***************************//

  
 BrancherHandle
  branch(Home home, const IntVarArgs& x,
         IntVarBranch vars, IntValBranch vals,
         int n,int m,int dosize0, IntBranchFilter bf=NULL) {
		using namespace Int;
		if (home.failed()) return BrancherHandle();
		ViewArray<IntView> xv(home,x);
		ViewSel<IntView>* vs[1] = { 
		  Branch::viewselint(home,vars) 
		};
		int a=n*(n-1)/2;
		int b=m*(m-1)/2;

		nSym=n*m*n*m;
		nSymparcom=n-1+m-1;
		std::cout<<nSym<<"\t"<<n<<","<<m<<"\n";
		nNum=m;
		 Aa=IntArgs(nSym);
		 Bb=IntArgs(nSym);
		 Cc=IntArgs(nSym);
		 Dd=IntArgs(nSym);
		 reF=IntArgs(nSym);
		 
		 
	 
		
		 
		IntArgs _A(a+b);
		IntArgs _B(a+b);
		
		 
		int r=0;
		int r1=0;
		subArray=IntArgs(nSym*x.size());
		re_index=IntArgs(nSym*x.size());
		for(int i=0;i<nSym*x.size();i++)
			subArray[i]=0;
		for(int i=0;i<nSym*x.size();i++)
			re_index[i]=-1;
		 
		for(int i=0;i<n-1;i++)
		{
			int j=i+1;
			for(int j=i+1;j<n;j++)
		  {
			_A[r1]=i+1;
			_B[r1++]=j+1;
		   }
			Aa[r]=i+1;
			Bb[r]=i+2;

			
			for(int s=i*m;s<(i+1)*m;s++)
				subArray[r*x.size()+s]=1;
			for(int s=(i+1)*m;s<(i+2)*m;s++)
				subArray[r*x.size()+s]=1;
			reF[r++]=1;
		}
		for(int i=0;i<m-1;i++)
		{  
			int j=i+1;
		   for(int j=i+1;j<m;j++)
		   {
			_A[r1]=i+1;
			_B[r1++]=j+1;
			 
		   } 
			 Aa[r]=i+1;
			 Bb[r]=i+2;
			 
			for(int s=i;s<n*m;s+=m)
				subArray[r*x.size()+s]=1;
			for(int s=i+1;s<n*m;s+=m)
				subArray[r*x.size()+s]=1;	
			 reF[r++]=2;
		}   
	 
 
	 
	 
		for(int i=0;i<n*(n-1)/2;i++)
			for(int j=0;j<m*(m-1)/2;j++)
			//int i=1,j=1;
			{ 	Aa[r]=_A[i];
				Bb[r]=_B[i];
				Cc[r]=_A[n*(n-1)/2+j];
				Dd[r]=_B[n*(n-1)/2+j];
				int i1=_A[i]-1;
				int j1=_B[i]-1;
				int i2=_A[n*(n-1)/2+j]-1;
				int j2=_B[n*(n-1)/2+j]-1;
				for(int s=i1*m;s<(i1+1)*m;s++)
				{ subArray[r*x.size()+s]=1;}
				for(int s=j1*m;s<(j1+1)*m;s++)
				{  subArray[r*x.size()+s]=1;}
				for(int s=i2;s<n*m;s+=m)
				{  	subArray[r*x.size()+s]=1;}
				for(int s=j2;s<n*m;s+=m)
				{  	subArray[r*x.size()+s]=1;}
				 
				reF[r++]=3;
				 
			}
		   
		///////////////////////////////////////////////////////////////
		nSym=r;  
	 
		std::cout<<"------------nSym is "<<nSym<<"\n";
		s_index=x.size()*dosize0;
		dosize=dosize0;
		IntVar z(home,0,s_index+1);
		 
		//post global constraints
		ReA=IntArgs(x.size()*dosize0*2);
		size1=6;
		//create MyA,MyB
		x_size=x.size();	 
		MyB=IntArgs(s_index*3);
			 
		for(int r=0;r<nSym;r++)	
		{ 
			int k=0; 
			for(int i=0;i<x_size;i++)
			if(subArray[r*x_size+i]==1)
			{	 
				re_index[r*x_size+i]=k; 
				 
				k++;
			}
			IntVarArgs y(home,k,0,dosize0-1);
			k=0;
			for(int i=0;i<x_size;i++)
			if(subArray[r*x_size+i]==1)
			{	y[k]=x[i];
				 
				k++;       
			}
			
			gwic(home,y,z,r);
		}	
		//build a columwise heuristic
		MVH=IntArgs(x.size());
		for(int i=0;i<m;i++)
			for(int j=0;j<n;j++)
				MVH[i*n+j]=j*m+i;
	 
			 
		return LReSBDSBrancher<IntView,1,int,2>::post
			(home,xv,z,vs,Branch::valselcommitint(home,x.size(),vals),bf,NULL);
    }
}


  

