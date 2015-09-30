//============================================================================
//                                  I B E X                                   
// File        : ibex_SolverQInter.cpp
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 13, 2014
// Last Update : Dec 26, 2014
//============================================================================


#include "ibex_SolverQInter.h"

#include "ibex_NoBisectableVariableException.h"
#include <algorithm>


using namespace std;

namespace ibex {


 SolverQInter::SolverQInter (Ctc& ctc, Bsc& bsc, CellBuffer& buffer, CtcQInter& ctcq)  : Solver (ctc,bsc,buffer), ctcq(ctcq) 
 { possiblesols=0;
   qvalid=ctcq.q;
   /*
   cout << " nb points " << ctcq.ctc_list.size() << endl;
   for(int i=0;i<ctcq.ctc_list.size();i++)
     nb_points_sol.push_back(0);
   */
 }

  SolverQInter::SolverQInter (Ctc& ctc, Bsc& bsc, CellBuffer& buffer, CtcQInter& ctcq, int qvalid)  :
    Solver (ctc,bsc,buffer), ctcq(ctcq) , qvalid(qvalid)
 { possiblesols=0;
 }

 
  /* the backtrackable  list of active points is put from the cell into the qinter constraint : shared pointer */
  void SolverQInter::manage_cell_info(Cell& c) {
    ctcq.points=c.get<QInterPoints>().points;
    ctcq.qmax= ctcq.points->size();
  }

  void SolverQInter::precontract(Cell& c) {
    manage_cell_info(c);
    ctc.enable_side_effects();
    ctcq.enable_side_effects();   // allows points list to be updated during the contraction 
  }

  void SolverQInter::postcontract(Cell& c) {
    ctcq.disable_side_effects();
    ctc.disable_side_effects();
  }

  void SolverQInter::prebisect(Cell& c)
  {  q1=ctcq.points->size();
    if (ctcq.qmax < q1) {
      //cout << " q1 reduit " << q1  << " " << ctcq.qmax << endl ; 
      q1=ctcq.qmax;
    }
     q2 = ctcq.midbox_activepoints_number(c.box); 
     
     ValidPoint* validpoint=&(c.get<ValidPoint>());

     if (! c.box.contains(*(validpoint->point)) || q2 >= validpoint->validpoints_number)
	   { 
	     Vector* mid  = new Vector (c.box.mid());
	     validpoint->validpoints_number=q2;
	     delete validpoint->point;
	     validpoint->point= mid;
	   }
	 
     if (q2== q1)  // stop condition  : all the compatible points are valid : no need of further bisections.
       {//cout << " q1 = q2 " << q1 << " " << q2 << endl;  
	 throw NoBisectableVariableException();
       }

    
  }

  Cell* SolverQInter::root_cell(const IntervalVector& init_box) 
  {Cell* root= new Cell(init_box) ; 
    root->add<QInterPoints>();
    root->add<ValidPoint>();
    QInterPoints* qinterpoints=&root->get<QInterPoints>();
    ValidPoint* validpoint=&root->get<ValidPoint>();
    qinterpoints->points= ctcq.points; 
    Vector* mid = new Vector(init_box.mid());
    validpoint->point= mid;
    validpoint->validpoints_number=0;
    ctcq.points_to_delete= false;  // ctcq.points will be deleted by the QInterPoints destructor.
    return root;}
 




  // compute the compatible measurements (stored in compatible_sol_points)
  // and the valid measurments at the middle of the box (stored in valid_sol_points) :
  // only the "currently maximal" solutions have their valid measurements stored.
  // each solution has a maximal boolean maximal_sol
  // the number of compatible measurements is the size of compatible_sol_points of the solution
  // the number of valid measurements is the size of valid_sol_points  of the solution
  // returns true if the number of compatible measurements is greater than ctcq.q and 
  // the if the number of valid measurements is greater than qvalid false otherwise
  // only these solutions will be used in the post-treatment


  bool SolverQInter::solution_test(Cell&  c)  {

    q1= ctcq.activepoints_contract_count(c.box);  // contract the measures for checking it is a solution (useful when the contraction is not with a fixpoint)
    // ctcq.points is updated
    //    cout <<  c.box <<  " q1 " << q1 << endl;

    if (q1 >= ctcq.q)
      // the compatible points are computed  here (because the need the current value of ctcq.points that will 
      // modified for valid points
      {set<int>* compatiblepoints=  new set<int>;
	list<int>::iterator iter = ctcq.points->begin() ;
	while (iter != ctcq.points->end())
	  {compatiblepoints->insert(*iter); iter++;}
	

	

	ValidPoint* validpoint=&(c.get<ValidPoint>());
	IntervalVector  vec (*(validpoint->point));
	// we use the same structure ctcq.points that is updated to contain the valid points
	q2=ctcq.activepoints_contract_count(vec);
	//	cout << " q " << ctcq.q << " q1 " << q1 << " q2 " << q2 << endl;
        if (q2 >= qvalid) {

	  
        
	set<int>* validpoints = new set<int>;
	iter = ctcq.points->begin() ;
	while (iter != ctcq.points->end())
	  {validpoints->insert(*iter); iter++;}
	//	double err_sol1= compute_err_sol(c.box);
	double err_sol1= compute_err_sol(vec);

	err_sol.push_back(err_sol1);
        if (trace)	cout << " solution " << err_sol.size() << "  " << vec << "  "  << err_sol1 << endl;	    
	list<set<int>*  >::iterator itersol= valid_sol_points.begin();
	list<set<int>*  >::iterator compatitersol= compatible_sol_points.begin();
	bool maximal=true;
	for (int i=0; i< maximal_sol.size(); i++)  // remove the previous solutions with a set of valid points included in the current one
	  { 
	    if (maximal_sol[i])
	      {if (
		   //(*itersol)->size() + ind <= validpoints->size()  && 
		   ((*itersol)->size() < validpoints->size() || 
		    ((*itersol)->size() ==validpoints->size() && err_sol[i] > err_sol1))
		   && 
		   includes (validpoints->begin(), validpoints->end(), (*itersol)->begin(), (*itersol)->end()))
		  { 
		    maximal_sol[i]=false;
		    delete *itersol; itersol=valid_sol_points.erase(itersol);
		    delete *compatitersol; compatitersol=compatible_sol_points.erase(compatitersol);
		    continue;}
		else if (
			 // (*itersol)->size() >= validpoints->size() + ind
			 ((*itersol)->size() > validpoints->size() || 
			  ((*itersol)->size() ==validpoints->size() && err_sol[i] < err_sol1))
                         &&
			 includes ( (*itersol)->begin(), (*itersol)->end(), validpoints->begin(), validpoints->end()))
		  {maximal=false; break;}
		
		itersol++; compatitersol++;}
	  }
	
	if (maximal) {  // the current solution is maximal (for the moment)
	  valid_sol_points.push_back(validpoints);
	  compatible_sol_points.push_back(compatiblepoints);
	  maximal_sol.push_back(true);
	}
	else
	  {maximal_sol.push_back(false);
	    delete compatiblepoints;
	    delete validpoints;}

	return true;
	}
	{delete compatiblepoints;  possiblesols++;return false;}
      }
    return false;
  }


  // destructor : the sets  valid_sol_points and compatible_sol_points are deleted
  SolverQInter::~SolverQInter() {
    list<set<int>* >::iterator itersol= valid_sol_points.begin(); 
    while (itersol != valid_sol_points.end())  {delete *itersol; itersol++;}
    itersol= compatible_sol_points.begin(); 
    while (itersol != compatible_sol_points.end())  {delete *itersol; itersol++;}
  }



  // report the maximal solutions 
  void SolverQInter::report_maximal_solutions(vector<IntervalVector> & res) {
    int kk=0;
    cout << "  " << res.size() << " " <<  compatible_sol_points.size() << " " << valid_sol_points.size() << endl;
    list<set<int>* >::iterator itersol= valid_sol_points.begin();
    list<set<int>* >::iterator compatitersol= compatible_sol_points.begin();
    for (int i=0; i<res.size(); i++)
      {
        if (maximal_sol[i]){
	  //	  if ((*itersol)->size() >=q)
	  kk++; cout << "Solution  " << kk << " " << i+1 << " " << res[i] << " " <<   (*compatitersol)->size() << " compatible measurements " 
			<< (*itersol)->size() << " valid measurements with error " <<  err_sol[i] << endl;

	    set<int>::iterator iter= (*itersol)->begin();
	    while (iter !=(*itersol)->end())
	      {cout << " " << *iter ;

		iter++; }
	    cout << endl;
	  
	  iter= (*compatitersol)->begin();
	    while (iter !=(*compatitersol)->end())
	      {cout << " " << *iter ;

		iter++; }
	    cout << endl;
	  itersol++;
	  compatitersol++;
	}
      }
  }



  // compare pairwise the maximal sols and if they share more than the half of measurements, keep the one that has more measurements
  //   in case of ties , depending on ind , keeps only one or both
  // old version , without taking account of the errors 
  /*
  void SolverQInter::keep_one_solution_pergroup (vector<IntervalVector> & res, int q,int ind) {
    


    list<set<int>* >::iterator itersol= valid_sol_points.begin();
    list<set<int>* >::iterator itersol2= valid_sol_points.begin();
    vector<int> max_sol;
    
    for (int i=0; i<res.size(); i++)
      max_sol.push_back(maximal_sol[i]);
    for (int i=0; i<res.size(); i++)
      if (maximal_sol[i])
          {  
	    if ((*itersol)->size() < q){  max_sol[i]=0;}
	    else{
             itersol2=itersol;
	     for (int j=i+1; j<res.size(); j++)
	       
	       if (maximal_sol[j])
		{itersol2++;
		  if (max_sol[j])
		    if(// neighbors (res[i] , res[j]) && 
		       sol_intersection (*(*itersol), *(*itersol2)))
		      if // ((*itersol)->size() > (*itersol2)->size() ||(*itersol2)->size()<q )
			((*itersol)->size() >= (*itersol2)->size() + ind||(*itersol2)->size()<q )
			  { max_sol[j]=0;}
			else if ((*itersol)->size() < (*itersol2)->size())
			  { max_sol[i]=0; break;}
		}
	    }
	    itersol++;
	  }
  
    int nbsol=0;
    for (int i=0; i<res.size(); i++)
      if (max_sol[i])
	{ nbsol++;
	  cout << " Sol " << nbsol << " " << i+1 << " " << res[i] << endl;
	}
  }
  */

  // sorting the solutions : first criterion : the size of validpoints, second criterion the error

  bool comparvalidsol (const pair <set<int>* , pair <int, double> > &p1, const pair <set<int>* , pair <int, double> > &p2) 
  {return ((p1.first)->size() > (p2.first)->size()
	   ||
	   ((p1.first)->size() == (p2.first)->size()
	    && p1.second.second < p2.second.second
	    ));}


 // compare pairwise the maximal sols and if they share more than the half of measurements, keep the one that has more measurements
  //   in case of ties ,
 void SolverQInter::keep_one_solution_pergroup (vector<IntervalVector> & res) {
    

    list<set<int>* >::iterator itersol= valid_sol_points.begin();

    vector<int> max_sol;

    vector<pair<set<int>* , pair <int , double> > > valid_sol_index;
    
    
    for (int i=0; i<res.size(); i++)
      max_sol.push_back(maximal_sol[i]);


    for (int i=0; i<res.size(); i++)
      if (maximal_sol[i])
	{  

	  valid_sol_index.push_back(make_pair (*itersol, make_pair (i, err_sol[i])));

	  itersol++;
	}


    sort(valid_sol_index.begin(),valid_sol_index.end(), comparvalidsol);

    for (int i=0; i< valid_sol_index.size() ; i++)
      if (max_sol[valid_sol_index[i].second.first])
	for (int j=i+1; j< valid_sol_index.size(); j++)
	  if (max_sol[valid_sol_index[j].second.first])
	    if(
	       sol_intersection (*(valid_sol_index[i].first), *(valid_sol_index[j].first)))
		  
		  { max_sol[valid_sol_index[j].second.first]=0;}
  
    int nbsol=0;
    cout << "sorted solutions " << endl;
    for (int i=0; i<valid_sol_index.size(); i++)
      if (max_sol[valid_sol_index[i].second.first])
	{ nbsol++;
	  cout << " Sol " << nbsol <<  " " << valid_sol_index[i].second.first+1 << res[valid_sol_index[i].second.first]
	       << " error " << err_sol [valid_sol_index[i].second.first] 
	       << "  " << (valid_sol_index[i].first)->size() << " valid points " <<endl;
	  set<int>::iterator iter= (valid_sol_index[i].first)->begin();
	  while (iter !=valid_sol_index[i].first->end())
	    {cout << " " << *iter ;

	      iter++; }
	  cout << endl;


	}
    // the number of solutions that have not reached the qvalid limit.
    cout << " possible sols " << possiblesols;

  }


  // test for 2 valid solutions to be compared ; sharing more than the half of the measurements



  bool  SolverQInter::sol_intersection( set<int>& solution1,  set<int>& solution2)
  {set <int> intersol;
    set_intersection (solution1.begin(), solution1.end(), solution2.begin() , solution2.end(), inserter (intersol, intersol.end()));
    /*
    cout << " solution 1 " << solution1.size() <<  endl;
    set<int>::iterator iter= (solution1.begin());
    while (iter !=solution1.end())
      {cout << " " << *iter ;
	iter++; }
    cout << endl;
    cout << " solution 2 " << solution2.size() << endl;
     iter= (solution2.begin());
    while (iter !=solution2.end())
      {cout << " " << *iter ;
	iter++; }
    cout << endl;
    cout << " intersection " <<  intersol.size() << endl;
    iter= (intersol.begin());
    while (iter !=intersol.end())
      {cout << " " << *iter ;
	iter++; }
    */
   
    if (intersol.size() > solution1.size() /2 || intersol.size() > solution2.size() /2) return true;
    else return false;
  }

  /* unused 
  bool  SolverQInter::neighbors (IntervalVector& v1, IntervalVector& v2)
{double d;
  double fact=10;
  for (int i=0; i< v1.size(); i++)
    { if (v1[i].lb() > v2[i].ub()) 
	d= v1[i].lb() - v2[i].ub();
      else
	if (v2[i].lb() > v1[i].ub())
	  d= v2[i].lb() - v1[i].ub();
      if (d > fact* bsc.prec(i)) return false;
    }
  
  return true;
}
  */			     

  double SolverQInter::compute_err_sol(IntervalVector& box){
    return ctcq.compute_err_sol(box);
  }

    

} // end namespace ibex



