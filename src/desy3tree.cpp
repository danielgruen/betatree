#include <iostream>
using namespace std;

#include <CCfits/CCfits>
using namespace CCfits;

#include "filter.h"

const int nmin=100; // minimum number of reference galaxies per leaf
const double mindeltam=0.1; // minimum width of leafs (mag)
const double brightcut=19.; // objects brighter than i=19 ignored from the tree
const double zp=30.;


void growTree(ObjectCollection *refcat, const vector<string> &criteria, const vector<double> detectionlimit,
              const vector<double> &minwidth, vector<double> min, vector<double> max, int nmin, string target, string target2, ObjectCollection *leaves)
{
  // target is something like DdsDs, target2 something like DdsDssq

  // stop if below size limit
  if(refcat->size()<2*nmin) {
     //cerr << "cell with n=" << refcat->size() << " not populated enough to split, saving" << endl; 
     Object *tmp = new Object(leaves->prototype);
     for(int i=0; i<criteria.size(); i++)
  	 {
  	 	tmp->setDoubleProperty(criteria[i]+"_min",min[i]);
  	 	tmp->setDoubleProperty(criteria[i]+"_max",max[i]);
  	 	
  	 }
  	 leaves->appendObject(tmp);
  	 return;
  }
  
  ObjectCollection *tsmallmax=0;
  ObjectCollection *tlargemax=0;
  
  // max snr cut
  int cmax=-1; double snrmax=0.; double medmax=0.; double shearsnrmax=0.;
  
  
  for(int i=0; i<criteria.size(); i++)
  {
     if(max[i]-min[i]<2.*minwidth[i]) // don't split any further
     {
       continue;
     }
  
     // get median
     refcat->sort(criteria[i],true);
     double med=(refcat->objects[int((refcat->size()-1)/2)]->doublePropertyValue(criteria[i])+
                refcat->objects[int((refcat->size()+1)/2)]->doublePropertyValue(criteria[i]))/2.;
 
     if((i==0 && med>detectionlimit[0]) || (i>0 && med+max[0]>detectionlimit[i]) ) 
     // don't split if any band is above detection limit at median split
     {
       //cerr << "# band " << criteria[i] << " is above detection limit, not splitting" << endl;
       continue;
     }
     
     // split in two
     ObjectCollection *tsmall = refcat->filter(Filter(criteria[i],min[i],med));
     ObjectCollection *tlarge = refcat->filter(Filter(criteria[i],med,max[i]));
     
     double asmall   = tsmall->average(target);
     double alarge   = tlarge->average(target);
     double ssmall   = tsmall->avgstdv(target);
     double slarge   = tlarge->avgstdv(target);
     double diffsnr  = fabs(asmall-alarge)/sqrt(ssmall*ssmall+slarge*slarge);
     double shearsnr = asmall*asmall+alarge*alarge;
     if(shearsnr>shearsnrmax && diffsnr>2.) {
       cmax=i; snrmax=diffsnr; medmax=med; shearsnrmax=shearsnr; 
       delete tsmallmax; tsmallmax=tsmall; 
       delete tlargemax; tlargemax=tlarge;
     } else {
       delete tsmall; delete tlarge;
     }
  }
  
  if(snrmax>2.) {
    //cerr << "found a significant difference when splitting by crit " << cmax << "=" << criteria[cmax] << "; branching" << endl;
    double dmax=max[cmax];
    max[cmax]=medmax;
    //cerr << "growing first branch at " << min[cmax] << "<" << criteria[cmax] << "<" << max[cmax] << endl;
    ObjectCollection *refcatlow = refcat->filter(Filter(criteria[cmax],min[cmax],max[cmax]));
    growTree(refcatlow, criteria, detectionlimit, minwidth, min, max, nmin, target, target2, leaves);
    delete refcatlow;

    max[cmax]=dmax; min[cmax]=medmax;
    //cerr << "growing second branch at " << min[cmax] << "<" << criteria[cmax] << "<" << max[cmax] << endl;
    ObjectCollection *refcathigh = refcat->filter(Filter(criteria[cmax],min[cmax],max[cmax]));
    growTree(refcathigh, criteria, detectionlimit, minwidth, min, max, nmin, target, target2, leaves);
    delete refcathigh;
  } else {
     //cerr << "difference insignificant; turning into leaf" << endl;
     Object *tmp = new Object(leaves->prototype);
     for(int i=0; i<criteria.size(); i++)
  	 {
  	 	tmp->setDoubleProperty(criteria[i]+"_min",min[i]);
  	 	tmp->setDoubleProperty(criteria[i]+"_max",max[i]);
  	 }
  	 leaves->appendObject(tmp);
  }
  
  delete tsmallmax; delete tlargemax;
  
}


ObjectCollection* tree(ObjectCollection *refcat, const vector<string> &criteria, const vector<double> &detectionlimit, 
                       const vector<double> &minwidth, int nmin, string target, string target2, double imin=20., double imax=25.)
// build a tree such that in each leaf there are at least nmin members (but more than 2*nmin if no significant difference in target)
{
  const int nkeys=criteria.size()*2;
  string *binnames =new string[nkeys];
  int    *binkeys  =new int[nkeys];
  Type   *bintypes =new Type[nkeys];
  
  vector<double> min(criteria.size(),-99);
  vector<double> max(criteria.size(),99);
  min[0]=imin; max[0]=imax;
  
  for(int i=0; i<criteria.size(); i++)
  {
      binnames[2*i]=criteria[i]+"_min";
      binnames[2*i+1]=criteria[i]+"_max";
      binkeys[2*i]=binkeys[2*i+1]=-1;
      bintypes[2*i]=bintypes[2*i+1]=doubleType;
  }
  
  ObjectPrototype *leavesPrototype = new ObjectPrototype(binkeys, binnames, bintypes, criteria.size()*2);
  ObjectCollection *leaves = new ObjectCollection(leavesPrototype);
  growTree(refcat, criteria, detectionlimit, minwidth, min, max, nmin, target, target2, leaves);
  delete []binnames; delete []binkeys; delete []bintypes; delete leavesPrototype;
  return leaves;
}



double flux_to_croppedmag(double flux)
{
  if(flux<1.) flux=1;
  return zp-2.5*log10(flux);
}


double mean_uncertainty(double mean, double meansq, int n)
{
  return sqrt((meansq-mean*mean)/(n-1.));
}

string ZeroPadNumber(int num, int length=3)
{
	stringstream ss;
	
	ss << num; 
	string ret;
	ss >> ret;
	
	// Append zero chars
	int str_length = ret.length();
	for (int i = 0; i < length - str_length; i++)
		ret = "0" + ret;
	return ret;
}

void leafBin(Object *leaf, ObjectCollection *pop, vector<string> criteria, 
	     string property, double pmin, double pmax, double psteps, double *a, double totalweight)
// bin property of objects from pop within leaf into a; make them have a total weight totalweight
{
     vector<Filter> vf;
     for(int j=0; j<criteria.size(); j++) {
       vf.push_back(Filter(
         criteria[j],
         leaf->doublePropertyValue(criteria[j]+"_min"),
         leaf->doublePropertyValue(criteria[j]+"_max")
       ));
     }
     vf.push_back(Filter(property,pmin,pmax));
     ObjectCollection *tmp = pop->filter(vf);    

     double relweight = totalweight/tmp->size();
     
     for(int j=0; j<tmp->size(); j++)
     {
       double t = (*tmp)[j]->doublePropertyValue(property);
       a[int((t-pmin)/psteps)] += relweight; 
     }
}


int main(int argc, char **argv)
{

  cerr << "# (1) reading catalog" << endl;
 
  ifstream f; f.open("../../george/spec_data.txt");
  int keys[] = {2,3,4,5,6,7,8,9,10,11};
  string names[] = {"redshift", "deep_r", "deep_i", "deep_z", "r", "i", "z", "err_r", "err_i", "err_z"};
  Type types[] = {doubleType, doubleType, doubleType, doubleType, doubleType, doubleType, doubleType, doubleType, doubleType, doubleType}; 


  ObjectCollection cat(f, keys, names, types, 10);

  vector<string> filters;
  filters.push_back("i");
  filters.push_back("r");
  filters.push_back("z");
  int nfilters=filters.size();

  vector<double> limits;
  for(int i=0; i<nfilters; i++)
	  limits.push_back(zp);

  vector<double> minwidth;
  vector<string> bands;  
  vector<string> criteria;
  
  for(int i=0; i<nfilters; i++)
  {
    bands.push_back(filters[i]);
    minwidth.push_back(mindeltam);
    if(i==0)
      criteria.push_back(filters[i]);
    else
      criteria.push_back(filters[i]+filters[0]);
    cat.transformColumn(filters[i],flux_to_croppedmag);
  }

  cat.transformColumnNew("redshiftsq",FilterFunctions::pow2,"redshift");

  cout << cat[0] << endl;


  // (4) build tree
  cerr << "# (2) building tree" << endl;

  ObjectCollection *t = tree(&cat, criteria, limits, minwidth, nmin, 
                             "redshift", "redshiftsq", 
                             brightcut, zp);
  
  
  // (5) getting averages in leaves
  cerr << "# (3) getting mean values in leaves" << endl;
  
  cat.treeAverage(t,criteria,"redshift");
  cat.treeAverage(t,criteria,"redshiftsq");
  
  cerr << "(6) writing output" << endl;

  // first extension: general information about the tree
  ObjectCollection *fields = new ObjectCollection();
  fields->createIntPropertyIfNecessary("BETATREE_VERSION",3,true);
  fields->createIntPropertyIfNecessary("BETATREE_NMIN",nmin,true);
  fields->createIntPropertyIfNecessary("BETATREE_NLEAVES",t->size(),true);
  fields->createDoublePropertyIfNecessary("BETATREE_MAXDEPTH",zp,true);
  
  for(int i=0; i<limits.size(); i++)
  	fields->createDoublePropertyIfNecessary("BETATREE_DEPTH_"+filters[i],limits[i],true);
  	
  t->createIntPropertyIfNecessary("leafid");
  for(int i=0; i<t->size(); i++)
    (*t)[i]->setIntProperty("leafid",i);
  
  
  // prepare p(z)
  const double pzmin=0.005;
  const double pzmax=4.005;
  const int npzbins=400;
  const double dpz=(pzmax-pzmin)/npzbins; // 0.01
  double *pz = new double[npzbins];
  for(int j=0; j<npzbins; j++) pz[j]=0.;
  
  for(int i=0; i<t->size(); i++)
  {    
      leafBin((*t)[i], ref, criteria, "z_phot", pzmin, pzmax, dpz, pz, nobj*beta);
  }
 


  cerr << "(6) writing output" << endl;

  // first extension: general information about the tree
  ObjectPrototype *prototype = new ObjectPrototype;
  ObjectCollection *fields = new ObjectCollection(prototype);
  Object *o = new Object(prototype);
  fields->appendObject(o);
  fields->createIntPropertyIfNecessary("BETATREE_VERSION",3,true);
  fields->createIntPropertyIfNecessary("BETATREE_NMIN",nmin,true);
  fields->createIntPropertyIfNecessary("BETATREE_NLEAVES",t->size(),true);
  fields->createDoublePropertyIfNecessary("BETATREE_BRIGHTCUT",brightcut,true);
  fields->createDoublePropertyIfNecessary("BETATREE_MAXDEPTH",zp,true);
  
  for(int i=0; i<limits.size(); i++)
  	fields->createDoublePropertyIfNecessary("BETATREE_DEPTH_"+filters[i],limits[i],true);
  	
  t->createIntPropertyIfNecessary("leafid");
  for(int i=0; i<t->size(); i++)
    (*t)[i]->setIntProperty("leafid",i);
  
  // prepare p(z)
  const double pzmin=0.005;
  const double pzmax=4.005;
  const int npzbins=400;
  const double dpz=(pzmax-pzmin)/npzbins; // 0.01
  double *pz = new double[npzbins];
  
  valarray<double> pofz(npzbins*t->size());

  vector<double> zmean(t->size());
  vector<double> sigmaz(t->size());

  for(int i=0; i<npzbins*t->size(); i++) pofz[i]=0.;
  
  for(int i=0; i<t->size(); i++)
  {    
      int nobj = (*t)[i]->intPropertyValue("nobj");

      if(nobj<1) continue;
      
      double lredshift = (*t)[i]->doublePropertyValue("redshift_mean");
      double lredshiftsq = (*t)[i]->doublePropertyValue("redshiftsq_mean");
     
      zmean[i]  = lredshift;
      sigmaz[i] = lredshiftsq-lredshift*lredshift;

      for(int j=0; j<npzbins; j++) pz[j]=0.;
      
      leafBin((*t)[i], &cat, criteria, "redshift", pzmin, pzmax, dpz, pz, 1./dpz);
      
      for(int j=0; j<npzbins; j++) {
        pofz[i*npzbins+j]=pz[j];
      }
  }
  
  
  fields->createDoublePropertyIfNecessary("BETATREE_POFZ_FIRSTBIN_ZMIN",pzmin,true);  
  fields->createDoublePropertyIfNecessary("BETATREE_POFZ_FIRSTBIN_ZMAX",pzmin+dpz,true);
  fields->createDoublePropertyIfNecessary("BETATREE_POFZ_NBINS",npzbins,true); 
  
  fields->writeToFITS(outputfile, fields->prototype->doublePropertyName, fields->prototype->intPropertyName, "TREE");
  t->writeToFITS(outputfile, t->prototype->doublePropertyName, t->prototype->intPropertyName, "LEAVES");
  
  auto_ptr<FITS> pFits(0);
    
  try
  {    
      pFits.reset( new FITS(outputfile,Write) );
  }
  catch (CCfits::FITS::CantOpen)
  {
      cerr << "cannot open" << outputfile << " for writing!" << endl;
      return 1;       
  }
  
  vector<string> colName(4,""); colName[0]="leafid"; colName[1]="zmean"; colName[2]="sigmaz"; colName[3]="pofz";
  vector<string> colForm(4,""); colForm[0]="1J";  colForm[1]="1D";  colForm[2]="1D"; colForm[3]=FilterFunctions::NumberToString(npzbins)+"D";
  vector<string> colUnit(4,"");
  Table* newTable = pFits->addTable("POFZ_LEAVES",t->size(),colName,colForm,colUnit);
  
  vector<int> leafid(t->size());
  for(int i=0; i<t->size(); i++) leafid[i]=i;
  newTable->column(colName[0]).write(leafid,1);
  newTable->column(colName[1]).write(zmean,1);
  newTable->column(colName[2]).write(sigmaz,1);
  newTable->column(colName[3]).write(pofz,t->size(),1);
 
  return 0;
}
