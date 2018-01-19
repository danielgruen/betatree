// betatree: build tree in color-magnitude space with variable bands
// mice variant that does not do cosmic variance and stuff like that

#include <iostream>
using namespace std;

#include "cfhtlib.h"
using namespace cfhtlib;
#include "buzzardlib.h"
using namespace buzzardlib;

#include <CCfits/CCfits>
using namespace CCfits;

const int nmin=1000; // minimum number of reference galaxies per leaf
const double mindeltam=0.1; // minimum width of leafs (mag)
const double brightcut=20.; // objects brighter than i=20 ignored from the tree

double betaratio(double beta1, double beta2) { return 2.*(beta1-beta2)/(beta1+beta2); }


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

double jackknife_sigma_4(double j1, double j2, double j3, double j4)
{
  double jmean = (j1+j2+j3+j4)/4.;
  double ssq = FilterFunctions::pow2(j1-jmean)+
               FilterFunctions::pow2(j2-jmean)+
               FilterFunctions::pow2(j3-jmean)+
               FilterFunctions::pow2(j4-jmean);
  return sqrt(ssq*3./4.);
}

double correct_beta_for_variance(double beta, double sigmabeta)
{
   if(beta<=0) return beta;
   return beta/(1.+FilterFunctions::pow2(sigmabeta/beta));
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
  if(argc!=2) {
    cerr << "syntax: " << argv[0] << " [outputfile]" << endl;
    return 0;
  }

  string outputfile=argv[1];

  const string detectionband = "r";
  const int nfilters = 4;
  const double maxdepth=50.; 
  string filters[nfilters];
  filters[0]=detectionband;
  
  vector<double> limits;

  filters[1]='g';
  filters[2]='i';
  filters[3]='z';

  for(int i=0; i<nfilters; i++) {
    limits.push_back(maxdepth);
  }
 
  // (2) read refcat
  cerr << "# (1) reading reference catalog" << endl;
  
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
  }
  ObjectCollection *ref    = read_refcat_prepared_buzzard(bands, "", vector<string>(0), "/nfs/slac/g/ki/ki23/des/dgruen/betatree/short101.fits");
  ref = ref->filter(Filter(filters[0],brightcut,maxdepth));


  // (4) build tree
  cerr << "# (2) building tree" << endl;

  ObjectCollection *t = tree(ref, criteria, limits, minwidth, nmin, 
                             "Z", "Z2", 
                             brightcut, maxdepth);
  
  
  // (5) getting averages in leaves
  cerr << "# (3) getting mean values in leaves" << endl;
  
  ref->treeAverage(t,criteria,"Z");
  ref->treeAverage(t,criteria,"Z2");
  
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
  fields->createDoublePropertyIfNecessary("BETATREE_MAXDEPTH",maxdepth,true);
  
  for(int i=0; i<limits.size(); i++)
  	fields->createDoublePropertyIfNecessary("BETATREE_DEPTH_"+filters[i],limits[i],true);
  	
  t->createIntPropertyIfNecessary("leafid");
  for(int i=0; i<t->size(); i++)
    (*t)[i]->setIntProperty("leafid",i);
  
  int nsum=0;
  
  // prepare p(z)
  const double pzmin=0.005;
  const double pzmax=4.005;
  const int npzbins=400;
  const double dpz=(pzmax-pzmin)/npzbins; // 0.01
  double *pz = new double[npzbins];
  
  valarray<double> pofz(npzbins*t->size());
  for(int i=0; i<npzbins*t->size(); i++) pofz[i]=0.;
  
  for(int i=0; i<t->size(); i++)
  {    
      int nobj = (*t)[i]->intPropertyValue("nobj");

      if(nobj<1) continue;
      
      nsum += nobj;
      
      for(int j=0; j<npzbins; j++) pz[j]=0.;
      
      leafBin((*t)[i], ref, criteria, "Z", pzmin, pzmax, dpz, pz, 1./dpz);
      
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
  
  vector<string> colName(2,""); colName[0]="leafid"; colName[1]="pofz";
  vector<string> colForm(2,""); colForm[0]="1J"; colForm[1]=FilterFunctions::NumberToString(npzbins)+"D";
  vector<string> colUnit(2,"");
  Table* newTable = pFits->addTable("POFZ_LEAVES",t->size(),colName,colForm,colUnit);
  
  vector<int> leafid(t->size());
  for(int i=0; i<t->size(); i++) leafid[i]=i;
  newTable->column(colName[0]).write(leafid,1);
  newTable->column(colName[1]).write(pofz,t->size(),1);
  
  return 0;
}
