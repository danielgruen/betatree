// betatree: build tree in color-magnitude space with variable bands
// mice variant that does not do cosmic variance and stuff like that

#include <iostream>
using namespace std;

#include "cfhtlib.h"
using namespace cfhtlib;

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
  if(argc < 5 || (argc%2)) {
    cerr << "syntax: " << argv[0] << " [zcluster] [output tree file] [maximum i band depth]"
         << " [bands_keys band_depth other than i...]" << endl;
    return 1;
  }

  const string detectionband = "i";
  const int nfilters = (argc-4)/2+1;
  const string aperture = "2";

  double zcluster = atof(argv[1]);
  assert(zcluster>0); assert(zcluster<1);
 
  string outputfile=argv[2];
  double maxdepth=atof(argv[3]);
  assert(maxdepth>brightcut);
  assert(maxdepth<=25);
  
  string filters[nfilters];
  filters[0]=detectionband;
  vector<double> limits;
  limits.push_back(maxdepth); // limit of detection band
  
  
  for(int i=4; i<argc; i++) {
    filters[(i-4)/2+1] = argv[i];
    i++;
    limits.push_back(atof(argv[i]));
    // S/N=10 at least at u=26, g=26.5, r=26.5, i=26, y=26, z=25.5, J=19, H=19, K=19
    int n=(i-5)/2+1;
    switch(filters[n][0]) {
      case 'u':
      case 'y':
        assert(limits[n]<=26.);
        break;
      case 'g':
      case 'r':
        assert(limits[n]<=26.5);
        break;
      case 'z':
        assert(limits[n]<=25.5);
        break;
      case 'J':
      case 'H':
      case 'K':
        assert(limits[n]<=19);
    	break;
      default:
        cerr << "unknown filter band " << filters[n] << ". exiting" << endl;
        return 1;
    }
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

    cout << "criterion " << i+1 << ": " << criteria[i] << endl;
  }
  ObjectCollection *ref    = read_refcat_prepared(bands, ZeroPadNumber(iz));
  ref = ref->filter(Filter(filters[0],brightcut,maxdepth));

  if(ref->prototype->doubleVkey("DdsDs"+ZeroPadNumber(iz))<0) { // create beta-columns if missing
    string zstring[] = {ZeroPadNumber(iz)};
    assign_DdsDs(ref, vector<double>(1,zcluster), zstring, "z_phot", "");
  }

  // (4) build tree
  cerr << "# (2) building tree" << endl;

  ObjectCollection *t = tree(ref, criteria, limits, minwidth, nmin, 
                             "DdsDs"+ZeroPadNumber(iz), "DdsDssq"+ZeroPadNumber(iz), 
                             brightcut, maxdepth);
  
  
  // (5) getting averages in leaves
  cerr << "# (3) getting mean values in leaves" << endl;
  
  ref->treeAverage(t,criteria,"DdsDs"+ZeroPadNumber(iz));
  ref->treeAverage(t,criteria,"DdsDssq"+ZeroPadNumber(iz));
  ref->treeAverage(t,criteria,"background"+ZeroPadNumber(iz));
  ref->treeAverage(t,criteria,"cluster"+ZeroPadNumber(iz));
  
  // (7) estimating photo-z uncertainty from comparison of D2 with C2015
  //     only if maxdepth <= 24.7
  if(maxdepth<=24.7) {
  	cerr << "# (3b) estimating photo-z uncertainty from comparison of D2 with C2015" << endl;
  	
  	vector<string> extracol; 
  	extracol.push_back("z_phot_C2015"); 
  	extracol.push_back("DdsDs_D2_"+ZeroPadNumber(iz));
  	extracol.push_back("DdsDs_C2015_"+ZeroPadNumber(iz));
  	extracol.push_back("DdsDssq_D2_"+ZeroPadNumber(iz));
  	extracol.push_back("DdsDssq_C2015_"+ZeroPadNumber(iz));
  	ObjectCollection *refC2015 = read_refcat_prepared(bands, ZeroPadNumber(iz), extracol, "D.C2015.i.cat");
  	
  	if(refC2015->prototype->doubleVkey("DdsDs_D2_"+ZeroPadNumber(iz))<0) { // create beta-columns if missing
  	  string zstring[] = {ZeroPadNumber(iz)};
  	  assign_DdsDs(refC2015, vector<double>(1,zcluster), zstring, "z_phot", "_D2_");
  	}  
  	if(refC2015->prototype->doubleVkey("DdsDs_C2015_"+ZeroPadNumber(iz))<0) { // create beta-columns if missing
  	  string zstring[] = {ZeroPadNumber(iz)};
  	  assign_DdsDs(refC2015, vector<double>(1,zcluster), zstring, "z_phot_C2015", "_C2015_");
  	}
  	
  	ObjectCollection *tC2015_C2015  = new ObjectCollection(t->prototype);
  	ObjectCollection *tC2015_D2  = new ObjectCollection(t->prototype);
  	tC2015_C2015->appendCollectionDeep(t);
  	tC2015_D2->appendCollectionDeep(t);
  	refC2015->treeAverage(tC2015_C2015,criteria,"DdsDs_C2015_"+ZeroPadNumber(iz));
  	refC2015->treeAverage(tC2015_D2,criteria,"DdsDs_D2_"+ZeroPadNumber(iz));
  	refC2015->treeAverage(tC2015_C2015,criteria,"DdsDssq_C2015_"+ZeroPadNumber(iz));
  	refC2015->treeAverage(tC2015_D2,criteria,"DdsDssq_D2_"+ZeroPadNumber(iz));
  	
  	tC2015_C2015->transformColumnNew("DdsDs_C2015_"+ZeroPadNumber(iz)+"_mean_uncertainty",mean_uncertainty,
  	                            "DdsDs_C2015_"+ZeroPadNumber(iz)+"_mean","DdsDssq_C2015_"+ZeroPadNumber(iz)+"_mean","nobj");
  	tC2015_D2->transformColumnNew("DdsDs_D2_"+ZeroPadNumber(iz)+"_mean_uncertainty",mean_uncertainty,
  	                            "DdsDs_D2_"+ZeroPadNumber(iz)+"_mean","DdsDssq_D2_"+ZeroPadNumber(iz)+"_mean","nobj");
	
  	t->extendCollection(*tC2015_C2015,"DdsDs_C2015_"+ZeroPadNumber(iz)+"_mean","DdsDs"+ZeroPadNumber(iz)+"_mean_C2015");
  	t->extendCollection(*tC2015_C2015,"DdsDs_C2015_"+ZeroPadNumber(iz)+"_mean_uncertainty","DdsDs"+ZeroPadNumber(iz)+"_mean_uncertainty_C2015");
  	t->extendCollection(*tC2015_D2,"DdsDs_D2_"+ZeroPadNumber(iz)+"_mean","DdsDs"+ZeroPadNumber(iz)+"_mean_D2");
  	t->extendCollection(*tC2015_D2,"DdsDs_D2_"+ZeroPadNumber(iz)+"_mean_uncertainty","DdsDs"+ZeroPadNumber(iz)+"_mean_uncertainty_D2");
  }

  // (4) estimating cosmic variance on weighted mean beta
  cerr << "# (4) estimating cosmic variance" << endl;
  ObjectCollection *train[4];  // single field reference catalog
  ObjectCollection *traint[4]; // tree filled with objects from just a single field
  
#pragma omp parallel for
  for(int f=0; f<4; f++) {
       train[f]=ref->filter(Filter("FIELD",f+0.99,f+1.01));
       traint[f] = new ObjectCollection(t->prototype);
       traint[f]->appendCollectionDeep(t);                      
       train[f]->treeAverage(traint[f],criteria,"DdsDs"+ZeroPadNumber(iz));
  }

  double betamean[4];
  
  t->createDoublePropertyIfNecessary("beta_noD1",-1,true);
  t->createDoublePropertyIfNecessary("beta_noD2",-1,true);
  t->createDoublePropertyIfNecessary("beta_noD3",-1,true);
  t->createDoublePropertyIfNecessary("beta_noD4",-1,true);

#pragma omp parallel for
  for(int f=0; f<4; f++) {
       // mean beta = sum_j(ntest_j betamall_j betamtrain_j)/sum_j(ntest_j betamall_j)
       //   where betamtrain_j = sum_{not f}(ntrain_j betatrain_j)/sum_{not f}(ntrain_j)
       //     and betamall_j   = (nall_j betaall_j)/(nall_j)

       // this assumes we're weighting sources by the same (all-based) weight in each jackknife case
       //  (otherwise we'll get different observables as well)
       
       double nte_bmall_bmtr=0.;
       double nte_bmall=0.;
       
       for(int j=0; j<t->size(); j++)
       {
         double ntr_btr=0.;
         double ntr=0.;
         for(int g=0; g<4; g++) {
           if(g==f) continue;
           if((*(traint[g]))[j]->intPropertyValue("nobj")==0) continue;
           ntr_btr += (*(traint[g]))[j]->intPropertyValue("nobj") *
                      (*(traint[g]))[j]->doublePropertyValue("DdsDs"+ZeroPadNumber(iz)+"_mean");
           ntr     += (*(traint[g]))[j]->intPropertyValue("nobj");
         }
         assert(ntr>0);
         double bmtr = ntr_btr/ntr; // mean value we think objects in box j have, based on jackknifed training sample
         
         double wmall = (*t)[j]->intPropertyValue("nobj") *
                        (*t)[j]->doublePropertyValue("DdsDs"+ZeroPadNumber(iz)+"_mean"); 
                         // density*mean value they do have, used for weighting
         
         nte_bmall_bmtr += wmall*bmtr;
         nte_bmall      += wmall;
         
         (*t)[j]->setDoubleProperty("beta_noD"+ZeroPadNumber(f+1,1),bmtr);
       }
       betamean[f]=nte_bmall_bmtr/nte_bmall;
  }
  
  double sigmabeta = 0;
  for(int f=0; f<4; f++) {
     sigmabeta += FilterFunctions::pow2(betamean[f]-(betamean[0]+betamean[1]+betamean[2]+betamean[3])/4.);
  }
  cerr << "4 beta estimates: " << betamean[0] << " " << betamean[1] << " " << betamean[2] << " " << betamean[3] << endl;
  sigmabeta = sqrt(3./4.*sigmabeta)/((betamean[0]+betamean[1]+betamean[2]+betamean[3])/4.);
  cerr << "relative cosmic variance: " << sigmabeta << endl;
  
  
  
  // (5) correcting beta noise bias  
  cerr << "# (5) correcting beta noise bias" << endl;
  
  t->transformColumnNew("sigmabeta",jackknife_sigma_4,"beta_noD1","beta_noD2","beta_noD3","beta_noD4");
  t->transformColumnNew("beta",correct_beta_for_variance, "DdsDs"+ZeroPadNumber(iz)+"_mean", "sigmabeta");

  
 
  
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
  fields->createDoublePropertyIfNecessary("BETATREE_RELATIVE_COSMIC_STDEV",sigmabeta,true);
  fields->createDoublePropertyIfNecessary("BETATREE_Z_LENS",zcluster,true);
  
  for(int i=0; i<limits.size(); i++)
  	fields->createDoublePropertyIfNecessary("BETATREE_DEPTH_"+filters[i],limits[i],true);
  	
  t->createIntPropertyIfNecessary("leafid");
  for(int i=0; i<t->size(); i++)
    (*t)[i]->setIntProperty("leafid",i);
  
  double betamean_C2015=0.;
  double betamean_DEEP=0.;
  double betamean_D2=0.;
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
      
	  if(maxdepth<24.7) {
          betamean_C2015 += nobj*(*t)[i]->doublePropertyValue("DdsDs"+ZeroPadNumber(iz)+"_mean_C2015");
          betamean_D2    += nobj*(*t)[i]->doublePropertyValue("DdsDs"+ZeroPadNumber(iz)+"_mean_D2");
	  }
	  
      double beta = (*t)[i]->doublePropertyValue("DdsDs"+ZeroPadNumber(iz)+"_mean");
      betamean_DEEP  += nobj*beta;
      nsum += nobj;
      
      for(int j=0; j<npzbins; j++) pz[j]=0.;
      
      leafBin((*t)[i], ref, criteria, "z_phot", pzmin, pzmax, dpz, pz, 1./dpz);
      
      for(int j=0; j<npzbins; j++) {
        pofz[i*npzbins+j]=pz[j];
      }
  }
  
  
  betamean_C2015 /= nsum;
  betamean_D2 /= nsum;
  betamean_DEEP /= nsum;
  fields->createDoublePropertyIfNecessary("BETATREE_BETAMEAN_C2015",betamean_C2015,true);
  fields->createDoublePropertyIfNecessary("BETATREE_BETAMEAN_D2",betamean_D2,true);
  fields->createDoublePropertyIfNecessary("BETATREE_BETAMEAN_DEEP",betamean_DEEP,true); 
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
