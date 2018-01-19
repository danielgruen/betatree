// buzzardlib.h: code snippets for reading Buzzard catalogs

#ifndef _BUZZARDLIB_
#define _BUZZARDLIB_

#include <cmath>
#include <cassert>
#include <time.h>
#include "filter.h"
#include "cosmo/cosmogeometry.h"
#include "cfhtlib.h"

namespace buzzardlib {


ObjectCollection* read_refcat_prepared_buzzard(vector<string> bands, 
                                       string zstring, 
                                       vector<string> extracolumns=vector<string>(0), 
                                       string filename="/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd/Chinchilla-0/addgalspostprocess/truth_rotated_Y3/Chinchilla-0Y3_v1.6_truth.101.fits") { // a tile with truth information

  srand (time(NULL));

  vector<string> photozcolumns;
  photozcolumns.push_back("Z");
  photozcolumns.push_back("TMAG");
  photozcolumns.push_back("SIZE");
  
  for(int i=0; i<extracolumns.size(); i++)
  {
    photozcolumns.push_back(extracolumns[i]);
  }
  
  ObjectCollection *ref;
  ref = new ObjectCollection(filename,1,photozcolumns);
 
  ref->prototype->doublePropertyName[ref->prototype->doubleVkey("TMAG1")]="gr";
  ref->prototype->doublePropertyName[ref->prototype->doubleVkey("TMAG2")]='r';
  ref->prototype->doublePropertyName[ref->prototype->doubleVkey("TMAG3")]="ir";
  ref->prototype->doublePropertyName[ref->prototype->doubleVkey("TMAG4")]="zr";

  // CPD cuts
  double maglim_r = 26.; // deep-field sample
  double psf_fwhm_r = 0.9 / 0.26;

  ref = ref->filter(Filter("r",0,-2.5*log10(0.5)+maglim_r));
  ref = ref->filter(Filter("SIZE",0.25*0.26*0.5*psf_fwhm_r,1.e30));

  ref->transformColumnNew("Z2",FilterFunctions::pow2,"Z");

  ref->transformColumn("gr",cfhtlib::mags_to_color,"r");
  ref->transformColumn("ir",cfhtlib::mags_to_color,"r");
  ref->transformColumn("zr",cfhtlib::mags_to_color,"r");

  return ref;
}

};

#endif
