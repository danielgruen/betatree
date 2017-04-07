// cfhtlib.h: code snippets for reading CFHT reference catalogs

#include <cmath>
#include <cassert>
#include <time.h>
#include "filter.h"
#include "cosmo/cosmogeometry.h"

namespace cfhtlib {

const string refcatpath=REFCATPATH;

const string numbers[]={"0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9", 
                       "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"};
// transformation functions

double offset=0.;
double seeing_flux_radius_pix=0.;
double zlens=0.;
double zwidth=0.;

double mag_offset(double mag) { 
// shift magnitude by offset
  if(50>mag && mag>-50) return mag+offset; // don't shift mag=99
  return mag;
}

double mags_to_color(double mag1, double mag2) { 
  if(50<mag1 || mag1<-50) return mag1;
  if(50<mag2 || mag2<-50) return mag2;
  
  return mag1-mag2; // don't subtract mag=99
}

double mag_auto(double mag_aper, double mag_auto_ref, double mag_aper_ref) {
// transform aper to auto magnitude using detection auto and measurement aper magnitude offset of object
  if(50>mag_aper && mag_aper>-50)  return mag_aper+(mag_auto_ref-mag_aper_ref); // don't shift mag=99
  return mag_aper;
}

double flux_radius_preseeing(double flux_radius) {
// subtract PSF flux radius in quadrature
  if(flux_radius<=seeing_flux_radius_pix) return 0.;
  return sqrt(flux_radius*flux_radius-seeing_flux_radius_pix*seeing_flux_radius_pix);
}

double cluster_slice(double zs) {
  if(fabs(zs-zlens)<zwidth*(1.+zlens)) return 1;
  return 0;
}

double foreground_slice(double zs) {
  if(zs<zlens && !cluster_slice(zs)) return 1;
  return 0;
}

double background_slice(double zs) {
  if(zs>zlens && !cluster_slice(zs)) return 1;
  return 0;
}

double angularDiameterDistance(double z1, double z2, int num=1000) {
   const double omegaM=0.3;
   const double omegaL=1.-omegaM;
   const double ckms = 2.99792458e5; // c in units of km/s
   const double h = 0.7;

   double sum = 0;
   double R2 = 1.0/(1.0+z1);
   double R1 = 1.0/(1.0+z2);
   double dR = (R2-R1)/((double)num);
   double R = R1;
   for (int i = 0 ; i < num ; i++, R += dR)
   {
      double term1 = omegaL*(R*R*R*R); // constant
      double term2 = omegaM*R;   // diluted matter

      double val1 = 1.0/sqrt(term1+term2);

      term1 = omegaL*(R+dR)*(R+dR)*(R+dR)*(R+dR);
      term2 = omegaM*(R+dR); 

      double val2 = 1.0/sqrt(term1+term2);

      sum += ((val1+val2)/2.0)*dR; // trapezium rule
      // note that da=a^2 dz; we multiply the a^2 into the square root terms
   }

   double result = sum*ckms/100./h/(1.0+z2); // 3000 MPc h^-1
   return result; // Mpc / rad  
}

double DdsDs(double zs) {
  if(zs<=zlens) return 0;
  return angularDiameterDistance(zlens,zs,1000)/angularDiameterDistance(0,zs,1000);
}

// methods you would actually want to use

ObjectCollection* read_refcat_prepared(vector<string> bands, 
                                       string zstring, 
                                       vector<string> extracolumns=vector<string>(0), 
                                       string filename="D.cat") {

  srand (time(NULL));

  vector<string> photozcolumns;
  photozcolumns.push_back("z_phot");
  photozcolumns.push_back("SED");
  //photozcolumns.push_back("FLUX_RADIUS_PRESEEING");
  photozcolumns.push_back("FIELD");
  for(int i=0; i<bands.size(); i++)
  {
    if(i==0)
      photozcolumns.push_back(bands[i]);
    else
      photozcolumns.push_back(bands[i]+bands[0]);
    photozcolumns.push_back(bands[i]+"err");
  }
  
  photozcolumns.push_back("DdsDs"+zstring);
  photozcolumns.push_back("DdsDssq"+zstring);
  photozcolumns.push_back("background"+zstring);
  photozcolumns.push_back("cluster"+zstring);
  
  for(int i=0; i<extracolumns.size(); i++)
  {
    photozcolumns.push_back(extracolumns[i]);
  }
  
  ObjectCollection *ref;
  ref = new ObjectCollection(refcatpath+filename,"OBJECTS",photozcolumns);
  
  return ref;
}

double _gain=0.;
double subtractshotnoise(double fluxerr, double flux)
{
   // fluxerr = sqrt(bgerr^2+flux/_gain)
   // fluxerr' = sqrt(max(0,fluxerr^2-flux/_gain))
   return sqrt(FilterFunctions::max(0,fluxerr*fluxerr-flux/_gain));
}
double addshotnoise(double fluxerr, double flux)
{
   // fluxerr' = sqrt(fluxerr^2+flux/_gain)
   if(fluxerr<=0) return fluxerr;
   return sqrt(fluxerr*fluxerr+flux/_gain);
}

ObjectCollection* filter_fluxerr(ObjectCollection *c, string aperture, vector<string> bands, 
vector<double> limits, const double *gain_at_zp30)
// demand that flux error is small enough to give S/N of 10
// at u=26, g=26.5, r=26.5, i=26, y=26, z=25.5, J=19, H=19, K=19
// note: Fabrice has re-scaled fluxes and flux errors from the original to match 
// a ZP of 31.4 instead of 30, which means they are not actual counts in the FITS images
// but re-scaled by a factor of 10^(0.4*1.4)=3.63
{
  vector<double> gaineff;
  double epsilon=1.e-9; 
  
  //FLUXERR_APER2=sqrt(A*sig^2+FLUX_APER2/GAIN)
  //GAIN = GAIN_AT_ZP_30/(10**(0.4*(ZP-30)))
  
  vector<Filter> vf;

  for (int i=0; i<bands.size(); i++)
  {
     // (1) pick an object with good flux measurement to determine ZP
     int o=0;
     while((*c)[o]->doublePropertyValue("FLUX_APER_"+bands[i]+aperture)<=0)
       o++;
       
     // (2) determine ZP
     double f=(*c)[o]->doublePropertyValue("FLUX_APER_"+bands[i]+aperture);
     double m=(*c)[o]->doublePropertyValue("MAG_APER_"+bands[i]+aperture);
     // m=-2.5lg(f)+ZP -> ZP=m+2.5lg(f)
     double zp=m+2.5*log10(f);
     cerr << "# ZP in " << bands[i] << " is " << zp << " because f=" << f << " and m=" << m << endl;
     _gain=gain_at_zp30[i]/(pow(10,0.4*(zp-30.)));
     gaineff.push_back(_gain);
     cerr << "# ZP corrected gain is " << _gain << endl;
     
     // (3) determine flux of object at given limiting magnitude
     // f=10**[(ZP-m)/2.5]
     double flim= pow(10,(zp-limits[i])/2.5);
     double ferrmax=flim/10.; // S/N of 10
     
     // (4) correct FLUXERR_APER for shot noise
     c->transformColumn("FLUXERR_APER_"+bands[i]+aperture,subtractshotnoise,"FLUX_APER_"+bands[i]+aperture);
     
     cerr << "# flim/10 is " << ferrmax << " when median ferr is "
          << c->median("FLUXERR_APER_"+bands[i]+aperture) << " and upper 10th percentile ferr is " 
          << c->percentile("FLUXERR_APER_"+bands[i]+aperture,0.9) << endl;
          
     vf.push_back(Filter("FLUXERR_APER_"+bands[i]+aperture,epsilon,ferrmax));
  }
  
  ObjectCollection *tmp = c->filter(vf);
  
  // add shot noise back on
  for (int i=0; i<bands.size(); i++)
  {
     _gain = gaineff[i];
     c->transformColumn("FLUXERR_APER_"+bands[i]+aperture,addshotnoise,"FLUX_APER_"+bands[i]+aperture);
  }
  
  return tmp;

}


ObjectCollection* read_refcat_C30(vector<string> bands, string aperture="2") {
    
  vector<string> refbands; vector<double> limits;

  const double gain1[] = {1026., 7095., 6919., 6807., 1882., 2148., 1116., 2379., 2134.};
  const double gain2[] = {2288., 5244., 5418., 5310., 2926., 1455.,  193.,  377.,  175.};
  const double gain3[] = { 972., 5830., 5759., 6359., 2445., 1824., 1517., 1532., 1378.};
  const double gain4[] = { 899., 6999., 6329., 6378., 1500., 1805., 1577., 2413., 1702.};  
  refbands.push_back("u"); limits.push_back(26.5);
  refbands.push_back("g"); limits.push_back(26.5);
  refbands.push_back("r"); limits.push_back(26.5);
  refbands.push_back("i"); limits.push_back(26.);
  refbands.push_back("y"); limits.push_back(25.5);
  refbands.push_back("z"); limits.push_back(25.5);
  
  vector<string> photozcolumns;
  photozcolumns.push_back("flag_photoz");
  photozcolumns.push_back("z_phot");
  photozcolumns.push_back("z_phot_C30");
  photozcolumns.push_back("MASK");
  photozcolumns.push_back("MAG_AUTO");    // i band in CFHT deep
  for(int i=0; i<refbands.size(); i++)
  {
    photozcolumns.push_back("MAG_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("FLUX_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("MAGERR_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("FLUXERR_APER_"+refbands[i]+aperture);
  }
  photozcolumns.push_back("Xpos"); // for matching later
  photozcolumns.push_back("Ypos");
  
  vector<Filter> goodzfilter;

  ObjectCollection *refD2 = new ObjectCollection(refcatpath+"D2.C30."+bands[0]+".cat", "PSSC", photozcolumns);
  double extinction_i=0.032;
    
  offset= -0.07-extinction_i;
  refD2->transformColumn("MAG_APER_u"+aperture,mag_offset);
  offset= -0.06-extinction_i;
  refD2->transformColumn("MAG_APER_g"+aperture,mag_offset);
  offset= -0.04-extinction_i;
  refD2->transformColumn("MAG_APER_r"+aperture,mag_offset);
  offset=  0.00-extinction_i;
  refD2->transformColumn("MAG_APER_i"+aperture,mag_offset);
  if(bands[0]=="i") refD2->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.01-extinction_i;
  refD2->transformColumn("MAG_APER_y"+aperture,mag_offset);
  if(bands[0]=="y") refD2->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.00-extinction_i;
  refD2->transformColumn("MAG_APER_z"+aperture,mag_offset);
    
  ObjectCollection *refD2f = filter_fluxerr(refD2,aperture,refbands,limits,gain2);
  refD2f->createIntPropertyIfNecessary("FIELD",2);
  
  cerr << "# " << refD2f->size() << " of " << refD2->size() << " objects left in D2 after fluxerr clipping" << endl;
      
  goodzfilter.push_back(Filter("flag_photoz",0,3));
  goodzfilter.push_back(Filter("MASK",0,0));
  goodzfilter.push_back(Filter("MAG_AUTO",20,24.7)); // CODEX cut means no objects fainter than 24.7 will be used
  goodzfilter.push_back(Filter("MAG_APER_"+bands[0]+aperture,-50,50));

  ObjectCollection *goodref = refD2f->filter(goodzfilter);
  
  cerr << "# good reference objects: " << goodref->size() << endl;
  
  // correct AUTO magnitudes from aperture magnitudes
  for(int i=0; i<refbands.size(); i++) {
    goodref->transformColumnNew(refbands[i],mag_auto,"MAG_APER_"+refbands[i]+aperture,"MAG_AUTO","MAG_APER_"+bands[0]+aperture);
    goodref->prototype->doublePropertyName[goodref->prototype->doubleVkey("MAGERR_APER_"+refbands[i]+aperture)]=refbands[i]+"err"; 
  }
  
  for(int i=1; i<bands.size(); i++)
  {
    goodref->transformColumnNew(bands[i]+bands[0],mags_to_color,bands[i],bands[0]);
  }
  
  return goodref;

}


ObjectCollection* read_refcat_C2015(vector<string> bands, string aperture="2") {
    
  vector<string> refbands; vector<double> limits;

  const double gain1[] = {1026., 7095., 6919., 6807., 1882., 2148., 1116., 2379., 2134.};
  const double gain2[] = {2288., 5244., 5418., 5310., 2926., 1455.,  193.,  377.,  175.};
  const double gain3[] = { 972., 5830., 5759., 6359., 2445., 1824., 1517., 1532., 1378.};
  const double gain4[] = { 899., 6999., 6329., 6378., 1500., 1805., 1577., 2413., 1702.};  
  refbands.push_back("u"); limits.push_back(26.5);
  refbands.push_back("g"); limits.push_back(26.5);
  refbands.push_back("r"); limits.push_back(26.5);
  refbands.push_back("i"); limits.push_back(26.);
  refbands.push_back("y"); limits.push_back(25.5);
  refbands.push_back("z"); limits.push_back(25.5);
  
  vector<string> photozcolumns;
  photozcolumns.push_back("flag_photoz");
  photozcolumns.push_back("z_phot");
  photozcolumns.push_back("z_phot_C2015");
  photozcolumns.push_back("MASK");
  photozcolumns.push_back("MAG_AUTO");    // i band in CFHT deep
  for(int i=0; i<refbands.size(); i++)
  {
    photozcolumns.push_back("MAG_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("FLUX_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("MAGERR_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("FLUXERR_APER_"+refbands[i]+aperture);
  }
  photozcolumns.push_back("Xpos"); // for matching later
  photozcolumns.push_back("Ypos");
  
  vector<Filter> goodzfilter;

  ObjectCollection *refD2 = new ObjectCollection(refcatpath+"D2.C2015."+bands[0]+".cat", "PSSC", photozcolumns);
  double extinction_i=0.032;
    
  offset= -0.07-extinction_i;
  refD2->transformColumn("MAG_APER_u"+aperture,mag_offset);
  offset= -0.06-extinction_i;
  refD2->transformColumn("MAG_APER_g"+aperture,mag_offset);
  offset= -0.04-extinction_i;
  refD2->transformColumn("MAG_APER_r"+aperture,mag_offset);
  offset=  0.00-extinction_i;
  refD2->transformColumn("MAG_APER_i"+aperture,mag_offset);
  if(bands[0]=="i") refD2->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.01-extinction_i;
  refD2->transformColumn("MAG_APER_y"+aperture,mag_offset);
  if(bands[0]=="y") refD2->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.00-extinction_i;
  refD2->transformColumn("MAG_APER_z"+aperture,mag_offset);
    
  ObjectCollection *refD2f = filter_fluxerr(refD2,aperture,refbands,limits,gain2);
  refD2f->createIntPropertyIfNecessary("FIELD",2);
  
  cerr << "# " << refD2f->size() << " of " << refD2->size() << " objects left in D2 after fluxerr clipping" << endl;
      
  goodzfilter.push_back(Filter("flag_photoz",0,3));
  goodzfilter.push_back(Filter("MASK",0,0));
  goodzfilter.push_back(Filter("MAG_AUTO",20,24.7)); // CODEX cut means no objects fainter than 24.7 will be used
  goodzfilter.push_back(Filter("MAG_APER_"+bands[0]+aperture,-50,50));

  ObjectCollection *goodref = refD2f->filter(goodzfilter);
  
  cerr << "# good reference objects: " << goodref->size() << endl;
  
  // correct AUTO magnitudes from aperture magnitudes
  for(int i=0; i<refbands.size(); i++) {
    goodref->transformColumnNew(refbands[i],mag_auto,"MAG_APER_"+refbands[i]+aperture,"MAG_AUTO","MAG_APER_"+bands[0]+aperture);
    goodref->prototype->doublePropertyName[goodref->prototype->doubleVkey("MAGERR_APER_"+refbands[i]+aperture)]=refbands[i]+"err"; 
  }
  
  for(int i=1; i<bands.size(); i++)
  {
    goodref->transformColumnNew(bands[i]+bands[0],mags_to_color,bands[i],bands[0]);
  }
  
  return goodref;

}



ObjectCollection* read_refcat(vector<string> bands, string aperture="2", string midfix="") {
// bands:    filter bands to be used for colors later, first one indicates detection bands
// aperture: aperture to be used
    
  vector<string> refbands; vector<double> limits;

  const double gain1[] = {1026., 7095., 6919., 6807., 1882., 2148., 1116., 2379., 2134.};
  const double gain2[] = {2288., 5244., 5418., 5310., 2926., 1455.,  193.,  377.,  175.};
  const double gain3[] = { 972., 5830., 5759., 6359., 2445., 1824., 1517., 1532., 1378.};
  const double gain4[] = { 899., 6999., 6329., 6378., 1500., 1805., 1577., 2413., 1702.};  
  refbands.push_back("u"); limits.push_back(26.5);
  refbands.push_back("g"); limits.push_back(26.5);
  refbands.push_back("r"); limits.push_back(26.5);
  refbands.push_back("i"); limits.push_back(26.);
  refbands.push_back("y"); limits.push_back(25.5);
  refbands.push_back("z"); limits.push_back(25.5);
  refbands.push_back("J"); limits.push_back(20.);
  refbands.push_back("H"); limits.push_back(19.5);
  refbands.push_back("K"); limits.push_back(19.);
  
  vector<string> photozcolumns;
  photozcolumns.push_back("flag_photoz");
  photozcolumns.push_back("z_phot");
  photozcolumns.push_back("SED");
  photozcolumns.push_back("MASK");
  photozcolumns.push_back("MASK_IR");
  photozcolumns.push_back("MAG_AUTO");    // i band in CFHT deep
  photozcolumns.push_back("FLUX_RADIUS");
  photozcolumns.push_back("NPIX");
  for(int i=0; i<refbands.size(); i++)
  {
    photozcolumns.push_back("MAG_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("FLUX_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("MAGERR_APER_"+refbands[i]+aperture);
    photozcolumns.push_back("FLUXERR_APER_"+refbands[i]+aperture);
  }
  photozcolumns.push_back("Xpos"); // for matching later
  photozcolumns.push_back("Ypos");
  photozcolumns.push_back("ALPHA_J2000"); 
  photozcolumns.push_back("DELTA_J2000"); 
 
  vector<Filter> goodzfilter;

  cerr << "# reading CFHT DEEP catalogs one by one and applying photometric offsets" << endl;

  
  ObjectCollection *refD1 = new ObjectCollection(refcatpath+"D1."+bands[0]+midfix+".photoz.cat", "OBJECTS", photozcolumns); // reference catalog, no size cut

  double extinction_i=0.045;
  offset= -0.07-extinction_i;
  refD1->transformColumn("MAG_APER_u"+aperture,mag_offset);
  offset= -0.11-extinction_i;
  refD1->transformColumn("MAG_APER_g"+aperture,mag_offset);
  offset= -0.07-extinction_i;
  refD1->transformColumn("MAG_APER_r"+aperture,mag_offset);
  offset=  0.00-extinction_i;
  refD1->transformColumn("MAG_APER_i"+aperture,mag_offset);
  if(bands[0]=="i") refD1->transformColumn("MAG_AUTO",mag_offset);
  offset= -0.01-extinction_i;
  refD1->transformColumn("MAG_APER_y"+aperture,mag_offset);
  if(bands[0]=="y") refD1->transformColumn("MAG_AUTO",mag_offset);
  offset= -0.02-extinction_i;
  refD1->transformColumn("MAG_APER_z"+aperture,mag_offset);
  offset=  0.17-extinction_i;
  refD1->transformColumn("MAG_APER_J"+aperture,mag_offset);
  offset=  0.28-extinction_i;
  refD1->transformColumn("MAG_APER_H"+aperture,mag_offset);
  offset=  0.25-extinction_i;
  refD1->transformColumn("MAG_APER_K"+aperture,mag_offset);
  
  seeing_flux_radius_pix=2.27;
  refD1->transformColumnNew("FLUX_RADIUS_PRESEEING",flux_radius_preseeing,"FLUX_RADIUS");
  
  ObjectCollection *ref = filter_fluxerr(refD1,aperture,refbands,limits,gain1);
  ref->createIntPropertyIfNecessary("FIELD",1);
  cerr << "# " << ref->size() << " objects left in D1 after fluxerr clipping" << endl;
  
  
  

  ObjectCollection *refD2 = new ObjectCollection(refcatpath+"D2."+bands[0]+midfix+".photoz.cat", "OBJECTS", photozcolumns);
  extinction_i=0.032;
    
  offset= -0.07-extinction_i;
  refD2->transformColumn("MAG_APER_u"+aperture,mag_offset);
  offset= -0.06-extinction_i;
  refD2->transformColumn("MAG_APER_g"+aperture,mag_offset);
  offset= -0.04-extinction_i;
  refD2->transformColumn("MAG_APER_r"+aperture,mag_offset);
  offset=  0.00-extinction_i;
  refD2->transformColumn("MAG_APER_i"+aperture,mag_offset);
  if(bands[0]=="i") refD2->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.01-extinction_i;
  refD2->transformColumn("MAG_APER_y"+aperture,mag_offset);
  if(bands[0]=="y") refD2->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.00-extinction_i;
  refD2->transformColumn("MAG_APER_z"+aperture,mag_offset);
  offset=  0.16-extinction_i;
  refD2->transformColumn("MAG_APER_J"+aperture,mag_offset);
  offset=  0.23-extinction_i;
  refD2->transformColumn("MAG_APER_H"+aperture,mag_offset);
  offset=  0.25-extinction_i;
  refD2->transformColumn("MAG_APER_K"+aperture,mag_offset);
    
  seeing_flux_radius_pix=2.37;
  refD2->transformColumnNew("FLUX_RADIUS_PRESEEING",flux_radius_preseeing,"FLUX_RADIUS");
  ObjectCollection *refD2f = filter_fluxerr(refD2,aperture,refbands,limits,gain2);
  refD2f->createIntPropertyIfNecessary("FIELD",2);
  ref->appendCollectionDeep(refD2f);
  
  cerr << "# " << refD2f->size() << " of " << refD2->size() << " objects left in D2 after fluxerr clipping" << endl;

  delete refD2f; delete refD2;
  
  
  
  
  
  ObjectCollection *refD3 = new ObjectCollection(refcatpath+"D3."+bands[0]+midfix+".photoz.cat", "OBJECTS", photozcolumns);
  extinction_i=0.014;
  
  offset= -0.03-extinction_i;
  refD3->transformColumn("MAG_APER_u"+aperture,mag_offset);
  offset= -0.05-extinction_i;
  refD3->transformColumn("MAG_APER_g"+aperture,mag_offset);
  offset= -0.04-extinction_i;
  refD3->transformColumn("MAG_APER_r"+aperture,mag_offset);
  offset=  0.00-extinction_i;
  refD3->transformColumn("MAG_APER_i"+aperture,mag_offset);
  if(bands[0]=="i") refD3->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.00-extinction_i;
  refD3->transformColumn("MAG_APER_y"+aperture,mag_offset);
  if(bands[0]=="y") refD3->transformColumn("MAG_AUTO",mag_offset);
  offset=  0.01-extinction_i;
  refD3->transformColumn("MAG_APER_z"+aperture,mag_offset);
  offset=  0.18-extinction_i;
  refD3->transformColumn("MAG_APER_J"+aperture,mag_offset);
  offset=  0.26-extinction_i;
  refD3->transformColumn("MAG_APER_H"+aperture,mag_offset);
  offset=  0.31-extinction_i;
  refD3->transformColumn("MAG_APER_K"+aperture,mag_offset);
  
  seeing_flux_radius_pix=2.37;
  refD3->transformColumnNew("FLUX_RADIUS_PRESEEING",flux_radius_preseeing,"FLUX_RADIUS");

  ObjectCollection *refD3f = filter_fluxerr(refD3,aperture,refbands,limits,gain3);
  refD3f->createIntPropertyIfNecessary("FIELD",3);
  ref->appendCollectionDeep(refD3f);
  
  cerr << "# " << refD3f->size() << " of " << refD3->size() << " objects left in D3 after fluxerr clipping" << endl;

  delete refD3f; delete refD3;
  
  
  ObjectCollection *refD4 = new ObjectCollection(refcatpath+"D4."+bands[0]+midfix+".photoz.cat", "OBJECTS", photozcolumns);
  extinction_i=0.042;
  
  offset= -0.02-extinction_i;
  refD4->transformColumn("MAG_APER_u"+aperture,mag_offset);
  offset= -0.09-extinction_i;
  refD4->transformColumn("MAG_APER_g"+aperture,mag_offset);
  offset= -0.05-extinction_i;
  refD4->transformColumn("MAG_APER_r"+aperture,mag_offset);
  offset=  0.00-extinction_i;
  refD4->transformColumn("MAG_APER_i"+aperture,mag_offset);
  if(bands[0]=="i") refD4->transformColumn("MAG_AUTO",mag_offset);
  offset= -0.02-extinction_i;
  refD4->transformColumn("MAG_APER_y"+aperture,mag_offset);
  if(bands[0]=="y") refD4->transformColumn("MAG_AUTO",mag_offset);
  offset= -0.01-extinction_i; 
  refD4->transformColumn("MAG_APER_z"+aperture,mag_offset);
  offset=  0.16-extinction_i;
  refD4->transformColumn("MAG_APER_J"+aperture,mag_offset);
  offset=  0.27-extinction_i;
  refD4->transformColumn("MAG_APER_H"+aperture,mag_offset);
  offset=  0.24-extinction_i;
  refD4->transformColumn("MAG_APER_K"+aperture,mag_offset);
  
  seeing_flux_radius_pix=2.37;
  refD4->transformColumnNew("FLUX_RADIUS_PRESEEING",flux_radius_preseeing,"FLUX_RADIUS");  

  ObjectCollection *refD4f = filter_fluxerr(refD4,aperture,refbands,limits,gain4);
  refD4f->createIntPropertyIfNecessary("FIELD",4);
  ref->appendCollectionDeep(refD4f);
  
  
  cerr << "# " << refD4f->size() << " of " << refD4->size() << " objects left in D4 after fluxerr clipping" << endl;

  delete refD4f; delete refD4;
  
    
    
    
  goodzfilter.push_back(Filter("flag_photoz",0,3));
  goodzfilter.push_back(Filter("MASK",0,0));
  goodzfilter.push_back(Filter("MASK_IR",0,0));
  goodzfilter.push_back(Filter("MAG_AUTO",20,24.7)); // CODEX cut means no objects fainter than 24.7 will be used
  goodzfilter.push_back(Filter("MAG_APER_"+bands[0]+aperture,-50,50));

  ObjectCollection *goodref = ref->filter(goodzfilter);
  
  cerr << "# good reference objects: " << goodref->size() << endl;
  
  // correct AUTO magnitudes from aperture magnitudes
  for(int i=0; i<refbands.size(); i++) {
    goodref->transformColumnNew(refbands[i],mag_auto,"MAG_APER_"+refbands[i]+aperture,"MAG_AUTO","MAG_APER_"+bands[0]+aperture);
    goodref->prototype->doublePropertyName[goodref->prototype->doubleVkey("MAGERR_APER_"+refbands[i]+aperture)]=refbands[i]+"err"; 
  }
  
  for(int i=1; i<bands.size(); i++)
  {
    goodref->transformColumnNew(bands[i]+bands[0],mags_to_color,bands[i],bands[0]);
  }
  
  return goodref;
}  


void assign_cluster_slice(ObjectCollection *refcat, double dzcluster, double dzwidth=0.06) 
{
   zlens=dzcluster; zwidth=dzwidth;
   refcat->transformColumnNew("cluster",cluster_slice,"z_phot");
}

void assign_DdsDs(ObjectCollection *refcat, vector<double> dzlens, const string *zstring, const string zphotcol="z_phot", const string midfix="") {

  for(int i=0; i<dzlens.size(); i++)
  {
    cerr << "# assigning DdsDs for lens redshift " << dzlens[i] << " to column " << "DdsDs"+midfix+zstring[i] << endl;
    zlens  = dzlens[i];
    zwidth = 0.06;
    refcat->transformColumnNew("DdsDs"+midfix+zstring[i],DdsDs,zphotcol);
    refcat->transformColumnNew("DdsDssq"+midfix+zstring[i],FilterFunctions::square,"DdsDs"+midfix+zstring[i]); // <(Dds/Ds)^2>
    refcat->transformColumnNew("background"+midfix+zstring[i],background_slice,zphotcol);
    refcat->transformColumnNew("cluster"+midfix+zstring[i],cluster_slice,zphotcol);
  }
  
}

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




};
