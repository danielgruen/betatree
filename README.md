# betatree
Decision tree pipeline for p(z) estimation from CFHT/Megacam photometry [arXiv:1610.01160]

## Installation

Clone this repository with
`git clone git@github.com:danielgruen/betatree.git`
then
`cd betatree`

You likely have to edit the `Makefile` to match the path to your installation of the 
CCfits library.

Then call
`make`
to build the `betatree` executable and download the required reference catalogs.

## Usage

Call
`./betatree [zcluster] [output tree file] [maximum i band depth] [bands_keys band_depths]`
to build a color-magnitude tree and write it to file, e.g.
`./betatree        0.5          tree.fits                     25    g 25.5 r 25.5 z 26` 
to generate a tree optimized for lenses at z=0.5, for a magnitude limited sample of 
galaxies with 20<i<25, using color information from i band data and g, r, z band data
down to magnitudes 25.5, 25.5, 26. We recommend choosing the depths such that your
photometric signal-to-noise ratio is at least 10 for galaxy brighter than these limits.
The resulting tree can be used for galaxies with 20<i<25 and any g, r, z flux.

Photometric information from the reference catalogs is complete and valid down to the 
magnitude limits described in arXiv:1610.01160. All magnitudes are MAG_AUTO (approximately
equal to the total magnitude) in the AB system for CFHT/Megacam filters, with the old i
band. The new i band is available and called 'y'.

## Output

The FITS file output contains three extensions.

### Tree information

* `BETATREE_BRIGHTCUT`,`BETATREE_MAXDEPTH`: lower and upper limit in detection band (=i) magnitude for which tree is valid
* `BETATREE_DEPTH_*`: depth in other bands up to which splits are made in the tree; should be at a SNR>=10 in your data
* `BETATREE_Z_LENS`: lens redshift for which the tree was optimized
* `BETATREE_BETAMEAN_DEEP`: lensing-weighted mean `\beta=D_ds/D_s` for the CFHTLS-DEEP magnitude-limited sample, estimated from the tree with 8-band template-fit photo-z
* `BETATREE_BETAMEAN_D2`,`BETATREE_BETAMEAN_C2015`: lensing-weighted mean `\beta` from D2 (COSMOS) only and from the Laigle+2015 COSMOS photo-z, indicative of overall biases; 0 if i band limit is above 24.7, where C2015 is not quite complete
* `BETATREE_RELATIVE_COSMIC_STDEV`: relative standard deviation `\sigma_\beta/\beta` of the lensing-weighted mean `\beta` estimated from variance of the four DEEP fields
* `BETATREE_POFZ_*`: description of the grid of the p(z) (in equal steps in z) saved in the third extension
* `BETATREE_NLEAVES`: number of leaves in tree
* `BETATREE_NMIN`: minimum number of reference galaxies per leaf
* `BETATREE_VERSION`: version information of the `betatree` code used 

### Leaf information

For each leaf, a line with
* `leafid`: running index
* `nobj`: number of objects in the CFHTLS-DEEP magnitude limited sample in leaf
* `i_min/max`: i band range of leaf
* `?i_min/max`: range of all colors used in this tree of leaf
* `beta`,`sigmabeta`: mean value and uncertainty of mean of `\beta` of all CFHTLS galaxies in leaf
 

### p(z)
