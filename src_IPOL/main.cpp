/**
 * @file main.cpp
 *
 * @brief Main function for Landscape Evolution Modelling.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

//! Global includes
#include <omp.h>
#include <time.h>


//! Local includes
#include "LibImages/LibImagesThreads.h" // to remove
#include "LibImages/LibImages.h" // to remove
#include "Utilities/Utilities.h" // to remove
#include "TWI/TWI.h" // to remove
#include "Utilities/Parameters.h"
#include "LEM/LEM.h"


using namespace std;


/**
 * @brief Main function.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/
int main(
  int argc,
  char** argv){
/*
  Image im0, im1, im2, im3;
  //im0.read("/home/ibal/Images/input/DEM/srtm_france/srtm_36_03_f32_downsampled_6.tif");
  //im1.read("/home/ibal/Images/input/DEM/srtm_france/srtm_36_04_f32_downsampled_6.tif");
  //im1.read("/home/ibal/Images/input/DEM/srtm_france/srtm_37_03_f32_downsampled_6.tif");
  im1.read("/home/ibal/Images/input/DEM/srtm_france/srtm_37_03_f32.tif");
  im3.read("/home/ibal/Images/input/DEM/srtm_france/srtm_37_04_f32_downsampled_6.tif");
  Image im(2 * im0.width(), 2 * im0.height(), 1);
  float min1, min2, min3, min0, max1, max2, max3, max0, mean1, mean2, mean3, mean0;
  im0.getStats(min0, max0, mean0);
  im1.getStats(min1, max1, mean1);
  im2.getStats(min2, max2, mean2);
  im3.getStats(min3, max3, mean3);

  float maxV = std::max(std::max(max0, max1), std::max(max2, max3));
  cout << "maxV = " << maxV << endl;
  const float limit = 1.01 * maxV;
  float minV = maxV - limit;
  cout << "minV = " << minV << endl;

  for (size_t i = 0; i < im0.height(); i++) {
    const float* iI0 = im0.getPtr(0, i);
    const float* iI1 = im1.getPtr(0, i);
    const float* iI2 = im2.getPtr(0, i);
    const float* iI3 = im3.getPtr(0, i);
    float* oI1 = im.getPtr(0, i);
    float* oI2 = im.getPtr(0, im0.height() + i);
    for (size_t j = 0; j < im0.width(); j++) {
      oI1[j] = std::max(iI0[j], minV) - minV;
      //oI1[j] = iI0[j];
      oI1[im0.width() + j] = std::max(iI2[j], minV) - minV;
      //oI1[im0.width() + j] = iI2[j];
      oI2[j] = std::max(iI1[j], minV) - minV;
      //oI2[j] = iI1[j];
      oI2[im0.width() + j] = std::max(iI3[j], minV) - minV;
      //oI2[im0.width() + j] = iI3[j];
    }
  }

  //Image im(im1.width(), 3 * im1.height() / 4, im1.channels());
  Image im(im1.width() / 2, im1.height() / 2, im1.channels());
  for (size_t c = 0; c < im.channels(); c++) {
    for (size_t i = 0; i < im.height(); i++) {
      const float* iI = im1.getPtr(c, i + im1.height() - im.height()) + im1.width() - im.width();
      float* oI = im.getPtr(c, i);
      for (size_t j = 0; j < im.width(); j++) {
        oI[j] = iI[j];
      }
    }
  }

  float minV, maxV, meanV;
  im.getStats(minV, maxV, meanV);
  for (size_t c = 0; c < im.channels(); c++) {
    for (size_t i = 0; i < im.height(); i++) {
      float* iI = im.getPtr(c, i);
      for (size_t j = 0; j < im.width(); j++) {
        const double val = (double) iI[j];
        iI[j] = (float) 255.0 * sqrt((val - minV) / (maxV - minV));
      }
    }
  }

  //im.write("/home/ibal/Images/input/DEM/srtm_france/srtm_france_downsampled_min.tif");
  //im.normalize(0, 255);
  im.write("/home/ibal/Images/input/DEM/srtm_france/dem_france_north_crop.png");
  return EXIT_SUCCESS;
*/
  //! To avoid use of denormalized numbers
  //_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
  //fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);

  //! Retrieve the arguments
  Parameters params;
  if (params.checkArgs(argc, argv) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }
  if (params.verbose()) {
#ifdef __AVX__
    cout << "AVX is used." << endl;
#endif // __AVX__
#ifdef __SSE__
    cout << "SSE is used." << endl;
#endif // __SSE__
    //params.print();
  }

  //! Choose the number of threads to use
#ifdef _OPENMP
  const int nb = getOptimalThreadsNumber(params.nbThreads(), params.inpName());
  //omp_set_dynamic(0);
  omp_set_num_threads(nb);
#endif // _OPENMP

  //! Initialization of time
  struct timespec time;
  clock_gettime(CLOCK_MONOTONIC, &time);

  //! Create the LEM class
  LEM lem;

  //! Initialize the environment
  lem.init(params);

  //! Finish the water initialization (only the water image is modified)
  lem.waterInitialization();

  //! Perform the main iteration
  if (params.autoParams()) {
    lem.autoSelect();
  }
  else {

    //! Perform the main iteration
    lem.runMainIterations();
  }

  //! Save the images
  lem.saveImages();

  //! Total amount of time
  if (params.verbose()) {
    cout << "-----------------------------------------------------------\n";
    getTime(time, "Total                                      ");
    cout << "-----------------------------------------------------------\n";
  }

  //! Exit
  return EXIT_SUCCESS;
}

