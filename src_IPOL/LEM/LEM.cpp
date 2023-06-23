/**
 * @file LEM.cpp
 *
 * @brief Class to handle the Landscape Evolution Model.
 *
 * @author Jérôme Darbon <jerome@math.ucla.edu> (original version).
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
   @author Miguel Colom - http://mcolom.info
 **/


//! Global includes
#include <omp.h>
#include <iostream>
#include <sstream>
#include <cfloat>
#include <cassert>


//! Local includes
#include "LEM.h"
#include "../Utilities/Utilities.h"
#include "../Utilities/Memory.h"
#include "../TWI/TWI.h"


using namespace std;


//! Default constructor
LEM::LEM() :

  //! Miscellaneous
  m_time            (timespec()),
#ifdef _OPENMP
  m_nbThreads       (omp_get_max_threads()),
#else
  m_nbThreads       (1),
#endif // _OPENMP
  m_waterHasBeenInitialized(false),

  //! Parameters
  m_params(NULL),

  //! Sizes
  m_width (0),
  m_height(0),
  m_border(1),

  //! Information values
  m_massOriginal(0.f),
  m_massTarget  (0.f),
  m_massCurrent (0.f),
  m_meanOriginal(0.0),
  m_uplift      (0.f),
  m_upliftAll   (0.f),
  m_iterCurrent (0),
  m_iterMax     (0),
  m_ratio       (0.f),
  m_nbCheck     (0),
  m_masses      ((float*) memalloc(16, m_nbThreads * sizeof(float))),
  m_norms       (NULL),
  m_nbs         (NULL),

  //! Images
  m_imLandscape    (NULL),
  m_imWaterLvl     (NULL),
  m_imSediment     (NULL),
  m_imLandscapeInit(NULL),
  m_imWaterLvlInit (NULL),
  m_imLandscapeRaw (NULL),
  m_imUplift       (NULL),

  //! Evolution schems
  m_evolutionWater                (NULL),
  m_evolutionCreep                (NULL),
  m_evolutionLandscape            (NULL),
  m_evolutionWaterAndSedimentation(NULL) {

  //! floatime initialization
  clock_gettime(CLOCK_MONOTONIC, &m_time);
}


//! Copy constructor.
LEM::LEM(
  const LEM& i_lem) :

  //! Miscellaneous
  m_time            (i_lem.m_time),
  m_nbThreads       (i_lem.m_nbThreads),
  m_waterHasBeenInitialized(i_lem.m_waterHasBeenInitialized),

  //! Parameters
  m_params(i_lem.m_params),

  //! Sizes
  m_width (i_lem.m_width),
  m_height(i_lem.m_height),
  m_border(i_lem.m_border),

  //! Information values
  m_massOriginal(i_lem.m_massOriginal),
  m_massTarget  (i_lem.m_massTarget),
  m_massCurrent (i_lem.m_massCurrent),
  m_meanOriginal(i_lem.m_meanOriginal),
  m_uplift      (i_lem.m_uplift),
  m_upliftAll   (i_lem.m_upliftAll),
  m_iterCurrent (i_lem.m_iterCurrent),
  m_iterMax     (i_lem.m_iterMax),
  m_ratio       (i_lem.m_ratio),
  m_nbCheck     (i_lem.m_nbCheck),
  m_masses      ((float*) memalloc(16, m_nbThreads * sizeof(float))),
  m_norms       ((float*) memalloc(16, m_nbThreads * (m_params->nbIter() + 1) * sizeof(float))),
  m_nbs         ((float*) memalloc(16, m_nbThreads * (m_params->nbIter() + 1) * sizeof(float))),

  //! Images
  m_imLandscape    (new Image(*i_lem.m_imLandscape)),
  m_imWaterLvl     (new Image(*i_lem.m_imWaterLvl)),
  m_imSediment     (new Image(*i_lem.m_imSediment)),
  m_imLandscapeInit(new Image(*i_lem.m_imLandscapeInit)),
  m_imWaterLvlInit (new Image(*i_lem.m_imWaterLvlInit)),
  m_imLandscapeRaw (new Image(*i_lem.m_imLandscapeRaw)),
  m_imUplift       (new Image(*i_lem.m_imUplift)),

  //! Evolution schems
  m_evolutionWater                (new WaterEvolution(*i_lem.m_evolutionWater)),
  m_evolutionCreep                (new CreepEvolution(*i_lem.m_evolutionCreep)),
  m_evolutionLandscape            (new LandscapeEvolution(*i_lem.m_evolutionLandscape)),
  m_evolutionWaterAndSedimentation(new WaterAndSedimentationEvolution(*i_lem.m_evolutionWaterAndSedimentation)) {

  //! Copy the data
  for (size_t n = 0; n < m_nbThreads; n++) {
    m_masses[n] = i_lem.m_masses[n];
  }
  for (size_t n = 0; n < m_nbThreads * (m_params->nbIter() + 1); n++) {
    m_norms[n] = i_lem.m_norms[n];
    m_nbs  [n] = i_lem.m_nbs  [n];
  }
}


//! Default destructor
LEM::~LEM() {
  releaseMemory();
  m_params = NULL;
}


//! Operator overload.
LEM& LEM::operator=(
  const LEM& i_lem) {

  if (&i_lem == this) {
    return *this;
  }

  releaseMemory();

  new (this) LEM(i_lem);
  return *this;
}


//! Initialize the environment. Must be call first.
void LEM::init(
  Parameters& p_params) {

  //! Release memory
  releaseMemory();

  //! Miscallaneous
  m_params = new Parameters(p_params);
  const bool verbose = m_params->verbose();
  m_masses = (float*) memalloc(16, m_nbThreads * sizeof(float));
  m_norms  = (float*) memalloc(16, m_nbThreads * (m_params->nbIter() + 1) * sizeof(float));
  m_nbs    = (float*) memalloc(16, m_nbThreads * (m_params->nbIter() + 1) * sizeof(float));

  //! For convenience
  const bool doRaw = m_params->noise() > 0.f || m_params->blur() > 0.f;

  //! Read the image
  Image imLandscape;
  imLandscape.read(p_params.inpName());
  //imLandscape.makeCone(false, false);
  //imLandscape.makePlan(85);
  
  printf("1) m_imLandscape: %d x %d\n", imLandscape.width(), imLandscape.height());

  //! Check that the image is gray, otherwise average the channels
  if (imLandscape.channels() > 1) {

    //! Special case
    if (imLandscape.channels() != 3) {
      cout << "The number of channels (" << imLandscape.channels();
      cout << ") is not handled." << endl;
      exit(EXIT_FAILURE);
    }

    //! RGB case
    if (verbose) cout << "Convert to gray level image." << endl;
    for (size_t i = 0; i < imLandscape.height(); i++) {
      float* iR = imLandscape.getPtr(0, i);
      const float* iG = imLandscape.getPtr(1, i);
      const float* iB = imLandscape.getPtr(2, i);
      for (size_t j = 0; j < imLandscape.width(); j++) {
        iR[j] = (iR[j] + iG[j] + iB[j]) / 3.f;
      }
    }
  }

  //! Fill in the rain matrix
  // Add the rain matrix to the params
  this->m_rainMatrix = new float[imLandscape.width() * imLandscape.height()];
  m_params->setRMatrix(this->m_rainMatrix);

  if (p_params.rainName()) {
      // The image is a matrix where the value of the pixel in
      // [0, 255] corresponds to a probability.
      Image imRain;
      imRain.read(p_params.rainName());

      // The size of the image should be the same as the landscape
      assert(imLandscape.width() == imRain.width());
      assert(imLandscape.height() == imRain.height());

      float *imRainPtr = imRain.getPtr();
      for (unsigned i = 0; i < imRain.width() * imRain.height(); i++) {
        this->m_rainMatrix[i] = (float)imRainPtr[i] / 255.0;
      }
  } else {
    // If no rain matrix is given, we assume that all entries of the image are equal to 1.0
    for (unsigned i = 0; i < imLandscape.width() * imLandscape.height(); i++) {
      this->m_rainMatrix[i] = 1.0;
    }
  }

  //! Add an uplift (only for automatic selection of the parameters)
  if (m_params->autoParams()) {
    const float uplift = 10.f;
    for (size_t i = 0; i < imLandscape.height(); i++) {
      float* iL = imLandscape.getPtr(0, i);
      for (size_t j = 0; j < imLandscape.width(); j++) {
        iL[j] += uplift;
      }
    }
  }

  //! Apply the zoom-in
  if (m_params->zoom() > 1) {
    Image im = imLandscape;
    im.sample(imLandscape, im.width() * im.height() * m_params->zoom() * m_params->zoom());
  }

  //! Keep a size of image compatible with the multiscale approach (3 scales)
  m_height = 4 * (imLandscape.height() / 4);
  m_width  = 4 * (imLandscape.width () / 4);

  //! Initialization of the images
  m_imWaterLvl  = new Image(m_width, m_height, 1, m_border);
  m_imLandscape = new Image(m_width, m_height, 1, m_border);
  m_imSediment  = new Image(m_width, m_height, 1, m_border);
  m_imUplift    = new Image(m_width, m_height, 1, m_border);
  if (doRaw) {
    m_imLandscapeRaw = new Image(m_width, m_height, 1, m_border);
    for (size_t i = 0; i < m_height; i++) {
      const float* iL = imLandscape.getPtr(0, i);
      float* oL = m_imLandscapeRaw->getPtr(0, i);

      for (size_t j = 0; j < m_width; j++) {
        oL[j] = iL[j];
      }
    }
  }

  //! Apply a Gaussian blur to the original image
  if (m_params->blur() > 0.f) {
    Image imG;
    imLandscape.applyGaussianBlur(imG, m_params->blur());
    imLandscape = imG;
  }

  //! Apply a Gaussian noise to the original image
  if (m_params->noise() > 0.f) {
    imLandscape.applyGaussianNoise(std::max(m_params->noise(), 0.f), m_params->oceanLevel() * 255);
  }

  //! Get min, max values of the landscape image
  float minVal, maxVal, mean, minValRaw, maxValRaw, meanRaw;
  if (doRaw) {
    m_imLandscapeRaw->getStats(minValRaw, maxValRaw, meanRaw);
  }
  imLandscape.getStats(minVal, maxVal, mean);

  //!Normalize the landscape between [0, 1]
  const float norm = 1.f / (maxVal - minVal);
  const float normRaw = 1.f / (maxValRaw - minValRaw);
  m_meanOriginal = 0.0;
  for (size_t i = 0; i < m_height; i++) {
    const float* iL = imLandscape.getPtr(0, i);
    float* oL = m_imLandscape->getPtr(0, i);
    float* oR = doRaw ? m_imLandscapeRaw->getPtr(0, i) : NULL;

    for (size_t j = 0; j < m_width; j++) {
      oL[j] = std::max(0.f, std::min(1.f, norm * (iL[j] - minVal)));
      m_meanOriginal += double(oL[j]);
      if (doRaw) {
        oR[j] = std::max(0.f, std::min(1.f, normRaw * (oR[j] - minVal)));
      }
    }
  }
  m_meanOriginal /= double(m_height * m_width);

  //! Initialize the water image
  const float value = m_params->initWaterLvl();
  for (size_t i = 0; i < m_height; i++) {
    float* oW = m_imWaterLvl->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      oW[j] = value;
    }
  }

  //! Evolution schemes
  m_evolutionWater                 = new WaterEvolution                ();
  m_evolutionCreep                 = new CreepEvolution                ();
  m_evolutionLandscape             = new LandscapeEvolution            ();
  m_evolutionWaterAndSedimentation = new WaterAndSedimentationEvolution();
  m_evolutionWater                ->init(m_width, m_height, m_border, p_params);
  
  // Print sizes
  printf("m_imLandscape: %d x %d\n", m_imLandscape->width(), m_imLandscape->height());

  //! Add the images to save
  m_imLandscapeInit = new Image(*m_imLandscape);
  m_imLandscapeInit->fillBorder(-1, m_params->borderLand());

  //! Mass original
  m_massOriginal = integrateLandscapeOverOceanLevel();

  //! Print some informations
  if (verbose) getTime(m_time, "LEM initialization                         ");
}


//! Finish the water initialization
void LEM::waterInitialization() {

  //! For convenience
  const bool verbose   = m_params->verbose();
  const size_t nbWater = m_params->nbWater();

  //! Border conditions for the landscape
  m_imLandscape->fillBorder(-1, m_params->borderLand());

  //! Multiscale
  if (m_params->useMultiScale()) {

    //! Convenience images
    Image imLandscapeRef(*m_imLandscape);
    Image imLandscape;
    Image imWaterLvl;

    //! Scale 2
    imLandscapeRef.zoomOut(imLandscape);
    imLandscape.zoomOut(*m_imLandscape);
    m_imWaterLvl->zoomOut(imWaterLvl);
    imWaterLvl.zoomOut(*m_imWaterLvl);
    waterInitializationCore(nbWater / 2);

    //! Scale 1
    imLandscapeRef.zoomOut(*m_imLandscape);
    imWaterLvl = *m_imWaterLvl;
    imWaterLvl.zoomIn(*m_imWaterLvl);
    waterInitializationCore(nbWater / 4);

    //! Scale 0
    *m_imLandscape = imLandscapeRef;
    imWaterLvl = *m_imWaterLvl;
    imWaterLvl.zoomIn(*m_imWaterLvl);
    waterInitializationCore(nbWater / 32);
  }
  else {
    waterInitializationCore(nbWater);
  }

  //! Check the values according to the ocean value
  const float oceanLvl = m_params->oceanLevel();
  for (size_t i = 0; i < m_height; i++) {
    float* iW = m_imWaterLvl->getPtr(0, i);
    const float* iL = m_imLandscape->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {

      //! If the landscape is above the ocean, keep it that way
      if (iL[j] >= oceanLvl) {
        iW[j] = std::max(iW[j], float(0)); // to ensure positive values
      }

      //! Otherwise, landscape + water should be equal to the ocean level
      else {
        iW[j] = std::max(oceanLvl - iL[j], float(0));
      }
    }
  }

  //! Add image to save
  m_imWaterLvlInit = new Image(*m_imWaterLvl);

  //! Keep track of the water initialization
  m_waterHasBeenInitialized = true;

  //! Print some informations
  if (verbose) getTime(m_time, "LEM finish water initialization            ");
}


//! Perform the main iteration after the water initialization.
void LEM::runMainIterations() {

  //! Check that the water has been initialized
  if (!m_waterHasBeenInitialized) {
    this->waterInitialization();
  }

  //! If automatic selection of E
  if (m_params->autoE()) {
    this->runAutoErosion();
  }

  //! Launch the core program
  this->runMainIterationsCore();
}

bool file_exists(const char *filename) {
    FILE *fp = NULL;
    if (NULL != (fp = fopen(filename, "rb"))) {
      fclose(fp);
      return true;
    }
    return false;
}

//! Perform the main iteration after the water initialization.
void LEM::runMainIterationsCore() {
  printf("runMainIterationsCore() *****************************\n");
    
  //! For convenience
  const bool verbose       = m_params->verbose();
  const bool addNoise      = m_params->addNoise() > 0.f;
  const float valNoise     = std::min(std::max(m_params->addNoise(), 0.f), 1.f);
  const size_t nbNoise     = m_params->iterNoise();
  const size_t borderLand  = m_params->borderLand();
  const size_t borderWater = m_params->borderWater();
  const size_t nbWater     = m_params->iterWater();
  m_evolutionWater->setRain(m_params->rainWater());
  m_evolutionWater->setRainMatrix(m_params->rainMatrix());

  //! Initialize the full evolution stopping criteria
  initializeInfos();

  //! Update the parameters values
  m_evolutionCreep                ->init(m_width, m_height, m_border, *m_params);
  m_evolutionLandscape            ->init(m_width, m_height, m_border, *m_params);
  m_evolutionWaterAndSedimentation->init(m_width, m_height, m_border, *m_params);

  //! Initialize the images
  for (size_t i = 0; i < m_height; i++) {
    const float* iL = m_imLandscapeInit->getPtr(0, i);
    const float* iW = m_imWaterLvlInit ->getPtr(0, i);
    float* oL = m_imLandscape->getPtr(0, i);
    float* oW = m_imWaterLvl ->getPtr(0, i);
    float* oC = m_imSediment ->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      oL[j] = iL[j];
      oW[j] = iW[j];
      oC[j] = 0.f;
    }
  }

  bool debug_load_images = true;

  // DEBUG: load m_imWaterLvl and m_imSediment from TIFF images
  // m_imSediment, m_imWaterLvl
  if (debug_load_images) {
    const char *debug_water_filename = "debug_water.tif";
    if (file_exists(debug_water_filename)) {
      printf("*** WARNING ***: loading %s\n", debug_water_filename);
      m_imWaterLvl->read(debug_water_filename, borderWater);
    }
    const char *debug_sediment_filename = "debug_sediment.tif";
    if (file_exists(debug_sediment_filename)) {
      printf("*** WARNING ***: loading %s\n", debug_sediment_filename);
      m_imSediment->read(debug_sediment_filename, borderWater);
    }
  }

  //! To keep track of the means
  double* meansCurrent = (double*) memalloc(16, m_nbThreads * sizeof(double));

  //! Initialize the progress bar
  ProgressBar bar;
  size_t nbBar = 0;
  bar.init(100, nbBar);

#ifdef _OPENMP
#pragma omp parallel
{
#endif

  //! Perform full evolution
  while (m_massCurrent >= m_massTarget && m_iterCurrent < m_iterMax) {

    //! Get the current thread id
#ifdef _OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif

    //! Add the noise every nbNoise iterations
    if (addNoise && (m_iterCurrent % nbNoise == 0)) {
#ifdef _OPENMP
#pragma omp single
#endif // _OPENMP
      m_imLandscape->applyGaussianNoise(valNoise, m_params->oceanLevel());
    }

#ifdef _OPENMP
#pragma omp barrier
#endif

    //! Perform the water and sedimentation evolution
    m_imLandscape->fillBorder(tid, borderLand);
    m_imSediment ->fillBorder(tid, borderWater);
    m_imWaterLvl ->fillBorder(tid, borderWater);
    m_evolutionWaterAndSedimentation->performOneStepEvolution(*m_imLandscape,
      *m_imSediment, *m_imWaterLvl, tid);

#ifdef _OPENMP
#pragma omp barrier
#endif

    //! Perform the landscape evolution
    m_imLandscape->fillBorder(tid, borderLand);
    m_imSediment ->fillBorder(tid, borderWater);
    m_imWaterLvl ->fillBorder(tid, borderWater);
    if (m_params->autoUplift() > 0) {

      //! Compute the transfert and update the mean
      m_evolutionLandscape->computeTransfert(*m_imWaterLvl, *m_imLandscape,
        *m_imSediment, meansCurrent[tid], tid);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

#ifdef _OPENMP
#pragma omp single
#endif // _OPENMP
      this->computeUplift(meansCurrent);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

      //! Update the landscape and apply the uplift
      m_evolutionLandscape->applyTransfert(*m_imLandscape,
        *m_imSediment, m_uplift, *m_imUplift, tid);
    }
    else {
      m_evolutionLandscape->performOneStepEvolution(*m_imWaterLvl,
      *m_imLandscape, *m_imSediment, tid);
    }

#ifdef _OPENMP
#pragma omp barrier
#endif

    //! Perform the creep evolution
    m_imLandscape->fillBorder(tid, borderLand);
    m_evolutionCreep->performOneStepEvolution(*m_imLandscape, tid);

#ifdef _OPENMP
#pragma omp barrier
#endif

    //! Perform a given number of water-only equation
    m_imLandscape->fillBorder(tid, borderLand);
    m_evolutionWaterAndSedimentation->setRain(m_params->rainWater());
    m_evolutionWaterAndSedimentation->setRainMatrix(m_params->rainMatrix());
    for (size_t nw = 0; nw < nbWater; nw++) {

      //! Set border condition
      m_imWaterLvl->fillBorder(tid, borderWater);
      m_imSediment->fillBorder(tid, borderWater);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

      //! Perform one step of water evolution
      m_evolutionWaterAndSedimentation->performOneStepEvolution(*m_imLandscape,
        *m_imSediment, *m_imWaterLvl, tid);
    }
    m_evolutionWaterAndSedimentation->setRain(m_params->rain());
    m_evolutionWaterAndSedimentation->setRainMatrix(m_params->rainMatrix());

#ifdef _OPENMP
#pragma omp barrier
#endif

    //! Update the landscape mass
    if (m_iterCurrent % m_nbCheck == 0) {
      integrateLandscapeOverOceanLevel(tid);
    }

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP
    //! Update the stability norm
    if (m_params->stabCurv()) {
      computeStabilityNorm_tid(tid);
    }

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

    //! Update info
#ifdef _OPENMP
#pragma omp single
#endif
    updateInfos(bar, nbBar);
  }
#ifdef _OPENMP
}
#endif

  //! Release memory
  memfree(meansCurrent);

  //! Print some informations
  if (verbose) {
    bar.update(100);
    getTime(m_time, "LEM perform full evolution                 ");
  }
}


//! Core part of the water initialization.
void LEM::waterInitializationCore(
  const size_t p_nbIter) {

  //! Initialization
  m_evolutionWater->init(m_imLandscape->width(), m_imLandscape->height(),
                         m_imLandscape->border(), *m_params);
  const size_t borderWater = m_params->borderWater();

  //! Boost
  if (m_params->boostWater() > 0.5f) {
    m_evolutionWater->setTimeStep(m_params->boostWater());
  }
  const size_t nbIterMax = (m_params->boostWater() > 0.5f ?
    p_nbIter / (2 * m_params->boostWater()) : p_nbIter);

  //! Initialization
  ProgressBar bar;
  bar.init(nbIterMax, 0);
  const size_t nbCheck = std::max(nbIterMax / 100, (size_t) 1);

#ifdef _OPENMP
#pragma omp parallel
{
#endif // _OPENMP

  //! Main loop
  for (size_t nbIter = 0; nbIter < nbIterMax; nbIter++) {

#ifdef _OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif

    //! Set border condition
    m_imWaterLvl->fillBorder(tid, borderWater);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

    //! Perform one step of water evolution
    m_evolutionWater->performOneStepEvolution(*m_imLandscape, *m_imWaterLvl, tid);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

    //! Display the progression every 100 iterations
#ifdef _OPENMP
#pragma omp single
#endif
    if (m_params->verbose() && nbIter % nbCheck == 0) {
      bar.update(std::min(nbIter, nbIterMax - 2));
    }
  }
#ifdef _OPENMP
}
#endif

  //! Display the final progression
  if (m_params->verbose()) {
    bar.update(std::max(int(nbIterMax) - 1, 0));
  }
}


//! Automatically estimate the value of erosion.
void LEM::runAutoErosion() {

  //! Verbose
  const bool verbose = m_params->verbose();
  if (verbose) {
    cout << "Estimate E" << endl;
  }

  //! Keep track of the wanted values
  const int   nbIter = m_params->nbIter();
  const float perc   = m_params->percLand();

  //! Boundaries of the erosion parameters
  m_params->setE(8.f);
  float ePrevious = 0.f;
  float eNext = m_params->e();

  //! Adjust the values for small iterations
  m_params->setI(std::min(nbIter, std::max(nbIter / 5, 50)));
  m_params->setP(1.f);
  m_params->setV(false);

  //! What percentage we should have reached with this amount of iterations
  const float pTarget = float(m_params->nbIter()) * perc / float(nbIter);

  //! Get the maximum boundary
  bool keepGoing = true;
  while (keepGoing) {

    //! Launch the core program with few iterations
    this->runMainIterationsCore();

    //! Check the current percentage of land removed
    const float pCurr = 1.f - m_massCurrent / m_massOriginal;

    //! Verbose
    if (verbose) {
      cout << "mO = " << m_massOriginal << ", mC = " << m_massCurrent << endl;
      cout << "E (" << m_params->e() << "): P = " << pCurr;
      cout << " (T = " << pTarget << ")" << endl;
    }

    //! Increase the upperbound if needed
    if (pCurr < pTarget) {
      eNext *= 2.f;
      m_params->setE(eNext);
    }
    else {
      keepGoing = false;
      m_params->setE((eNext + ePrevious) / 2.f);
    }
  }

  //! While E doesn't fit it's purpose
  keepGoing = true;
  int nbLoop = 0;
  while (keepGoing && nbLoop < 25) {

    //! Launch the core program with few iterations
    this->runMainIterationsCore();

    //! Check the current percentage of land removed
    const float pCurr = 1.f - m_massCurrent / m_massOriginal;

    //! Verbose
    if (verbose) {
      cout << "mO = " << m_massOriginal << ", mC = " << m_massCurrent << endl;
      cout << "E (" << m_params->e() << "): P = " << pCurr;
      cout << " (T = " << pTarget << ")" << endl;
    }

    //! Stop the loop if close to the target
    if (fabsf((pCurr - pTarget) / pTarget) < 1e-5) {
      keepGoing = false;
    }

    //! Otherwise adjust the erosion parameter
    else {
      if (pCurr > pTarget) {
        eNext = m_params->e();
        m_params->setE((eNext + ePrevious) / 2.f);
      }
      else {
        ePrevious = m_params->e();
        m_params->setE((eNext + ePrevious) / 2.f);
      }
    }

    //! Keep track of the number of iterations
    nbLoop++;
  }

  //! Put back the wanted values
  m_params->setI(nbIter);
  m_params->setP(perc);
  m_params->setV(verbose);
}


//! Integrate landscape over ocean level.
float LEM::integrateLandscapeOverOceanLevel() const {

  //! Initialization
  float result         = 0.f;
  const float oceanLvl = m_params->oceanLevel();

  for (size_t i = 0; i < m_height; i++) {
    const float* iL = m_imLandscape->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      result += std::max(iL[j] - oceanLvl, (float) 0);
    }
  }

  //! Return the result
  return result;
}


//! Integrate landscape over ocean level.
void LEM::integrateLandscapeOverOceanLevel(
  const size_t p_tid) {

  //! Initialization
  float result         = 0.f;
  const float oceanLvl = m_params->oceanLevel();
  const size_t hBegin  = m_imLandscape->nHeight(p_tid);
  const size_t hEnd    = m_imLandscape->nHeight(p_tid + 1);

#ifdef __AVX__
  __m256 zMean = _mm256_setzero_ps();
  const __m256 zOceanLvl = _mm256_set1_ps(oceanLvl);
  const __m256 z0 = _mm256_setzero_ps();
  const __m256 zU = _mm256_set1_ps(m_upliftAll);
#endif // __AVX__
#ifdef __SSE__
  __m128 xMean = _mm_setzero_ps();
  const __m128 xOceanLvl = _mm_set1_ps(oceanLvl);
  const __m128 x0 = _mm_setzero_ps();
  const __m128 xU = _mm_set1_ps(m_upliftAll);
#endif // __SSE__

  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iL = m_imLandscape->getPtr(0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {
      zMean += _mm256_max_ps(_mm256_loadu_ps(iL + j) - zOceanLvl - zU, z0);
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {
      xMean += _mm_max_ps(_mm_loadu_ps(iL + j) - xOceanLvl - xU, x0);
    }
#endif // __SSE__

    //! Normal version
    for (; j < m_width; j++) {
      result += std::max(iL[j] - oceanLvl - m_upliftAll, 0.f);
    }
  }

  //! Get the value
#ifdef __AVX__
  float zVal[8];
  _mm256_storeu_ps(zVal, zMean);
  for (size_t n = 0; n < 8; n++) {
    result += zVal[n];
  }
#endif // __AVX__
#ifdef __SSE__
  float xVal[4];
  _mm_storeu_ps(xVal, xMean);
  for (size_t n = 0; n < 4; n++) {
    result += xVal[n];
  }
#endif // __SSE__

  //! Return the result
  m_masses[p_tid] = result;
}


//! Integrate water over ocean level.
void LEM::integrateWaterOverOceanLevel(
  float& o_mass,
  float& o_max,
  float& o_norm) const {

  //! Initialization
  o_mass   =  0.f;
  float nb =  0.f;
  o_max    = -1.f;
  o_norm   =  0.f;
  const float oceanLvl = m_params->oceanLevel();

  //! Get the result
  for (size_t i = 0; i < m_height; i++) {
    const float* iL = m_imLandscape->getPtr(0, i);
    const float* iW = m_imWaterLvl ->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      if (iL[j] > oceanLvl) {
        const float value = iW[j];
        nb     += 1.f;
        o_mass += value;
        o_max   = std::max(o_max, value);
        o_norm += value * value;
      }
    }
  }

  //! Get the result normalized
  o_mass /= nb;
  o_norm = sqrt(o_norm / nb);
}

//! Save images after the main evolution.
void LEM::saveImages() {

  //! For convenience
  const string outFolder  = m_params->outName();
  Image im, imTwiO, imTwiF;
  std::vector<std::string> names;
  std::vector<Image *> imToSave;
  const bool verbose = m_params->verbose();
  const bool doRaw   = m_params->blur() > 0.f || m_params->noise() > 0.f;

  //! Un-apply the zoom if needed
  if (m_params->zoom() > 1 && !m_params->noZoomOut()) {
    applyZoomOut();
    if (verbose) getTime(m_time, "LEM apply zoom-out                         ");
  }

  //! Compute stats and save them
  computeStats();
  if (verbose) getTime(m_time, "LEM compute statistics                     ");

  //! Apply quad visualization only if the input image is too small
  const bool applyQuad = false;//m_width * m_height < 1024 * 1024;

  //! Compute some statistics for the images (min, max and mean values)
  float minW, maxW, meanW, minC, maxC, meanC, minI, maxI, meanI, minL, maxL, meanL;
  m_imWaterLvlInit->getStats(minI, maxI, meanI);
  m_imWaterLvl    ->getStats(minW, maxW, meanW);
  m_imSediment    ->getStats(minC, maxC, meanC);
  m_imLandscape   ->getStats(minL, maxL, meanL);

  //! To avoid problem of saturation
  maxW = std::min(maxW, 1.f);
  maxI = std::min(maxI, 1.f);
  maxC = std::min(maxC, 1.f);

  //! Initialize the images and names
  //! Mandatory images
  size_t nbImages = 5;
  imToSave.push_back(new Image(*m_imLandscapeInit));
  names.push_back(outFolder + "originalLandscape.tif");
  imToSave.push_back(new Image(m_width, m_height,
    m_imWaterLvlInit->channels()));
  names.push_back(outFolder + "initWater.tif");
  imToSave.push_back(new Image(*m_imLandscape));
  names.push_back(outFolder + "finalLandscape.tif");
  imToSave.push_back(new Image(m_width, m_height,
    m_imWaterLvl->channels()));
  names.push_back(outFolder + "finalWater.tif");
  imToSave.push_back(new Image(m_width, m_height,
    m_imSediment->channels()));
  names.push_back(outFolder + "finalSediment.tif");
  if (doRaw) {
    imToSave.push_back(new Image(*m_imLandscapeRaw));
    names.push_back(outFolder + "rawLandscape.tif");
    nbImages++;
  }

  //! If TWI wanted
  size_t indTwiO = 0, indTwiF = 0;
  float minTwiO, maxTwiO, minTwiF, maxTwiF, mean;
  if (m_params->computeTWI()) {
    TWI twi;
    nbImages += 2;
    imToSave.push_back(new Image(m_width, m_height, 1));
    names.push_back(outFolder + "originalLandscape_twi.tif");
    indTwiO = imToSave.size() - 1;
    twi.computeTWI(*m_imLandscapeInit, imTwiO);
    imTwiO.getStats(minTwiO, maxTwiO, mean);

    imToSave.push_back(new Image(m_width, m_height, 1));
    names.push_back(outFolder + "finalLandscape_twi.tif");
    indTwiF = imToSave.size() - 1;
    twi.computeTWI(*m_imLandscape, imTwiF);
    imTwiF.getStats(minTwiF, maxTwiF, mean);
  }

  //! If topography wanted
  size_t indTopoO = 0, indTopoF = 0;
  if (m_params->visualizeTopo() != NULL) {
    nbImages += 2;
    m_imLandscapeInit->getTopo(im, *m_imWaterLvlInit, m_params->oceanLevel(),
      0, 1, m_params->visualizeTopo(), m_params->waterColor());
    imToSave.push_back(new Image(im));
    names.push_back(outFolder + "originalLandscape_topo.png");
    indTopoO = imToSave.size() - 1;

    m_imLandscape->getTopo(im, *m_imWaterLvl, m_params->oceanLevel(), 0, 1,
      m_params->visualizeTopo(), m_params->waterColor());
    imToSave.push_back(new Image(im));
    names.push_back(outFolder + "finalLandscape_topo.png");
    indTopoF = imToSave.size() - 1;
    if (doRaw) {
      nbImages++;
      Image imWater(m_width, m_height, 1);
      m_imLandscapeRaw->getTopo(im, imWater, m_params->oceanLevel(), 0, 1,
        m_params->visualizeTopo(), m_params->waterColor());
      imToSave.push_back(new Image(im));
      names.push_back(outFolder + "rawLandscape_topo.png");
    }
  }

  //! Save the 3D models (must be done at the end)
  const size_t nbPoints = 300 * 300;
  size_t indMeshO = 0, indMeshF = 0;
  if (m_params->visualizeMesh()) {

    //! Compute the size of the output image
    const float r = sqrt(float(nbPoints) / float(m_height * m_width));
    const size_t w = m_height * m_width <= nbPoints ? m_width : m_width  * r;
    const size_t h = m_height * m_width <= nbPoints ? m_width : m_height * r;
    imToSave.push_back(new Image(w, h, m_imLandscapeInit->channels()));
    names.push_back(outFolder + "originalLandscape.obj");
    indMeshO = imToSave.size() - 1;

    imToSave.push_back(new Image(w, h, m_imLandscape->channels()));
    names.push_back(outFolder + "finalLandscape.obj");
    indMeshF = imToSave.size() - 1;
  }

  //! Because of the automatic uplift, the highest value of landscape may be > 1
  const float factor = 1.f;  // 255.f / std::max(1.f, maxL);

  //! Pre-process the images
#ifdef _OPENMP
#pragma omp parallel
{
  const size_t tid = omp_get_thread_num();
#else
  const size_t tid = 0;
#endif // _OPENMP

  //! First, the original landscape
  imToSave[0]->multiply(factor, tid);

  //! Second, the original water level
  m_imWaterLvlInit->applySqrtNormalization(*imToSave[1], tid, minI, maxI, minI,
    maxI);    
//  m_imWaterLvlInit->*imToSave[1];
//  imToSave[1]->multiply(factor, tid);
//  m_imWaterLvlInit->multiply(factor, *imToSave[1] );
    
  //! Third, the final landscape
  imToSave[2]->multiply(factor, tid);

  //! Fourth, the final water level
  m_imWaterLvl->applySqrtNormalization(*imToSave[3], tid, minW, maxW, minW,
    maxW);


  //! Fifth, the final water times concentration
  m_imSediment->applySqrtNormalization(*imToSave[4], tid, minC,
    maxC, minC, maxC);

  //! Sixth, the raw landscape
  if (doRaw) {
    imToSave[5]->multiply(factor, tid);
  }

  //! If meshes wanted
  if (m_params->visualizeMesh()) {
    m_imLandscapeInit ->sample(*imToSave[indMeshO], nbPoints, tid);
    imToSave[indMeshO]->multiply(1, tid);
    m_imLandscape     ->sample(*imToSave[indMeshF], nbPoints, tid);
    imToSave[indMeshF]->multiply(1, tid);
  }

  //! If twi wanted
  if (m_params->computeTWI()) {
    imTwiO.applySqrtNormalization(*imToSave[indTwiO], tid, minTwiO, maxTwiO, minTwiO, maxTwiO);
    imTwiF.applySqrtNormalization(*imToSave[indTwiF], tid, minTwiF, maxTwiF, minTwiF, maxTwiF);
  }

#ifdef _OPENMP
}
#endif // _OPENMP

  //! Save all the images
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t n = 0; n < nbImages; n++) {
    imToSave[n]->write(names[n].c_str(), applyQuad);
  }

  //! saveAsObj is parallelized itself
  if (m_params->visualizeMesh()) {
    imToSave[indMeshO]->saveAsObj(names[indMeshO]);
    imToSave[indMeshF]->saveAsObj(names[indMeshF]);
    if (m_params->visualizeTopo() != NULL) {
      imToSave[0]->saveAsObj(outFolder + "originalLandscape_color.obj",
        *m_imWaterLvlInit, *imToSave[indTopoO], m_params->oceanLevel(), 0.3);
      imToSave[2]->saveAsObj(outFolder + "finalLandscape_color.obj",
        *m_imWaterLvl, *imToSave[indTopoF], m_params->oceanLevel(), 0.3);
    }
  }

  //! Release memory
  for (size_t n = 0; n < imToSave.size(); n++) {
    delete imToSave[n];
    imToSave[n] = NULL;
  }
  imToSave.clear();
  names.clear();

  //! Print time informations
  if (verbose) getTime(m_time, "LEM save final images                      ");
}


//! Compute some stats.
void LEM::computeStats() {

  //! For convenience
  const bool verbose = m_params->verbose();

  //! Iteration reached
  writeValue("nbIter", m_params->outName(), m_iterCurrent);
  if (verbose) {
    cout << "Iter reached  : " << m_iterCurrent << " / " << m_iterMax << endl;
  }

  //! Percentage of mass left
  const float percReached = 1.f - m_massCurrent / m_massOriginal;
  writeValue("Percentage", m_params->outName(), percReached);
  if (verbose) {
    cout << "Mass reached  : " << 100.f * percReached << " %" << endl;
  }

  //! Current mass of water
  float waterMass, waterMax, waterNorm;
  this->integrateWaterOverOceanLevel(waterMass, waterMax, waterNorm);
  writeValue("Water", m_params->outName(), waterMass);
  writeValue("WaterMax", m_params->outName(), waterMax);
  writeValue("WaterNorm", m_params->outName(), waterNorm);
  if (verbose) {
    cout << "Water mass    : " << waterMass << endl;
  }

  //! Stability norm
  float stabilityNorm;
  this->computeStabilityNorm(stabilityNorm);
  writeValue("StabilityNorm", m_params->outName(), stabilityNorm);
  if (verbose) {
    cout << "Stability norm: " << stabilityNorm * 1000 << " x 10-3" << endl;
  }

  //! Second stability norm
  float newNorm;
  this->computeNorm(newNorm);
  writeValue("Norm", m_params->outName(), newNorm);
  if (verbose) {
    cout << "New norm      : " << newNorm * 1000 << " x 10-3" << endl;
  }

  //! Stability curve
  if (m_params->stabCurv()) {
    cout << "iter current = " << m_iterCurrent << endl;
    std::vector<float> norms(m_iterCurrent);
    for (size_t n = 0; n < m_iterCurrent; n++) {
      float norm = 0.f;
      float nb   = 0.f;
      for (size_t k = 0; k < m_nbThreads; k++) {
        norm += m_norms[n * m_nbThreads + k];
        nb   += m_nbs  [n * m_nbThreads + k];
      }
      norms[n] = sqrt(norm / nb);
    }
    writeCurve("StabilityCurve", m_params->outName(), norms);
  }

  //! Automatic selection of Erosion parameter
  if (m_params->autoE()) {
    if (verbose) {
      cout << "E estimated to: " << m_params->e() << endl;
    }
    writeValue("Erosion", m_params->outName(), m_params->e());
  }

  //! Total uplift
  writeValue("Uplift", m_params->outName(), m_upliftAll);
  if (verbose) {
    cout << "Uplift total  : " << m_upliftAll << endl;
  }
}


//! Adjust the time step parameter.
void LEM::adjustTimeStep() {

  //! Compute the advancement
  const float ratioIter = 1.f - (float(m_iterMax) - float(m_iterCurrent)) / float(m_iterMax);
  const float ratioMass = (m_massCurrent - m_massOriginal) / (m_massTarget - m_massOriginal);
  m_ratio = std::max(ratioIter, ratioMass);

  //! Use the boost
  if (m_params->boostIter() > 0.5f) {

    //! Initialization
    float dt = m_params->boostIter();
    float ratioNb = 1.f;
    std::vector<std::pair<float, float> > values;
    values.push_back(std::make_pair(0.f, dt));

    //! Pre-compute values
    while (dt > 0.5f) {
      dt = std::max(0.5f, dt / 2);
      ratioNb /= 2.f;
      values.push_back(std::make_pair(1.f - ratioNb, dt));
    }
    values.push_back(std::make_pair(1.f, 0.5f));

    //! Adjust timeStep according to the advancement
    bool doContinue = true;
    size_t ind = 0;
    while (doContinue && ind < values.size()) {
      if (m_ratio >= values[ind].first) {
        ind++;
      }
      else {
        doContinue = false;
      }
    }
    dt = values[ind - 1].second;

    //! Update the new time step
    m_evolutionCreep                ->setTimeStep(dt);
    m_evolutionLandscape            ->setTimeStep(dt);
    m_evolutionWaterAndSedimentation->setTimeStep(dt);
  }
}


//! Initialize the informations parameters
void LEM::initializeInfos() {

  //! Check the original EmergedMass
  m_massCurrent = m_massOriginal;

  //! Deduce the stopping criteria, aka the target of mass removed
  m_massTarget  = m_massOriginal * (float(1) - m_params->percLand());

  //! Iteration parameters
  m_iterCurrent = 0;
  m_iterMax     = m_params->nbIter();
  m_nbCheck     = 1;

  //! Uplift values
  m_uplift    = 0.f;
  m_upliftAll = 0.f;
}


//! Update the informations parameters
void LEM::updateInfos(
  ProgressBar& i_bar,
  size_t& io_nbBar) {

  //! Get the current mass
  m_massCurrent = 0.f;
  for (size_t n = 0; n < m_nbThreads; n++) {
    m_massCurrent += m_masses[n];
  }

  //! Adjust the parameters values
  adjustTimeStep();

  //! Print informations
  if (m_params->verbose()) {
    const size_t nb = std::min(98, (int) floor(m_ratio * 100.f));
    if (io_nbBar < nb) {
      io_nbBar = nb;
      i_bar.update(io_nbBar);
    }
  }

  //! Update the counter
  m_iterCurrent++;
}


//! Compute the stability norm over the original landscape z0 and the evolved
//! landscape z1.
void LEM::computeStabilityNorm(
  float& o_norm) {

  //! Initialization
  double sumS = 0.0, sumN = 0.0;

  //! Call the parallelized version of this function
#ifdef _OPENMP
#pragma omp parallel
{
  const size_t tid = omp_get_thread_num();
#else
  const size_t tid = 0;
#endif // _OPENMP

  this->computeStabilityNorm_tid(tid);

#ifdef _OPENMP
}
#endif // _OPENMP

  //! Get the result
  for (size_t n = 0; n < m_nbThreads; n++) {
    sumS += m_norms[m_iterCurrent * m_nbThreads + n];
    sumN += m_nbs  [m_iterCurrent * m_nbThreads + n];
  }
  o_norm = sqrt(sumS / sumN);
}


//! Compute the stability norm over the original landscape z0 and the evolved
//! landscape z1.
void LEM::computeStabilityNorm_tid(
  const size_t p_tid) {

  //! Initialization
  double sumS = 0.0, sumN = 0.0, sumH = 0.0;
  const float oceanLvl = m_params->oceanLevel();
  const int hBegin     = m_imLandscape->nHeight(p_tid);
  const int hEnd       = m_imLandscape->nHeight(p_tid + 1);

#ifdef __AVX__
  __m256 zsum = _mm256_setzero_ps();
  __m256 znb  = _mm256_setzero_ps();
  __m256 zh   = _mm256_setzero_ps();
  const __m256 z4 = _mm256_set1_ps(0.25f);
  const __m256 z0 = _mm256_setzero_ps();
  const __m256 z1 = _mm256_set1_ps(1.f);
  const __m256 zlvl = _mm256_set1_ps(oceanLvl);
#endif // __AVX__
#ifdef __SSE__
  __m128 xsum = _mm_setzero_ps();
  __m128 xnb  = _mm_setzero_ps();
  __m128 xh   = _mm_setzero_ps();
  const __m128 x4 = _mm_set1_ps(0.25f);
  const __m128 x0 = _mm_setzero_ps();
  const __m128 x1 = _mm_set1_ps(1.f);
  const __m128 xlvl = _mm_set1_ps(oceanLvl);
#endif // __SSE__

  //! Loop over the pixels
  for (int i = hBegin; i < hEnd; i++) {
    const float* iL0T = m_imLandscapeInit->getPtr(0, i - 1);
    const float* iL0C = m_imLandscapeInit->getPtr(0, i    );
    const float* iL0B = m_imLandscapeInit->getPtr(0, i + 1);
    const float* iL1T = m_imLandscape    ->getPtr(0, i - 1);
    const float* iL1C = m_imLandscape    ->getPtr(0, i    );
    const float* iL1B = m_imLandscape    ->getPtr(0, i + 1);
    int j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < int(m_width) - 8; j += 8) {
      const __m256 h0x = _mm256_loadu_ps(iL0C + j + 1)
                       - _mm256_loadu_ps(iL0C + j - 1);
      const __m256 h1x = _mm256_loadu_ps(iL1C + j + 1)
                       - _mm256_loadu_ps(iL1C + j - 1);
      const __m256 h0y = _mm256_loadu_ps(iL0B + j)
                       - _mm256_loadu_ps(iL0T + j);
      const __m256 h1y = _mm256_loadu_ps(iL1B + j)
                       - _mm256_loadu_ps(iL1T + j);
      const __m256 mask = _mm256_cmp_ps(_mm256_loadu_ps(iL1C + j), zlvl, _CMP_GT_OS);
      zsum += applyMask256_ps(mask, ((h0x - h1x) * (h0x - h1x) + (h0y - h1y) * (h0y - h1y)) * z4, z0);
      znb  += applyMask256_ps(mask, z1, z0);
      zh   += applyMask256_ps(mask, (h0x * h0x + h1x * h1x + h0y * h0y + h1y * h1y) * z4, z0);
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < int(m_width) - 4; j += 4) {
      const __m128 h0x = _mm_loadu_ps(iL0C + j + 1)
                       - _mm_loadu_ps(iL0C + j - 1);
      const __m128 h1x = _mm_loadu_ps(iL1C + j + 1)
                       - _mm_loadu_ps(iL1C + j - 1);
      const __m128 h0y = _mm_loadu_ps(iL0B + j)
                       - _mm_loadu_ps(iL0T + j);
      const __m128 h1y = _mm_loadu_ps(iL1B + j)
                       - _mm_loadu_ps(iL1T + j);
      const __m128 mask = _mm_cmpgt_ps(_mm_loadu_ps(iL1C + j), xlvl);
      xsum += applyMask_ps(mask, ((h0x - h1x) * (h0x - h1x) + (h0y - h1y) * (h0y - h1y)) * x4, x0);
      xnb  += applyMask_ps(mask, x1, x0);
      xh   += applyMask_ps(mask, (h0x * h0x + h1x * h1x + h0y * h0y + h1y * h1y) * x4, x0);
    }
#endif // __SSE__

    //! Normal version
    for (; j < int(m_width); j++) {
      if (iL1C[j] > oceanLvl) {

        //! Horizontal gradient
        const double h0x = iL0C[j + 1] - iL0C[j - 1];
        const double h1x = iL1C[j + 1] - iL1C[j - 1];

        //! Vertical gradient
        const double h0y = iL0B[j] - iL0T[j];
        const double h1y = iL1B[j] - iL1T[j];

        //! Compute the Sobolev norm
        sumS += ((h0x - h1x) * (h0x - h1x) + (h0y - h1y) * (h0y - h1y)) * 0.25;
        sumH += (h0x * h0x + h1x * h1x + h0y * h0y + h1y * h1y) * 0.25;
        sumN += 1.0;
      }
    }
  }

  //! Obtain the norm
#ifdef __AVX__
  float zval[8]; _mm256_storeu_ps(zval, zsum);
  float zvn [8]; _mm256_storeu_ps(zvn, znb);
  float zvh [8]; _mm256_storeu_ps(zvh, zh);
  sumS += zval[0] + zval[1] + zval[2] + zval[3] + zval[4] + zval[5] + zval[6] + zval[7];
  sumN += zvn[0] + zvn[1] + zvn[2] + zvn[3] + zvn[4] + zvn[5] + zvn[6] + zvn[7];
  sumH += zvh[0] + zvh[1] + zvh[2] + zvh[3] + zvh[4] + zvh[5] + zvh[6] + zvh[7];
#endif // __AVX__
#ifdef __SSE__
  float xval[4]; _mm_storeu_ps(xval, xsum);
  float xvn [4]; _mm_storeu_ps(xvn, xnb);
  float xvh [4]; _mm_storeu_ps(xvh, xh);
  sumS += xval[0] + xval[1] + xval[2] + xval[3];
  sumN += xvn[0] + xvn[1] + xvn[2] + xvn[3];
  sumH += xvh[0] + xvh[1] + xvh[2] + xvh[3];
#endif // __SSE__
  m_norms[m_iterCurrent * m_nbThreads + p_tid] = 2 * sumS;
  //m_nbs  [m_iterCurrent * m_nbThreads + p_tid] = sumN;
  m_nbs  [m_iterCurrent * m_nbThreads + p_tid] = sumH;
}


//! Compute another stability norm.
void LEM::computeNorm(
  float& o_norm) const {

  //! Initialization
  double sumS = 0.0, sumH = 0.0;
  const float oceanLvl = m_params->oceanLevel();

  //! Loop over the pixels
  for (size_t i = 0; i < m_height; i++) {
    const float* iL0 = m_imLandscapeInit->getPtr(0, i);
    const float* iL1 = m_imLandscape->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      if (iL1[j] > oceanLvl) {
        const double v0 = iL0[j];
        const double v1 = iL1[j];
        sumS += (v0 - v1) * (v0 - v1);
        sumH += v0 * v0 + v1 * v1;
      }
    }
  }

  //! Return the norm
  o_norm = sqrt(2 * sumS / sumH);
}


//! Compute the current uplift term.
void LEM::computeUplift(
  const double* i_means) {

  //! Compute the average uplift to perform
  double meanCurrent = 0.0;
  for (size_t k = 0; k < m_nbThreads; k++) {
    meanCurrent += i_means[k];
  }
  meanCurrent /= double(m_width * m_height);

  //! Update the values
  m_uplift = m_meanOriginal - meanCurrent;
  m_upliftAll += m_uplift;

  //! Compute the uniform uplift
  if (m_params->autoUplift() == 1) {
    for (size_t i = 0; i < m_height; i++) {
      float* iU = m_imUplift->getPtr(0, i);

      for (size_t j = 0; j < m_width; j++) {
        iU[j] = m_uplift;
      }
    }
  }
  else {
    double mean = 0.0;
    for (size_t i = 0; i < m_height; i++) {
      float* iU = m_imUplift->getPtr(0, i);
      const float* iL = m_imLandscape->getPtr(0, i);

      for (size_t j = 0; j < m_width; j++) {
        const double value = sqrt(sqrt((double(std::max(iL[j], 0.f)))));
        mean += value;
        iU[j] = value;
      }
    }
    mean /= double(m_width * m_height);

    const double norm = m_uplift / mean;
    for (size_t i = 0; i < m_height; i++) {
      float* iU = m_imUplift->getPtr(0, i);

      for (size_t j = 0; j < m_width; j++) {
        iU[j] *= norm;
      }
    }
  }
}


//! Auto selection of the parameters
void LEM::autoSelect() {

  //! Be sure that the water has been initialized
  if (!m_waterHasBeenInitialized) {
    this->waterInitialization();
  }

  //! Save the current state
  Image imLand             = *m_imLandscape;
  Image imWater            = *m_imWaterLvl;
  Image imSediment         = *m_imSediment;
  const bool verbose       = m_params->verbose();

  //! Compute the range of estimation
  const float vs = m_params->s();
  const float vm = m_params->expM();
  const float vn = m_params->expN();
  const float ds = m_params->deltaS();
  const float dm = m_params->deltaM();
  const float dn = m_params->deltaN();
  std::vector<float> valS, valM, valN;
  /*getValues(vs * 1e6, valS);
  getValues(vm, valM);
  getValues(vn, valN);
  for (size_t ns = 0; ns < valS.size(); ns++) {
    valS[ns] *= 1e-6;
  }*/
  valM.push_back(vm - dm); valM.push_back(vm); valM.push_back(vm + dm);
  valN.push_back(vn - dn); valN.push_back(vn); valN.push_back(vn + dn);
  valS.push_back(vs - ds); valS.push_back(vs); valS.push_back(vs + ds);

  //! Initialization
  float bestNorm = DBL_MAX;
  float bestS = vs, bestM = vm, bestN = vn, bestE = m_params->e();
  m_params->setV(false);
  m_params->setW(false);

  //! Progress bar
  ProgressBar bar;
  bar.init(valM.size() * valN.size() * valS.size(), 0);

  //! Loops over all the values to look at
  for (size_t im = 0, k = 0; im < valM.size(); im++) {
    m_params->setM(valM[im]);

    for (size_t in = 0; in < valN.size(); in++) {
      m_params->setN(valN[in]);

      for (size_t is = 0; is < valS.size(); is++, k++) {
        m_params->setS(valS[is]);

        //! Set the image
        *m_imLandscape = imLand;
        *m_imWaterLvl  = imWater;
        *m_imSediment  = imSediment;

        //! Set the parameters
        m_evolutionCreep                ->init(m_width, m_height, m_border, *m_params);
        m_evolutionLandscape            ->init(m_width, m_height, m_border, *m_params);
        m_evolutionWaterAndSedimentation->init(m_width, m_height, m_border, *m_params);

        //! Adjust the erosion coefficient
        this->runAutoErosion();

        //! Run the main iterations
        this->runMainIterationsCore();

        //! Compute the norm
        float norm;
        this->computeStabilityNorm(norm);

        //! Keep best parameters
        if (bestNorm > norm) {
          bestNorm = norm;
          bestS = m_params->s();
          bestM = m_params->expM();
          bestN = m_params->expN();
          bestE = m_params->e();
        }

        //! Display the progress
        if (verbose) bar.update(k);
      }
    }
  }

  //! For convenience
  m_params->setV(verbose);

  //! Save the estimated values
  m_params->setS(bestS);
  writeValue("BestS", m_params->outName(), m_params->s() * 1e6);
  m_params->setE(bestE);
  writeValue("BestE", m_params->outName(), m_params->e());
  m_params->setM(bestM);
  writeValue("BestM", m_params->outName(), m_params->expM());
  m_params->setN(bestN);
  writeValue("BestN", m_params->outName(), m_params->expN());

  //! Set the image
  *m_imLandscape = imLand;
  *m_imWaterLvl  = imWater;
  *m_imSediment  = imSediment;

  //! Set the parameters
  m_evolutionCreep                ->init(m_width, m_height, m_border, *m_params);
  m_evolutionLandscape            ->init(m_width, m_height, m_border, *m_params);
  m_evolutionWaterAndSedimentation->init(m_width, m_height, m_border, *m_params);

  //! Run the main iterations
  this->runMainIterationsCore();
}


//! Un-zoom the images if needed.
void LEM::applyZoomOut() {

  //! For convenience
  const bool doRaw = m_params->blur() > 0.f || m_params->noise() > 0.f;

  //! Initialization of all the images
  Image imLandscape     = *m_imLandscape;
  Image imLandscapeInit = *m_imLandscapeInit;
  Image imWaterLvl      = *m_imWaterLvl;
  Image imWaterLvlInit  = *m_imWaterLvlInit;
  Image imSediment      = *m_imSediment;
  Image imLandscapeRaw;
  if (doRaw) {
    imLandscapeRaw = *m_imLandscapeRaw;
  }

  //! Get the factor of sampling
  const size_t sizeOut = m_width * m_height / m_params->zoom() / m_params->zoom();

  //! Get the new size (must be computed as in sample function from LibImages)
  const size_t width  = sqrt(float(sizeOut * m_width ) / float(m_height));
  const size_t height = sqrt(float(sizeOut * m_height) / float(m_width ));
  m_width  = width;
  m_height = height;

  //! Initialize the images before the parallelization
  m_imLandscape    ->init(m_width, m_height, 1, 1);
  m_imLandscapeInit->init(m_width, m_height, 1, 1);
  m_imWaterLvl     ->init(m_width, m_height, 1, 1);
  m_imWaterLvlInit ->init(m_width, m_height, 1, 1);
  m_imSediment     ->init(m_width, m_height, 1, 1);
  if (doRaw) {
    m_imLandscapeRaw->init(m_width, m_height, 1, 1);
  }

#ifdef _OPENMP
#pragma omp parallel
{
  const size_t tid = omp_get_thread_num();
#else
  const size_t tid = 0;
#endif // _OPENMP

  //! Apply the zoom-out
  imLandscape    .sample(*m_imLandscape              , sizeOut, tid);
  imLandscapeInit.sample(*m_imLandscapeInit          , sizeOut, tid);
  imWaterLvl     .sample(*m_imWaterLvl               , sizeOut, tid);
  imWaterLvlInit .sample(*m_imWaterLvlInit           , sizeOut, tid);
  imSediment     .sample(*m_imSediment, sizeOut, tid);
  if (doRaw) {
    imLandscapeRaw.sample(*m_imLandscapeRaw, sizeOut, tid);
  }

#ifdef _OPENMP
}
#endif // _OPENMP
}


//! Release memory.
void LEM::releaseMemory() {

  //! Release memory
  memfree(m_masses);
  memfree(m_norms);
  memfree(m_nbs);

  //! Release images
  if (m_imWaterLvl != NULL) {
    delete m_imWaterLvl;
    m_imWaterLvl = NULL;
  }
  if (m_imLandscape != NULL) {
    delete m_imLandscape;
    m_imLandscape = NULL;
  }
  if (m_imSediment != NULL) {
    delete m_imSediment;
    m_imSediment = NULL;
  }
  if (m_imWaterLvlInit != NULL) {
    delete m_imWaterLvlInit;
    m_imWaterLvlInit = NULL;
  }
  if (m_imLandscapeInit != NULL) {
    delete m_imLandscapeInit;
    m_imLandscapeInit = NULL;
  }
  if (m_imLandscapeRaw != NULL) {
    delete m_imLandscapeRaw;
    m_imLandscapeRaw = NULL;
  }
  if (m_imUplift != NULL) {
    delete m_imUplift;
    m_imUplift = NULL;
  }

  //! Release memory for Evolution schemes
  if (m_evolutionWater != NULL) {
    delete m_evolutionWater;
    m_evolutionWater = NULL;
  }
  if (m_evolutionCreep != NULL) {
    delete m_evolutionCreep;
    m_evolutionCreep = NULL;
  }
  if (m_evolutionLandscape != NULL) {
    delete m_evolutionLandscape;
    m_evolutionLandscape = NULL;
  }
  if (m_evolutionWaterAndSedimentation != NULL) {
    delete m_evolutionWaterAndSedimentation;
    m_evolutionWaterAndSedimentation = NULL;
  }
  
  //! Release rain matrix
  if (this->m_rainMatrix != NULL)
    delete[] this->m_rainMatrix;

  //! Release memory for the parameters
  if (m_params != NULL) {
    delete m_params;
    m_params = NULL;
  }
}





