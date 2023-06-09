/**
 * @file TWI.cpp
 *
 * @brief Functions to compute the Topographic Wetness Index (Beven index).
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <malloc.h>
#include <omp.h>
#include <queue>
#include <iostream>
#include <sstream>
#include <algorithm>


//! Local includes
#include "TWI.h"
#include "../Utilities/Utilities.h"
#include "ConnectedComponent.h"


using namespace std;



//! Default constructor
TWI::TWI() :

  //! Miscellaneous
  m_width (0),
  m_height(0),
  m_border(0),

  //! Parameters
  m_useDiag(true),

  //! Images
  m_imOriginal     (NULL),
  m_im             (NULL),
  m_imRegularized  (NULL),
  m_imCatchmentArea(NULL) {
}


//! Copy constructor.
TWI::TWI(
  const TWI& i_twi) :

  //! Miscellaneous
  m_width (i_twi.m_width),
  m_height(i_twi.m_height),
  m_border(i_twi.m_border),

  //! Parameters
  m_useDiag(i_twi.m_useDiag),

  //! Images
  m_imOriginal     (new Image(*i_twi.m_imOriginal)),
  m_im             (new Image(*i_twi.m_im)),
  m_imRegularized  (new Image(*i_twi.m_imRegularized)),
  m_imCatchmentArea(new Image(*i_twi.m_imCatchmentArea)) {
}


//! Default destructor
TWI::~TWI() {
  releaseMemory();
}


//! Operator overload.
TWI& TWI::operator=(
  const TWI& i_twi) {

  if (&i_twi == this) {
    return *this;
  }

  releaseMemory();

  new (this) TWI(i_twi);
  return *this;
}


//! Perform the water initialization.
void TWI::computeTWI(
  Image const& i_im,
  Image &o_im) {

  //! Initialize the images
  this->init(i_im);

  //! Be carefull, result will be in marker and thus will erase the data in the
  //! marker
  this->applyGeodesicReconstruction();

  //! Perform the landscape regularization
  this->performRegularization();
  //m_imRegularized = new Image(*m_imOriginal);

  //! Compute the catchment area
  this->computeCatchmentArea();

  //! Compute the final index
  this->computeIndex(o_im);
}


//! Initialize the images.
void TWI::init(
  Image const& i_im) {

  //! Release the memory before initialization
  releaseMemory();

  //! Size
  m_width  = i_im.width ();
  m_height = i_im.height();
  m_border = 1;

  //! Images initialization
  m_imOriginal      = new Image(m_width, m_height, 1, m_border);
  m_im              = new Image(m_width, m_height, 1, m_border);
  m_imCatchmentArea = new Image(m_width, m_height, 1, m_border);

  //! Fill the images
  for (size_t i = 0; i < m_height; i++) {
    const float* iI = i_im.getPtr(0, i);
    float* oO = m_imOriginal->getPtr(0, i);
    float* oI = m_im->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      oO[j] = iI[j];

      //! Set pixel to 0 if you want it to be minimal
      //! Otherwise set it to the max value
      if ((i == 0) || (j == 0) || (i == m_height - 1) || (j == m_width - 1)) {
        oI[j] = iI[j];
      }
      else {
        oI[j] = 1.f;
      }
    }
  }
}


//! Apply a geodesic reconstruction by in-place erosion.
void TWI::applyGeodesicReconstruction() {

  //! Initialization
  size_t nbIter = 0;
  const size_t maxIter = 3;
  bool isStable = false;
  std::queue< std::pair<int, int> > queue;

  //! Sequential version
  while (!isStable && nbIter < maxIter) {

    //! Variable update
    nbIter++;
    isStable = true;

    //! Raster scan
    for (int i = 0; i < (int) m_height; i++) {
      const float* iK = m_imOriginal->getPtr(0, i);
      float* iR0 = m_im->getPtr(0, i);
      const float* iR1 = m_im->getPtr(0, i - 1);

      for (int j = 0; j < (int) m_width; j++) {
        const float value = std::max(iK[j], std::min(iR0[j],
                                            std::min(iR0[j - 1], iR1[j])));

        if (value != iR0[j]) {
          isStable = false;
          iR0[j] = value;
        }
      }
    }

    //! Anti-raster scan
    for (int i = int(m_height) - 1; i >= 0; i--) {
      const float* iK = m_imOriginal->getPtr(0, i);
      float* iR0 = m_im->getPtr(0, i);
      const float* iR1 = m_im->getPtr(0, i + 1);

      for (int j = int(m_width) - 1; j >= 0; j--) {
        const float value = std::max(iK[j], std::min(iR0[j],
                                            std::min(iR0[j + 1], iR1[j])));

        if (value != iR0[j]) {
          isStable = false;
          iR0[j] = value;
        }
      }
    }
  }

  //! Update the queue with the pixel to treat
  for (int i = 0; i < int(m_height); i++) {
    const float* iK0 = m_imOriginal->getPtr(0, i);
    const float* iK1 = m_imOriginal->getPtr(0, i + 1);
    float* iR0 = m_im->getPtr(0, i);
    const float* iR1 = m_im->getPtr(0, i + 1);

    for (int j = 0; j < int(m_width); j++) {
      if (iR0[j] < iR1[j] && iK1[j] < iR1[j]) {
        queue.push(std::make_pair(i, j));
      }
      if (iR0[j] < iR0[j + 1] && iK0[j + 1] < iR0[j + 1]) {
        queue.push(std::make_pair(i, j));
      }
    }
  }

  //! We now perform the queue-based version
  while (!queue.empty()) {
    std::pair<int,int> p = queue.front();
    queue.pop();
    const int i = p.first;
    const int j = p.second;
    const float* iT = m_imOriginal->getPtr(0, i - 1);
    const float* iC = m_imOriginal->getPtr(0, i);
    const float* iB = m_imOriginal->getPtr(0, i + 1);
    float* oT = m_im->getPtr(0, i - 1);
    float* oC = m_im->getPtr(0, i);
    float* oB = m_im->getPtr(0, i + 1);
    const float val = oC[j];

    if (val < oC[j - 1] && oC[j - 1] != iC[j - 1]) {
      oC[j - 1] = std::max(val, iC[j - 1]);
      queue.push(std::make_pair(i, j - 1));
    }

    if (val < oC[j + 1] && oC[j + 1] != iC[j + 1]) {
      oC[j + 1] = std::max(val, iC[j + 1]);
      queue.push(std::make_pair(i, j + 1));
    }

    if (val < oT[j] && oT[j] != iT[j]) {
      oT[j] = std::max(val, iT[j]);
      queue.push(std::make_pair(i - 1, j));
    }

    if (val < oB[j] && oB[j] != iB[j]) {
      oB[j] = std::max(val, iB[j]);
      queue.push(std::make_pair(i + 1, j));
    }

    if (m_useDiag) {
      if (val < oT[j - 1] && oT[j - 1] != iT[j - 1]) {
        oT[j - 1] = std::max(val, iT[j - 1]);
        queue.push(std::make_pair(i - 1, j - 1));
      }

      if (val < oT[j + 1] && oT[j + 1] != iT[j + 1]) {
        oT[j + 1] = std::max(val, iT[j + 1]);
        queue.push(std::make_pair(i - 1, j + 1));
      }

      if (val < oB[j - 1] && oB[j - 1] != iB[j - 1]) {
        oB[j - 1] = std::max(val, iB[j - 1]);
        queue.push(std::make_pair(i + 1, j - 1));
      }

      if (val < oB[j + 1] && oB[j + 1] != iB[j + 1]) {
        oB[j + 1] = std::max(val, iB[j + 1]);
        queue.push(std::make_pair(i + 1, j + 1));
      }
    }
/*
    const float maskT = m_imOriginal->getPtr(0, i - 1)[j];
    const float maskB = m_imOriginal->getPtr(0, i + 1)[j];
    const float maskL = m_imOriginal->getPtr(0, i)[j - 1];
    const float maskR = m_imOriginal->getPtr(0, i)[j + 1];

    float* iRT = m_im->getPtr(0, i - 1);
    float* iRC = m_im->getPtr(0, i);
    float* iRB = m_im->getPtr(0, i + 1);

    const float marker       = iRC[j];
    const float markerTop    = iRT[j];
    const float markerBottom = iRB[j];
    const float markerLeft   = iRC[j - 1];
    const float markerRight  = iRC[j + 1];

    if (marker < markerLeft && markerLeft != maskLeft) {
      iRC[j - 1] = std::max(marker, maskLeft);
      queue.push(std::make_pair(i, j - 1));
    }
    if (marker < markerTop && markerTop != maskTop) {
      iRT[j] = std::max(marker, maskTop);
      queue.push(std::make_pair(i - 1, j));
    }
    if (marker < markerRight && markerRight != maskRight) {
      iRC[j + 1] = std::max(marker, maskRight);
      queue.push(std::make_pair(i, j + 1));
    }
    if (marker < markerBottom && markerBottom != maskBottom) {
      iRB[j] = std::max(marker, maskBottom);
      queue.push(std::make_pair(i + 1, j));
    }*/
  }
}


//! Perform a regularization of the landscape.
void TWI::performRegularization() {

  //! For connected components extraction the border is set to -1.
  m_im->setBorder(-1);

  //! Initialization of the regularized version
  m_imRegularized = new Image(*m_im);

  //! Initialization of the list of connected components
  std::vector<ConnectedComponent> allCC;
  allCC.reserve(m_height * m_width);

  //! Image to keep track of already-processed pixels
  Image isProcessed(m_width, m_height, 1, m_border);

  //! Extract the connected components
  for (size_t i = 0; i < m_height; i++) {
    const float* iP = isProcessed.getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      if (iP[j] == 0) {
        ConnectedComponent cc;
        cc.extract(*m_im, isProcessed, i, j, m_useDiag);
        allCC.push_back(cc);
      }
    }
  }

  //! For each connected component, update its flag
  for (size_t n = 0; n < allCC.size(); n++) {
    for (size_t l = 0; l < allCC[n].getNbElements(); l++) {

      //! Get the current position
      const std::pair<int, int> s = allCC[n].getElement(l);

      //! Check for lower neighbor
      if (!m_im->hasLowerNeighbor(s.first, s.second, 0, m_useDiag)) {
        allCC[n].setProcess(true);
      }
    }
  }

  //! We process all connected components that are flagged true
  //! i.e, there is at least one point in the connected component that does not
  //! have a lower neighbor
  for (size_t n = 0; n < allCC.size(); n++) {
    if (allCC[n].getProcess()) {
      allCC[n].process(*m_im, *m_imRegularized, m_useDiag);
    }
  }
}


//! Compute the catchment area for each pixel.
void TWI::computeCatchmentArea() {

  //! For convenience
  const float eps  = 1e-15f;
  const float rain = 1e-7;
  const size_t nb = m_useDiag ? 8 : 4;

  //! Initialize the borders of both landscape
  m_imRegularized->setBorder(-1);

  //! Sort the sites / pixel positions
  std::vector<std::pair<int, int> > sites(m_width * m_height);
  for (size_t i = 0; i < m_height; i++) {
    for (size_t j = 0; j < m_width; j++) {
      sites[i * m_width + j] = make_pair(i, j);
    }
  }
  std::sort(sites.begin(), sites.end(), topologicalSort(*m_imRegularized));

  //! Flow passing through edges
  //! Each site contains the set of incoming edges with flow on it
  const size_t nh = m_height + 2;
  const size_t nw = m_width  + 2;
  std::vector<Flow  > phi(nw * nh * nb, Flow(0, 0, 0));
  std::vector<size_t> phiNb(nw * nh, 0);

  //! We process sites from higher to lower
  //! During the process, the set of incoming arcs are updated on the fly for
  //! lower pixels
  for (size_t k = 0; k < m_width * m_height; k++) {

    //! Current position
    const int i   = sites[k].first;
    const int j   = sites[k].second;
    const int ind = (i + 1) * nw + j + 1;

    //! We first update the incoming flow
    float incomingFlow = 0;
    for (size_t n = 0; n < phiNb[ind]; n++) {
      incomingFlow += phi[nb * ind + n].flow;
    }
    incomingFlow += rain;

    //! We need to get the number of lower sites
    float nbLowerSites = 0;
    const float* iRT = m_imRegularized->getPtr(0, i - 1);
    const float* iRC = m_imRegularized->getPtr(0, i);
    const float* iRB = m_imRegularized->getPtr(0, i + 1);
    const float val = iRC[j];
    nbLowerSites += val > iRC[j - 1];
    nbLowerSites += val > iRC[j + 1];
    nbLowerSites += val > iRT[j];
    nbLowerSites += val > iRB[j];
    if (m_useDiag) {
      nbLowerSites += val > iRT[j - 1];
      nbLowerSites += val > iRT[j + 1];
      nbLowerSites += val > iRB[j - 1];
      nbLowerSites += val > iRB[j + 1];
    }

    //! Remark: nb_lower_sites cannot be 0 since the landscape have been regularized
    //! Now we can update phi
    //! more precisely, each site contains the set of incoming edges
    const Flow flow(i, j, incomingFlow / nbLowerSites);
    if (val > iRC[j - 1]) {
      size_t* iN = &phiNb[ind - 1];
      phi[nb * (ind - 1) + ((*iN)++)] = flow;
    }
    if (val > iRC[j + 1]) {
      size_t* iN = &phiNb[ind + 1];
      phi[nb * (ind + 1) + ((*iN)++)] = flow;
    }
    if (val > iRT[j]) {
      size_t* iN = &phiNb[ind - nw];
      phi[nb * (ind - nw) + ((*iN)++)] = flow;
    }
    if (val > iRB[j]) {
      size_t* iN = &phiNb[ind + nw];
      phi[nb * (ind + nw) + ((*iN)++)] = flow;
    }
    if (m_useDiag) {
      if (val > iRT[j - 1]) {
        size_t* iN = &phiNb[ind - 1 - nw];
        phi[nb * (ind - 1 - nw) + ((*iN)++)] = flow;
      }
      if (val > iRT[j + 1]) {
        size_t* iN = &phiNb[ind + 1 - nw];
        phi[nb * (ind + 1 - nw) + ((*iN)++)] = flow;
      }
      if (val > iRB[j - 1]) {
        size_t* iN = &phiNb[ind - 1 + nw];
        phi[nb * (ind - 1 + nw) + ((*iN)++)] = flow;
      }
      if (val > iRB[j + 1]) {
        size_t* iN = &phiNb[ind + 1 + nw];
        phi[nb * (ind + 1 + nw) + ((*iN)++)] = flow;
      }
    }
  }

  //! Release memory
  sites.clear();

  //! We now define theta(i,j)
  //! recall that phi(i,j) is stored at site j
  //! later, we will need the set of outgoing edges at site i ; thus, we create
  //! theta by storing outgoing flows.
  //! Also, note we need to consider points that are in the border of the image
  Image theta  (m_width, m_height, 1, m_border);
  Image thetaNb(m_width, m_height, 1, m_border);
  for (int i = -1; i < int(m_height) + 1; i++) {
    const float* iR = m_imRegularized->getPtr(0, i);

    for (int j = -1; j < int(m_width) + 1; j++) {

      //! We should discard the 4 corners
      if ((i == -1            && j == -1) ||
          (i == -1            && j == int(m_width)) ||
          (i == int(m_height) && j == -1) ||
          (i == int(m_height) && j == int(m_width))) {
        continue;
      }

      //! Local value of the regularized landscape
      const float value = iR[j];

      //! We iterate on incoming edges toward site (i,j)
      const size_t ind = (i + 1) * nw + j + 1;
      for (size_t n = 0; n < phiNb[ind]; n++) {
        const int p   = phi[nb * ind + n].posI;
        const int q   = phi[nb * ind + n].posJ;
        const float flow  = phi[nb * ind + n].flow;
        const float steep = m_imRegularized->getPtr(0, p)[q] - value;
        theta  .getPtr(0, p)[q] += flow / (eps + steep);
        thetaNb.getPtr(0, p)[q]++;
      }
    }
  }

  //! Now we compute the level of water for each pixel and normalize it.
  for (size_t i = 0; i < m_height; i++) {
    float* oW = m_imCatchmentArea->getPtr(0, i);
    const float* iT = theta  .getPtr(0, i);
    const float* iN = thetaNb.getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      oW[j] = iT[j] / (iN[j] * 10);
    }
  }
}


//! Compute the index.
void TWI::computeIndex(
  Image& o_im) {

  //! Size initialization
  o_im.init(m_width, m_height, 1);

  //! For convenience
  const float eps = 1e-15f;
  const float norm = 1.f / sqrtf(2.f);

  for (size_t i = 0; i < m_height; i++) {
    const float* iT = m_imRegularized->getPtr(0, i - 1);
    const float* iC = m_imRegularized->getPtr(0, i);
    const float* iB = m_imRegularized->getPtr(0, i + 1);
    const float* iA = m_imCatchmentArea->getPtr(0, i);
    float* oI = o_im.getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {

      //! Compute the maximum slope
      float maxSlope = eps;
      maxSlope = std::max(maxSlope, fabsf(iT[j - 1] - iC[j]) * norm);
      maxSlope = std::max(maxSlope, fabsf(iT[j] - iC[j]));
      maxSlope = std::max(maxSlope, fabsf(iT[j + 1] - iC[j]) * norm);
      maxSlope = std::max(maxSlope, fabsf(iC[j - 1] - iC[j]));
      maxSlope = std::max(maxSlope, fabsf(iC[j + 1] - iC[j]));
      maxSlope = std::max(maxSlope, fabsf(iB[j - 1] - iC[j]) * norm);
      maxSlope = std::max(maxSlope, fabsf(iB[j] - iC[j]));
      maxSlope = std::max(maxSlope, fabsf(iB[j + 1] - iC[j]) * norm);

      //! Get the index
      oI[j] = log(std::max(eps, iA[j] / maxSlope));
      //oI[j] = iA[j];
      //oI[j] = iC[j];
    }
  }
}


//! Release memory
void TWI::releaseMemory() {

  if (m_im != NULL) {
    delete m_im;
    m_im = NULL;
  }

  if (m_imCatchmentArea != NULL) {
    delete m_imCatchmentArea;
    m_imCatchmentArea = NULL;
  }

  if (m_imOriginal != NULL) {
    delete m_imOriginal;
    m_imOriginal = NULL;
  }

  if (m_imRegularized != NULL) {
    delete m_imRegularized;
    m_imRegularized = NULL;
  }
}












