/**
 * @file WaterAndSedimentationEvolution.cpp
 *
 * @brief Class to handle the Water and Sedimentation Evolution scheme.
 *
 * @author Jérôme Darbon <jerome@math.ucla.edu> (original version).
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <malloc.h>
#include <omp.h>
#include <math.h>
#include <assert.h> // to remove
#include <cmath> // to remove


//! Local includes
#include "WaterAndSedimentationEvolution.h"
#include "../Utilities/Utilities.h"


using namespace std;


//! Default constructor
WaterAndSedimentationEvolution::WaterAndSedimentationEvolution() :

  //! Parameters
  m_params(NULL),
  m_dt(0),
  m_rain(0),
  m_rainMatrix(NULL),

  //! Size
  m_width (0),
  m_height(0),

  //! Images
  m_imWaterTmp   (NULL),
  m_imSedimentTmp(NULL) {

}


//! Copy constructor
WaterAndSedimentationEvolution::WaterAndSedimentationEvolution(
  const WaterAndSedimentationEvolution& i_wse) :

  //! Parameters
  m_params(i_wse.m_params),
  m_dt(i_wse.m_dt),
  m_rain(i_wse.m_rain),
  m_rainMatrix(i_wse.m_rainMatrix),

  //! Size
  m_width (i_wse.m_width),
  m_height(i_wse.m_height),

  //! Images
  m_imWaterTmp   (new ImageX<float>(*i_wse.m_imWaterTmp   )),
  m_imSedimentTmp(new ImageX<float>(*i_wse.m_imSedimentTmp)) {

}


//! Default destructor
WaterAndSedimentationEvolution::~WaterAndSedimentationEvolution() {

  //! Release memory
  releaseMemory();
  m_params = NULL;
}

//! Operator overload.
WaterAndSedimentationEvolution& WaterAndSedimentationEvolution::operator=(
  const WaterAndSedimentationEvolution& i_wse) {

  if (&i_wse == this) {
    return *this;
  }

  releaseMemory();

  new (this) WaterAndSedimentationEvolution(i_wse);
  return *this;
}


//! Initialize the images and the parameters
void WaterAndSedimentationEvolution::init(
  const size_t p_width,
  const size_t p_height,
  const size_t p_border,
  Parameters const& p_params) {

  //! Release memory
  releaseMemory();

  //! Parameters
  m_params = &p_params;
  m_dt     = m_params->timeStep();
  m_rain   = m_params->rain();
  m_rainMatrix   = m_params->rainMatrix();
  m_alpha  = p_params.alpha();
  
  //! Sizes
  m_height = p_height;
  m_width  = p_width;

  //! Images
  m_imWaterTmp                   = new ImageX<float>(m_width, m_height, 1, p_border);
  m_imSedimentTmp = new ImageX<float>(m_width, m_height, 1, p_border);
}


//! Perform a one step evolution of the water and sedimentation scheme.
void WaterAndSedimentationEvolution::performOneStepEvolution(
  const Image& i_imLandscape,
  Image& io_imSediment,
  Image& io_imWaterLvl,
  const size_t p_tid) {

  //! For convenience
  const float timeStep = m_dt;
  const float rain     = m_rain;
  const float* rainMatrix     = (float*)m_rainMatrix;
  const size_t hBegin  = i_imLandscape.nHeight(p_tid);
  const size_t hEnd    = i_imLandscape.nHeight(p_tid + 1);
  const float weight   = 1.f / sqrt(2.);
  const float oceanLvl = m_params->oceanLevel();
  const float alpha     = m_alpha;
  
#if defined(__AVX__)
  //! AVX version
  const __m256 zTimeStep = _mm256_set1_ps(timeStep);
  const __m256 zRain     = _mm256_set1_ps(rain);
  const __m256 zOceanLvl = _mm256_set1_ps(oceanLvl);
  const __m256 zWeight   = _mm256_set1_ps(weight);
  const __m256 z0        = _mm256_setzero_ps();
  const __m256 zalpha = _mm256_set1_ps(alpha);
  
#endif
#if defined(__SSE__)
  //! SSE version
  const __m128 xTimeStep = _mm_set1_ps(timeStep);
  const __m128 xRain     = _mm_set1_ps(rain);
  const __m128 xOceanLvl = _mm_set1_ps(oceanLvl);
  const __m128 xWeight   = _mm_set1_ps(weight);
  const __m128 x0        = _mm_setzero_ps();
  const __m128 xalpha   = _mm_set1_ps(alpha);
  
#endif


  //! Compute the transfer
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iWT = io_imWaterLvl.getPtr(0, i - 1);
    const float* iWC = io_imWaterLvl.getPtr(0, i);
    const float* iWB = io_imWaterLvl.getPtr(0, i + 1);
    const float* iLT = i_imLandscape.getPtr(0, i - 1);
    const float* iLC = i_imLandscape.getPtr(0, i);
    const float* iLB = i_imLandscape.getPtr(0, i + 1);
    const float* iCT = io_imSediment.getPtr(0, i - 1);
    const float* iCC = io_imSediment.getPtr(0, i);
    const float* iCB = io_imSediment.getPtr(0, i + 1);
    float* oTw = m_imWaterTmp->getPtr(p_tid, 0, i);
    float* oTc = m_imSedimentTmp->getPtr(p_tid, 0, i);
    size_t j = 0;

#if defined(__AVX__)
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! For convenience: central value
      const __m256 wc = _mm256_loadu_ps(iWC + j);
      const __m256 lc = _mm256_loadu_ps(iLC + j);
      const __m256 cc = _mm256_loadu_ps(iCC + j);

      //! Compute the steep
      const __m256 steepL  = _mm256_loadu_ps(iWC + j - 1) + _mm256_loadu_ps(iLC + j - 1) - wc - lc;
      const __m256 steepR  = _mm256_loadu_ps(iWC + j + 1) + _mm256_loadu_ps(iLC + j + 1) - wc - lc;
      const __m256 steepT  = _mm256_loadu_ps(iWT + j    ) + _mm256_loadu_ps(iLT + j    ) - wc - lc;
      const __m256 steepB  = _mm256_loadu_ps(iWB + j    ) + _mm256_loadu_ps(iLB + j    ) - wc - lc;
      const __m256 steepTL = _mm256_loadu_ps(iWT + j - 1) + _mm256_loadu_ps(iLT + j - 1) - wc - lc;
      const __m256 steepTR = _mm256_loadu_ps(iWT + j + 1) + _mm256_loadu_ps(iLT + j + 1) - wc - lc;
      const __m256 steepBL = _mm256_loadu_ps(iWB + j - 1) + _mm256_loadu_ps(iLB + j - 1) - wc - lc;
      const __m256 steepBR = _mm256_loadu_ps(iWB + j + 1) + _mm256_loadu_ps(iLB + j + 1) - wc - lc;

      //! Update the transfer
      __m256 wValue = applyMask256_ps(_mm256_cmp_ps(steepL , z0, _CMP_GT_OS), _mm256_loadu_ps(iWC + j - 1), wc) * steepL  +
                      applyMask256_ps(_mm256_cmp_ps(steepR , z0, _CMP_GT_OS), _mm256_loadu_ps(iWC + j + 1), wc) * steepR  +
                      applyMask256_ps(_mm256_cmp_ps(steepT , z0, _CMP_GT_OS), _mm256_loadu_ps(iWT + j    ), wc) * steepT  +
                      applyMask256_ps(_mm256_cmp_ps(steepB , z0, _CMP_GT_OS), _mm256_loadu_ps(iWB + j    ), wc) * steepB  +
           zWeight * (applyMask256_ps(_mm256_cmp_ps(steepTL, z0, _CMP_GT_OS), _mm256_loadu_ps(iWT + j - 1), wc) * steepTL +
                      applyMask256_ps(_mm256_cmp_ps(steepTR, z0, _CMP_GT_OS), _mm256_loadu_ps(iWT + j + 1), wc) * steepTR +
                      applyMask256_ps(_mm256_cmp_ps(steepBL, z0, _CMP_GT_OS), _mm256_loadu_ps(iWB + j - 1), wc) * steepBL +
                      applyMask256_ps(_mm256_cmp_ps(steepBR, z0, _CMP_GT_OS), _mm256_loadu_ps(iWB + j + 1), wc) * steepBR);

      __m256 cValue = applyMask256_ps(_mm256_cmp_ps(steepL , z0, _CMP_GT_OS), _mm256_loadu_ps(iCC + j - 1), cc) * steepL  +
                      applyMask256_ps(_mm256_cmp_ps(steepR , z0, _CMP_GT_OS), _mm256_loadu_ps(iCC + j + 1), cc) * steepR  +
                      applyMask256_ps(_mm256_cmp_ps(steepT , z0, _CMP_GT_OS), _mm256_loadu_ps(iCT + j    ), cc) * steepT  +
                      applyMask256_ps(_mm256_cmp_ps(steepB , z0, _CMP_GT_OS), _mm256_loadu_ps(iCB + j    ), cc) * steepB  +
           zWeight * (applyMask256_ps(_mm256_cmp_ps(steepTL, z0, _CMP_GT_OS), _mm256_loadu_ps(iCT + j - 1), cc) * steepTL +
                      applyMask256_ps(_mm256_cmp_ps(steepTR, z0, _CMP_GT_OS), _mm256_loadu_ps(iCT + j + 1), cc) * steepTR +
                      applyMask256_ps(_mm256_cmp_ps(steepBL, z0, _CMP_GT_OS), _mm256_loadu_ps(iCB + j - 1), cc) * steepBL +
                      applyMask256_ps(_mm256_cmp_ps(steepBR, z0, _CMP_GT_OS), _mm256_loadu_ps(iCB + j + 1), cc) * steepBR);

      //! Store the values
      _mm256_storeu_ps(oTw + j, zTimeStep * wValue);
      _mm256_storeu_ps(oTc + j, zTimeStep * cValue);
    }
#endif
#if defined(__SSE__)
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! For convenience: central value
      const __m128 wc = _mm_loadu_ps(iWC + j);
      const __m128 lc = _mm_loadu_ps(iLC + j);
      const __m128 cc = _mm_loadu_ps(iCC + j);

      //! Compute the steep
      const __m128 steepL  = _mm_loadu_ps(iWC + j - 1) + _mm_loadu_ps(iLC + j - 1) - wc - lc;
      const __m128 steepR  = _mm_loadu_ps(iWC + j + 1) + _mm_loadu_ps(iLC + j + 1) - wc - lc;
      const __m128 steepT  = _mm_loadu_ps(iWT + j    ) + _mm_loadu_ps(iLT + j    ) - wc - lc;
      const __m128 steepB  = _mm_loadu_ps(iWB + j    ) + _mm_loadu_ps(iLB + j    ) - wc - lc;
      const __m128 steepTL = _mm_loadu_ps(iWT + j - 1) + _mm_loadu_ps(iLT + j - 1) - wc - lc;
      const __m128 steepTR = _mm_loadu_ps(iWT + j + 1) + _mm_loadu_ps(iLT + j + 1) - wc - lc;
      const __m128 steepBL = _mm_loadu_ps(iWB + j - 1) + _mm_loadu_ps(iLB + j - 1) - wc - lc;
      const __m128 steepBR = _mm_loadu_ps(iWB + j + 1) + _mm_loadu_ps(iLB + j + 1) - wc - lc;

      //! Update the transfer
      __m128 wValue = applyMask_ps(_mm_cmpgt_ps(steepL , x0), _mm_loadu_ps(iWC + j - 1), wc) * steepL  +
                      applyMask_ps(_mm_cmpgt_ps(steepR , x0), _mm_loadu_ps(iWC + j + 1), wc) * steepR  +
                      applyMask_ps(_mm_cmpgt_ps(steepT , x0), _mm_loadu_ps(iWT + j    ), wc) * steepT  +
                      applyMask_ps(_mm_cmpgt_ps(steepB , x0), _mm_loadu_ps(iWB + j    ), wc) * steepB  +
           xWeight * (applyMask_ps(_mm_cmpgt_ps(steepTL, x0), _mm_loadu_ps(iWT + j - 1), wc) * steepTL +
                      applyMask_ps(_mm_cmpgt_ps(steepTR, x0), _mm_loadu_ps(iWT + j + 1), wc) * steepTR +
                      applyMask_ps(_mm_cmpgt_ps(steepBL, x0), _mm_loadu_ps(iWB + j - 1), wc) * steepBL +
                      applyMask_ps(_mm_cmpgt_ps(steepBR, x0), _mm_loadu_ps(iWB + j + 1), wc) * steepBR);

      __m128 cValue = applyMask_ps(_mm_cmpgt_ps(steepL , x0), _mm_loadu_ps(iCC + j - 1), cc) * steepL  +
                      applyMask_ps(_mm_cmpgt_ps(steepR , x0), _mm_loadu_ps(iCC + j + 1), cc) * steepR  +
                      applyMask_ps(_mm_cmpgt_ps(steepT , x0), _mm_loadu_ps(iCT + j    ), cc) * steepT  +
                      applyMask_ps(_mm_cmpgt_ps(steepB , x0), _mm_loadu_ps(iCB + j    ), cc) * steepB  +
           xWeight * (applyMask_ps(_mm_cmpgt_ps(steepTL, x0), _mm_loadu_ps(iCT + j - 1), cc) * steepTL +
                      applyMask_ps(_mm_cmpgt_ps(steepTR, x0), _mm_loadu_ps(iCT + j + 1), cc) * steepTR +
                      applyMask_ps(_mm_cmpgt_ps(steepBL, x0), _mm_loadu_ps(iCB + j - 1), cc) * steepBL +
                      applyMask_ps(_mm_cmpgt_ps(steepBR, x0), _mm_loadu_ps(iCB + j + 1), cc) * steepBR);

      //! Store the values
      _mm_storeu_ps(oTw + j, xTimeStep * wValue);
      _mm_storeu_ps(oTc + j, xTimeStep * cValue);
    }
#endif
    //! Normal version to finish the loop
    for (; j < m_width; j++) {

      //! For convenience: central value
      const float wc = iWC[j];
      const float lc = iLC[j];
      const float cc = iCC[j];

      //! Compute the steep
      const float steepL  = iWC[j - 1] + iLC[j - 1] - wc - lc;
      const float steepR  = iWC[j + 1] + iLC[j + 1] - wc - lc;
      const float steepT  = iWT[j    ] + iLT[j    ] - wc - lc;
      const float steepB  = iWB[j    ] + iLB[j    ] - wc - lc;
      const float steepTL = iWT[j - 1] + iLT[j - 1] - wc - lc;
      const float steepTR = iWT[j + 1] + iLT[j + 1] - wc - lc;
      const float steepBL = iWB[j - 1] + iLB[j - 1] - wc - lc;
      const float steepBR = iWB[j + 1] + iLB[j + 1] - wc - lc;

      //! Update the transfer
      float wValue = (steepL  > 0 ? iWC[j - 1] : wc) * steepL  +
                     (steepR  > 0 ? iWC[j + 1] : wc) * steepR  +
                     (steepT  > 0 ? iWT[j    ] : wc) * steepT  +
                     (steepB  > 0 ? iWB[j    ] : wc) * steepB  +
           weight * ((steepTL > 0 ? iWT[j - 1] : wc) * steepTL +
                     (steepTR > 0 ? iWT[j + 1] : wc) * steepTR +
                     (steepBL > 0 ? iWB[j - 1] : wc) * steepBL +
                     (steepBR > 0 ? iWB[j + 1] : wc) * steepBR);

      float cValue = (steepL  > 0 ? iCC[j - 1] : cc) * steepL  +
                     (steepR  > 0 ? iCC[j + 1] : cc) * steepR  +
                     (steepT  > 0 ? iCT[j    ] : cc) * steepT  +
                     (steepB  > 0 ? iCB[j    ] : cc) * steepB  +
           weight * ((steepTL > 0 ? iCT[j - 1] : cc) * steepTL +
                     (steepTR > 0 ? iCT[j + 1] : cc) * steepTR +
                     (steepBL > 0 ? iCB[j - 1] : cc) * steepBL +
                     (steepBR > 0 ? iCB[j + 1] : cc) * steepBR);

      //! Store the values
      oTw[j] = timeStep * wValue;
      oTc[j] = timeStep * cValue;
    }
  }

#ifdef _OPENMP
#pragma omp barrier
#endif

  //! Update the water
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iTw = m_imWaterTmp->getPtr(p_tid, 0, i);
    const float* iTc = m_imSedimentTmp->getPtr(p_tid, 0, i);
    const float* iL = i_imLandscape.getPtr(0, i);
    const float* iR = (float*)&rainMatrix[i*m_width];
    float* oW = io_imWaterLvl.getPtr(0, i);
    float* oC = io_imSediment.getPtr(0, i);
    
    size_t j = 0;

#if defined(__AVX__)
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! We don't want to have negative water
      _mm256_storeu_ps(oW + j, applyMask256_ps(_mm256_cmp_ps(_mm256_loadu_ps(iL + j), zOceanLvl, _CMP_LT_OS), zOceanLvl - _mm256_loadu_ps(iL + j),
        _mm256_max_ps(_mm256_loadu_ps(oW + j) + zalpha * _mm256_loadu_ps(iTw + j) + zRain * _mm256_loadu_ps(iR + j) * zTimeStep, z0)));

      //! We don't want to have negative sedimentation AND
      //! We don't want to have more sedimentation than landscape
      _mm256_storeu_ps(oC + j, _mm256_max_ps(_mm256_loadu_ps(oC + j) + zalpha * _mm256_loadu_ps(iTc + j), z0));
    }
#endif
#if defined(__SSE__)
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! We don't want to have negative water
      _mm_storeu_ps(oW + j, applyMask_ps(_mm_cmplt_ps(_mm_loadu_ps(iL + j), xOceanLvl), xOceanLvl - _mm_loadu_ps(iL + j),
        _mm_max_ps(_mm_loadu_ps(oW + j) + xalpha * _mm_loadu_ps(iTw + j) + xRain *_mm_loadu_ps(iR + j) * xTimeStep, x0)));

      //! We don't want to have negative sedimentation AND
      //! We don't want to have more sedimentation than landscape
      _mm_storeu_ps(oC + j, _mm_max_ps(_mm_loadu_ps(oC + j) + xalpha * _mm_loadu_ps(iTc + j), x0));
    }
#endif
    //! Normal version to finish the loop
    for (; j < m_width; j++) {

      //! We don't want to have negative water AND
      //! The ocean can't move
      if (iL[j] < oceanLvl) {
        oW[j] = oceanLvl - iL[j];
      }
      else {
        oW[j] = std::max(oW[j] + m_alpha * iTw[j] + timeStep * rain * iR[j], 0.f);
      }

      //! We don't want to have negative sedimentation
      oC[j] = std::max(oC[j] + m_alpha * iTc[j], 0.f);
    }
  }
}


//! Release memory.
void WaterAndSedimentationEvolution::releaseMemory() {

  //! Release memory
  if (m_imWaterTmp != NULL) {
    delete m_imWaterTmp;
    m_imWaterTmp = NULL;
  }

  if (m_imSedimentTmp != NULL) {
    delete m_imSedimentTmp;
    m_imSedimentTmp = NULL;
  }
}
