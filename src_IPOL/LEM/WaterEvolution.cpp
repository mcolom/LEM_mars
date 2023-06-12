/**
 * @file WaterEvolution.cpp
 *
 * @brief Class to handle the Water Evolution scheme.
 *
 * @author Jérôme Darbon <jerome@math.ucla.edu> (original version).
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <malloc.h>
#include <omp.h>
#include <math.h>
#include <xmmintrin.h>
#include <x86intrin.h>


//! Local includes
#include "WaterEvolution.h"
#include "../Utilities/Utilities.h"


using namespace std;


//! Default constructor
WaterEvolution::WaterEvolution() :

  //! Parameters
  m_dt(0),
  m_rain(0),
  m_rainMatrix(NULL),
  m_oceanLvl(0),

  //! Size
  m_width (0),
  m_height(0),

  //! Images
  m_imTmp(NULL) {

}


//! Copy constructor
WaterEvolution::WaterEvolution(
  const WaterEvolution& i_we) :

  //! Parameters
  m_dt(i_we.m_dt),
  m_rain(i_we.m_rain),
  m_rainMatrix(i_we.m_rainMatrix),
  m_oceanLvl(i_we.m_oceanLvl),

  //! Size
  m_width (i_we.m_width),
  m_height(i_we.m_height),

  //! Images
  m_imTmp(new ImageX<float>(*i_we.m_imTmp)) {

}


//! Default destructor
WaterEvolution::~WaterEvolution() {

  //! Release memory
  releaseMemory();
}

//! Operator overload.
WaterEvolution& WaterEvolution::operator=(
  const WaterEvolution& i_we) {

  if (&i_we == this) {
    return *this;
  }

  releaseMemory();

  new (this) WaterEvolution(i_we);
  return *this;
}


//! Initialize the images and the parameters
void WaterEvolution::init(
  const size_t p_width,
  const size_t p_height,
  const size_t p_border,
  Parameters const& p_params) {

  //! Release memory
  releaseMemory();

  //! Parameters
  m_dt       = p_params.timeStep();
  m_rain     = p_params.rain();
  m_rainMatrix     = p_params.rainMatrix();
  m_oceanLvl = p_params.oceanLevel();
  m_mars_a = p_params.mars_a();

  //! Sizes
  m_height = p_height;
  m_width  = p_width;

  //! Images
  m_imTmp  = new ImageX<float>(m_width, m_height, 1, p_border);
}


//! Perform a one step evolution of the water scheme.
void WaterEvolution::performOneStepEvolution(
  const Image& i_imLandscape,
  Image& io_imWaterLvl,
  const size_t p_tid) {

  //! For convenience
  const float timeStep = m_dt;
  const float rain     = m_rain;
  const float* rainMatrix     = (float*)m_rainMatrix;
  const size_t hBegin  = i_imLandscape.nHeight(p_tid);
  const size_t hEnd    = i_imLandscape.nHeight(p_tid + 1);
  const float oceanLvl = m_oceanLvl;
  const float mars_a     = m_mars_a;

#if defined(__AVX__)
  //! AVX version
  const __m256 zTimeStep = _mm256_set1_ps(timeStep);
  const __m256 zRain     = _mm256_set1_ps(rain);
  const __m256 zWeight   = _mm256_set1_ps(1.f / sqrt(2.));
  const __m256 z0        = _mm256_setzero_ps();
  const __m256 zOceanLvl = _mm256_set1_ps(oceanLvl);
  const __m256 zMars_a = _mm256_set1_ps(mars_a);
#endif
#if defined(__SSE__)
  //! SSE version
  const __m128 xTimeStep = _mm_set1_ps(timeStep);
  const __m128 xRain     = _mm_set1_ps(rain);
  const __m128 xWeight   = _mm_set1_ps(1.f / sqrt(2.));
  const __m128 x0        = _mm_setzero_ps();
  const __m128 xOceanLvl = _mm_set1_ps(oceanLvl);
  const __m128 xMars_a   = _mm_set1_ps(mars_a);
#endif
  const float weight   = 1.f / sqrt(2.f);


  //! Compute the transfer
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iWT = io_imWaterLvl.getPtr(0, i - 1);
    const float* iWC = io_imWaterLvl.getPtr(0, i);
    const float* iWB = io_imWaterLvl.getPtr(0, i + 1);
    const float* iLT = i_imLandscape.getPtr(0, i - 1);
    const float* iLC = i_imLandscape.getPtr(0, i);
    const float* iLB = i_imLandscape.getPtr(0, i + 1);
    float* oT = m_imTmp->getPtr(p_tid, 0, i);
    size_t j = 0;

#if defined(__AVX__)
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! Get the 8 neighbors and central value
      const __m256 wc = _mm256_loadu_ps(iWC + j);
      const __m256 lc = _mm256_loadu_ps(iLC + j);

      //! Compute the steep
      const __m256 steepL  = _mm256_loadu_ps(iWC + j - 1) - wc + _mm256_loadu_ps(iLC + j - 1) - lc;
      const __m256 steepR  = _mm256_loadu_ps(iWC + j + 1) - wc + _mm256_loadu_ps(iLC + j + 1) - lc;
      const __m256 steepT  = _mm256_loadu_ps(iWT + j    ) - wc + _mm256_loadu_ps(iLT + j    ) - lc;
      const __m256 steepB  = _mm256_loadu_ps(iWB + j    ) - wc + _mm256_loadu_ps(iLB + j    ) - lc;
      const __m256 steepTL = _mm256_loadu_ps(iWT + j - 1) - wc + _mm256_loadu_ps(iLT + j - 1) - lc;
      const __m256 steepTR = _mm256_loadu_ps(iWT + j + 1) - wc + _mm256_loadu_ps(iLT + j + 1) - lc;
      const __m256 steepBL = _mm256_loadu_ps(iWB + j - 1) - wc + _mm256_loadu_ps(iLB + j - 1) - lc;
      const __m256 steepBR = _mm256_loadu_ps(iWB + j + 1) - wc + _mm256_loadu_ps(iLB + j + 1) - lc;

      //! Update the transfer
      __m256 value = applyMask256_ps(_mm256_cmp_ps(steepL , z0, _CMP_GT_OS), _mm256_loadu_ps(iWC + j - 1), wc) * steepL  +
                     applyMask256_ps(_mm256_cmp_ps(steepR , z0, _CMP_GT_OS), _mm256_loadu_ps(iWC + j + 1), wc) * steepR  +
                     applyMask256_ps(_mm256_cmp_ps(steepT , z0, _CMP_GT_OS), _mm256_loadu_ps(iWT + j    ), wc) * steepT  +
                     applyMask256_ps(_mm256_cmp_ps(steepB , z0, _CMP_GT_OS), _mm256_loadu_ps(iWB + j    ), wc) * steepB  +
          zWeight * (applyMask256_ps(_mm256_cmp_ps(steepTL, z0, _CMP_GT_OS), _mm256_loadu_ps(iWT + j - 1), wc) * steepTL +
                     applyMask256_ps(_mm256_cmp_ps(steepTR, z0, _CMP_GT_OS), _mm256_loadu_ps(iWT + j + 1), wc) * steepTR +
                     applyMask256_ps(_mm256_cmp_ps(steepBL, z0, _CMP_GT_OS), _mm256_loadu_ps(iWB + j - 1), wc) * steepBL +
                     applyMask256_ps(_mm256_cmp_ps(steepBR, z0, _CMP_GT_OS), _mm256_loadu_ps(iWB + j + 1), wc) * steepBR);

      //! Store the value
      _mm256_storeu_ps(oT + j, value * zTimeStep);
    }
#endif
#if defined(__SSE__)
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! Get the 8 neighbors and central value
      const __m128 wc = _mm_loadu_ps(iWC + j);
      const __m128 lc = _mm_loadu_ps(iLC + j);

      //! Compute the steep
      const __m128 steepL  = _mm_loadu_ps(iWC + j - 1) - wc + _mm_loadu_ps(iLC + j - 1) - lc;
      const __m128 steepR  = _mm_loadu_ps(iWC + j + 1) - wc + _mm_loadu_ps(iLC + j + 1) - lc;
      const __m128 steepT  = _mm_loadu_ps(iWT + j    ) - wc + _mm_loadu_ps(iLT + j    ) - lc;
      const __m128 steepB  = _mm_loadu_ps(iWB + j    ) - wc + _mm_loadu_ps(iLB + j    ) - lc;
      const __m128 steepTL = _mm_loadu_ps(iWT + j - 1) - wc + _mm_loadu_ps(iLT + j - 1) - lc;
      const __m128 steepTR = _mm_loadu_ps(iWT + j + 1) - wc + _mm_loadu_ps(iLT + j + 1) - lc;
      const __m128 steepBL = _mm_loadu_ps(iWB + j - 1) - wc + _mm_loadu_ps(iLB + j - 1) - lc;
      const __m128 steepBR = _mm_loadu_ps(iWB + j + 1) - wc + _mm_loadu_ps(iLB + j + 1) - lc;

      //! Update the transfer
      __m128 value = applyMask_ps(_mm_cmpgt_ps(steepL , x0), _mm_loadu_ps(iWC + j - 1), wc) * steepL  +
                     applyMask_ps(_mm_cmpgt_ps(steepR , x0), _mm_loadu_ps(iWC + j + 1), wc) * steepR  +
                     applyMask_ps(_mm_cmpgt_ps(steepT , x0), _mm_loadu_ps(iWT + j    ), wc) * steepT  +
                     applyMask_ps(_mm_cmpgt_ps(steepB , x0), _mm_loadu_ps(iWB + j    ), wc) * steepB  +
          xWeight * (applyMask_ps(_mm_cmpgt_ps(steepTL, x0), _mm_loadu_ps(iWT + j - 1), wc) * steepTL +
                     applyMask_ps(_mm_cmpgt_ps(steepTR, x0), _mm_loadu_ps(iWT + j + 1), wc) * steepTR +
                     applyMask_ps(_mm_cmpgt_ps(steepBL, x0), _mm_loadu_ps(iWB + j - 1), wc) * steepBL +
                     applyMask_ps(_mm_cmpgt_ps(steepBR, x0), _mm_loadu_ps(iWB + j + 1), wc) * steepBR);

      //! Store the value
      _mm_storeu_ps(oT + j, value * xTimeStep);
    }
#endif
    //! Normal version because the width is not mandatory to be a multiple of 4
    for (; j < m_width; j++) {

      //! For convenience: central value
      const float wc = iWC[j];
      const float lc = iLC[j];

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
      // This is tau_rho in the pseudocode
      oT[j] = timeStep * ((steepL  > 0 ? iWC[j - 1] : wc) * steepL  +
                          (steepR  > 0 ? iWC[j + 1] : wc) * steepR  +
                          (steepT  > 0 ? iWT[j    ] : wc) * steepT  +
                          (steepB  > 0 ? iWB[j    ] : wc) * steepB  +
                weight * ((steepTL > 0 ? iWT[j - 1] : wc) * steepTL +
                          (steepTR > 0 ? iWT[j + 1] : wc) * steepTR +
                          (steepBL > 0 ? iWB[j - 1] : wc) * steepBL +
                          (steepBR > 0 ? iWB[j + 1] : wc) * steepBR));
    }
  }


#ifdef _OPENMP
#pragma omp barrier
#endif

  //! Update the water
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iT = m_imTmp->getPtr(p_tid, 0, i);
    const float* iL = i_imLandscape.getPtr(0, i);
    float* oW = io_imWaterLvl.getPtr(0, i);
    const float* iR = (float*)&rainMatrix[i*m_width];
    size_t j = 0;

// With AVX and SSE we store in [oW + j] the result of applying a mask according to the max operator.
#if defined(__AVX__)
    //! AVX version
    for (; j < m_width - 8; j+=8) {

      //! We don't want to have negative water
      // _mm256_storeu_ps(.);
      // _mm256_storeu_ps(oW + j, applyMask256_ps(.));
      // Argument of applyMask256_ps:
      //    _mm256_cmp_ps(_mm256_loadu_ps(iL + j), zOceanLvl), --> Comparison
      //    zOceanLvl - _mm256_loadu_ps(iL + j), --> when h < theta
      //    _mm256_max_ps(_mm256_loadu_ps(oW + j) + _mm256_loadu_ps(iT + j) + zRain * _mm256_loadu_ps(iR + j) * zTimeStep, z0) --> when h >= theta
      _mm256_storeu_ps(oW + j, applyMask256_ps(_mm256_cmp_ps(_mm256_loadu_ps(iL + j), zOceanLvl, _CMP_LT_OS), zOceanLvl - _mm256_loadu_ps(iL + j),
        _mm256_max_ps(_mm256_loadu_ps(oW + j) + zMars_a * _mm256_loadu_ps(iT + j) + zRain * _mm256_loadu_ps(iR + j) * zTimeStep, z0)));
    }
#endif
#if defined(__SSE__)
    //! SSE version
    for (; j < m_width - 4; j+=4) {
      //! We don't want to have negative water
      _mm_storeu_ps(oW + j, applyMask_ps(_mm_cmplt_ps(_mm_loadu_ps(iL + j), xOceanLvl), xOceanLvl - _mm_loadu_ps(iL + j),
        _mm_max_ps(_mm_loadu_ps(oW + j) + xMars_a * _mm_loadu_ps(iT + j) + xRain * _mm_loadu_ps(iR + j) * xTimeStep, x0)));
    }
#endif
    //! Normal version
    for (; j < m_width; j++) {

      //! We don't want to have negative water AND
      //! The ocean can't move
      if (iL[j] < oceanLvl) {
        oW[j] = oceanLvl - iL[j];
      }
      else {
        // oW: output water
        // iT: temp image. It's tau_rho in the pseudocode
        // iR: rain matrix
        oW[j] = std::max(oW[j] + m_mars_a * iT[j] + timeStep * rain * iR[j], 0.f);
      }
    }
  }
}


//! Release memory.
void WaterEvolution::releaseMemory() {

  //! Release memory
  if (m_imTmp != NULL) {
    delete m_imTmp;
    m_imTmp = NULL;
  }
}
