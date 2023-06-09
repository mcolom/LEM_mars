/**
 * @file CreepEvolution.cpp
 *
 * @brief Class to handle the Creep Evolution scheme.
 *
 * @author Jérôme Darbon <jerome@math.ucla.edu> (original version).
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 * @author Miguel Colom - http://mcolom.info
 **/


//! Global includes
#include <malloc.h>
#include <omp.h>
#include <math.h>


//! Local includes
#include "LandscapeEvolution.h"


using namespace std;


//! Default constructor
LandscapeEvolution::LandscapeEvolution() :

  //! Parameters
  m_params(NULL),
  m_dt(0),

  //! Size
  m_width (0),
  m_height(0),

  //! Images
  m_imLandscapeVariation(NULL),

  //! Function
  m_powM(NULL),
  m_powN(NULL)
#ifdef __AVX__
  , m_powM256(NULL)
  , m_powN256(NULL)
#endif // __AVX__
#ifdef __SSE__
  , m_powM128(NULL)
  , m_powN128(NULL)
#endif // __SSE__
   {
}


//! Copy constructor.
LandscapeEvolution::LandscapeEvolution(
  const LandscapeEvolution& i_le) :

  //! Parameters
  m_params(i_le.m_params),
  m_dt(m_params->timeStep()),

  //! Miscellaneous
  m_width (i_le.m_width),
  m_height(i_le.m_height),

  //! Images
  m_imLandscapeVariation(new ImageX<float>(*i_le.m_imLandscapeVariation)),

  //! Function
  m_powM(i_le.m_powM),
  m_powN(i_le.m_powN)
#ifdef __AVX__
  , m_powM256(i_le.m_powM256)
  , m_powN256(i_le.m_powN256)
#endif // __AVX__
#ifdef __SSE__
  , m_powM128(i_le.m_powM128)
  , m_powN128(i_le.m_powN128)
#endif // __SSE__
  {

}


//! Default destructor
LandscapeEvolution::~LandscapeEvolution() {

  //! Release memory
  releaseMemory();
  m_params = NULL;
}

//! Operator overload.
LandscapeEvolution& LandscapeEvolution::operator=(
  const LandscapeEvolution& i_le) {

  if (&i_le == this) {
    return *this;
  }

  releaseMemory();

  new (this) LandscapeEvolution(i_le);
  return *this;
}


//! Initialize the images and the parameters
void LandscapeEvolution::init(
  const size_t p_width,
  const size_t p_height,
  const size_t p_border,
  Parameters const& p_params) {

  //! Release memory
  releaseMemory();

  //! Parameters
  m_params = &p_params;
  m_dt     = m_params->timeStep();

  //! Sizes
  m_height = p_height;
  m_width  = p_width;

  //! Images
  m_imLandscapeVariation = new ImageX<float>(m_width, m_height, 1, p_border);

  //! Function
  m_powM = getPowPtr<float>(m_params->expM());
  m_powN = getPowPtr<float>(m_params->expN() / 2.f);
#ifdef __AVX__
  m_powM256 = getPowPtr256<__m256>(m_params->expM());
  m_powN256 = getPowPtr256<__m256>(m_params->expN() / 2.f);
#endif // __AVX__
#ifdef __SSE__
  m_powM128 = getPowPtr128<__m128>(m_params->expM());
  m_powN128 = getPowPtr128<__m128>(m_params->expN() / 2.f);
#endif // __SSE__
}


//! Perform a one step evolution of the landscape scheme.
void LandscapeEvolution::performOneStepEvolution(
  const Image &i_imWaterLvl,
  Image &io_imLandscape,
  Image &io_imWaterTimesConcentration,
  const size_t p_tid) {

  //! For convenience
  const float timeStep    = m_dt;
  const float expM        = m_params->expM      ();
  const float expN        = m_params->expN      () / 2.f;
  const float erosionE    = m_params->e         ();
  const float erosionS    = m_params->s         ();
  const float oceanLvl    = m_params->oceanLevel();
  const size_t hBegin = i_imWaterLvl.nHeight(p_tid);
  const size_t hEnd   = i_imWaterLvl.nHeight(p_tid + 1);

#ifdef __AVX__
  //! AVX version
  const __m256 zTimeStep = _mm256_set1_ps(timeStep);
  const __m256 zErosionE = _mm256_set1_ps(erosionE);
  const __m256 zErosionS = _mm256_set1_ps(erosionS);
  const __m256 z0        = _mm256_setzero_ps();
  const __m256 zWeight   = _mm256_set1_ps(0.5);
  const __m256 zOceanLvl = _mm256_set1_ps(oceanLvl);
#endif // __AVX__
#ifdef __SSE__
  //! SSE version
  const __m128 xTimeStep = _mm_set1_ps(timeStep);
  const __m128 xErosionE = _mm_set1_ps(erosionE);
  const __m128 xErosionS = _mm_set1_ps(erosionS);
  const __m128 x0        = _mm_setzero_ps();
  const __m128 xWeight   = _mm_set1_ps(0.5);
  const __m128 xOceanLvl = _mm_set1_ps(oceanLvl);
#endif // __SSE__

  //! Step 8: compute erosion
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iLC = io_imLandscape.getPtr(0, i    );
    const float* iLT = io_imLandscape.getPtr(0, i - 1);
    const float* iLB = io_imLandscape.getPtr(0, i + 1);
    const float* iWC = i_imWaterLvl  .getPtr(0, i    );
    const float* iWT = i_imWaterLvl  .getPtr(0, i - 1);
    const float* iWB = i_imWaterLvl  .getPtr(0, i + 1);
    const float* iC  = io_imWaterTimesConcentration.getPtr(0, i);
    float* oV = m_imLandscapeVariation->getPtr(p_tid, 0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! Initialization
      const __m256 wc = _mm256_loadu_ps(iWC + j);
      const __m256 lc = _mm256_loadu_ps(iLC + j);
      __m256 amountToErode = _mm256_setzero_ps();
      const __m256 hc = lc + wc;

      //! Horizontal case
      {
        const __m256 hl = _mm256_loadu_ps(iLC + j - 1) + _mm256_loadu_ps(iWC + j - 1);
        const __m256 hr = _mm256_loadu_ps(iLC + j + 1) + _mm256_loadu_ps(iWC + j + 1);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - hl, hc - hr), z0);
        amountToErode += tmp * tmp;
      }

      //! Vertical case
      {
        const __m256 ht = _mm256_loadu_ps(iLT + j) + _mm256_loadu_ps(iWT + j);
        const __m256 hb = _mm256_loadu_ps(iLB + j) + _mm256_loadu_ps(iWB + j);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - ht, hc - hb), z0);
        amountToErode += tmp * tmp;
      }

      //! Diagonal case
      {
        const __m256 htl = _mm256_loadu_ps(iLT + j - 1) + _mm256_loadu_ps(iWT + j - 1);
        const __m256 hbr = _mm256_loadu_ps(iLB + j + 1) + _mm256_loadu_ps(iWB + j + 1);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - htl, hc - hbr), z0);
        amountToErode += tmp * tmp * zWeight;
      }
      {
        const __m256 htr = _mm256_loadu_ps(iLT + j + 1) + _mm256_loadu_ps(iWT + j + 1);
        const __m256 hbl = _mm256_loadu_ps(iLB + j - 1) + _mm256_loadu_ps(iWB + j - 1);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - htr, hc - hbl), z0);
        amountToErode += tmp * tmp * zWeight;
      }

      //! Get the final amount to erode
      const __m256 delta = expN == expM ?  m_powM256(wc * amountToErode, expM) :
        m_powM256(wc, expM) * m_powN256(amountToErode, expN);
      amountToErode = applyMask256_ps(_mm256_cmp_ps(wc, z0, _CMP_GT_OS),
        - zErosionE * delta, z0);

      //! Check that the erosion won't be submarine (but the sedimentation can)
      const __m256 mask = _mm256_cmp_ps(lc, zOceanLvl, _CMP_LE_OS);
      amountToErode = applyMask256_ps(mask, z0, amountToErode);

      //! Step 9: computation of sedimentation^{n+1}
      const __m256 sValue = applyMask256_ps(_mm256_cmp_ps(wc, z0, _CMP_GT_OS),
        zErosionS * _mm256_loadu_ps(iC + j) / wc, z0);

      //! Step 10. update the landscape variation
      //! Note that the creep is performed elsewhere/indepently of the
      //! erosion/sedimentation
      _mm256_storeu_ps(oV + j, zTimeStep * (amountToErode + sValue));
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! Initialization
      const __m128 wc = _mm_loadu_ps(iWC + j);
      const __m128 lc = _mm_loadu_ps(iLC + j);
      __m128 amountToErode = _mm_setzero_ps();
      const __m128 hc = lc + wc;

      //! Horizontal case
      {
        const __m128 hl = _mm_loadu_ps(iLC + j - 1) + _mm_loadu_ps(iWC + j - 1);
        const __m128 hr = _mm_loadu_ps(iLC + j + 1) + _mm_loadu_ps(iWC + j + 1);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - hl, hc - hr), x0);
        amountToErode += tmp * tmp;
      }

      //! Vertical case
      {
        const __m128 ht = _mm_loadu_ps(iLT + j) + _mm_loadu_ps(iWT + j);
        const __m128 hb = _mm_loadu_ps(iLB + j) + _mm_loadu_ps(iWB + j);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - ht, hc - hb), x0);
        amountToErode += tmp * tmp;
      }

      //! Diagonal case
      {
        const __m128 htl = _mm_loadu_ps(iLT + j - 1) + _mm_loadu_ps(iWT + j - 1);
        const __m128 hbr = _mm_loadu_ps(iLB + j + 1) + _mm_loadu_ps(iWB + j + 1);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - htl, hc - hbr), x0);
        amountToErode += tmp * tmp * xWeight;
      }
      {
        const __m128 htr = _mm_loadu_ps(iLT + j + 1) + _mm_loadu_ps(iWT + j + 1);
        const __m128 hbl = _mm_loadu_ps(iLB + j - 1) + _mm_loadu_ps(iWB + j - 1);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - htr, hc - hbl), x0);
        amountToErode += tmp * tmp * xWeight;
      }

      //! Get the final amount to erode
      const __m128 delta = expN == expM ? m_powM128(wc * amountToErode, expM) :
        m_powM128(wc, expM) * m_powN128(amountToErode, expN);
      amountToErode = applyMask_ps(_mm_cmpgt_ps(wc, x0),
        -xErosionE * delta, x0);

      //! Check that the erosion won't be submarine (but the sedimentation can)
      const __m128 mask = _mm_cmplt_ps(lc, xOceanLvl);
      amountToErode = applyMask_ps(mask, x0, amountToErode);

      //! Step 9: computation of sedimentation^{n+1}
      const __m128 sValue = applyMask_ps(_mm_cmpgt_ps(wc, x0),
        xErosionS * _mm_loadu_ps(iC + j) / wc, x0);

      //! Step 10. update the landscape variation
      //! Note that the creep is performed elsewhere/indepently of the
      //! erosion/sedimentation
      _mm_storeu_ps(oV + j, xTimeStep * (amountToErode + sValue));
    }
#endif // __SSE__

    //! Normal version to finish the loop
    for (; j < m_width; j++) {

      //! Initialization
      float amountToErode = 0.f, sValue = 0.f;
      const float h_ij = iLC[j] + iWC[j];

      //! Horizontal case
      {
        const float h_m1 = iLC[j - 1] + iWC[j - 1];
        const float h_p1 = iLC[j + 1] + iWC[j + 1];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp;
      }

      //! Vertical case
      {
        const float h_m1 = iLT[j] + iWT[j];
        const float h_p1 = iLB[j] + iWB[j];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp;
      }

      //! Diagonales
      {
        const float h_m1 = iLT[j - 1] + iWT[j - 1];
        const float h_p1 = iLB[j + 1] + iWB[j + 1];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp * 0.5;
      }
      {
        const float h_m1 = iLT[j + 1] + iWT[j + 1];
        const float h_p1 = iLB[j - 1] + iWB[j - 1];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp * 0.5;
      }

      //! Get the final amount to erode
      if (iWC[j] > 0) {
        amountToErode = -erosionE * m_powM(iWC[j]       , expM)
                                  * m_powN(amountToErode, expN);
        sValue        = erosionS * iC[j] / iWC[j];
      }
      else {
        amountToErode = 0.f;
        sValue        = 0.f;
      }

      //! Check that the erosion won't be submarine (but the sedimentation can)
      if (iLC[j] < oceanLvl) {
        amountToErode = 0.f;
      }

      //! Step 10. update the landscape variation
      //! Note that the creep is performed elsewhere/indepently of the
      //! erosion/sedimentation
      oV[j] = timeStep * (amountToErode + sValue);
    }
  }

#ifdef _OPENMP
#pragma omp barrier
#endif

  //! Update the landscape
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iV = m_imLandscapeVariation->getPtr(p_tid, 0, i);
    float* oL = io_imLandscape.getPtr(0, i);
    float* oC = io_imWaterTimesConcentration.getPtr(0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const __m256 land = _mm256_loadu_ps(oL + j);
      __m256 val = _mm256_max_ps(_mm256_min_ps(_mm256_loadu_ps(iV + j), _mm256_loadu_ps(oC + j)), -land);

      //! Update the landscape and sediment
      _mm256_storeu_ps(oL + j, _mm256_loadu_ps(oL + j) + val);
      _mm256_storeu_ps(oC + j, _mm256_loadu_ps(oC + j) - val);
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const __m128 land = _mm_loadu_ps(oL + j);
      __m128 val = _mm_max_ps(_mm_min_ps(_mm_loadu_ps(iV + j), _mm_loadu_ps(oC + j)), -land);

      //! Update the landscape and sediment
      _mm_storeu_ps(oL + j, _mm_loadu_ps(oL + j) + val);
      _mm_storeu_ps(oC + j, _mm_loadu_ps(oC + j) - val);
    }
#endif // __SSE__

    //! Normal version to finish the loop
    for (; j < m_width; j++) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const float value = std::max(std::min(iV[j], oC[j]), -oL[j]);

      //! Update the landscape and sediment
      oL[j] += value;
      oC[j] -= value;
    }
  }
}


//! Compute the transfert and update the new mean of the landscape.
void LandscapeEvolution::computeTransfert(
  const Image& i_imWaterLvl,
  const Image &i_imLandscape,
  const Image &i_imWaterTimesConcentration,
  double &o_mean,
  const size_t p_tid) {

  //! For convenience
  const float timeStep    = m_dt;
  const float expM        = m_params->expM      ();
  const float expN        = m_params->expN      () / 2.f;
  const float erosionE    = m_params->e         ();
  const float erosionS    = m_params->s         ();
  const float oceanLvl    = m_params->oceanLevel();
  const size_t hBegin = i_imWaterLvl.nHeight(p_tid);
  const size_t hEnd   = i_imWaterLvl.nHeight(p_tid + 1);

#ifdef __AVX__
  //! AVX version
  const __m256 zTimeStep = _mm256_set1_ps(timeStep);
  const __m256 zErosionE = _mm256_set1_ps(erosionE);
  const __m256 zErosionS = _mm256_set1_ps(erosionS);
  const __m256 z0        = _mm256_setzero_ps();
  const __m256 zWeight   = _mm256_set1_ps(0.5);
  const __m256 zOceanLvl = _mm256_set1_ps(oceanLvl);
#endif // __AVX__
#ifdef __SSE__
  //! SSE version
  const __m128 xTimeStep = _mm_set1_ps(timeStep);
  const __m128 xErosionE = _mm_set1_ps(erosionE);
  const __m128 xErosionS = _mm_set1_ps(erosionS);
  const __m128 x0        = _mm_setzero_ps();
  const __m128 xWeight   = _mm_set1_ps(0.5);
  const __m128 xOceanLvl = _mm_set1_ps(oceanLvl);
#endif // __SSE__

  //! Step 8: compute erosion
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iLC = i_imLandscape.getPtr(0, i    );
    const float* iLT = i_imLandscape.getPtr(0, i - 1);
    const float* iLB = i_imLandscape.getPtr(0, i + 1);
    const float* iWC = i_imWaterLvl  .getPtr(0, i    );
    const float* iWT = i_imWaterLvl  .getPtr(0, i - 1);
    const float* iWB = i_imWaterLvl  .getPtr(0, i + 1);
    const float* iC  = i_imWaterTimesConcentration.getPtr(0, i);
    float* oV = m_imLandscapeVariation->getPtr(p_tid, 0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! Initialization
      const __m256 wc = _mm256_loadu_ps(iWC + j);
      const __m256 lc = _mm256_loadu_ps(iLC + j);
      __m256 amountToErode = _mm256_setzero_ps();
      const __m256 hc = lc + wc;

      //! Horizontal case
      {
        const __m256 hl = _mm256_loadu_ps(iLC + j - 1) + _mm256_loadu_ps(iWC + j - 1);
        const __m256 hr = _mm256_loadu_ps(iLC + j + 1) + _mm256_loadu_ps(iWC + j + 1);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - hl, hc - hr), z0);
        amountToErode += tmp * tmp;
      }

      //! Vertical case
      {
        const __m256 ht = _mm256_loadu_ps(iLT + j) + _mm256_loadu_ps(iWT + j);
        const __m256 hb = _mm256_loadu_ps(iLB + j) + _mm256_loadu_ps(iWB + j);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - ht, hc - hb), z0);
        amountToErode += tmp * tmp;
      }

      //! Diagonal case
      {
        const __m256 htl = _mm256_loadu_ps(iLT + j - 1) + _mm256_loadu_ps(iWT + j - 1);
        const __m256 hbr = _mm256_loadu_ps(iLB + j + 1) + _mm256_loadu_ps(iWB + j + 1);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - htl, hc - hbr), z0);
        amountToErode += tmp * tmp * zWeight;
      }
      {
        const __m256 htr = _mm256_loadu_ps(iLT + j + 1) + _mm256_loadu_ps(iWT + j + 1);
        const __m256 hbl = _mm256_loadu_ps(iLB + j - 1) + _mm256_loadu_ps(iWB + j - 1);
        const __m256 tmp = _mm256_max_ps(_mm256_max_ps(hc - htr, hc - hbl), z0);
        amountToErode += tmp * tmp * zWeight;
      }

      //! Get the final amount to erode
      const __m256 delta = expN == expM ?  m_powM256(wc * amountToErode, expM) :
        m_powM256(wc, expM) * m_powN256(amountToErode, expN);
      amountToErode = applyMask256_ps(_mm256_cmp_ps(wc, z0, _CMP_GT_OS),
        - zErosionE * delta, z0);

      //! Check that the erosion won't be submarine (but the sedimentation can)
      const __m256 mask = _mm256_cmp_ps(lc, zOceanLvl, _CMP_LE_OS);
      amountToErode = applyMask256_ps(mask, z0, amountToErode);

      //! Step 9: computation of sedimentation^{n+1}
      const __m256 sValue = applyMask256_ps(_mm256_cmp_ps(wc, z0, _CMP_GT_OS),
        zErosionS * _mm256_loadu_ps(iC + j) / wc, z0);

      //! Step 10. update the landscape variation
      //! Note that the creep is performed elsewhere/indepently of the
      //! erosion/sedimentation
      _mm256_storeu_ps(oV + j, zTimeStep * (amountToErode + sValue));
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! Initialization
      const __m128 wc = _mm_loadu_ps(iWC + j);
      const __m128 lc = _mm_loadu_ps(iLC + j);
      __m128 amountToErode = _mm_setzero_ps();
      const __m128 hc = lc + wc;

      //! Horizontal case
      {
        const __m128 hl = _mm_loadu_ps(iLC + j - 1) + _mm_loadu_ps(iWC + j - 1);
        const __m128 hr = _mm_loadu_ps(iLC + j + 1) + _mm_loadu_ps(iWC + j + 1);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - hl, hc - hr), x0);
        amountToErode += tmp * tmp;
      }

      //! Vertical case
      {
        const __m128 ht = _mm_loadu_ps(iLT + j) + _mm_loadu_ps(iWT + j);
        const __m128 hb = _mm_loadu_ps(iLB + j) + _mm_loadu_ps(iWB + j);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - ht, hc - hb), x0);
        amountToErode += tmp * tmp;
      }

      //! Diagonal case
      {
        const __m128 htl = _mm_loadu_ps(iLT + j - 1) + _mm_loadu_ps(iWT + j - 1);
        const __m128 hbr = _mm_loadu_ps(iLB + j + 1) + _mm_loadu_ps(iWB + j + 1);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - htl, hc - hbr), x0);
        amountToErode += tmp * tmp * xWeight;
      }
      {
        const __m128 htr = _mm_loadu_ps(iLT + j + 1) + _mm_loadu_ps(iWT + j + 1);
        const __m128 hbl = _mm_loadu_ps(iLB + j - 1) + _mm_loadu_ps(iWB + j - 1);
        const __m128 tmp = _mm_max_ps(_mm_max_ps(hc - htr, hc - hbl), x0);
        amountToErode += tmp * tmp * xWeight;
      }

      //! Get the final amount to erode
      const __m128 delta = expN == expM ? m_powM128(wc * amountToErode, expM) :
        m_powM128(wc, expM) * m_powN128(amountToErode, expN);
      amountToErode = applyMask_ps(_mm_cmpgt_ps(wc, x0),
        -xErosionE * delta, x0);

      //! Check that the erosion won't be submarine (but the sedimentation can)
      const __m128 mask = _mm_cmplt_ps(lc, xOceanLvl);
      amountToErode = applyMask_ps(mask, x0, amountToErode);

      //! Step 9: computation of sedimentation^{n+1}
      const __m128 sValue = applyMask_ps(_mm_cmpgt_ps(wc, x0),
        xErosionS * _mm_loadu_ps(iC + j) / wc, x0);

      //! Step 10. update the landscape variation
      //! Note that the creep is performed elsewhere/indepently of the
      //! erosion/sedimentation
      _mm_storeu_ps(oV + j, xTimeStep * (amountToErode + sValue));
    }
#endif // __SSE__

    //! Normal version to finish the loop
    for (; j < m_width; j++) {

      //! Initialization
      float amountToErode = 0.f, sValue = 0.f;
      const float h_ij = iLC[j] + iWC[j];

      //! Horizontal case
      {
        const float h_m1 = iLC[j - 1] + iWC[j - 1];
        const float h_p1 = iLC[j + 1] + iWC[j + 1];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp;
      }

      //! Vertical case
      {
        const float h_m1 = iLT[j] + iWT[j];
        const float h_p1 = iLB[j] + iWB[j];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp;
      }

      //! Diagonales
      {
        const float h_m1 = iLT[j - 1] + iWT[j - 1];
        const float h_p1 = iLB[j + 1] + iWB[j + 1];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp * 0.5;
      }
      {
        const float h_m1 = iLT[j + 1] + iWT[j + 1];
        const float h_p1 = iLB[j - 1] + iWB[j - 1];
        const float tmp  = std::max(std::max(h_ij - h_m1, h_ij - h_p1), 0.f);
        amountToErode += tmp * tmp * 0.5;
      }

      //! Get the final amount to erode
      if (iWC[j] > 0) {
        amountToErode = -erosionE * m_powM(iWC[j]       , expM)
                                  * m_powN(amountToErode, expN);
        sValue        = erosionS * iC[j] / iWC[j];
      }
      else {
        amountToErode = 0.f;
        sValue        = 0.f;
      }

      //! Check that the erosion won't be submarine (but the sedimentation can)
      if (iLC[j] < oceanLvl) {
        amountToErode = 0.f;
      }

      //! Step 10. update the landscape variation
      //! Note that the creep is performed elsewhere/indepently of the
      //! erosion/sedimentation
      oV[j] = timeStep * (amountToErode + sValue);
    }
  }

#ifdef _OPENMP
#pragma omp barrier
#endif

  //! Update the mean according to the variation
  o_mean = 0.0;
#ifdef __AVX__
  __m256 zMean = _mm256_setzero_ps();
#endif // __AVX__
#ifdef __SSE__
  __m128 xMean = _mm_setzero_ps();
#endif // __SSE__

  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iV = m_imLandscapeVariation->getPtr(p_tid, 0, i);
    const float* iL = i_imLandscape.getPtr(0, i);
    const float* iC = i_imWaterTimesConcentration.getPtr(0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const __m256 land = _mm256_loadu_ps(iL + j);
      __m256 val = _mm256_min_ps(_mm256_loadu_ps(iV + j), _mm256_loadu_ps(iC + j));

      //! Update the mean
      zMean += land + val;
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const __m128 land = _mm_loadu_ps(iL + j);
      __m128 val = _mm_min_ps(_mm_loadu_ps(iV + j), _mm_loadu_ps(iC + j));

      //! Update the mean
      xMean += land + val;
    }
#endif // __SSE__

    //! Normal version to finish the loop
    for (; j < m_width; j++) {

      //! we don't want to remove more sedimentation that already is
      const float value = std::min(iV[j], iC[j]);
      o_mean += iL[j] + value;
    }
  }

  //! Update the mean value
#ifdef __AVX__
  float zmeans[8];
  _mm256_storeu_ps(zmeans, zMean);
  o_mean += double(zmeans[0]) + double(zmeans[1]) + double(zmeans[2])
          + double(zmeans[3]) + double(zmeans[4]) + double(zmeans[5])
          + double(zmeans[6]) + double(zmeans[7]);
#endif // __AVX__
#ifdef __SSE__
  float xmeans[4];
  _mm_storeu_ps(xmeans, xMean);
  o_mean += double(xmeans[0]) + double(xmeans[1]) + double(xmeans[2])
          + double(xmeans[3]);
#endif // __SSE__
}


//! Update both landscape height and water times concentration and apply an
//! uplift.
void LandscapeEvolution::applyTransfert(
  Image& io_imLandscape,
  Image& io_imWaterTimesConcentration,
  const float p_uplift,
  const Image& i_imUplift,
  const size_t p_tid) const {

  //! For convenience
  const size_t hBegin = io_imLandscape.nHeight(p_tid);
  const size_t hEnd   = io_imLandscape.nHeight(p_tid + 1);

#ifdef __AVX__
  //! AVX version
  const __m256 zUplift = _mm256_set1_ps(p_uplift);
#endif // __AVX__
#ifdef __SSE__
  //! SSE version
  const __m128 xUplift = _mm_set1_ps(p_uplift);
#endif // __SSE__

  //! Update the landscape
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iV = m_imLandscapeVariation->getPtr(p_tid, 0, i);
    const float* iU = i_imUplift.getPtr(0, i);
    float* oL = io_imLandscape.getPtr(0, i);
    float* oC = io_imWaterTimesConcentration.getPtr(0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const __m256 land = _mm256_loadu_ps(oL + j);
      __m256 val = _mm256_max_ps(_mm256_min_ps(_mm256_loadu_ps(iV + j), _mm256_loadu_ps(oC + j)), -(land + zUplift));

      //! Update the landscape and sediment
      _mm256_storeu_ps(oL + j, _mm256_loadu_ps(oL + j) + val + _mm256_loadu_ps(iU + j));
      _mm256_storeu_ps(oC + j, _mm256_loadu_ps(oC + j) - val);
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const __m128 land = _mm_loadu_ps(oL + j);
      __m128 val = _mm_max_ps(_mm_min_ps(_mm_loadu_ps(iV + j), _mm_loadu_ps(oC + j)), -(land + xUplift));

      //! Update the landscape and sediment
      _mm_storeu_ps(oL + j, _mm_loadu_ps(oL + j) + val + _mm_loadu_ps(iU + j));
      _mm_storeu_ps(oC + j, _mm_loadu_ps(oC + j) - val);
    }
#endif // __SSE__

    //! Normal version to finish the loop
    for (; j < m_width; j++) {

      //! We don't want the landscape to be negative AND
      //! we don't want to remove more sedimentation that already is
      const float value = std::max(std::min(iV[j], oC[j]), -(oL[j] + p_uplift));

      //! Update the landscape and sediment
      oL[j] += value + iU[j];
      oC[j] -= value;
    }
  }
}


//! Release memory.
void LandscapeEvolution::releaseMemory() {

  //! Release memory
  if (m_imLandscapeVariation != NULL) {
    delete m_imLandscapeVariation;
    m_imLandscapeVariation = NULL;
  }
}









