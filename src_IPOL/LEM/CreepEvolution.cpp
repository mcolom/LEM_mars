/**
 * @file CreepEvolution.cpp
 *
 * @brief Class to handle the Creep Evolution scheme.
 *
 * @author Jérôme Darbon <jerome@math.ucla.edu> (original version).
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
   @author Miguel Colom - http://mcolom.info
 **/


//! Global includes
#include <malloc.h>
#include <omp.h>
#include <math.h>
#include <xmmintrin.h>
#include <x86intrin.h>


//! Local includes
#include "CreepEvolution.h"


using namespace std;


//! Default constructor
CreepEvolution::CreepEvolution() :
  //! Parameters
  m_params(NULL),
  m_dt(0),

  //! Size
  m_width (0),
  m_height(0),

  //! Images
  m_imLandscapeVariation(NULL) {

}


//! Copy constructor.
CreepEvolution::CreepEvolution(
  const CreepEvolution& i_ce) :

  //! Parameters
  m_params(i_ce.m_params),
  m_dt(m_params->timeStep()),

  //! Miscellaneous
  m_width (i_ce.m_width),
  m_height(i_ce.m_height),

  //! Images
  m_imLandscapeVariation(new ImageX<float>(*i_ce.m_imLandscapeVariation)) {

}


//! Default destructor
CreepEvolution::~CreepEvolution() {
  releaseMemory();
  m_params = NULL;
}

//! Operator overload.
CreepEvolution& CreepEvolution::operator=(
  const CreepEvolution& i_ce) {

  if (&i_ce == this) {
    return *this;
  }

  releaseMemory();

  new (this) CreepEvolution(i_ce);
  return *this;
}


//! Initialize the images and the parameters
void CreepEvolution::init(
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
}


//! Perform a one step evolution of the water scheme.
void CreepEvolution::performOneStepEvolution(
  Image &io_imLandscape,
  const size_t p_tid) {

  //! For convenience
  const float timeStep = m_dt;
  const size_t hBegin  = io_imLandscape.nHeight(p_tid);
  const size_t hEnd    = io_imLandscape.nHeight(p_tid + 1);
  const float weight   = 1.f / sqrt(2);
  const float oceanLvl = m_params->oceanLevel();

  //! Check CFL conditions for the creep value: c * dt <= 1/4
  const float creep    = std::min(m_params->creep(), 1.f / (4.f * timeStep));

#ifdef __AVX__
  //! AVX version
  const __m256 zWeight   = _mm256_set1_ps(weight);
  const __m256 zCreep    = _mm256_set1_ps(creep);
  const __m256 zTimeStep = _mm256_set1_ps(timeStep);
  const __m256 z4        = _mm256_set1_ps(4.f);
  const __m256 zOceanLvl = _mm256_set1_ps(oceanLvl);
#endif // __AVX__
#ifdef __SSE__
  //! SSE version
  const __m128 xWeight   = _mm_set1_ps(weight);
  const __m128 xCreep    = _mm_set1_ps(creep);
  const __m128 xTimeStep = _mm_set1_ps(timeStep);
  const __m128 x4        = _mm_set1_ps(4.f);
  const __m128 xOceanLvl = _mm_set1_ps(oceanLvl);
#endif // __SSE__

  //! InternalVariationLandscape is an image used internally
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iLT = io_imLandscape.getPtr(0, i - 1);
    const float* iLC = io_imLandscape.getPtr(0, i    );
    const float* iLB = io_imLandscape.getPtr(0, i + 1);
    float* oV = m_imLandscapeVariation->getPtr(p_tid, 0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! Compute the creep with 8 neigbors
      __m256 value = _mm256_loadu_ps(iLC + j - 1) + _mm256_loadu_ps(iLC + j + 1) +
                     _mm256_loadu_ps(iLT + j    ) + _mm256_loadu_ps(iLB + j    ) -
                z4 * _mm256_loadu_ps(iLC + j    ) +

           zWeight * (_mm256_loadu_ps(iLT + j - 1) + _mm256_loadu_ps(iLT + j + 1) +
                      _mm256_loadu_ps(iLB + j - 1) + _mm256_loadu_ps(iLB + j + 1) -
                 z4 * _mm256_loadu_ps(iLC + j    ));

      //! Store the result
      _mm256_storeu_ps(oV + j, value * zCreep);
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! Compute the creep with 8 neigbors
      __m128 value = _mm_loadu_ps(iLC + j - 1) + _mm_loadu_ps(iLC + j + 1) +
                     _mm_loadu_ps(iLT + j    ) + _mm_loadu_ps(iLB + j    ) -
                x4 * _mm_loadu_ps(iLC + j    ) +

           xWeight * (_mm_loadu_ps(iLT + j - 1) + _mm_loadu_ps(iLT + j + 1) +
                      _mm_loadu_ps(iLB + j - 1) + _mm_loadu_ps(iLB + j + 1) -
                 x4 * _mm_loadu_ps(iLC + j    ));

      //! Store the result
      _mm_storeu_ps(oV + j, value * xCreep);
    }
#endif // __SSE__

    //! Normal version
    for (; j < m_width; j++) {
      oV[j] = creep * (iLC[j - 1] + iLC[j + 1] + iLT[j    ] + iLB[j    ] - 4 * iLC[j]
           + weight * (iLT[j - 1] + iLT[j + 1] + iLB[j - 1] + iLB[j + 1] - 4 * iLC[j]));
    }
  }

#ifdef _OPENMP
#pragma omp barrier
#endif

  //! Update the landscape
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iV = m_imLandscapeVariation->getPtr(p_tid, 0, i);
    float* oL = io_imLandscape.getPtr(0, i);
    size_t j = 0;

#ifdef __AVX__
    //! AVX version
    for (; j < m_width - 8; j += 8) {

      //! We don't want the landscape to be negative
      const __m256 land  = _mm256_loadu_ps(oL + j);
      __m256 delta = _mm256_max_ps(zTimeStep * _mm256_loadu_ps(iV + j), zOceanLvl - land);

      //! Update the landscape value
      _mm256_storeu_ps(oL + j, land + delta);
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < m_width - 4; j += 4) {

      //! We don't want the landscape to be negative
      const __m128 land = _mm_loadu_ps(oL + j);
      __m128 delta = _mm_max_ps(xTimeStep * _mm_loadu_ps(iV + j), xOceanLvl - land);

      //! Update the landscape value
      _mm_storeu_ps(oL + j, land + delta);
    }
#endif // __SSE__

    //! Normal version
    for (; j < m_width; j++) {

      //! We don't want the landscape to be negative
      const float delta = std::max(timeStep * iV[j], oceanLvl - oL[j]);

      //! Update the landscape value
      oL[j] += delta;
    }
  }
}


//! Release memory.
void CreepEvolution::releaseMemory() {

  //! Release memory
  if (m_imLandscapeVariation != NULL) {
    delete m_imLandscapeVariation;
    m_imLandscapeVariation = NULL;
  }
}









