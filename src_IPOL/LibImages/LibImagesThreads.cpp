/**
 * @brief Class to deal with images.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

//! Global includes
#include <iostream>
#include <omp.h>
#include <xmmintrin.h>
#include <x86intrin.h>
#include <math.h>
#include <algorithm>
#include <string.h>

//! Local includes
#include "LibImagesThreads.h"
#include "../Utilities/Utilities.h"
#include "../Utilities/Memory.h"


using namespace std;


//! Default constructor
template <class T>
ImageX<T>::ImageX() :

  //! Size
  m_width   (0),
  m_height  (0),
  m_channels(0),
  m_border  (0),

  //! Miscalleneous
#ifdef _OPENMP
  m_nbThreads(omp_get_max_threads()),
#else
  m_nbThreads(1),
#endif // _OPENMP

  //! Pointers
  m_ptrs((T**) memalloc(16, m_nbThreads * sizeof(T*))),
  m_heights((size_t*) memalloc(16, (m_nbThreads + 1) * sizeof(size_t))) {

  //! To handle the different sizes
  for (size_t n = 0; n < m_nbThreads; n++) {
    m_ptrs[n] = NULL;
  }
}


//! Surcharged constructor with a specified size.
template <class T>
ImageX<T>::ImageX(
  const size_t i_width,
  const size_t i_height,
  const size_t i_channels,
  const size_t i_border) :

  //! Size
  m_width   (i_width),
  m_height  (i_height),
  m_channels(i_channels),
  m_border  (i_border),

  //! Miscalleneous
#ifdef _OPENMP
  m_nbThreads(omp_get_max_threads()),
#else
  m_nbThreads(1),
#endif // _OPENMP

  //! Pointers
  m_ptrs((T**) memalloc(16, m_nbThreads * sizeof(T*))),
  m_heights((size_t*) memalloc(16, (m_nbThreads + 1) * sizeof(size_t))) {

  //! Initialize the beginning and end for each slices of the image
  initializeHeights(m_height, m_heights, m_nbThreads);

  //! The image is sliced into m_nbThreads different lines. Those lines are
  //! keep in memory in different location.
  for (size_t n = 0; n < m_nbThreads; n++) {
    const size_t height = m_heights[n + 1] - m_heights[n] + 2 * m_border;
    const size_t width  = m_width + 2 * m_border;

    //! Allocate the memory
    m_ptrs[n] = (T*) memalloc(16, (width * height * m_channels + 128) * sizeof(T));

    //! Fill it with zeros
    for (size_t k = 0; k < height * width * m_channels; k++) {
      m_ptrs[n][k] = T(0);
    }
  }
}


//! Copy constructor.
template <class T>
ImageX<T>::ImageX(
  const ImageX<T>& i_im) :

  //! Size
  m_width   (i_im.m_width),
  m_height  (i_im.m_height),
  m_channels(i_im.m_channels),
  m_border  (i_im.m_border),

  //! Miscalleneous
  m_nbThreads(i_im.m_nbThreads),

  //! Pointers
  m_ptrs((T**) memalloc(16, m_nbThreads * sizeof(T*))),
  m_heights((size_t*) memalloc(16, (m_nbThreads + 1) * sizeof(size_t))) {

  //! Initialize the beginning and end for each slices of the image
  initializeHeights(m_height, m_heights, m_nbThreads);

  //! The image is sliced into m_nbThreads different lines. Those lines are
  //! keep in memory in different location.
  for (size_t n = 0; n < m_nbThreads; n++) {
    const size_t height = m_heights[n + 1] - m_heights[n] + 2 * m_border;
    const size_t width  = m_width + 2 * m_border;

    //! Allocate the memory
    m_ptrs[n] = (T*) memalloc(16, (width * height * m_channels + 128) * sizeof(T));

    //! Fill it with zeros
    for (size_t k = 0; k < height * width * m_channels; k++) {
      m_ptrs[n][k] = T(0);
    }
  }

  //! Copy the data
#ifdef _OPENMP
#pragma omp parallel
{
  const size_t n = omp_get_thread_num();
#else
  const size_t n = 0;
#endif
  const T* iI = i_im.m_ptrs[n];
  T* oI = m_ptrs[n];
  const size_t height = m_heights[n + 1] - m_heights[n] + 2 * m_border;
  const size_t width  = m_width + 2 * m_border;

  //! Fill it with the same values
  for (size_t k = 0; k < height * width * m_channels; k++) {
    oI[k] = iI[k];
  }
#ifdef _OPENMP
}
#endif
}


//! Default destructor
template <class T>
ImageX<T>::~ImageX() {

  //! Release the memory
  for (size_t n = 0; n < m_nbThreads; n++) {
    memfree(m_ptrs[n]);
  }

  memfree(m_ptrs);
  memfree(m_heights);
}


//! Get the pointer to the channel c, row i, according to the current active
//! thread n.
template <class T>
T* ImageX<T>::getPtr(
  const int p_n,
  const int p_c,
  const int p_i) const {

  //! Get the current number of lines
  const int height = m_heights[p_n + 1] - m_heights[p_n] + 2 * m_border;
  const int width  = m_width + 2 * m_border;

  return m_ptrs[p_n] + (p_c * height + p_i + int(m_border)
                        - int(m_heights[p_n])) * width + int(m_border);
}


//! Instance that can be used by this class
template class ImageX<float>;
template class ImageX<double>;
template class ImageX<int>;
template class ImageX<long double>;














