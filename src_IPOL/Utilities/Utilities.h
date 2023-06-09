#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED


//! Global includes
#include <string.h>
#include <time.h>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <xmmintrin.h>
#include <x86intrin.h>
#include <vector>


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief Compute an amount of time elapsed.
 *
 * @param io_start : when the clock started; Will be
 *        reinitialized just after;
 * @param p_name : name of the section to print.
 *
 * @return none.
 **/
void getTime(
  struct timespec &io_start,
  const char* p_name);


/**
 * @brief Open a file with the given name in the given path, write a string then
 *        close the file. Will exit the program if the file can't be opened.
 *
 * @return none.
 **/
void writeValue(
  const std::string &i_fileName,
  const std::string &i_pathName,
  const std::string &p_value);


/**
 * @brief Surcharged function to write a float value.
 *
 * @return none.
 **/
void writeValue(
  const std::string &i_fileName,
  const std::string &i_pathName,
  const float p_value);


/**
 * @brief Write a curve.
 **/
void writeCurve(
  const std::string &i_fileName,
  const std::string &i_pathName,
  std::vector<float> const& i_curve);


/**
 * @brief Write the estimation of automatic selection of parameters.
 **/
void writeParameters(
  const std::string &i_fileName,
  const std::string &i_pathName,
  const std::vector<float> &i_valNorm,
  const std::vector<float> &i_valCreep,
  const std::vector<float> &i_valSedimentation,
  const std::vector<float> &i_valExpM,
  const std::vector<float> &i_valExpN);


/**
 * @brief Write the automatic selection of the best parameters with the range of
 *        the tested parameters.
 **/
void writeBestParameters(
  const std::string &i_inpPath,
  const std::string &i_outPath,
  const float p_bestNorm,
  const std::vector<float> p_bestValues,
  const std::vector<float> p_rangeMin,
  const std::vector<float> p_rangeMax,
  const std::vector<std::string> p_paramNames);


/**
 * @brief Get left and right values from the center one: 0.5 -> 0.4, 0.5, 0.6.
 **/
void getValues(
  const float p_value,
  std::vector<float>& o_values);


/**
 * @brief Structure to sort an array according to an image.
 **/
struct topologicalSort
{
  const Image & imRef;
  topologicalSort(const Image& i_im) : imRef(i_im) {}

  bool operator()(std::pair<int, int> s1, std::pair<int, int> s2) {
    const int i1 = s1.first;
    const int j1 = s1.second;
    const int i2 = s2.first;
    const int j2 = s2.second;
    const float val1 = imRef.getPtr(0, i1)[j1];
    const float val2 = imRef.getPtr(0, i2)[j2];

    if (val1 < val2) {
      return false;
    }

    if (val1 > val2) {
      return true;
    }

    if (i1 < i2) {
      return false;
    }
    else if ((i1 == i2) && (j1 < j2)) {
      return false;
    }

    return true;
  }
};


/**
 * @brief Return the optimal cut of the lines of an image according to
 *        the number of threads available.
 *
 * @param i_height: the height of the image;
 * @param o_heights: will contain the lines position (must be already allocated)
 * @param p_nbThreads: the number of threads available.
 **/
void initializeHeights(
  const size_t i_height,
  size_t* o_heights,
  const size_t p_nbThreads);


/**
 * @brief For convenience.
 **/
//! pow(value, 0.5)
template <class T>
T mySqrt(const T p_value, const T) {
  return sqrt(p_value);
}

//! pow(value, exp)
template <class T>
T myPow(const T p_value, const T p_exp) {
  return pow(p_value, p_exp);
}

//! pow(value, 0.25)
template <class T>
T myPow25(const T p_value, const T) {
  return sqrt(sqrt(p_value));
}

//! pow(value, 1)
template <class T>
T myId(const T p_value, const T) {
  return p_value;
}

//! pow(value, 2)
template <class T>
T mySqr2(const T p_value, const T) {
  return p_value * p_value;
}

//! pow(value, 1.5)
template <class T>
T myPow15(const T p_value, const T) {
  return sqrt(p_value) * p_value;
}


/**
 * @brief Choose which pow function to use according to expM/expN parameters.
 **/
template <class T>
T (*getPowPtr(const T p_value))(const T, const T) {
  if (p_value == 0.25) {
    return &myPow25<T>;
  }
  else if (p_value == 0.5) {
    return &mySqrt<T>;
  }
  else if (p_value == 1.0) {
    return &myId<T>;
  }
  else if (p_value == 1.5) {
    return &myPow15<T>;
  }
  else if (p_value == 2.0) {
    return &mySqr2<T>;
  }
  else {
    return &myPow<T>;
  }
}

#ifdef __SSE__
/**
 * @brief If else then SSE implementation.
 **/
inline __m128 applyMask_ps(
  const __m128 i_mask,
  const __m128 i_then,
  const __m128 i_else){

  // return mask ? then : else
  return _mm_or_ps(_mm_and_ps(i_mask, i_then), _mm_andnot_ps(i_mask, i_else));
}

//! pow(value, 0.5)
template <class T>
T mySqrt_128(const T p_value, const float) {
  return _mm_sqrt_ps(p_value);
}

//! pow(value, exp)
template <class T>
T myPow_128(const T p_value, const float p_exp) {
  float values[4];
  _mm_storeu_ps(values, p_value);
  values[0] = pow(values[0], p_exp);
  values[1] = pow(values[1], p_exp);
  values[2] = pow(values[2], p_exp);
  values[3] = pow(values[3], p_exp);
  return _mm_loadu_ps(values);
}

//! pow(value, 0.25)
template <class T>
T myPow25_128(const T p_value, const float) {
  return _mm_sqrt_ps(_mm_sqrt_ps(p_value));
}

//! pow(value, 1)
template <class T>
T myId_128(const T p_value, const float) {
  return p_value;
}

//! pow(value, 2)
template <class T>
T mySqr2_128(const T p_value, const float) {
  return p_value * p_value;
}

//! pow(value, 1.5)
template <class T>
T myPow15_128(const T p_value, const float) {
  return _mm_sqrt_ps(p_value) * p_value;
}

/**
 * @brief Choose which pow function to use according to expM/expN parameters.
 **/
template <class T>
T (*getPowPtr128(const float p_value))(const T, const float) {
  if (p_value == 0.25) {
    return &myPow25_128<T>;
  }
  else if (p_value == 0.5) {
    return &mySqrt_128<T>;
  }
  else if (p_value == 1.0) {
    return &myId_128<T>;
  }
  else if (p_value == 1.5) {
    return &myPow15_128<T>;
  }
  else if (p_value == 2.0) {
    return &mySqr2_128<T>;
  }
  else {
    return &myPow_128<T>;
  }
}
#endif // __SSE__

#ifdef __AVX__
/**
 * @brief If else then AVX implementation.
 **/
inline __m256 applyMask256_ps(
  const __m256 i_mask,
  const __m256 i_then,
  const __m256 i_else){

  // return mask ? then : else
  return _mm256_or_ps(_mm256_and_ps(i_mask, i_then), _mm256_andnot_ps(i_mask, i_else));
}

//! pow(value, 0.5)
template <class T>
__m256 mySqrt_256(const __m256 p_value, const float) {
  return _mm256_sqrt_ps(p_value);
}

//! pow(value, exp)
template <class T>
__m256 myPow_256(const __m256 p_value, const float p_exp) {
  float values[8];
  _mm256_storeu_ps(values, p_value);
  values[0] = pow(values[0], p_exp);
  values[1] = pow(values[1], p_exp);
  values[2] = pow(values[2], p_exp);
  values[3] = pow(values[3], p_exp);
  values[4] = pow(values[4], p_exp);
  values[5] = pow(values[5], p_exp);
  values[6] = pow(values[6], p_exp);
  values[7] = pow(values[7], p_exp);
  return _mm256_loadu_ps(values);
}

//! pow(value, 0.25)
template <class T>
T myPow25_256(const T p_value, const float) {
  return _mm256_sqrt_ps(_mm256_sqrt_ps(p_value));
}

//! pow(value, 1)
template <class T>
T myId_256(const T p_value, const float) {
  return p_value;
}

//! pow(value, 2)
template <class T>
T mySqr2_256(const T p_value, const float) {
  return p_value * p_value;
}

//! pow(value, 1.5)
template <class T>
T myPow15_256(const T p_value, const float) {
  return _mm256_sqrt_ps(p_value) * p_value;
}

/**
 * @brief Choose which pow function to use according to expM/expN parameters.
 **/
template <class T>
T (*getPowPtr256(const float p_value))(const T, const float) {
  if (p_value == 0.25) {
    return &myPow25_256<T>;
  }
  else if (p_value == 0.5) {
    return &mySqrt_256<T>;
  }
  else if (p_value == 1) {
    return &myId_256<T>;
  }
  else if (p_value == 1.5) {
    return &myPow15_256<T>;
  }
  else if (p_value == 2) {
    return &mySqr2_256<T>;
  }
  else {
    return &myPow_256<T>;
  }
}
#endif // __AVX__

/**
 * @brief Get the size of a png image without reading it entirely.
 **/
void getSize(
  const char* i_name,
  size_t &o_width,
  size_t &o_height);


/**
 * @brief Return the initial number of threads to use with OPENMP
 *
 * @param p_nbThreads: maximal number of threads wanted by the user;
 * @param p_name: name of the image (to have an idea of the size).
 *
 * @return the optimal number of threads to use.
 **/
int getOptimalThreadsNumber(
  const int p_nbThreads,
  const char* p_name);


/**
 * @brief Convert a float to a string with a given precision.
 **/
template <class T>
inline std::string float2str(
  const T p_val,
  const size_t p_precision) {
  std::stringstream stream;
  stream << std::fixed << std::setprecision(p_precision) << p_val;
  return stream.str();
}


/**
 * @brief Read nodes of a palette of color for topography visualization.
 **/
void readPalette(
  const char* i_fileName,
  std::vector<float> &o_nodes);


#endif // UTILITIES_H_INCLUDED
