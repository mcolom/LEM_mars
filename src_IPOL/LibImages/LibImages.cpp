/**
 * @brief Class to deal with images, and do conversion from/to cv::Mat.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 * @author Miguel Colom - http://mcolom.info
 **/

//! Global includes
#include <malloc.h>
#include <complex.h> // Must be included BEFORE fftw
#include <fftw3.h>
#include <tiffio.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <xmmintrin.h>
#include <x86intrin.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <vector>
#include <sys/types.h>
#include <unistd.h>

//! Local includes
#include "LibImages.h"
#include "mt19937ar.h"
#include "../Utilities/Utilities.h"
#include "../Utilities/Memory.h"


using namespace std;


//! Default constructor
Image::Image() :

  //! Size
  m_width   (0),
  m_height  (0),
  m_channels(0),
  m_border  (0),
  m_size    (0),

  //! Miscalleneous
#ifdef _OPENMP
  m_nbThreads(omp_get_max_threads()),
#else
  m_nbThreads(1),
#endif // _OPENMP

  //! Pointers
  m_ptr     (NULL),
  m_heights (NULL) {

}


//! Surcharged constructor with a specified size.
Image::Image(
  const size_t i_width,
  const size_t i_height,
  const size_t i_channels,
  const size_t i_border) :

  //! Size
  m_width   (i_width),
  m_height  (i_height),
  m_channels(i_channels),
  m_border  (i_border),
  m_size    (i_channels * (i_width + 2 * i_border) * (i_height + 2 * i_border)),

  //! Miscalleneous
#ifdef _OPENMP
  m_nbThreads(omp_get_max_threads()),
#else
  m_nbThreads(1),
#endif // _OPENM

  //! Pointers
  m_ptr     ((float*) memalloc(16, m_size * sizeof(float))),
  m_heights ((size_t*) memalloc(16, (m_nbThreads + 1) * sizeof(size_t))) {

  //! Fill the image with zeros
  size_t k = 0;

#if defined(__AVX__)
  //! AVX version
  for (; k < m_size - 8; k += 8) {
    _mm256_storeu_ps(m_ptr + k, _mm256_setzero_ps());
  }
#endif
#if defined(__SSE__)
  //! SSE version
  for (; k < m_size - 4; k += 4) {
    _mm_storeu_ps(m_ptr + k, _mm_setzero_ps());
  }
#endif
  //! Normal version
  for (; k < m_size; k++) {
    m_ptr[k] = float(0);
  }

  //! Initialize the beginning and end for each slices of the image
  initializeHeights(m_height, m_heights, m_nbThreads);
}


//! Copy constructor.
Image::Image(
  const Image& i_im) :

  //! Size
  m_width   (i_im.m_width),
  m_height  (i_im.m_height),
  m_channels(i_im.m_channels),
  m_border  (i_im.m_border),
  m_size    (i_im.m_size),

  //! Miscalleneous
  m_nbThreads(i_im.m_nbThreads),

  //! Pointers
  m_ptr     ((float*) memalloc(16, m_size * sizeof(float))),
  m_heights ((size_t*) memalloc(16, (m_nbThreads + 1) * sizeof(size_t))) {

  //! Copy the data
  size_t k = 0;
#if defined(__AVX__)
  //! AVX version
  for (; k < m_size - 8; k += 8) {
    _mm256_storeu_ps(m_ptr + k, _mm256_loadu_ps(i_im.m_ptr + k));
  }
#endif
#if defined(__SSE__)
  //! SSE version
  for (; k < m_size - 4; k += 4) {
    _mm_storeu_ps(m_ptr + k, _mm_loadu_ps(i_im.m_ptr + k));
  }
#endif
  //! Normal version
  for (; k < m_size; k++) {
    m_ptr[k] = i_im.m_ptr[k];
  }

  //! Copy the slices
  for (size_t k = 0; k < m_nbThreads + 1; k++) {
    m_heights[k] = i_im.m_heights[k];
  }
}


//! Default destructor
Image::~Image() {

  //! Release the memory
  releaseMemory();
}


//! Operator overload.
Image& Image::operator=(
  const Image& i_im) {

  if (&i_im == this) {
    return *this;
  }

  releaseMemory();

  new (this) Image(i_im);
  return *this;
}


//! Generic read of an image.
void Image::read(
  const char* p_name,
  const size_t i_border) {

  //! Parameters check
  if (NULL == p_name) {
    cout << "Null name" << endl;
    exit(EXIT_FAILURE);
  }
  //! Check the extension
  string ext = p_name;
  const int pos = ext.find_last_of(".");
  ext = ext.substr(pos + 1, ext.size());

  //! Call the right function
  if (ext == "png") {
    this->readPng(p_name, i_border);
  }
  else if (ext == "tif" || ext == "tiff") {
    this->readTiff(p_name, i_border);
  }
  else {
    cout << "Extension " << ext << " not known. Abort." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Read png-only image.
void Image::readPng(
  const char* p_name,
  const size_t i_border) {

  //! Initialization
  const size_t png_sig_len = 4;
  int transform = PNG_TRANSFORM_IDENTITY;

  //! Open the PNG input file
  FILE *volatile fp = NULL;
  if (0 == strcmp(p_name, "-")) {
    fp = stdin;
  }
  else if (NULL == (fp = fopen(p_name, "rb"))) {
    cout << "Can't open the file " << p_name << "." << endl;
    exit(EXIT_FAILURE);
  }

  //! Read in some of the signature bytes and check this signature
  png_byte png_sig[4];
  if ((png_sig_len != fread(png_sig, 1, png_sig_len, fp))
    || 0 != png_sig_cmp(png_sig, (png_size_t) 0, png_sig_len)) {
    cout << "Bad signature." << endl;
    (void) fclose(fp);
    exit(EXIT_FAILURE);
  }

  //! Create and initialize the png_struct with the default stderr and error
  //! handling
  png_structp png_ptr;
  if (NULL == (png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL))) {
    cout << "Can't initialize the png_struct." << endl;
    (void) fclose(fp);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    exit(EXIT_FAILURE);
  }

  //! Allocate/initialize the memory for image information
  png_infop info_ptr;
  if (NULL == (info_ptr = png_create_info_struct(png_ptr))) {
    cout << "Can't allocate the memory for image information." << endl;
    (void) fclose(fp);
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    exit(EXIT_FAILURE);
  }

  //! Set error handling
  if (0 != setjmp(png_jmpbuf(png_ptr))) {
    cout << "Can't set the error handling." << endl;
    (void) fclose(fp);
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    exit(EXIT_FAILURE);
  }

  //! Set up the input control using standard C streams
  png_init_io(png_ptr, fp);

  //! let libpng know that some bytes have been read
  png_set_sig_bytes(png_ptr, png_sig_len);

  //! Set the read filter transforms, to get 8bit RGB whatever the original
  //! file may contain:
  //! PNG_TRANSFORM_STRIP_16      strip 16-bit samples to 8 bits
  //! PNG_TRANSFORM_PACKING       expand 1, 2 and 4-bit
  //!                             samples to bytes
  transform |= (PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING);

  //! Read in the entire image at once
  png_read_png(png_ptr, info_ptr, transform, NULL);

  //! Get image informations
  m_width    = (size_t) png_get_image_width (png_ptr, info_ptr);
  m_height   = (size_t) png_get_image_height(png_ptr, info_ptr);
  m_channels = (size_t) png_get_channels    (png_ptr, info_ptr);
  png_bytepp row_pointers;
  row_pointers = png_get_rows(png_ptr, info_ptr);

  //! Allocate the image
  new (this) Image(m_width, m_height, m_channels, i_border);

  //! Deinterlace and convert png RGB RGB RGB 8bit to RRR GGG BBB
  //! the image is deinterlaced layer after layer
  //! this generic loop also works for one single channel
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < m_height; i++) {
      float* oI = this->getPtr(c, i);
      png_bytep row_ptr = row_pointers[i] + c;

      for (size_t j = 0; j < m_width; j++) {
        oI[j] = (float) *row_ptr;
        row_ptr += m_channels;
      }
    }
  }

  //! Test if image is really a color image and exclude the alpha channel
  if (m_channels == 2) {
    cout << "Warning: this image has only 2 channels. Only one will be kept\n";
    m_channels = 1;
  }
  else if (m_channels >= 3) {
    const size_t wh = m_width * m_height;
    size_t k = 0;
    const float* iR = this->getPtr(0, 0);
    const float* iG = this->getPtr(1, 0);
    const float* iB = this->getPtr(2, 0);
    while (k < wh && iR[k] == iG[k] && iR[k] == iB[k]) {
      k++;
    }
    m_channels = (k == wh ? 1 : 3);
  }

  //! clean up and free any memory allocated, close the file
  (void) fclose(fp);
  png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
}


//! Read a tiff-only image.
void Image::readTiff(
  const char* p_name,
  const size_t i_border) {

  //! Open the tiff file
  TIFF *tif = TIFFOpen(p_name, "r");

  //! Check that the file has been open
  if (!tif) {
    cout << "Unable to read TIFF file " << p_name << ". Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Initialization
  uint32 w = 0, h = 0;
  uint16 spp = 0, bps = 0, fmt = 0;

  //! Get the metadata
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
  TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &fmt);

  //! Check the metadata
  if (spp != 1 || bps != (uint16) sizeof(float) * 8 ||
      fmt != SAMPLEFORMAT_IEEEFP) {
    cout << "readTiff: metadata non-conform. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Allocate the image
  new (this) Image(w, h, 1, i_border);

  //! Read the values
  for (size_t i = 0; i < m_height; i++) {
    float* oI = this->getPtr(0, i);
    if (TIFFReadScanline(tif, oI, i, 0) < 0) {
      cout << "readTiff: error reading row " << i << endl;
      exit(EXIT_FAILURE);
    }
  }

  //! Close the file
  TIFFClose(tif);
}


//! Generic write image.
void Image::write(
  const char* p_name,
  const bool p_quad) const {

  //! parameters check
  if (0 >= m_width || 0 >= m_height || 0 >= m_channels) {
    cout << "Inconsistent size." << endl;
    exit(EXIT_FAILURE);
  }

  //! Parameters check
  if (NULL == p_name) {
    cout << "Null name" << endl;
    exit(EXIT_FAILURE);
  }
  //! Check the extension
  string ext = p_name;
  const int pos = ext.find_last_of(".");
  ext = ext.substr(pos + 1, ext.size());

  //! Call the right function
  if (ext == "png") {
    this->writePng(p_name, p_quad);
  }
  else if (ext == "tif" || ext == "tiff") {
    this->writeTiff(p_name, p_quad);
  }
  else {
    cout << "Extension " << ext << " not known. Abort." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Write an image via the libpng write function
void Image::writePng(
  const std::string &p_name,
  const bool p_quad) const {

  //! open the PNG output file
  FILE *volatile fp;
  if (0 == strcmp(p_name.c_str(), "-")) {
    fp = stdout;
  }
  else if (NULL == (fp = fopen(p_name.c_str(), "wb"))) {
    exit(EXIT_FAILURE);
  }

  //! For convenience
  const size_t w = p_quad ? 2 * m_width  : m_width ;
  const size_t h = p_quad ? 2 * m_height : m_height;

  //! allocate the interlaced array and row pointers
  const size_t whc = w* h * m_channels;
  png_byte *idata = NULL;
  if (NULL == (idata = (png_byte *) malloc(whc * sizeof(png_byte)))) {
    (void) fclose(fp);
    exit(EXIT_FAILURE);
  }

  png_bytep *row_pointers = NULL;
  if (NULL == (row_pointers = (png_bytep *) malloc(h * sizeof(png_bytep)))) {
    (void) fclose(fp);
    memfree(idata);
    exit(EXIT_FAILURE);
  }

  //! create and initialize the png_struct with the default stderr and error
  //! handling
  png_structp png_ptr;
  if (NULL == (png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL))) {
    (void) fclose(fp);
    memfree(idata);
    memfree(row_pointers);
    exit(EXIT_FAILURE);
  }

  //! allocate/initialize the memory for image information
  png_infop info_ptr;
  if (NULL == (info_ptr = png_create_info_struct(png_ptr))) {
    (void) fclose(fp);
    memfree(idata);
    memfree(row_pointers);
    png_destroy_write_struct(&png_ptr, NULL);
    exit(EXIT_FAILURE);
  }

  //! set error handling
  if (0 != setjmp(png_jmpbuf(png_ptr))) {
    (void) fclose(fp);
    memfree(idata);
    memfree(row_pointers);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    exit(EXIT_FAILURE);
  }

  //! set up the input control using standard C streams
  png_init_io(png_ptr, fp);

  //! set image informations
  const png_byte bit_depth = 8;
  int color_type;
  switch (m_channels) {
    case 1:
      color_type = PNG_COLOR_TYPE_GRAY;
      break;
    case 2:
      color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
      break;
    case 3:
      color_type = PNG_COLOR_TYPE_RGB;
      break;
    case 4:
      color_type = PNG_COLOR_TYPE_RGB_ALPHA;
      break;
    default:
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      memfree(row_pointers);
      memfree(idata);
      (void) fclose(fp);
      exit(EXIT_FAILURE);
  }
  const int interlace = PNG_INTERLACE_ADAM7;
  const int compression = PNG_COMPRESSION_TYPE_BASE;
  const int filter = PNG_FILTER_TYPE_BASE;

  //! set image header
  png_set_IHDR(png_ptr, info_ptr, (png_uint_32) w, (png_uint_32) h,
    bit_depth, color_type, interlace, compression, filter);
  png_write_info(png_ptr, info_ptr);

  //!interlace and convert RRR GGG BBB to RGB RGB RGB
  //! the image is interlaced layer after layer
  //! this involves more memory exchange, but allows a generic loop
  for (size_t c = 0; c < m_channels; c++) {
    png_byte *idata_ptr = idata + c;

    for (size_t i = 0; i < h; i++) {
      const float* iI = this->getPtr(c, p_quad ? i / 2 : i);

      for (size_t j = 0; j < w; j++) {
        const int tmp = (int) floor((float) iI[p_quad ? j / 2 : j] + 0.5f);
        *idata_ptr = (png_byte) std::min(255, std::max(0, tmp));
        idata_ptr += m_channels;
      }
    }
  }

  //! set row pointers
  for (size_t i = 0; i < h; i++) {
    row_pointers[i] = idata + (size_t) (m_channels * w * i);
  }

  //! write out the entire image and end it
  png_write_image(png_ptr, row_pointers);
  png_write_end(png_ptr, info_ptr);

  //! clean up and memfree any memory allocated, close the file
  (void) fclose(fp);
  memfree(idata);
  memfree(row_pointers);
  png_destroy_write_struct(&png_ptr, &info_ptr);
}


//! Write an image via the libtiff write function
void Image::writeTiff(
  const std::string &p_name,
  const bool p_quad) const {

  float im_min, im_max, im_mean;
  this->getStats(im_min, im_max, im_mean);
  printf("Saving %s (min=%.2f, max=%.2f, mean=%.2f)\n", p_name.c_str(), im_min, im_max, im_mean);

  //! Open the file
  TIFF *tif = TIFFOpen(p_name.c_str(), "w");

  //! Check that the file has been created
  if (!tif) {
    cout << "Unable to write TIFF file " << p_name << endl;
    exit(EXIT_FAILURE);
  }

  //! For convenience
  const size_t h = m_height * (p_quad ? 2 : 1);
  const size_t w = m_width  * (p_quad ? 2 : 1);
  uint32 rowsperstrip;
  float* line = (float*) memalloc(16, (p_quad ? w : m_width) * sizeof(float));

  //! Header of the tiff file
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32) w);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32) h);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, (uint16) m_channels);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, (uint16) sizeof(float) * 8);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  rowsperstrip = TIFFDefaultStripSize(tif, (uint32) h);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

  //! Write the file data
  for (size_t c = 0, ok = 1; ok && c < m_channels; c++) {
    for (size_t i = 0; ok && i < h; i++) {
      const float* iI = this->getPtr(c, p_quad ? i / 2 : i);

      //! Copy the line
      if (p_quad) {
        for (size_t j = 0; j < m_width; j++) {
          line[2 * j + 0] = line[2 * j + 1] = iI[j];
        }
      }
      else {
        for (size_t j = 0; j < m_width; j++) {
          line[j] = iI[j];
        }
      }

      //! Write the line
      if (TIFFWriteScanline(tif, line, (uint32) i, (tsample_t) c) < 0) {
        cout << "WriteTIFF: error writing row " << i << endl;
        ok = 0;
      }
    }
  }

  //! Release memory
  memfree(line);

  //! Close the file
  TIFFClose(tif);
}


//! Add a border to the current image.
void Image::addBorder(
  const size_t i_border) {

  //! Temporary image with border
  Image tmp(m_width, m_height, m_channels, i_border);

  //! Copy the inner image into it
  this->copyInner(tmp);

  //! Retrieve the bigger image
  *this = tmp;
}


//! Set a value for the border of the current image.
void Image::setBorder(
  const float p_value) {

#if defined(__AVX__)
  const __m256 zValue = _mm256_set1_ps(p_value);
#endif
#if defined(__SSE__)
  const __m128 xValue = _mm_set1_ps(p_value);
#endif

  //! For every channels
  for (size_t c = 0; c < m_channels; c++) {

    //! Top and bottom borders
    for (int n = 0; n < int(m_border); n++) {
      float* iT = this->getPtr(c, - n - 1);
      float* iB = this->getPtr(c, m_height + n);
      int j = -int(m_border);

#if defined(__AVX__)
      //! AVX version
      for (; j < int(m_width + m_border) - 8; j += 8) {
        _mm256_storeu_ps(iT + j, zValue);
        _mm256_storeu_ps(iB + j, zValue);
      }
#endif
#if defined(__SSE__)
      //! SSE version
      for (; j < int(m_width + m_border) - 4; j += 4) {
        _mm_storeu_ps(iT + j, xValue);
        _mm_storeu_ps(iB + j, xValue);
      }
#endif
      //! Normal version
      for (; j < int(m_width + m_border); j++) {
        iT[j] = p_value;
        iB[j] = p_value;
      }
    }

    //! Left and right borders
    for (int i = -int(m_border); i < int(m_height + m_border); i++) {
      float* iT = this->getPtr(c, i);

      for (int n = 0; n < int(m_border); n++) {
        iT[-n - 1] = p_value;
        iT[m_width + n] = p_value;
      }
    }
  }
}


//! Copy the inner image of the current image to another image.
void Image::copyInner(
  Image& o_im,
  const int p_tid) const {

  //! Check the size of the output image.
  if (o_im.width   () != m_width   ||
      o_im.height  () != m_height  ||
      o_im.channels() != m_channels) {

    o_im.init(m_width, m_height, m_channels, o_im.border());
  }

  //! For convenience
  const size_t hBegin = p_tid >= 0 ? m_heights[p_tid] : 0;
  const size_t hEnd   = p_tid >= 0 ? m_heights[p_tid + 1] : m_height;

  //! Copy the values
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = hBegin; i < hEnd; i++) {
      const float* iI = this->getPtr(c, i);
      float* oI = o_im.getPtr(c, i);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        _mm256_storeu_ps(oI + j, _mm256_loadu_ps(iI + j));
      }
#endif
#if defined(__SSE__)
      //! AVX version
      for (; j < m_width - 4; j += 4) {
        _mm_storeu_ps(oI + j, _mm_loadu_ps(iI + j));
      }
#endif
      //! Normal version
      for (; j < m_width; j++) {
        oI[j] = iI[j];
      }
    }
  }
}


//! Apply a zoom-in of the input image as array and normalize it.
void Image::applySqrtNormalization(
  Image& o_im,
  const float p_min,
  const float p_max) const {

  //! Get the min and max values
  float minVal, maxVal, mean;
  this->getStats(minVal, maxVal, mean);
  minVal = std::max(0.f, sqrtf(minVal));
  maxVal = std::max(0.f, sqrtf(maxVal));

  //! To avoid problems of saturation
  maxVal = std::min(float(1), maxVal);

  //! Initialization of the return image
  o_im.init(m_width, m_height, m_channels, 0);

  //! For convenience
  const float norm = (p_max - p_min) / (maxVal - minVal);
#if defined(__AVX__)
  const __m256 zMin  = _mm256_set1_ps(minVal);
  const __m256 zNorm = _mm256_set1_ps(norm);
  const __m256 zpMin = _mm256_set1_ps(p_min);
  const __m256 z0    = _mm256_setzero_ps();
#endif
#if defined(__SSE__)
  const __m128 xMin  = _mm_set1_ps(minVal);
  const __m128 xNorm = _mm_set1_ps(norm);
  const __m128 xpMin = _mm_set1_ps(p_min);
  const __m128 x0    = _mm_setzero_ps();
#endif

  //! Normalization
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < m_height; i++) {
      const float* iI = this->getPtr(c, i);
      float* oI = o_im.getPtr(c, i);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        _mm256_storeu_ps(oI + j, zpMin + zNorm *
          (_mm256_sqrt_ps(_mm256_max_ps(z0, _mm256_loadu_ps(iI + j))) - zMin));
      }
#endif
#if defined(__SSE__)
      //! SSE version
      for (; j < m_width - 4; j += 4) {
        _mm_storeu_ps(oI + j, xpMin + xNorm *
          (_mm_sqrt_ps(_mm_max_ps(x0, _mm_loadu_ps(iI + j))) - xMin));
      }
#endif
      for (; j < m_width; j++) {
        oI[j] = p_min + norm * (sqrt(std::max(float(0), iI[j])) - minVal);
      }
    }
  }
}


//! Apply a zoom-in of the input image as array and normalize it.
void Image::applySqrtNormalization(
  Image& o_im,
  const size_t p_tid,
  const float p_min,
  const float p_max,
  const float i_min,
  const float i_max) const {

  //! The output image must have been initialized
  if (o_im.width   () != m_width  ||
      o_im.height  () != m_height ||
      o_im.channels() != m_channels) {
    cout << "applyQuadSqrtNormalization: the output image must have been ";
    cout << "initialized before." << endl;
    exit(EXIT_FAILURE);
  }

  //! For convenience
  const float minVal = sqrt(std::max(float(0), i_min));
  const float maxVal = sqrt(std::max(float(0), i_max));
  const float norm   = (p_max - p_min) / (maxVal - minVal);
#if defined(__AVX__)
  const __m256 zMin  = _mm256_set1_ps(minVal);
  const __m256 zNorm = _mm256_set1_ps(norm);
  const __m256 zpMin = _mm256_set1_ps(p_min);
  const __m256 z0    = _mm256_setzero_ps();
#endif
#if defined(__SSE__)
  const __m128 xMin  = _mm_set1_ps(minVal);
  const __m128 xNorm = _mm_set1_ps(norm);
  const __m128 xpMin = _mm_set1_ps(p_min);
  const __m128 x0    = _mm_setzero_ps();
#endif

  //! Normalization
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = m_heights[p_tid]; i < m_heights[p_tid + 1]; i++) {
      const float* iI = this->getPtr(c, i);
      float* oI = o_im.getPtr(c, i);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        _mm256_storeu_ps(oI + j, zpMin + zNorm *
          (_mm256_sqrt_ps(_mm256_max_ps(z0, _mm256_loadu_ps(iI + j))) - zMin));
      }
#endif
#if defined(__SSE__)
      //! SSE version
      for (; j < m_width - 4; j += 4) {
        _mm_storeu_ps(oI + j, xpMin + xNorm *
          (_mm_sqrt_ps(_mm_max_ps(x0, _mm_loadu_ps(iI + j))) - xMin));
      }
#endif
      for (; j < m_width; j++) {
        oI[j] = p_min + norm * (sqrt(std::max(float(0), iI[j])) - minVal);
      }
    }
  }
}


//! Multiply the image by a constant.
void Image::multiply(
  const float p_constant,
  const int p_tid) {

  //! For convenience
  const size_t hBegin = p_tid >= 0 ? m_heights[p_tid] : 0;
  const size_t hEnd   = p_tid >= 0 ? m_heights[p_tid + 1] : m_height;
#if defined(__AVX__)
  const __m256 zConstant = _mm256_set1_ps(p_constant);
#endif
#if defined(__SSE__)
  const __m128 xConstant = _mm_set1_ps(p_constant);
#endif

  //! Loop over channels
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = hBegin; i < hEnd; i++) {
      float* oI = this->getPtr(c, i);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        _mm256_storeu_ps(oI + j, _mm256_loadu_ps(oI + j) * zConstant);
      }
#endif
#if defined(__SSE__)
      //! SSE version
      for (; j < m_width - 4; j += 4) {
        _mm_storeu_ps(oI + j, _mm_loadu_ps(oI + j) * xConstant);
      }
#endif
      //! Normal version
      for (; j < m_width; j++) {
        oI[j] *= p_constant;
      }
    }
  }
}


//! Return a zoom-out of an image by simple extraction
void Image::zoomOut(
  Image& o_im) const {

  //! New image 4 times smaller
  o_im.init(m_width / 2, m_height / 2, m_channels, m_border);

  //! Apply a Gaussian to the image
  Image imBlur;
  this->applyGaussianBlur(imBlur, float(1.4));

  //! Get the values
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < o_im.height(); i++) {
      const float* iI = imBlur.getPtr(c, 2 * i);
      float* oI = o_im.getPtr(c, i);

      for (size_t j = 0; j < o_im.width(); j++) {
        oI[j] = iI[2 * j];
      }
    }
  }
}


//! Return a zoom-in of an image by simple extraction
void Image::zoomIn(
  Image& o_im) const {

  //! New image 4 times bigger
  o_im.init(m_width * 2, m_height * 2, m_channels, m_border);

  //! Get the values
  for (size_t chnl = 0; chnl < m_channels; chnl++) {

    //! Inner image
    for (size_t i = 0; i < m_height - 1; i++) {
      const float* iI1 = this->getPtr(chnl, i);
      const float* iI2 = this->getPtr(chnl, i + 1);
      float* oI1 = o_im.getPtr(chnl, 2 * i);
      float* oI2 = o_im.getPtr(chnl, 2 * i + 1);

      for (size_t j = 0; j < m_width - 1; j++) {

        //! Get the 4 corners
        const float a = iI1[j];
        const float b = iI1[j + 1];
        const float c = iI2[j];
        const float d = iI2[j + 1];

        //! Get the values
        oI1[2 * j    ] = a;
        oI1[2 * j + 1] = (a + b) / float(2);
        oI2[2 * j    ] = (a + c) / float(2);
        oI2[2 * j + 1] = (a + b + c + d) / float(4);
      }
    }

    //! Right border
    for (size_t i = 0; i < m_height - 1; i++) {
      const float* iI1 = this->getPtr(chnl, i);
      const float* iI2 = this->getPtr(chnl, i + 1);
      float* oI1 = o_im.getPtr(chnl, 2 * i);
      float* oI2 = o_im.getPtr(chnl, 2 * i + 1);

      //! Get the 4 corners
      const float a = iI1[m_width - 1];
      const float b = 2 * iI1[m_width - 1] - iI1[m_width - 2];
      const float c = iI2[m_width - 2];
      const float d = 2 * iI2[m_width - 1] - iI2[m_width - 2];

        //! Get the values
        oI1[2 * (m_width - 1)    ] = a;
        oI1[2 * (m_width - 1) + 1] = (a + b) / float(2);
        oI2[2 * (m_width - 1)    ] = (a + c) / float(2);
        oI2[2 * (m_width - 1) + 1] = (a + b + c + d) / float(4);
    }


    //! Bottom border
    const float* iI1 = this->getPtr(chnl, m_height - 1);
    const float* iI2 = this->getPtr(chnl, m_height - 2);
    float* oI1 = o_im.getPtr(chnl, 2 * (m_height - 1));
    float* oI2 = o_im.getPtr(chnl, 2 * (m_height - 1) + 1);

    for (size_t j = 0; j < m_width - 1; j++) {

      //! Get the 4 corners
      const float a = iI1[j];
      const float b = 2 * iI1[j] - iI2[j];
      const float c = iI1[j + 1];
      const float d = 2 * iI1[j + 1] - iI2[j + 1];

        //! Get the values
        oI1[2 * j    ] = a;
        oI1[2 * j + 1] = (a + b) / float(2);
        oI2[2 * j    ] = (a + c) / float(2);
        oI2[2 * j + 1] = (a + b + c + d) / float(4);
    }

    //! Bottom-right corner
    const float a = iI1[m_width - 1];
    const float b = 2 * iI1[m_width - 1] - iI1[m_width - 2];
    const float c = 2 * iI1[m_width - 1] - iI2[m_width - 1];
    const float d = 2 * iI1[m_width - 1] - iI2[m_width - 2];

    oI1[2 * (m_width - 1)    ] = a;
    oI1[2 * (m_width - 1) + 1] = (a + b) / float(2);
    oI2[2 * (m_width - 1)    ] = (a + c) / float(2);
    oI2[2 * (m_width - 1) + 1] = (a + b + c + d) / float(4);
  }
}


//! In-place normalization of an image.
void Image::normalize(
  const float p_min,
  const float p_max) {

  //! Get the min and max values
  float minVal, maxVal, mean;
  this->getStats(minVal, maxVal, mean);
  const float norm = (p_max - p_min) / (maxVal - minVal);

#if defined(__AVX__)
  const __m256 zpMin = _mm256_set1_ps(p_min);
  const __m256 zMin  = _mm256_set1_ps(minVal);
  const __m256 zNorm = _mm256_set1_ps(norm);
#endif
#if defined(__SSE__)
  const __m128 xpMin = _mm_set1_ps(p_min);
  const __m128 xMin  = _mm_set1_ps(minVal);
  const __m128 xNorm = _mm_set1_ps(norm);
#endif

  //! Normalization
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < m_height; i++) {
      float* oI = this->getPtr(c, i);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        _mm256_storeu_ps(oI + j, zpMin + zNorm * (_mm256_loadu_ps(oI + j) - zMin));
      }
#endif
#if defined(__SSE__)
      //! SSE version
      for (; j < m_width - 4; j+= 4) {
        _mm_storeu_ps(oI + j, xpMin + xNorm * (_mm_loadu_ps(oI + j) - xMin));
      }
#endif
      //! Normal version
      for (; j < m_width; j++) {
        oI[j] = p_min + norm * (oI[j] - minVal);
      }
    }
  }
}


//! Save the image as a 3D model in OBJ format.
void Image::saveAsObj(
  const std::string &p_fileName) const {

  //! Open the file
  ofstream file(p_fileName.c_str());

  //! Subsampling until the size is ok
  const size_t maxSize    = std::max(m_width, m_height);

  //! For convenience
  const float cw = float(100) / float(maxSize);
  const float ch = float(100) / float(maxSize);
  const float cl = float(10) / float(255);
  const int hh = m_height / 2;
  const int hw = m_width  / 2;

  //! Check that the file is open
  if (file.is_open()) {

    //! Buffer for each threads
    std::vector<std::string> points(m_nbThreads);
    std::vector<std::string> faces (m_nbThreads);

    //! Values (all have the same intensity between [0, 1]
#ifdef _OPENMP
#pragma omp parallel
{
    //! Get the current thread id
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // _OPENMP

    for (int i = (int) m_heights[tid]; i < (int) m_heights[tid + 1]; i++) {
      const float* iI = this->getPtr(0, i);

      for (int j = 0; j < (int) m_width; j++) {
        const float v  = iI[j];
        //const float a = std::max(v, float(10));

        //! Points cloud
        {
          ostringstream oss;
          oss << "v " << (j - hw) * cw << " " << (i - hh) * ch << " " << v * cl;
          //oss << " " << a << " " << a << " " << a << "\n";
          points[tid] += oss.str() + "\n";
        }

        //! Faces
        if (i < int(m_height - 1) && j < int(m_width - 1)) {
          ostringstream oss1;
          oss1 <<  i      * m_width + j + 1 << " ";
          oss1 << (i + 1) * m_width + j + 1 << " ";
          oss1 <<  i      * m_width + j + 2 << "\n";
          faces[tid] += "f " + oss1.str();
          ostringstream oss2;
          oss2 <<  i      * m_width + j + 2 << " ";
          oss2 << (i + 1) * m_width + j + 1 << " ";
          oss2 << (i + 1) * m_width + j + 2 << "\n";
          faces[tid] += "f " + oss2.str();
        }
      }
    }
#ifdef _OPENMP
}
#endif // _OPENMP

    //! Header
    file << "####\n";
    file << "#\n";
    file << "# OBJ File Generated by LEM\n";
    file << "#\n";
    file << "####\n";
    file << "# Object LEM.obj\n";
    file << "#\n";
    file << "# Vertices: "<< (m_width * m_height) << "\n";
    file << "# Faces: 0\n";
    file << "#\n";
    file << "####\n";

    //! Write the buffers
    for (size_t n = 0; n < m_nbThreads; n++) {
      file << points[n];
    }
    file << "\n";
    for (size_t n = 0; n < m_nbThreads; n++) {
      file << faces[n];
    }

    //! Bottom
    file << "# " << m_width * m_height << " vertices, 0 vertices normals\n\n";
    file << "# 0 faces, 0 coords texture\n\n";
    file << "# End of File\n";

    //! Close the file
    file.close();
  }
  else {
    cout << "Can't save the file " << p_fileName << "!" << endl;
    exit(EXIT_FAILURE);
  }
}


//! Save the image as a 3D model in OBJ format.
void Image::saveAsObj(
  const std::string &p_fileName,
  const Image& i_imWater,
  const float p_threshold,
  const float p_oceanLevel) const {

  //! Open the file
  ofstream file(p_fileName.c_str());

  //! For convenience
  const float cw = float(100) / float(std::max(m_width, m_height));
  const float ch = float(100) / float(std::max(m_width, m_height));
  const float cl = float(10) / float(255);
  const int hh = m_height / 2;
  const int hw = m_width  / 2;

  //! Check that the file is open
  if (file.is_open()) {

    //! Buffer for each threads
    std::vector<std::string> points(m_nbThreads);
    std::vector<std::string> faces (m_nbThreads);

    //! Values (all have the same intensity between [0, 1]
#ifdef _OPENMP
#pragma omp parallel
{
    //! Get the current thread id
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // _OPENMP

    for (int i = (int) m_heights[tid]; i < (int) m_heights[tid + 1]; i++) {
      const float* iI = this->getPtr(0, i);
      const float* iW = i_imWater.getPtr(0, i);

      for (int j = 0; j < (int) m_width; j++) {
        const float v  = iI[j];
        //const float a = std::max(v, float(10));
        const float w = iW[j];

        //! Points cloud
        if (v <= p_oceanLevel) {
          ostringstream oss;
          oss << "v " << (j - hw) * cw << " " << (i - hh) * ch << " " << p_oceanLevel * cl;
          //oss << " 0 0 192\n";
          points[tid] += oss.str() + "\n";
        }
        else if (w > p_threshold) {
          ostringstream oss;
          oss << "v " << (j - hw) * cw << " " << (i - hh) * ch << " " << (v + w) * cl;
          //oss << " 0 0 192\n";
          points[tid] += oss.str() + "\n";
        }
        else {
          ostringstream oss;
          oss << "v " << (j - hw) * cw << " " << (i - hh) * ch << " " << v * cl;
          //oss << " " << a << " " << a << " " << a << "\n";
          points[tid] += oss.str() + "\n";
        }

        //! Faces
        if (i < int(m_height - 1) && j < int(m_width - 1)) {
          ostringstream oss1;
          oss1 <<  i      * m_width + j + 1 << " ";
          oss1 << (i + 1) * m_width + j + 1 << " ";
          oss1 <<  i      * m_width + j + 2 << "\n";
          faces[tid] += "f " + oss1.str();
          ostringstream oss2;
          oss2 <<  i      * m_width + j + 2 << " ";
          oss2 << (i + 1) * m_width + j + 1 << " ";
          oss2 << (i + 1) * m_width + j + 2 << "\n";
          faces[tid] += "f " + oss2.str();
        }
      }
    }
#ifdef _OPENMP
}
#endif // _OPENMP

    //! Header
    file << "####\n";
    file << "#\n";
    file << "# OBJ File Generated by LEM\n";
    file << "#\n";
    file << "####\n";
    file << "# Object LEM.obj\n";
    file << "#\n";
    file << "# Vertices: "<< m_width * m_height << "\n";
    file << "# Faces: 0\n";
    file << "#\n";
    file << "####\n";

    //! Write the buffers
    for (size_t n = 0; n < m_nbThreads; n++) {
      file << points[n];
    }
    file << "\n";
    for (size_t n = 0; n < m_nbThreads; n++) {
      file << faces[n];
    }

    //! Bottom
    file << "# " << m_width * m_height << " vertices, 0 vertices normals\n\n";
    file << "# 0 faces, 0 coords texture\n\n";
    file << "# End of File\n";

    //! Close the file
    file.close();
  }
  else {
    cout << "Can't save the file " << p_fileName << "!" << endl;
    exit(EXIT_FAILURE);
  }
}


//! Save the image as a 3D model in OBJ format with false color.
void Image::saveAsObj(
  const std::string &p_fileName,
  const Image& i_imWater,
  const Image& i_imColor,
  const float p_oceanLevel,
  const float p_norm) const {

  //! Get some statistics from the water image
  float minVal, maxVal, mean;
  i_imWater.getStats(minVal, maxVal, mean);

  //! Open the file
  ofstream file(p_fileName.c_str());

  //! For convenience
  const float cw = 1;//float(100) / float(std::max(m_width, m_height));
  const float ch = 1;//float(100) / float(std::max(m_width, m_height));
  const float cl = p_norm;//float(10) / float(255);
  const int hh = m_height / 2;
  const int hw = m_width  / 2;

  //! Check that the file is open
  if (file.is_open()) {

    //! Buffer for each threads
    std::vector<std::string> points(m_nbThreads);
    std::vector<std::string> faces (m_nbThreads);

    //! Values (all have the same intensity between [0, 1])
#ifdef _OPENMP
#pragma omp parallel
{
    //! Get the current thread id
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // _OPENMP

    for (int i = (int) m_heights[tid]; i < (int) m_heights[tid + 1]; i++) {
      const float* iI = this->getPtr(0, i);
      const float* iW = i_imWater.getPtr(0, i);
      const float* iR = i_imColor.getPtr(0, i);
      const float* iG = i_imColor.getPtr(1, i);
      const float* iB = i_imColor.getPtr(2, i);

      for (int j = 0; j < (int) m_width; j++) {
        const float v = iI[j];
        const float w = iW[j];

        //! Points cloud
        if (v <= p_oceanLevel) {
          ostringstream oss;
          oss << "v " << (j - hw) * cw << " " << (i - hh) * ch << " " << p_oceanLevel * cl;
          oss << " " << iR[j] << " " << iG[j] << " " << iB[j];
          points[tid] += oss.str() + "\n";
        }
        else if (w > mean) {
          ostringstream oss;
          oss << "v " << (j - hw) * cw << " " << (i - hh) * ch << " " << (v + w) * cl;
          oss << " " << iR[j] << " " << iG[j] << " " << iB[j];
          points[tid] += oss.str() + "\n";
        }
        else {
          ostringstream oss;
          oss << "v " << (j - hw) * cw << " " << (i - hh) * ch << " " << v * cl;
          oss << " " << iR[j] << " " << iG[j] << " " << iB[j];
          points[tid] += oss.str() + "\n";
        }

        //! Faces
        if (i < int(m_height - 1) && j < int(m_width - 1)) {
          ostringstream oss1;
          oss1 <<  i      * m_width + j + 1 << " ";
          oss1 << (i + 1) * m_width + j + 1 << " ";
          oss1 <<  i      * m_width + j + 2 << "\n";
          faces[tid] += "f " + oss1.str();
          ostringstream oss2;
          oss2 <<  i      * m_width + j + 2 << " ";
          oss2 << (i + 1) * m_width + j + 1 << " ";
          oss2 << (i + 1) * m_width + j + 2 << "\n";
          faces[tid] += "f " + oss2.str();
        }
      }
    }
#ifdef _OPENMP
}
#endif // _OPENMP

    //! Header
    file << "####\n";
    file << "#\n";
    file << "# OBJ File Generated by LEM\n";
    file << "#\n";
    file << "####\n";
    file << "# Object LEM.obj\n";
    file << "#\n";
    file << "# Vertices: "<< m_width * m_height << "\n";
    file << "# Faces: 0\n";
    file << "#\n";
    file << "####\n";

    //! Write the buffers
    for (size_t n = 0; n < m_nbThreads; n++) {
      file << points[n];
    }
    file << "\n";
    for (size_t n = 0; n < m_nbThreads; n++) {
      file << faces[n];
    }

    //! Bottom
    file << "# " << m_width * m_height << " vertices, 0 vertices normals\n\n";
    file << "# 0 faces, 0 coords texture\n\n";
    file << "# End of File\n";

    //! Close the file
    file.close();
  }
  else {
    cout << "Can't save the file " << p_fileName << "!" << endl;
    exit(EXIT_FAILURE);
  }
}


//! Apply a Gaussian filter.
void Image::applyGaussianBlur(
  Image& o_im,
  const float p_sigma) const {

  //! Initialization
  int sP = ceil(float(6) * p_sigma);
  sP += 1 - (sP % 2);
  const int sP2 = (sP - 1) / 2;
  o_im = Image(m_width, m_height, m_channels, m_border);

  //! Check the size of the image
  if (sP2 > (int) std::min(m_height - 1, m_width - 1)) {
    const float sigma = float(std::min(m_height - 1, m_width - 1)) / float(6);
    cout << "gaussianBlur: sigma has been adjusted from " << p_sigma << " to ";
    cout << sigma << " due to the size of the image." << endl;
    this->applyGaussianBlur(o_im, sigma);
    return;
  }

  //! Compute the Kernel
  float* gKern  =(float*) memalloc(16, sP * sizeof(float));
  float norm = 0.f;
  for (int n = 0; n < sP; n++) {
    gKern[n] = exp(-float((n - sP2) * (n - sP2)) / (float(2) * p_sigma * p_sigma));
    norm += gKern[n];
  }
  norm = float(1) / norm;
  for (int n = 0; n < sP; n++) {
    gKern[n] *= norm;
  }

  //! Temporary image
  Image im1(m_width, m_height, m_channels, sP2);
  Image im2(m_width, m_height, m_channels, sP2);

#ifdef _OPENMP
#pragma omp parallel
{

  const size_t tid = omp_get_thread_num();
#else
  const size_t tid = 0;
#endif // _OPENMP

  //! Temporary image to avoid boundary effect
  this->copyInner(im1, tid);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

  im1.fillBorder(tid);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

  const int hBegin = m_heights[tid];
  const int hEnd   = m_heights[tid + 1];

  for (size_t c = 0; c < m_channels; c++) {

    //! Compute the blur in column
    for (int i = hBegin; i < hEnd; i++) {
      float* oI = im2.getPtr(c, i);

      for (int n = 0; n < sP; n++) {
        const float* iI = im1.getPtr(c, i + n - sP2);
        const float val = gKern[n];
        size_t j = 0;

#if defined(__AVX__)
        //! AVX version
        const __m256 zVal = _mm256_set1_ps(val);
        for (; j < m_width - 8; j += 8) {
          _mm256_storeu_ps(oI + j, _mm256_loadu_ps(oI + j) +
                                   _mm256_loadu_ps(iI + j) * zVal);
        }
#endif
#if defined(__SSE__)
        //! SSE version
        const __m128 xVal = _mm_set1_ps(val);
        for (; j < m_width - 4; j += 4) {
          _mm_storeu_ps(oI + j, _mm_loadu_ps(oI + j) +
                                _mm_loadu_ps(iI + j) * xVal);
        }
#endif
        //! Normal version
        for (; j < m_width; j++) {
          oI[j] += val * iI[j];
        }
      }
    }
  }

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

  im2.fillBorder(tid);

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP

  for (size_t c = 0; c < m_channels; c++) {

    //! Compute the blur in row
    for (int i = hBegin; i < hEnd; i++) {
      float* oI = o_im.getPtr(c, i);

      for (int n = 0; n < sP; n++) {
        const float* iI = im2.getPtr(c, i);
        const float val = gKern[n];
        int j = 0;

#if defined(__AVX__)
        //! AVX version
        const __m256 zVal = _mm256_set1_ps(val);
        for (; j < int(m_width) - 8; j += 8) {
          _mm256_storeu_ps(oI + j, _mm256_loadu_ps(oI + j) +
                                   _mm256_loadu_ps(iI + j + n - sP2) * zVal);
        }
#endif
#if defined(__SSE__)
        //! SSE version
        const __m128 xVal = _mm_set1_ps(val);
        for (; j < int(m_width) - 4; j += 4) {
          _mm_storeu_ps(oI + j, _mm_loadu_ps(oI + j) +
                                _mm_loadu_ps(iI + j + n - sP2) * xVal);
        }
#endif
        //! Normal version
        for (; j < (int) m_width; j++) {
          oI[j] += val * iI[j + n - sP2];
        }
      }
    }
  }
#ifdef _OPENMP
}
#endif

  //! Release memory
  memfree(gKern);
}


//! Apply in-place a Gaussian noise of standard deviation sigma.
void Image::applyGaussianNoise(
  const float p_sigma,
  const float p_oceanLvl) {

  //! Noise initialization
  mt_init_genrand((unsigned long int) time (NULL) +
                  (unsigned long int) getpid());

  //! For each pixels
  for (size_t k = 0; k < m_size; k++) {

    //! Add the noise
    const double a = mt_genrand_res53();
    const double b = mt_genrand_res53();

    //! Don't add noise under the ocean level
    if (m_ptr[k] > p_oceanLvl) {
      m_ptr[k] += p_sigma * (float) (sqrtl(-2.l * log(a)) * cos(2.l * M_PI * b));
    }
  }
}


//! Symmetrize the borders.
void Image::fillBorder(
  const int p_tid,
  const size_t p_sym) {

  //! Check the size
  const size_t hBegin = p_tid >= 0 ? m_heights[p_tid    ] : 0;
  const size_t hEnd   = p_tid >= 0 ? m_heights[p_tid + 1] : m_height;

  //! For convenience
  const int w = m_width;
  const int h = m_height;
  const int b = m_border;

  //! For each channel
  for (size_t c = 0; c < m_channels; c++) {

    //! Left and right
    for (size_t i = hBegin; i < hEnd; i++) {
      float* iI = this->getPtr(c, i);

      for (int j = 0; j < (int) m_border; j++) {

        switch(p_sym) {
          case 0: //! Neumann conditions
              iI[j - b] = iI[0];
              iI[w + j] = iI[w - 1];
            break;
          case 1: //! Symmetry cb abc ba
            iI[j - b] = iI[b - j];
            iI[w + j] = iI[w - j - 2];
            break;
          case 2: //! Affine interpolation
            iI[j - b] = std::max((iI[1] - iI[0]) * (j - b) + iI[0], 0.f);
            iI[w + j] = std::max((iI[w - 1] - iI[w - 2]) * (j + 1) + iI[w - 1], 0.f);
            break;
          case 3: //! Dirichlet conditions
            iI[j - b] = 0.f;
            iI[w + j] = 0.f;
            break;
          default: //! Unknown conditions
            cout << "Unknown value for fillBorder. Abort." << endl;
            exit(EXIT_FAILURE);
            break;
        }
      }
    }

    //! Top
    if (hBegin == 0) {
      for (int i = 0; i < (int) m_border; i++) {
        const float* iT = this->getPtr(c, b - i);
        const float* iT0 = this->getPtr(c, 0);
        const float* iT1 = this->getPtr(c, 1);
        float* oT = this->getPtr(c, i - b);

        for (int j = -int(m_border); j < int(m_width + m_border); j++) {

          switch (p_sym) {
            case 0: //! Neumann conditions
              oT[j] = iT0[j];
              break;
            case 1: //! Symmetry cb abc ba
              oT[j] = iT[j];
              break;
            case 2: //! Affine interpolation
              oT[j] = std::max((iT1[j] - iT0[j]) * (i - b) + iT0[j], 0.f);
              break;
            case 3: //! Dirichlet conditions
              oT[j] = 0.f;
              break;
            default: //! Unknown conditions
              cout << "Unknown value for fillBorder. Abort." << endl;
              exit(EXIT_FAILURE);
              break;
          }
        }
      }
    }

    //! Bottom
    if (hEnd == m_height) {
      for (int i = 0; i < (int) m_border; i++) {
        const float* iB = this->getPtr(c, h - i - 2);
        const float* iBw1 = this->getPtr(c, h - 1);
        const float* iBw2 = this->getPtr(c, h - 2);
        float* oB = this->getPtr(c, h + i);

        for (int j = -int(m_border); j < int(m_width + m_border); j++) {

          switch (p_sym) {
            case 0: //! Neumann conditions
              oB[j] = iBw1[j];
              break;
            case 1: //! Symmetry cb abc ba
              oB[j] = iB[j];
              break;
            case 2: //! Affine interpolation
              oB[j] = std::max((iBw1[j] - iBw2[j]) * (i + 1) + iBw1[j], 0.f);
              break;
            case 3: //! Dirichlet conditions
              oB[j] = 0.f;
              break;
            default: //! Unknown conditions
              cout << "Unknown value for fillBorder. Abort." << endl;
              exit(EXIT_FAILURE);
              break;
          }
        }
      }
    }
  }
}


//! Release the memory
void Image::releaseMemory() {

  //! Release the memory
  memfree(m_ptr);
  memfree(m_heights);
}


//! Upsampling until (width * height) is close to size.
void Image::sample(
  Image& o_im,
  const size_t p_size,
  const int p_tid) const {

  //! Initialize the image
  const size_t width  = sqrt(float(p_size * m_width ) / float(m_height));
  const size_t height = sqrt(float(p_size * m_height) / float(m_width ));
  if (p_tid >= 0) {
    if (o_im.width   () != width  ||
        o_im.height  () != height ||
        o_im.channels() != m_channels) {
      cout << "The output image must have been initialized." << endl;
      cout << "w = " << o_im.width() << " (" << width << ")" << endl;
      cout << "h = " << o_im.height() << " (" << height << ")" << endl;
      cout << "c = " << o_im.channels() << " (" << m_channels << ")" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {
    o_im.init(width, height, m_channels, m_border);
  }

  //! Get the step
  const float di = float(m_height - 1) / float(height - 1);
  const float dj = float(m_width  - 1) / float(width  - 1);

  //! For convenience
  const size_t hBegin = p_tid >= 0 ? o_im.nHeight(p_tid    ) : 0;
  const size_t hEnd   = p_tid >= 0 ? o_im.nHeight(p_tid + 1) : height;

  //! Loop over channels
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = hBegin; i < hEnd; i++) {
      float* oI = o_im.getPtr(c, i);

      //! Get the current row position
      const size_t i0  = (size_t) floor(float(i) * di);
      const size_t i1  = std::min(i0 + 1, m_height - 1);
      const float* iI0 = this->getPtr(c, i0);
      const float* iI1 = this->getPtr(c, i1);

      //! Get the weight
      const float wi = float(i) * di - float(i0);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      const __m256 zwi = _mm256_set1_ps(wi);
      const __m256 z1  = _mm256_set1_ps(1.f);
      for (; j < width - 8; j += 8) {
        const size_t j00 = (size_t) floor(float(j) * dj);
        const size_t j01 = std::min(m_width - 1, j00 + 1);
        const size_t j10 = (size_t) floor(float(j + 1) * dj);
        const size_t j11 = std::min(m_width - 1, j10 + 1);
        const size_t j20 = (size_t) floor(float(j + 2) * dj);
        const size_t j21 = std::min(m_width - 1, j20 + 1);
        const size_t j30 = (size_t) floor(float(j + 3) * dj);
        const size_t j31 = std::min(m_width - 1, j30 + 1);
        const size_t j40 = (size_t) floor(float(j + 4) * dj);
        const size_t j41 = std::min(m_width - 1, j40 + 1);
        const size_t j50 = (size_t) floor(float(j + 5) * dj);
        const size_t j51 = std::min(m_width - 1, j50 + 1);
        const size_t j60 = (size_t) floor(float(j + 6) * dj);
        const size_t j61 = std::min(m_width - 1, j60 + 1);
        const size_t j70 = (size_t) floor(float(j + 7) * dj);
        const size_t j71 = std::min(m_width - 1, j70 + 1);
        const __m256 zwj = _mm256_set_ps(float(j + 7) * dj - float(j70),
                                         float(j + 6) * dj - float(j60),
                                         float(j + 5) * dj - float(j50),
                                         float(j + 4) * dj - float(j40),
                                         float(j + 3) * dj - float(j30),
                                         float(j + 2) * dj - float(j20),
                                         float(j + 1) * dj - float(j10),
                                         float(j + 0) * dj - float(j00));
        const __m256 a = _mm256_set_ps(iI0[j70], iI0[j60], iI0[j50], iI0[j40], iI0[j30], iI0[j20], iI0[j10], iI0[j00]);
        const __m256 b = _mm256_set_ps(iI1[j70], iI1[j60], iI1[j50], iI1[j40], iI1[j30], iI1[j20], iI1[j10], iI1[j00]);
        const __m256 c = _mm256_set_ps(iI0[j71], iI0[j61], iI0[j51], iI0[j41], iI0[j31], iI0[j21], iI0[j11], iI0[j01]);
        const __m256 d = _mm256_set_ps(iI1[j71], iI1[j61], iI1[j51], iI1[j41], iI1[j31], iI1[j21], iI1[j11], iI1[j01]);
        _mm256_storeu_ps(oI + j, (z1 - zwj) * ((z1 - zwi) * a + zwi * b) +
                                       zwj  * ((z1 - zwi) * c + zwi * d));
      }
#endif
#if defined(__SSE__)
      //! SSE version
      const __m128 xwi = _mm_set1_ps(wi);
      const __m128 x1  = _mm_set1_ps(1.f);
      for (; j < width - 4; j += 4) {
        const size_t j00 = (size_t) floor(float(j) * dj);
        const size_t j01 = std::min(m_width - 1, j00 + 1);
        const size_t j10 = (size_t) floor(float(j + 1) * dj);
        const size_t j11 = std::min(m_width - 1, j10 + 1);
        const size_t j20 = (size_t) floor(float(j + 2) * dj);
        const size_t j21 = std::min(m_width - 1, j20 + 1);
        const size_t j30 = (size_t) floor(float(j + 3) * dj);
        const size_t j31 = std::min(m_width - 1, j30 + 1);
        const __m128 xwj = _mm_set_ps(float(j + 3) * dj - float(j30),
                                      float(j + 2) * dj - float(j20),
                                      float(j + 1) * dj - float(j10),
                                      float(j + 0) * dj - float(j00));
        const __m128 a = _mm_set_ps(iI0[j30], iI0[j20], iI0[j10], iI0[j00]);
        const __m128 b = _mm_set_ps(iI1[j30], iI1[j20], iI1[j10], iI1[j00]);
        const __m128 c = _mm_set_ps(iI0[j31], iI0[j21], iI0[j11], iI0[j01]);
        const __m128 d = _mm_set_ps(iI1[j31], iI1[j21], iI1[j11], iI1[j01]);
        _mm_storeu_ps(oI + j, (x1 - xwj) * ((x1 - xwi) * a + xwi * b) +
                                    xwj  * ((x1 - xwi) * c + xwi * d));
      }
#endif
      //! Normal version
      for (; j < width; j++) {
        //! Get the current column position
        const size_t j0 = (size_t) floor(float(j) * dj);
        const size_t j1 = std::min(j0 + 1, m_width - 1);

        //! Get the weight
        const float wj = float(j) * dj - float(j0);

        //! Bilinear interpolation
        oI[j] = (1.f - wj) * ((1.f - wi) * iI0[j0] + wi * iI1[j0]) +
                       wj  * ((1.f - wi) * iI0[j1] + wi * iI1[j1]);
      }
    }
  }
}


//! Get the topography in false color of a landscape, according to the water
//! level and the ocean level.
void Image::getTopo(
  Image& o_im,
  const Image& i_imWater,
  const float p_oceanLvl,
  const float p_min,
  const float p_max,
  const char* p_folder,
  const float p_waterThreshold) const {

  //! Check that the folder is not null
  if (p_folder == NULL) {
    cout << "Please provide a non-null folder path." << endl;
    exit(EXIT_FAILURE);
  }

  //! For convenience
#if defined(__AVX__)
  const __m256 z1 = _mm256_set1_ps(0.1);
  const __m256 z2 = _mm256_set1_ps(0.2);
  const __m256 z8 = _mm256_set1_ps(0.8);
#endif
#if defined(__SSE__)
  const __m128 x1 = _mm_set1_ps(0.1);
  const __m128 x2 = _mm_set1_ps(0.2);
  const __m128 x8 = _mm_set1_ps(0.8);
#endif // _SSE__

  //! Initialization of the images
  Image imColor, imShadow, imCauchy, imTmp, imLaplacien;
  o_im.init(m_width, m_height, 3, m_border);

  //! Get the colorized map
  this->getColorized(imColor, i_imWater, p_oceanLvl, p_min, p_max, p_folder,
    p_waterThreshold);

  //! Get the shadow image (80% of the shading)
  this->getShadow(imShadow);
  imShadow.blurCauchy(imTmp, 0.5);
  imTmp.stretch(imShadow, 0.01, 0, 0.9);

  //! Get the laplacian (20% of the shading)
  this->getLaplacien(imLaplacien, -1.f);
  imLaplacien.blurCauchy(imTmp, 0.5);
  imTmp.stretch(imLaplacien, 0.01, 0, 0.9);

  //! Get the blend between the colorized map and the shadow
  for (size_t c = 0; c < 3; c++) {
    for (size_t i = 0; i < m_height; i++) {
      const float* iC = imColor    .getPtr(c, i);
      const float* iS = imShadow   .getPtr(0, i);
      const float* iL = imLaplacien.getPtr(0, i);
      float* oI = o_im.getPtr(c, i);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        _mm256_storeu_ps(oI + j, _mm256_loadu_ps(iC + j) * (z1 + z8 *
                      _mm256_loadu_ps(iS + j) + z2 * _mm256_loadu_ps(iL + j)));
      }
#endif
#if defined(__SSE__)
      //! SSE version
      for (; j < m_width - 4; j += 4) {
        _mm_storeu_ps(oI + j, _mm_loadu_ps(iC + j) * (x1 + x8 *
                            _mm_loadu_ps(iS + j) + x2 * _mm_loadu_ps(iL + j)));
      }
#endif
      //! Normal version
      for (; j < m_width; j++) {
        oI[j] = iC[j] * (0.1 + 0.8 * iS[j] + 0.2 * iL[j]);
      }
    }
  }
}


//! Get a colorized version of a landscape image.
void Image::getColorized(
  Image& o_im,
  const Image& i_imWater,
  const float p_oceanLvl,
  const float p_min,
  const float p_max,
  const char* p_folder,
  const float p_waterThreshold) const {

  //! Check size
  if (i_imWater.height() != m_height || i_imWater.width() != m_width) {
    cout << "Size of water and im should be the same. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Open the file that contains the palette nodes for water
  std::vector<float> nTmp;
  readPalette((string(p_folder) + "topographic_water.gpl").c_str(), nTmp);
  const size_t nbWater = nTmp.size() / 3;

  //! Reverse the water color palette
  std::vector<float> nWater(nTmp.size());
  for (size_t n = 0; n < nbWater; n++) {
    nWater[3 * n + 0] = nTmp[3 * (nbWater - n - 1) + 0];
    nWater[3 * n + 1] = nTmp[3 * (nbWater - n - 1) + 1];
    nWater[3 * n + 2] = nTmp[3 * (nbWater - n - 1) + 2];
  }

  //! Open the file that contains the palette nodes for ground
  std::vector<float> nGround;
  readPalette((string(p_folder) + "topographic_ground.gpl").c_str(), nGround);
  const size_t nbGround = nGround.size() / 3;

  //! Get the mean, min and max value for the water
  float minVal, maxVal, mean;
  i_imWater.getStats(minVal, maxVal, mean);

  //! Initialization of the output image
  o_im.init(m_width, m_height, 3, m_border);
  const float normWater = (maxVal > minVal ? float(nbWater - 1) / (maxVal - minVal) : 1);
  const float normGround = float(nbGround - 1) / (p_max - p_min);

  //! Apply the colorization
#ifdef _OPENMP
#pragma omp parallel
{
  const size_t tid = omp_get_thread_num();
#else
  const size_t tid = 0;
#endif // _OPENMP

  const size_t hBegin = m_heights[tid];
  const size_t hEnd   = m_heights[tid + 1];

  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iI = this->getPtr(0, i);
    const float* iW = i_imWater.getPtr(0, i);
    float* oR = o_im.getPtr(0, i);
    float* oG = o_im.getPtr(1, i);
    float* oB = o_im.getPtr(2, i);

    for (size_t j = 0; j < m_width; j++) {

      //! Initialization
      const float valLand  = p_max * pow((iI[j] - p_oceanLvl) / p_max, 0.25);
      const float valWater = iW[j];

      //! If the landscape is above ocean level or the water is high enough
      if (iI[j] <= p_oceanLvl || valWater > p_waterThreshold * mean) {
        const size_t k0 = (size_t) floor((valWater - minVal) * normWater);
        const size_t k1 = std::min(nbWater - 1, k0 + 1);
        const float   x = (valWater - minVal) * normWater - float(k0);
        oR[j] = nWater[3 * k0 + 0] * (1.f - x) + nWater[3 * k1 + 0] * x;
        oG[j] = nWater[3 * k0 + 1] * (1.f - x) + nWater[3 * k1 + 1] * x;
        oB[j] = nWater[3 * k0 + 2] * (1.f - x) + nWater[3 * k1 + 2] * x;
      }

      //! Check for boundary conditions
      else if (valLand <= p_min) {
        oR[j] = nGround[0];
        oG[j] = nGround[1];
        oB[j] = nGround[2];
      }
      else if (valLand >= p_max) {
        oR[j] = nGround[3 * (nbGround - 1) + 0];
        oG[j] = nGround[3 * (nbGround - 1) + 1];
        oB[j] = nGround[3 * (nbGround - 1) + 2];
      }

      //! Generic value -> interpolation between the nodes
      else {
        const size_t k0 = (size_t) floor((valLand - p_min) * normGround);
        const size_t k1 = std::min(nbGround - 1, k0 + 1);
        const float   x = (valLand - p_min) * normGround - float(k0);
        oR[j] = nGround[3 * k0 + 0] * (1.f - x) + nGround[3 * k1 + 0] * x;
        oG[j] = nGround[3 * k0 + 1] * (1.f - x) + nGround[3 * k1 + 1] * x;
        oB[j] = nGround[3 * k0 + 2] * (1.f - x) + nGround[3 * k1 + 2] * x;
      }
    }
  }
#ifdef _OPENMP
}
#endif // _OPENMP
}


//! Get Shadow image.
void Image::getShadow(
  Image& o_im) const {

  //! Initialization
  o_im.init(m_width, m_height, m_channels, m_border);

#ifdef _OPENMP
#pragma omp parallel
{
  const size_t tid = omp_get_thread_num();
#else
  const size_t tid = 0;
#endif // _OPENMP

  //! For convenience
  const int hBegin = m_heights[tid];
  const int hEnd   = m_heights[tid + 1];

  //! Process each pixel
  for (size_t c = 0; c < m_channels; c++) {
    for (int i = hBegin; i < hEnd; i++) {
      const int it = m_border > 0 ? i - 1 : std::max(0, i - 1);
      const int ib = m_border > 0 ? i + 1 : std::min((int) m_height - 1, i + 1);
      const float* iT = this->getPtr(0, it);
      const float* iC = this->getPtr(0, i);
      const float* iB = this->getPtr(0, ib);
      float* oI = o_im.getPtr(0, i);

      for (int j = 0; j < (int) m_width; j++) {
        const int jl = m_border > 0 ? j - 1 : std::max(j - 1, 0);
        const int jr = m_border > 0 ? j + 1 : std::min(j + 1, (int) m_width - 1);

        //! Compute the vector for the gradient in X
        const float vdx[3] = {1, 0, (iC[jr] - iC[jl]) * 0.5f};

        //! Compute the vector for the gradient in Y
        const float vdy[3] = {0, 1, (iB[j] - iT[j]) * 0.5f};

        //! Vector for the sun (lighting)
        const float sun[3] = {-1, -1, 1};

        //! Get the normal of the surface
        const float nor[3] = {vdx[1] * vdy[2] - vdx[2] * vdy[1],
                          vdx[2] * vdy[0] - vdx[0] * vdy[2],
                          vdx[0] * vdy[1] - vdx[1] * vdy[0]};

        //! Compute the scalar product to get the shadow
        oI[j] = nor[0] * sun[0] + nor[1] * sun[1] + nor[2] * sun[2];
      }
    }
  }

#ifdef _OPENMP
}
#endif // _OPENMP
}


//! Stretch an image based on quantile.
void Image::stretch(
  Image& o_im,
  const float p_quantile,
  const float p_min,
  const float p_max) const {

  //! Initialization
  o_im.init(m_width, m_height, m_channels, m_border);

  //! Channel independency
  for (size_t c = 0; c < m_channels; c++) {

    //! Cast the values into a vector
    std::vector<float> values(m_width * m_height);
    for (size_t i = 0; i < m_height; i++) {
      const float* iI = this->getPtr(c, i);
      float* oV = &values[i * m_width];

      for (size_t j = 0; j < m_width; j++) {
        oV[j] = iI[j];
      }
    }

    //! Get the index of first and last quantile wanted
    const size_t nQuantile = (size_t) (p_quantile * float(m_width * m_height));

    //! Partial sort of the values (to avoid to sort the whole array)
    std::partial_sort(values.begin(), values.begin() + nQuantile + 1, values.end());
    const float qMin = values[nQuantile];
    std::partial_sort(values.begin(), values.begin() + nQuantile + 1, values.end(), std::greater<float>());
    const float qMax = values[nQuantile];

    //! For convenience
    const float norm = p_max / (qMax - qMin);

    //! Stretch the values
    for (size_t i = 0; i < m_height; i++) {
      const float* iI = this->getPtr(c, i);
      float* oI = o_im.getPtr(c, i);

      for (size_t j = 0; j < m_width; j++) {
        const float value = iI[j];

        //! Saturation
        if (value <= qMin) {
          oI[j] = p_min;
        }
        else if (value >= qMax) {
          oI[j] = p_max;
        }

        //! Stretching
        else {
          oI[j] = p_min + norm * (value - qMin);
        }
      }
    }
  }
}


//! Blur an image with a Cauchy kernel by using the fftw.
void Image::blurCauchy(
  Image& o_im,
  const float p_sigma) const {

  //! Multithread for fftw initialization
#ifdef _OPENMP
  fftwf_init_threads();
  fftwf_plan_with_nthreads(m_nbThreads);
#endif // _OPENMP

  //! Size initialization
  o_im.init(m_width, m_height, m_channels, m_border);
  const int h = m_height;
  const int w = m_width;
  const float normF = 1.f / float(w * h);

  //! Fill the kernel over the image
  float* kernel = (float*) malloc(w * h * sizeof(float));
  float normK = 0.f;
  const float sigma2 = p_sigma * p_sigma;
  for (int i = 0, k = 0; i < h; i++) {
    for (int j = 0; j < w; j++, k++) {

      //! Get the coordinates
      const float x = j < w / 2 ? j : j - w;
      const float y = i < h / 2 ? i : i - h;

      //! Get the kernel value
      const float value = sigma2 / (sigma2 + x * x + y * y);
      kernel[k] = value;
      normK += value;
    }
  }
  normK = 1.f / normK;

  //! Normalize the kernel
  fftwf_complex *vec = (fftwf_complex*) fftwf_malloc(w * h * sizeof(fftwf_complex));
  for (int k = 0; k < w * h; k++) {
    vec[k] = float(kernel[k] * normK);
  }
  memfree(kernel);

  //! DFT of the kernel
  fftwf_complex *dftK = (fftwf_complex*) fftwf_malloc(w * h * sizeof(fftwf_complex));
  fftwf_plan planK = fftwf_plan_dft_2d(h, w, vec, dftK, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(planK);
  fftwf_destroy_plan(planK);

  //! Channel independency
    for (size_t c = 0; c < m_channels; c++) {

    //! DFT of the image
    fftwf_complex *dftI = (fftwf_complex*) fftwf_malloc(w * h * sizeof(fftwf_complex));
    for (int i = 0; i < h; i++) {
      const float* iI = this->getPtr(c, i);

      for (int j = 0; j < w; j++) {
        vec[i * w + j] = (float) iI[j];
      }
    }
    fftwf_plan planI = fftwf_plan_dft_2d(h, w, vec, dftI, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(planI);
    fftwf_destroy_plan(planI);

    //! Point-wise multiplication
    for (int k = 0; k < w * h; k++) {
      dftI[k] *= dftK[k];
    }

    //! DFT inverse of the image
    fftwf_plan plan = fftwf_plan_dft_2d(h, w, dftI, vec, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_free(dftI);

    //! Normalization to get the final result
    for (int i = 0; i < h; i++) {
      float* oI = o_im.getPtr(c, i);

      for (int j = 0; j < w; j++) {
        oI[j] = (float) crealf(vec[i * w + j] * normF);
      }
    }
  }

  //! Cleanup
  fftwf_free(vec);
  fftwf_free(dftK);
#ifdef _OPENMP
  fftwf_cleanup_threads();
#else
  fftwf_cleanup();
#endif
}


//! Get stats of an image (min, max, and mean values).
void Image::getStats(
  float& o_min,
  float& o_max,
  float& o_mean) const {

  //! Initialization
  o_min  = this->getPtr(0, 0)[0];
  o_max  = o_min;
  double mean = 0;

#if defined(__AVX__)
  __m256 zMin  = _mm256_loadu_ps(this->getPtr(0, 0));
  __m256 zMax  = zMin;
  __m256 zMean = _mm256_setzero_ps();
#endif
#if defined(__SSE__)
  __m128 xMin  = _mm_loadu_ps(this->getPtr(0, 0));
  __m128 xMax  = xMin;
  __m128 xMean = _mm_setzero_ps();
#endif

  //! Update the values
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < m_height; i++) {
      const float* iI = this->getPtr(c, i);
      size_t j = 0;

#if defined(__AVX__)
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        const __m256 zVal = _mm256_loadu_ps(iI + j);
        zMin   = _mm256_min_ps(zMin, zVal);
        zMax   = _mm256_max_ps(zMax, zVal);
        zMean += zVal;
      }
#endif
#if defined(__SSE__)
      //! SSE version
      for (; j < m_width - 4; j += 4) {
        const __m128 xVal = _mm_loadu_ps(iI + j);
        xMin   = _mm_min_ps(xMin, xVal);
        xMax   = _mm_max_ps(xMax, xVal);
        xMean += xVal;
      }
#endif
      //! Normal version
      for (; j < m_width; j++) {
        const float value = iI[j];
        o_min = std::min(o_min, value);
        o_max = std::max(o_max, value);
        mean += (double) value;
      }
    }
  }

  //! Get the value
#if defined(__AVX__)
  float zminVal[8], zmaxVal[8], zmeanVal[8];
  _mm256_storeu_ps(zminVal , zMin );
  _mm256_storeu_ps(zmaxVal , zMax );
  _mm256_storeu_ps(zmeanVal, zMean);
  for (size_t n = 0; n < 8; n++) {
    o_min = std::min(o_min, zminVal[n]);
    o_max = std::max(o_max, zmaxVal[n]);
    mean += (double) zmeanVal[n];
  }
#endif
#if defined(__SSE__)
  float xminVal[4], xmaxVal[4], xmeanVal[4];
  _mm_storeu_ps(xminVal , xMin );
  _mm_storeu_ps(xmaxVal , xMax );
  _mm_storeu_ps(xmeanVal, xMean);
  for (size_t n = 0; n < 4; n++) {
    o_min = std::min(o_min, xminVal[n]);
    o_max = std::max(o_max, xmaxVal[n]);
    mean += (double) xmeanVal[n];
  }
#endif // __SSE__

  // [ToDo] What happens if neither __AVX__ nor __SSE__ are defined?

  //! Normalize the sum
  o_mean = (float) (mean / double(m_width * m_height * m_channels));
}


//! Get the HSV visualization of an image
void Image::getHSV(
  Image &o_im,
  const float p_min,
  const float p_max,
  const int p_tid) const {

  //! Check size
  if (p_tid >=0) {
    if (o_im.width() != m_width || o_im.height() != m_height || o_im.channels() != 3) {
      cout << "Image::getHSV : the output image must be initialized." << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {
    o_im.init(m_width, m_height, 3);
  }

  //! Get the min, max values
  float minVal = p_min, maxVal = p_max;
  if (p_tid < 0) {
    float mean;
    this->getStats(minVal, maxVal, mean);
  }

  //! For convenience
  const size_t hBegin = p_tid < 0 ? 0        : m_heights[p_tid];
  const size_t hEnd   = p_tid < 0 ? m_height : m_heights[p_tid + 1];
  const float norm    = 1.f / std::max(maxVal - minVal, 1.f);

  //! LUT to convert HSV to RGB
  float* lut = (float*) memalloc(16, 3 * 361 * sizeof(float));
  for (int k = 0; k < 360; k++) {
    const float C = 255.f;
    const float X = C * (1 - fabs(float(k%120) / 60.f - 1));
    if (k < 60) {
      lut[k * 3 + 0] = 0; lut[k * 3 + 1] = X; lut[k * 3 + 2] = C;
    }
    else if (k < 120) {
      lut[k * 3 + 0] = 0; lut[k * 3 + 1] = C; lut[k * 3 + 2] = X;
    }
    else if (k < 180) {
      lut[k * 3 + 0] = X; lut[k * 3 + 1] = C; lut[k * 3 + 2] = 0;
    }
    else if (k < 240) {
      lut[k * 3 + 0] = C; lut[k * 3 + 1] = X; lut[k * 3 + 2] = 0;
    }
    else if (k < 300) {
      lut[k * 3 + 0] = C; lut[k * 3 + 1] = 0; lut[k * 3 + 2] = X;
    }
    else {
      lut[k * 3 + 0] = X; lut[k * 3 + 1] = 0; lut[k * 3 + 2] = C;
    }
  }

  //! To avoid problem
  lut[360 * 3 + 0] = 4.25; lut[360 * 3 + 1] = 0; lut[360 * 3 + 2] = 255;

  //! Convert HSV to RGB (S = V = 1)
  for (size_t i = hBegin; i < hEnd; i++) {
    const float* iI = this->getPtr(0, i);
    float* oR = o_im.getPtr(0, i);
    float* oG = o_im.getPtr(1, i);
    float* oB = o_im.getPtr(2, i);

    for (size_t j = 0; j < m_width; j++) {
      const float x = 120.f + 239.f * norm * (iI[j] - minVal);
      const int   a  = std::max(std::min((int) x, 359), 0);
      const int   b  = a + 1;
      const float cb = lut[a * 3 + 0];
      const float db = lut[b * 3 + 0];
      const float cg = lut[a * 3 + 1];
      const float dg = lut[b * 3 + 1];
      const float cr = lut[a * 3 + 2];
      const float dr = lut[b * 3 + 2];
      oR[j] = x * (dr - cr) + cr * b - dr * a;
      oG[j] = x * (dg - cg) + cg * b - dg * a;
      oB[j] = x * (db - cb) + cb * b - db * a;
    }
  }

  //! Release memory
  memfree(lut);
}


//! Check if the current pixel as a lower neighbor.
bool Image::hasLowerNeighbor(
  const int p_i,
  const int p_j,
  const size_t p_c,
  const bool p_useDiag) const {

  //! For convenience
  const float* iIT = this->getPtr(p_c, p_i - 1);
  const float* iIC = this->getPtr(p_c, p_i);
  const float* iIB = this->getPtr(p_c, p_i + 1);
  const float value = iIC[p_j];

  //! Left neighbor
  if (iIC[p_j - 1] < value) {
    return true;
  }

  //! Right neighbor
  if (iIC[p_j + 1] < value) {
    return true;
  }

  //! Top neighbor
  if (iIT[p_j] < value) {
    return true;
  }

  //! Bottom neighbor
  if (iIB[p_j] < value) {
    return true;
  }

  //! Diagonal elements
  if (p_useDiag) {
    //! Top left
    if (iIT[p_j - 1] < value) {
      return true;
    }

    //! Top right
    if (iIT[p_j + 1] < value) {
      return true;
    }

    //! Bottom left
    if (iIB[p_j - 1] < value) {
      return true;
    }

    //! Bottom right
    if (iIB[p_j + 1] < value) {
      return true;
    }
  }
  return false;
}


//! Get the laplacien of the image multiplied by a factor.
void Image::getLaplacien(
  Image& o_im,
  const float p_factor) const {

  //! Initialize the size of the output image
  if (o_im.width() != m_width || o_im.height() != m_height ||
      o_im.channels() != m_channels) {
    o_im.init(m_width, m_height, m_channels, m_border);
  }

#ifdef _OPENMP
#pragma omp parallel
{
  const size_t tid = omp_get_thread_num();
#else
  const size_t tid = 0;
#endif // _OPENMP

  //! For convenience
  const int hBegin = m_heights[tid];
  const int hEnd   = m_heights[tid + 1];

  //! Process each pixel
  for (size_t c = 0; c < m_channels; c++) {
    for (int i = hBegin; i < hEnd; i++) {
      const int it = m_border > 0 ? i - 1 : std::max(0, i - 1);
      const int ib = m_border > 0 ? i + 1 : std::min((int) m_height - 1, i + 1);
      const float* iT = this->getPtr(c, it);
      const float* iC = this->getPtr(c, i);
      const float* iB = this->getPtr(c, ib);
      float* oI = o_im.getPtr(c, i);

      for (int j = 0; j < (int) m_width; j++) {
        const int jl = m_border > 0 ? j - 1 : std::max(j - 1, 0);
        const int jr = m_border > 0 ? j + 1 : std::min(j + 1, (int) m_width - 1);

        oI[j] = p_factor * (iT[j] + iB[j] + iC[jl] + iC[jr] - 4 * iC[j]);
      }
    }
  }

#ifdef _OPENMP
}
#endif // _OPENMP
}


//! Make a cone image.
void Image::makeCone(
  const bool p_descend,
  const bool p_half) {

  //! Size initialization
  const size_t w  = 512;
  const size_t h  = w / (p_half ? 2 : 1);
  const size_t h2 = p_half ? h : h / 2;
  const size_t w2 = w / 2;
  this->init(w, h, 1);

  //! For convenience
  const float norm = 255.f / sqrt(float((w2-1) * (w2-1) + (h2-1) * (h2-1)));

  //! First half of the cone
  for (size_t i = 0; i < h2; i++) {
    float* oI = this->getPtr(0, h2 - i - 1);

    //! First quarter of the cone
    for (size_t j = 0; j < w2; j++) {

      //! Compute the radius
      const float r = sqrt(float(j * j) + float(i * i));

      //! Get the value of the cone
      oI[w2 + j] = (p_descend ? r * norm : 255.f - r * norm);

      //! Apply the symmetry
      oI[w2 - j - 1] = oI[w2 + j];
    }
  }

  //! Second half of the cone
  if (!p_half) {
    for (size_t i = 0; i < h2; i++) {
      const float* iI = this->getPtr(0, i);
      float* oI = this->getPtr(0, h - i - 1);

      for(size_t j = 0; j < w; j++) {
        oI[j] = iI[j];
      }
    }
  }
}


//! Make a plan
void Image::makePlan(
  const size_t p_theta) {

  //! For convenience
  const double theta = (p_theta % 90 <= 45 ? p_theta % 90 : 90 - (p_theta % 90)) * M_PI / 180;

  //! Size initialization
  const int w  = 512;
  const int h  = w;
  const int w2 = w / 2;
  const int h2 = h / 2;
  this->init(w, h, 1);

  //! Get the coordinate of the point
  const double ax = 0;
  const double ay = 0;
  const double az = 128;

  //! Get the coordinates of two vectors
  const double ux = h2;
  const double uy = h2 * tan(theta);
  const double uz = 0;
  const double vx = -w2 * tan(theta);
  const double vy = w2;
  const double vz = 128;

  //! Get the cartesian equation of the plan
  const double a = uy * vz - uz * vy;
  const double b = uz * vx - ux * vz;
  const double c = ux * vy - uy * vx;
  const double d = -(a * ax + b * ay + c * az);

  //! Get the plan
  for (int i = 0; i < h; i++) {
    float* oI = this->getPtr(0, i);

    for (int j = 0; j < w; j++) {

      //! Get the coordinates of the current point
      const double x = h2 - i;
      const double y = j - w2;

      //! Get the value
      oI[j] = -(a * x + b * y + d) / c;
    }
  }
}


//! For debug only. Print an image with its border.
void Image::print() const {

  for (size_t c = 0; c < m_channels; c++) {
    for (int i = -int(m_border); i < int(m_height + m_border); i++) {
      const float* iI = this->getPtr(c, i);
      if (i == 0 || i == int(m_height)) {
        cout << endl;
      }

      for (int j = -int(m_border); j < int(m_width + m_border); j++) {
        if (j == 0 || j == int(m_width)) {
          cout << " ";
        }
        cout << iI[j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;
}












