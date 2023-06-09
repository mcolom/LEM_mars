/**
 * @file Utilites.cpp
 *
 * @brief Miscellaneous functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <png.h>
#include <tiffio.h>
#include <omp.h>


//! Local includes
#include "Utilities.h"


using namespace std;


//! Compute an amount of time elapsed.
void getTime(
  struct timespec &io_start,
  const char* p_name) {

  //! Check the current time
  struct timespec finish;
  clock_gettime(CLOCK_MONOTONIC, &finish);

  //! Compute the elapsed time
  double elapsed = (finish.tv_sec - io_start.tv_sec) * 1000;
  elapsed += (finish.tv_nsec - io_start.tv_nsec) / 1000000.0;

  //! Print the result
  cout << p_name << ": (ms) = " << elapsed << endl;

  //! Start a new timer
  clock_gettime(CLOCK_MONOTONIC, &io_start);
}


//! Open a file with the given name in the given path, write a string then
//! close the file. Will exit the program if the file can't be opened.
void writeValue(
  const std::string &i_fileName,
  const std::string &i_pathName,
  const std::string &p_value) {

  //! Open the file
  FILE *file = fopen((i_pathName + i_fileName + ".txt").c_str(), "w");

  //! Check that the file has been opened correctly
  if (file != NULL) {

    //! Write the value in it
    fputs((const char*) p_value.c_str(), file);

    //! Close the file
    fclose(file);
  }
  else {
    cout << "Error: the file " << i_fileName << " couldn't be opened in the ";
    cout << "folder " << i_pathName << "." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Surcharged function to write a float value.
void writeValue(
  const std::string &i_fileName,
  const std::string &i_pathName,
  const float p_value) {

  //! Cast the float value into string
  std::ostringstream ss;
  ss << p_value;

  //! Call the generic function
  writeValue(i_fileName, i_pathName, ss.str() + "\n");
}


//! Write a curve.
void writeCurve(
  const std::string &i_fileName,
  const std::string &i_pathName,
  std::vector<float> const& i_curve) {

  //! Open the file
  ofstream myfile;
  myfile.open((i_pathName + i_fileName + ".txt").c_str());

  //! Write the values
  for (size_t n = 0; n < i_curve.size(); n++) {
    myfile << i_curve[n] << endl;
  }

  //! Close the file
  myfile.close();
}


//! Write the estimation of automatic selection of parameters.
void writeParameters(
  const std::string &i_fileName,
  const std::string &i_pathName,
  const std::vector<float> &i_valNorm,
  const std::vector<float> &i_valCreep,
  const std::vector<float> &i_valSedimentation,
  const std::vector<float> &i_valExpM,
  const std::vector<float> &i_valExpN) {

  //! Check the size of the arrays
  const size_t N = i_valNorm.size();
  if (i_valCreep        .size() != N ||
      i_valSedimentation.size() != N ||
      i_valExpM         .size() != N ||
      i_valExpN         .size() != N) {
    cout << "writeParameters: the size of the arrays must identical." << endl;
    exit(EXIT_FAILURE);
  }

  //! Open the file
  ofstream file;
  file.open((i_pathName + i_fileName + ".txt").c_str());

  //! Check that the file has been opened correctly
  if (file.is_open()) {

    //! Write the number of parameters
    file << "5\n";

    //! Write the name of the tested parameters
    file << "Norm\t" << "Creep\t" << "Sedimentation\t" << "ExpM\t" << "ExpN\n";

    //! Write the value in it
    for (size_t n = 0; n < N; n++) {
      file << i_valNorm[n] << "\t" << i_valCreep[n] << "\t";
      file << i_valSedimentation[n] << "\t" << i_valExpM[n] << "\t";
      file << i_valExpN[n] << "\n";
    }

    //! Close the file
    file.close();
  }
  else {
    cout << "Error: the file " << i_fileName << " couldn't be opened in the ";
    cout << "folder " << i_pathName << "." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Write the automatic selection of the best parameters with the range of the
//! tested parameters.
void writeBestParameters(
  const std::string &i_inpPath,
  const std::string &i_outPath,
  const float p_bestNorm,
  const std::vector<float> p_bestValues,
  const std::vector<float> p_rangeMin,
  const std::vector<float> p_rangeMax,
  const std::vector<std::string> p_paramNames) {

  //! Check the size of the arrays
  const size_t N = p_bestValues.size();
  if (p_rangeMin.size() != N ||
      p_rangeMax.size() != N ||
      p_paramNames.size() != N) {
    cout << "writeBestParameters : all input arrays must have the same size.\n";
    exit(EXIT_FAILURE);
  }

  //! Extract the name of the input image
  const size_t pos1 = i_inpPath.find_last_of("/");
  const size_t pos2 = i_inpPath.find_last_of(".");
  const std::string imName = i_inpPath.substr(pos1 + 1, pos2 - pos1 - 1);
  const std::string fileName = i_outPath + imName + "_bestParams.txt";

  //! Open the file
  ofstream file;
  file.open(fileName.c_str());

  //! Check that the file has been opened correctly
  if (file.is_open()) {

    //! Norm reached
    file << "Norm reached: " << p_bestNorm << endl << endl;

    //! Write the set of parameters
    file << "Best parameters:" << endl;
    for (size_t k = 0; k < N; k++) {
      file << "\t-\t" << p_paramNames[k] << ":\t" << p_bestValues[k] << "\t";
      file << "[" << p_rangeMin[k] << ", " << p_rangeMax[k] << "]\n";
    }

    //! Close the file
    file.close();
  }
  else {
    cout << "Error: the file " << fileName << " couldn't be opened." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Get left and right values from the center one: 0.5 -> 0.4, 0.5, 0.6.
void getValues(
  const float p_value,
  std::vector<float>& o_values) {

  float val = p_value;
  float p = 1.0;
  while (val / 10.f - floor(val / 10.f) == 0.f) {
    val /= 10.f;
    p /= 10.f;
  }
  while (val > floor(val)) {
    val *= 10.f;
    p *= 10.f;
  }

  o_values.clear();

  if (val > 0.f) {
    o_values.push_back(val > 1.f ? (val - 1) / p : 0.9f / p);
  }
  o_values.push_back(val / p);
  o_values.push_back((val + 1) / p);
}


//! Return the optimal cut of the lines of an image according to the number of
//! threads available.
void initializeHeights(
  const size_t i_height,
  size_t* o_heights,
  const size_t p_nbThreads) {

  //! Get the initial number of lines for each thread to deal with.
  const size_t nb = p_nbThreads;
  const size_t h = ceil(float(i_height) / float(nb));

  //! Initialize the number of lines
  size_t* lines = (size_t*) malloc(nb * sizeof(size_t));
  for (size_t n = 0; n < nb; n++) {
    lines[n] = h;
  }

  //! Remove one by one the lines that are not needed
  int nCur = nb - 1;
  int k = int(h * nb) - int(i_height);
  while (k > 0) {
    lines[nCur]--;
    if (nCur > 0) {
      nCur--;
    }
    else {
      nCur = nb - 1;
    }
    k--;
  }

  //! Get the lines position
  o_heights[0] = 0;
  for (size_t n = 0; n < nb; n++) {
    o_heights[n + 1] = o_heights[n] + lines[n];
  }

  //! Release memory
  free(lines);
}


//! Get the size of a png image without reading it entirely.
void getSize(
  const char* i_name,
  size_t &o_width,
  size_t &o_height) {

  //! Get the extension of the image name.
  string ext = i_name;
  const int pos = ext.find_last_of(".");
  ext = ext.substr(pos + 1, ext.size());

  if (ext == "png") {
    const size_t png_sig_len = 4;

    png_byte png_sig[png_sig_len];
    png_structp png_ptr;
    png_infop info_ptr;
    FILE *volatile fp = NULL;
    int transform = PNG_TRANSFORM_IDENTITY;

    //! open the PNG input file
    if (0 == strcmp(i_name, "-")) {
      fp = stdin;
    }
    else if (NULL == (fp = fopen(i_name, "rb"))) {
      cout << "can't open the file " << i_name << endl;
      exit(EXIT_FAILURE);
    }

    //! read in some of the signature bytes and check this signature
    if ((png_sig_len != fread(png_sig, 1, png_sig_len, fp))
        || 0 != png_sig_cmp(png_sig, (png_size_t) 0, png_sig_len)) {
      cout << "bad signature." << endl;
      exit(EXIT_FAILURE);
    }

    //! create and initialize the png_struct with the default stderr and error
    //! handling
    if (NULL == (png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                                  NULL, NULL, NULL))) {
      cout << "Can't initialize the png_struct." << endl;
      exit(EXIT_FAILURE);
    }

    //! allocate/initialize the memory for image information
    if (NULL == (info_ptr = png_create_info_struct(png_ptr))) {
      cout << "Can't initialize the memory." << endl;
      exit(EXIT_FAILURE);
    }

    //! set error handling
    if (0 != setjmp(png_jmpbuf(png_ptr))) {
      //! if we get here, we had a problem reading the file.
      //! free all of the memory associated with the png_ptr and info_ptr.
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      if (NULL != fp && stdin != fp) {
        (void) fclose(fp);
      }
      cout << "Can't read the file." << endl;
      exit(EXIT_FAILURE);
    }

    //! set up the input control using standard C streams
    png_init_io(png_ptr, fp);

    //! let libpng know that some bytes have been read
    png_set_sig_bytes(png_ptr, png_sig_len);

    /*
     * set the read filter transforms, to get 8bit RGB whatever the
     * original file may contain:
     * PNG_TRANSFORM_STRIP_16      strip 16-bit samples to 8 bits
     * PNG_TRANSFORM_PACKING       expand 1, 2 and 4-bit
     *                             samples to bytes
     */
    transform |= (PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING);

    //! read in the entire image at once
    png_read_png(png_ptr, info_ptr, transform, NULL);

    //! Get image informations
    o_width  = (size_t) png_get_image_width (png_ptr, info_ptr);
    o_height = (size_t) png_get_image_height(png_ptr, info_ptr);

    //! clean up and free any memory allocated, close the file
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    if (NULL != fp && stdin != fp) {
      (void) fclose(fp);
    }
  }
  else if (ext == "tif" || ext == "tiff") {

    //! Open the tiff file
    TIFF *tif = TIFFOpen(i_name, "r");

    //! Check that the file has been open
    if (!tif) {
      cout << "Unable to read TIFF file " << i_name << ". Abort." << endl;
      exit(EXIT_FAILURE);
    }

    //! Initialization
    uint32 w = 0, h = 0;

    //! Get the metadata
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);

    o_width  = (size_t) w;
    o_height = (size_t) h;

    //! Close the file
    TIFFClose(tif);
  }
  else {
    cout << "getSize: The input image hasn't the right extension name. Abort.";
    cout << endl;
    exit(EXIT_FAILURE);
  }
}


//! Set initial number of threads to use with OPENMP
int getOptimalThreadsNumber(
  const int p_nbThreads,
  const char* p_name) {

  //! Get the size of the image.
  size_t width, height;
  getSize(p_name, width, height);

  //! Choose the optimal number of threads
  int nb;
  if (width * height <= 64 * 64) {
    nb = std::min(2, p_nbThreads);
  }
  else if (width * height <= 128 * 128) {
    nb = std::min(4, p_nbThreads);
  }
  else if (width * height <= 192 * 192) {
    nb = std::min(6, p_nbThreads);
  }
  else if (width * height <= 512 * 512) {
    nb = std::min(8, p_nbThreads);
  }
  else if (width * height <= 1024 * 1024) {
    nb = std::min(16, p_nbThreads);
  }
  else {
    nb = p_nbThreads;
  }

  return nb;
}


//! Read nodes of a palette of color for topography visualization.
void readPalette(
  const char* i_fileName,
  std::vector<float> &o_nodes) {

  //! Open the file that contains the palette nodes
  ifstream myFile(i_fileName);
  if (!myFile) {
    cout << "Can not open the file " << i_fileName << endl;
    exit(EXIT_FAILURE);
  }

  //! Initialization
  o_nodes.clear();

  //! Read the file
  string line;
  while (std::getline(myFile, line)) {
    istringstream iss(line);
    float r, g, b;
    string color;
    iss >> r >> g >> b >> color;
    if (color.find("#") != string::npos) {
      o_nodes.push_back(r);
      o_nodes.push_back(g);
      o_nodes.push_back(b);
    }
  }

  //! Close the file
  myFile.close();

}














