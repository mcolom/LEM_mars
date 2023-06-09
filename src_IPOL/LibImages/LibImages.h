#ifndef LIBIMAGES_H_INCLUDED
#define LIBIMAGES_H_INCLUDED

//! Global includes
#include <stdlib.h>
#include <string>
#include <png.h>
#include <vector>

//! Local includes


/**
 * @brief Small Image class, which contains basic manipulations.
 *
 **/
class Image {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    Image();


    /**
     * @brief Surcharged constructor with a specified size.
     *
     * @param i_width, i_height, i_channels : size of the wanted image;
     * @param i_border: size of the border to add around the image.
     **/
    Image(
      const size_t i_width,
      const size_t i_height,
      const size_t i_channels,
      const size_t i_border = 0);


    /**
     * @brief Copy constructor.
     **/
    Image(
      const Image& i_im);


    /**
     * @brief Operator overload.
     **/
    Image& operator=(
      const Image& i_im);


    /**
     * @brief Getter
     **/
    size_t width    () const {return m_width    ;}
    size_t height   () const {return m_height   ;}
    size_t channels () const {return m_channels ;}
    size_t border   () const {return m_border   ;}
    size_t nHeight  (const size_t p_n) const {return m_heights[p_n];}


    /**
     * @brief Initialize an image full of zeros.
     **/
    void init(
      const size_t i_width,
      const size_t i_height,
      const size_t i_channels = 1,
      const size_t i_border = 0) {
      this->releaseMemory();
      new (this) Image(i_width, i_height, i_channels, i_border);}

    void init(
      const Image& i_im) {
      this->init(i_im.m_width, i_im.m_height, i_im.m_channels, i_im.m_border);}


    //! Default destructor
    ~Image();


    /**
     * @brief Generic read of an image. Call readPng or readTiff according to
     *        the extension of the input name of the image to read.
     **/
    void read(
      const char* p_name,
      const size_t i_border = 0);


    /**
     * @brief Generic write of an image. Call writePng or writeTiff according to
     *        the extension of the input name of the image to write.
     *
     * @param p_name : path to the image to save;
     * @param p_quad: if true, apply a zoom-in by duplication of 1 pixel into 4.
     *
     * @return none.
     **/
     void write(
      const char* p_name,
      const bool p_quad = false) const;


    /**
     * @brief Get the pointer to the channel c, row i.
     *
     * @param p_c : current channel to adress,
     * @param p_i : current row to adress.
     *
     * @return the pointer of I(c, i, 0);
     **/
    float* getPtr(
      const int p_c = 0,
      const int p_i = 0) const {
      return m_ptr + (p_c * int(m_height + 2 * m_border) + p_i + (int) m_border)
                          * int(m_width  + 2 * m_border) + (int) m_border;
    }


    /**
     * @brief Add a border to the current image.
     *
     * @param i_border: size of the border to add.
     *
     * @return none.
     **/
    void addBorder(
      const size_t i_border);


    /**
     * @brief Set a value for the border of the current image.
     *
     * @param p_value : value for the initialization.
     *
     * @return none.
     **/
    void setBorder(
      const float p_value);


    /**
     * @brief Apply a zoom-in of the input image as array and normalize it.
     *
     * @param p_min, p_max: range of the normalized image.
     *
     * @return the zoomed-in normalized image.
     **/
    void applySqrtNormalization(
      Image& o_im,
      const float p_min,
      const float p_max) const;


    /**
     * @brief Apply a zoom-in of the input image as array and normalize it.
     *
     * @param o_im: MUST be initialized before;
     * @param p_min, p_max: range of the normalized image;
     * @param i_min, i_max: min/max values of the current image.
     *
     * @return the zoomed-in normalized image.
     **/
    void applySqrtNormalization(
      Image& o_im,
      const size_t p_tid,
      const float p_min,
      const float p_max,
      const float i_min,
      const float i_max) const;


    /**
     * @brief Multiply the image by a constant.
     **/
    void multiply(
      const float p_constant,
      const int p_tid = -1);


    /**
     * @brief Return a zoom-out of an image.
     **/
    void zoomOut(
      Image& o_im) const;


    /**
     * @brief Return a zoom-in of an image.
     **/
    void zoomIn(
      Image& o_im) const;

    /**
     * @brief In-place normalization of an image.
     **/
    void normalize(
      const float p_min = 0,
      const float p_max = 255);

    /**
     * @brief Apply a Gaussian blur.
     */
    void applyGaussianBlur(
      Image& o_im,
      const float p_sigma) const;

    /**
     * @brief Apply in-place a Gaussian noise of standard deviation sigma.
     *
     * @param p_sigma: standard deviation of the noise to add;
     * @param p_oceanLvl: don't add noise under the ocean level.
     **/
    void applyGaussianNoise(
      const float p_sigma,
      const float p_oceanLvl = 0.f);

    /**
     * @brief Symmetrize the borders.
     *
     * @param p_tid: current thread id. If -1, then apply it on the full image.
     * @param p_sym: kind of symmetry to apply:
     *               0:       44 456 66
     *               1:       65 456 54
     *               2:       23 456 78
     *               3:       00 456 00 (Dirichlet)
     *
     * @return none.
     **/
    void fillBorder(
      const int p_tid = -1,
      const size_t p_sym = 0);

    /**
     * @brief Save the image as a 3D model in OBJ format.
     **/
    void saveAsObj(
      const std::string &p_fileName) const;

    /**
     * @brief Save the image as a 3D model in OBJ format, and add a water.
     **/
    void saveAsObj(
      const std::string &p_fileName,
      const Image& i_imWater,
      const float p_threshold,
      const float p_oceanLevel) const;

    /**
     * @brief Save the image as a 3D model in OBJ format with false color.
     **/
    void saveAsObj(
      const std::string &p_fileName,
      const Image& i_imWater,
      const Image& i_imColor,
      const float p_oceanLevel,
      const float p_norm) const;

    /**
     * @brief Upsampling until (width * height) is close to size.
     **/
    void sample(
      Image& o_im,
      const size_t p_size,
      const int p_tid = -1) const;

    /**
     * @brief Get the topography in false color of a landscape, according to the
     *        water level and the ocean level.
     **/
    void getTopo(
      Image& o_im,
      const Image& i_imWater,
      const float p_oceanLvl,
      const float p_min,
      const float p_max,
      const char* p_folder,
      const float p_waterThreshold = 1.f) const;

    /**
     * @brief Get a colorized version of a landscape image.
     **/
    void getColorized(
      Image& o_im,
      const Image& i_imWater,
      const float p_oceanLvl,
      const float p_min,
      const float p_max,
      const char* p_folder,
      const float p_waterThreshold = 1.f) const;

    /**
     * @brief Get Shadow image.
     **/
    void getShadow(
      Image& o_im) const;

    /**
     *@brief Stretch an image based on quantile.
     **/
    void stretch(
      Image& o_im,
      const float p_quantile = 0.01,
      const float p_min = 0,
      const float p_max = 255) const;

    /**
     * @brief Blur an image with a Cauchy kernel by using the fftw.
     **/
    void blurCauchy(
      Image& o_im,
      const float p_sigma) const;

    /**
     * @brief Get stats of an image (min, max, and mean values).
     **/
    void getStats(
      float& o_min,
      float& o_max,
      float& o_mean) const;

    /**
     * @brief Apply a HSV visualization of an image.
     *
     * @param o_im: must be initialized with 3 channels.
     **/
    void getHSV(
      Image& o_im,
      const float p_min = 0,
      const float p_max = 1,
      const int p_tid = -1) const;

    /**
     * @brief Check if the current pixel as a lower neighbor.
     **/
    bool hasLowerNeighbor(
      const int p_i,
      const int p_j,
      const size_t p_c,
      const bool p_useDiag = false) const;

    /**
     * @brief Get the laplacien of the image multiplied by a factor.
     **/
    void getLaplacien(
      Image& o_im,
      const float p_factor) const;

    /**
     * @brief Make a cone image.
     **/
    void makeCone(
      const bool p_descend = true,
      const bool p_half = true);

    /**
     * @brief Make a plan.
     *
     * @param p_theta: belongs to [0, 360[. Define the orientation of the plan.
     **/
    void makePlan(
      const size_t p_theta);


    /**
     * @brief For debug only. Print an image with its border.
     **/
    void print() const;


  private:
    /**
     * @brief Copy the inner image of the current image to another image.
     *
     * @param o_im : will recieve the inner image of this, without modifying its
     *               borders.
     *
     * @return none.
     **/
    void copyInner(
      Image& o_im,
      const int p_tid =-1) const;


    /**
     * @brief Read an image via the Libpng library. Will exit the main program
     *        in case of problem.
     *
     * @param i_name : path to the image which will be filled into p_ptr;
     * @param i_border : size of the border to add around the image (will be
     *                   initialized full of zeros.
     *
     * @return none.
     **/
    void readPng(
      const char* p_name,
      const size_t i_border = 0);


    /**
     * @brief Read an image via the Libtiff library. Will exit the main program
     *        in case of problem.
     *
     * @param i_name : path to the image which will be filled into p_ptr;
     * @param i_border : size of the border to add around the image (will be
     *                   initialized full of zeros.
     *
     * @return none.
     **/
    void readTiff(
      const char* p_name,
      const size_t i_border = 0);


    /**
     * @brief Write an image via the Libpng library. Will exit the main problem
     *        in case of problem.
     *
     * @param i_name: path to the image to save;
     * @param p_quad: if true, apply a zoom-in by duplication of 1 pixel into 4.
     **/
    void writePng(
      const std::string &p_name,
      const bool p_quad = false) const;


    /**
     * @brief Write an image via the Libtiff library. Will exit the main problem
     *        in case of problem.
     *
     * @param i_name: path to the image to save;
     * @param p_quad: if true, apply a zoom-in by duplication of 1 pixel into 4.
     **/
    void writeTiff(
      const std::string &p_name,
      const bool p_quad = false) const;


    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


  //! Data members
  private:

    //! Size
    size_t m_width;
    size_t m_height;
    size_t m_channels;
    size_t m_border;
    size_t m_size;

    //! Miscalleneous
    size_t m_nbThreads;

    //! Pointers
    float* m_ptr;
    size_t* m_heights;
};
#else
class Image;

#endif // LIBIMAGES_H_INCLUDED
