#ifndef TWI_H_INCLUDED
#define TWI_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief Small structure that contains a flow and its position.
 **/
struct Flow {
  int posI;
  int posJ;
  float flow;

  Flow (int i, int j, float f) : posI(i), posJ(j), flow(f) {}
};


/**
 * @brief Class to perform the Topographic Wetness Index
 **/
class TWI {

  //! Public methods
  public:

    /**
     * @brief Default constructor.
     **/
    TWI();

    /**
     * @brief Copy constructor.
     **/
    TWI(
      const TWI& i_twi);

    /**
     * @brief Operator overload.
     **/
    TWI& operator=(
      const TWI& i_twi);

    /**
     * @brief Default destructor.
     **/
    ~TWI();

    /**
     * @brief Compute the Topographic Wetness Index.
     *
     * @param i_im: contains the original landscape that will be holes filled;
     * @param o_im: will contain the Topographic Wetness Index.
     *
     * @return none.
     **/
    void computeTWI(
      Image const& i_im,
      Image &o_im);


  //! Private functions
  private:

    /**
     * @brief Release memory.
     **/
    void releaseMemory();

    /**
     * @brief Initialize the images.
     **/
    void init(
      Image const& i_im);

    /**
     * @brief Apply a geodesic reconstruction by in-place erosion.
     **/
    void applyGeodesicReconstruction();

    /**
     * @brief Perform a regularization of the image.
     **/
    void performRegularization();

    /**
     * @brief Compute the catchment area for each pixel.
     **/
    void computeCatchmentArea();

    /**
     * @brief Compute the index.
     **/
    void computeIndex(
      Image& o_im);


  //! Data members
  private:

    //! Miscellaneous
    size_t m_width;
    size_t m_height;
    int m_border;

    //! Parameters
    bool m_useDiag;

    //! Images
    Image *m_imOriginal;
    Image *m_im;
    Image *m_imRegularized;
    Image *m_imCatchmentArea;
};
#else
class TWI;


#endif // TWI_H_INCLUDED
