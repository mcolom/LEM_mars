#ifndef WATEREVOLUTION_H_INCLUDED
#define WATEREVOLUTION_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"
#include "../Utilities/Parameters.h"
#include "../LibImages/LibImagesThreads.h" // to remove


/**
 * @brief Small class to describe the Water Evolution scheme.
 **/
class WaterEvolution {

  //! Public methods
  public:

    /**
     * @brief Default constructor.
     **/
    WaterEvolution();

    /**
     * @brief Copy constructor.
     **/
    WaterEvolution(
      const WaterEvolution& i_we);

    /**
     * @brief Operator overload.
     **/
    WaterEvolution& operator=(
      const WaterEvolution& i_we);

    /**
     * @brief Default destructor.
     **/
    ~WaterEvolution();

    /**
     * @brief Initialize the images and parameters of the Water Evolution.
     *
     * @param p_width, p_height: size of the images;
     * @param p_border: size of the border;
     * @param p_params: see Parameters().
     **/
    void init(
      const size_t p_width,
      const size_t p_height,
      const size_t p_border,
      Parameters const& p_params);

    /**
     * @brief Perform a one step evolution of the water scheme.
     *
     * @param io_env: see Environment.
     *
     * @return none.
     **/
    void performOneStepEvolution(
      const Image& i_imLandscape,
      Image& io_imWaterLvl,
      const size_t p_tid);

    /**
     * @brief Update time step.
     **/
    void setTimeStep(const float i_dt  ) {m_dt   = i_dt  ;}
    void setRain    (const float i_rain) {m_rain = i_rain;}
    void setRainMatrix(const float* i_rainMatrix) {m_rainMatrix = (float*)i_rainMatrix;}


  //! Private functions
    /**
     * @brief Release memory.
     **/
    void releaseMemory();


  //! Data members
  private:

    //! Parameters
    float m_dt;
    float m_rain;
    float* m_rainMatrix;
    float m_oceanLvl;
    float m_mars_a;

    //! Size
    size_t m_width;
    size_t m_height;

    //! Images
    ImageX<float> *m_imTmp;
};
#else
class WaterEvolution;

#endif // WATEREVOLUTION_H_INCLUDED
