#ifndef WAfloatERANDSEDIMENTATIONEVOLUTION_H_INCLUDED
#define WAfloatERANDSEDIMENTATIONEVOLUTION_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"
#include "../LibImages/LibImagesThreads.h"
#include "../Utilities/Parameters.h"


/**
 * @brief Small class to describe the Water Evolution scheme.
 **/
class WaterAndSedimentationEvolution {

  //! Public methods
  public:

    /**
     * @brief Default constructor.
     **/
    WaterAndSedimentationEvolution();

    /**
     * @brief Copy constructor.
     **/
    WaterAndSedimentationEvolution(
      const WaterAndSedimentationEvolution& i_wse);

    /**
     * @brief Operator overload.
     **/
    WaterAndSedimentationEvolution& operator=(
      const WaterAndSedimentationEvolution& i_wse);

    /**
     * @brief Default destructor.
     **/
    ~WaterAndSedimentationEvolution();

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
     * @brief Perform one step of the water and sedimentation evolution.
     *
     * @param io_env: see Environment.
     *
     * @return none.
     **/
    void performOneStepEvolution(
      const Image& i_imLandscape,
      Image& io_imSediment,
      Image& io_imWaterLvl,
      const size_t p_tid);

    /**
     * @brief Update time step.
     **/
    void setTimeStep(const float i_dt  ) {m_dt   = i_dt  ;}
    void setRain    (const float i_rain) {m_rain = i_rain;}
    void setRainMatrix    (const float* i_rainMatrix) {m_rainMatrix = (float*)i_rainMatrix;}


  //! Private functions
    /**
     * @brief Release memory.
     **/
    void releaseMemory();


  //! Data members
  private:

    //! Parameters
    const Parameters* m_params;
    float m_dt;
    float m_rain;
    float* m_rainMatrix;

    //! Size
    size_t m_width;
    size_t m_height;

    //! Images
    ImageX<float> *m_imWaterTmp;
    ImageX<float> *m_imSedimentTmp;
};
#else
class WaterAndSedimentationEvolution;

#endif // WAfloatERANDSEDIMENTATIONEVOLUTION_H_INCLUDED
