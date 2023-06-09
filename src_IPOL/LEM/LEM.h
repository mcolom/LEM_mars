#ifndef LEM_H_INCLUDED
#define LEM_H_INCLUDED


//! Global includes
#include <time.h>
#include <string.h>
#include <vector>


//! Local includes
#include "../Utilities/Parameters.h"
#include "../LibImages/LibImages.h"
#include "../Utilities/ProgressBar.h"
#include "WaterEvolution.h"
#include "CreepEvolution.h"
#include "LandscapeEvolution.h"
#include "WaterAndSedimentationEvolution.h"


/**
 * @brief Class to perform the Landscape Evolution Model.
 **/
class LEM {

  //! Public methods
  public:
    /**
     * @brief Default constructor.
     **/
    LEM();

    /**
     * @brief Copy constructor.
     **/
    LEM(
      const LEM& i_lem);

    /**
     * @brief Operator overload.
     **/
    LEM& operator=(
      const LEM& i_lem);

    /**
     * @brief Default destructor.
     **/
    ~LEM();

    /**
     * @brief Initialize the environment. Must be call first.
     **/
    void init(
      Parameters& p_params);

    /**
     * @brief Finish the water initialization.
     **/
    void waterInitialization();

    /**
     * @brief Perform the main iteration after the water initialization.
     **/
    void runMainIterations();

    /**
     * @brief Automatic selection of parameters.
     **/
    void autoSelect();

    /**
     * @brief Save images after the main evolution.
     **/
    void saveImages();


  //! Private interne functions
  private:

    /**
     * @brief Core part of the water initialization.
     **/
    void waterInitializationCore(
      const size_t p_nbIter);

    /**
     * @brief Core part of the main iteration.
     **/
    void runMainIterationsCore();

    /**
     * @brief Automatically estimate the value of erosion.
     **/
    void runAutoErosion();

    /**
     * @brief Integrate landscape over ocean level.
     *
     * @return the number of landscape higher than the ocean level.
     **/
    float integrateLandscapeOverOceanLevel() const;
    void integrateLandscapeOverOceanLevel(
      const size_t p_tid);

    /**
     * @brief Integrate water over ocean level.
     *
     * @return the normalize number of water for landscape higher than the ocean
     *         level.
     **/
    void integrateWaterOverOceanLevel(
      float& o_mass,
      float& o_max,
      float& o_norm) const;

    /**
     * @brief Compute some stats.
     **/
    void computeStats();

    /**
     * @brief Release memory.
     **/
    void releaseMemory();

    /**
     * @brief Adjust the time step parameter.
     **/
     void adjustTimeStep();

    /**
     * @brief Initialize the info parameters.
     **/
    void initializeInfos();

    /**
     * @brief Update the info parameters.
     **/
    void updateInfos(
      ProgressBar& i_bar,
      size_t &io_nbBar);

    /**
     * @brief Compute the stability norm over the original landscape z0 and the
     *        evolved landscape z1.
     **/
    void computeStabilityNorm(
      float& o_norm);

    /**
     * @brief Compute the stability norm over the original landscape z0 and the
     *        evolved landscape z1.
     **/
    void computeStabilityNorm_tid(
      const size_t p_tid);

    /**
     * @brief Compute another stability norm.
     **/
    void computeNorm(
      float& o_norm) const;

    /**
     * @brief Compute the current uplift term.
     **/
    void computeUplift(
      const double* i_means);

    /**
     * @brief Un-zoom the images if needed.
     **/
    void applyZoomOut();


  //! Data members
  private:

    //! Miscellaneous
    struct timespec m_time;
    size_t m_nbThreads;
    bool m_waterHasBeenInitialized;

    //! Parameters
    Parameters *m_params;

    //! Sizes
    size_t m_width;
    size_t m_height;
    size_t m_border;

    //! Information values
    float m_massOriginal;
    float m_massTarget;
    float m_massCurrent;
    double m_meanOriginal;
    float m_uplift;
    float m_upliftAll;
    size_t m_iterCurrent;
    size_t m_iterMax;
    float m_ratio;
    size_t m_nbCheck;
    float* m_masses;
    float* m_norms;
    float* m_nbs;

    //! Images
    Image* m_imLandscape;
    Image* m_imWaterLvl;
    Image* m_imSediment;
    Image* m_imLandscapeInit;
    Image* m_imWaterLvlInit;
    Image* m_imLandscapeRaw;
    Image* m_imUplift;
    
    //! Rain matrix
    float *m_rainMatrix = NULL;
    //! Rain parameter
    float m_rain;

    //! Evolution schemes
    WaterEvolution                 *m_evolutionWater;
    CreepEvolution                 *m_evolutionCreep;
    LandscapeEvolution             *m_evolutionLandscape;
    WaterAndSedimentationEvolution *m_evolutionWaterAndSedimentation;
};
#else
class LEM;

#endif // LEM_H_INCLUDED
