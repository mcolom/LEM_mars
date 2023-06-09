#ifndef CREEPEVOLUTION_H_INCLUDED
#define CREEPEVOLUTION_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"
#include "../LibImages/LibImagesThreads.h"
#include "../Utilities/Parameters.h"


/**
 * @brief Small class to describe the Creep Evolution scheme.
 **/
class CreepEvolution {

  //! Public methods
  public:

    /**
     * @brief Default constructor.
     **/
    CreepEvolution();

    /**
     * @brief Copy constructor.
     **/
    CreepEvolution(
      const CreepEvolution& i_ce);


    /**
     * @brief Operator overload.
     **/
    CreepEvolution& operator=(
      const CreepEvolution& i_ce);

    /**
     * @brief Default destructor.
     **/
    ~CreepEvolution();

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
     * @brief Perform a one step evolution of the creep scheme.
     *
     * @param io_env: see Environment.
     *
     * @return none.
     **/
    void performOneStepEvolution(
      Image &io_imLandscape,
      const size_t p_tid);

    /**
     * @brief Update time step.
     **/
    void setTimeStep(const float i_dt) {m_dt = i_dt;}


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

    //! Size
    size_t m_width;
    size_t m_height;

    //! Images
    ImageX<float> *m_imLandscapeVariation;
};
#else
class CreepEvolution;

#endif // CREEPEVOLUTION_H_INCLUDED
