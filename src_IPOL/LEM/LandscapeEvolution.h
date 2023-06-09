#ifndef LANDSCAPEEVOLUTION_H_INCLUDED
#define LANDSCAPEEVOLUTION_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"
#include "../LibImages/LibImagesThreads.h"
#include "../Utilities/Parameters.h"
#include "../Utilities/Utilities.h"


/**
 * @brief Small class to describe the Creep Evolution scheme.
 **/
class LandscapeEvolution {

  //! Public methods
  public:

    /**
     * @brief Default constructor.
     **/
    LandscapeEvolution();

    /**
     * @brief Copy constructor.
     **/
    LandscapeEvolution(
      const LandscapeEvolution& i_le);


    /**
     * @brief Operator overload.
     **/
    LandscapeEvolution& operator=(
      const LandscapeEvolution& i_le);

    /**
     * @brief Default destructor.
     **/
    ~LandscapeEvolution();

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
     * @brief Perform a one step evolution of the landscape scheme.
     *
     * @param io_env: see Environment.
     *
     * @return none.
     **/
    void performOneStepEvolution(
      const Image &i_imWaterLvl,
      Image &io_imLandscape,
      Image &io_imWaterTimesConcentration,
      const size_t p_tid);

    /**
     * @brief Update time step.
     **/
    void setTimeStep(const float i_dt) {m_dt = i_dt;}

    /**
     * @brief Compute the transfert and update the new mean of the landscape.
     **/
    void computeTransfert(
      const Image& i_imWaterLvl,
      const Image &i_imLandscape,
      const Image &i_imWaterTimesConcentration,
      double &o_mean,
      const size_t p_tid);

    /**
     * @brief Update both landscape height and water times concentration and
     *        apply an uplift.
     **/
    void applyTransfert(
      Image& io_imLandscape,
      Image& io_imWaterTimesConcentration,
      const float p_uplift,
      const Image& i_imUplift,
      const size_t p_tid) const;


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
    ImageX<float>* m_imLandscapeVariation;

    //! Pow function
    float (*m_powM)(const float, const float);
    float (*m_powN)(const float, const float);
#ifdef __AVX__
    __m256 (*m_powM256)(const __m256, const float);
    __m256 (*m_powN256)(const __m256, const float);
#endif // __AVX__
#ifdef __SSE__
    __m128 (*m_powM128)(const __m128, const float);
    __m128 (*m_powN128)(const __m128, const float);
#endif // __SSE__
};
#else
class LandscapeEvolution;

#endif // LANDSCAPEEVOLUTION_H_INCLUDED
