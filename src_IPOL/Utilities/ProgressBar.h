#ifndef PROGRESSBAR_H_INCLUDED
#define PROGRESSBAR_H_INCLUDED


//! Global includes
#include <stdlib.h>


//! Local includes


/**
 * @brief Progress bar class.
 **/
class ProgressBar {

  public:
    /**
     * @brief Default constructor.
     **/
    ProgressBar();

    /**
     * @brief Initialization constructor.
     **/
    ProgressBar(
      const size_t p_nbTotal,
      const size_t p_nbDisplay = 50);

    /**
     * @brief Default destructor.
     **/
    ~ProgressBar();

    /**
     * @brief Initialization.
     **/
    void init(
      const size_t p_nbTotal,
      const size_t p_current);

    /**
     * @brief Update and show.
     **/
    void update(
      const size_t p_current);

  private:
    /**
     * @brief Display the current state.
     **/
    void display();

  //! Data members
  private:
    size_t m_current;
    size_t m_total;
    size_t m_nbDisplay;
};
#else
class ProgressBar;

#endif // PROGRESSBAR_H_INCLUDED
