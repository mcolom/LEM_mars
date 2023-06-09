#ifndef CONNECTEDCOMPONENT_H_INCLUDED
#define CONNECTEDCOMPONENT_H_INCLUDED


//! Global includes
#include <vector>
#include <stdlib.h>


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief Short class describing a connected component.
 **/
class ConnectedComponent {

  //! Public methods
  public:

    /**
     * @brief Default constructor.
     **/
     ConnectedComponent();

    /**
     * @brief Copy constructor.
     **/
    ConnectedComponent(
      const ConnectedComponent& i_cc);

    /**
     * @brief Operator overload.
     **/
    ConnectedComponent& operator=(
      const ConnectedComponent& i_cc);

     /**
      * @brief Default destructor.
      **/
    ~ConnectedComponent();

    /**
     * @brief Setters.
     **/
    void setLevel  (const float i_level  ) {m_level       = i_level  ;}
    void setProcess(const bool  i_process) {m_toBeProcess = i_process;}
    void addElement(const std::pair<int, int> &i_element) {
      m_elements.push_back(i_element);}

    /**
     * @brief Getters.
     **/
    float  getLevel     () const {return m_level          ;}
    bool   getProcess   () const {return m_toBeProcess    ;}
    size_t getNbElements() const {return m_elements.size();}
    std::pair<int, int> getElement(const size_t p_n) const {
      return m_elements[p_n];}

    /**
     * @brief Extract connected components of the current pixel.
     *
     * @param i_im : image on which to perform the extraction;
     * @param io_mask : mask of already processed pixels. Will be updated.
     * @param p_i, p_j : current position of the pixel.
     *
     * @return none.
     **/
    void extract(
      const Image& i_im,
      Image& io_mask,
      const int p_i,
      const int p_j,
      const bool p_useDiag = false);

    /**
     * @brief Process the regularization of the image according to the current
     *        set of connected pixels.
     *
     * @param i_im : input reference image;
     * @param o_im : will have a regularized version of i_im.
     *
     * @return none.
     **/
    void process(
      const Image& i_im,
      Image &o_im,
      const bool p_useDiag = false);

  //! Private methods

    /**
     * @brief Release memory.
     **/
    void releaseMemory();


  //! Data members
  private:
    float m_level;
    bool m_toBeProcess;
    std::vector<std::pair<int, int> > m_elements;
};
#else
class ConnectedComponent;

#endif // CONNECTEDCOMPONENT_H_INCLUDED
