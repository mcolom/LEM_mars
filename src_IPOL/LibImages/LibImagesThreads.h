#ifndef LIBIMAGESTHREADS_H_INCLUDED
#define LIBIMAGESTHREADS_H_INCLUDED

//! Global includes
#include <stdlib.h>
#include <string>
#include <png.h>

//! Local includes


/**
 * @brief Small Image class, which contains basic manipulations.
 *
 **/
template <class T>
class ImageX {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    ImageX();


    /**
     * @brief Surcharged constructor with a specified size.
     *
     * @param i_width, i_height, i_channels : size of the wanted image;
     * @param i_border: size of the border to add around the image.
     **/
    ImageX(
      const size_t i_width,
      const size_t i_height,
      const size_t i_channels,
      const size_t i_border = 0);


    /**
     * @brief Copy constructor.
     **/
    ImageX(
      const ImageX<T>& i_im);


    /**
     * @brief Operator overload.
     **/
    ImageX<T>& operator=(
      const ImageX<T>& i_im) {
      new (this) ImageX<T>(i_im);
      return *this;}


    /**
     * @brief Getter
     **/
    size_t width   () const {return m_width   ;}
    size_t height  () const {return m_height  ;}
    size_t channels() const {return m_channels;}
    size_t border  () const {return m_border  ;}
    size_t nHeight (const size_t p_n) const {return m_heights[p_n];}


    //! Default destructor
    ~ImageX();


    /**
     * @brief Get the pointer to the channel c, row i, according to the current
     *        active thread n.
     *
     * @param p_n : current thread;
     * @param p_c : current channel to adress;
     * @param p_i : current row to adress.
     *
     * @return the pointer of I(n, c, i, 0);
     **/
    T* getPtr(
      const int p_n = 0,
      const int p_c = 0,
      const int p_i = 0) const;


  //! Data members
  private:

    //! Size
    size_t m_width;
    size_t m_height;
    size_t m_channels;
    size_t m_border;

    //! Miscalleneous
    size_t m_nbThreads;

    //! Pointers
    T** m_ptrs;
    size_t* m_heights;
};
#else
template <class T>
class ImageX;

#endif // LIBIMAGES_H_INCLUDED
