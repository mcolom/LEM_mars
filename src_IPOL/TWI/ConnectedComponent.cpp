/**
 * @file ConnectedComponent.cpp
 *
 * @brief Class to handle connected components.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <queue>


//! Local includes
#include "ConnectedComponent.h"


using namespace std;


//! Default constructor
ConnectedComponent::ConnectedComponent() :

  m_level      (0),
  m_toBeProcess(false),
  m_elements   (std::vector<std::pair<int, int> >()) {
}


//! Copy constructor.
ConnectedComponent::ConnectedComponent(
  const ConnectedComponent& i_cc) :

  m_level      (i_cc.m_level),
  m_toBeProcess(i_cc.m_toBeProcess),
  m_elements   (i_cc.m_elements) {
}


//! Default destructor
ConnectedComponent::~ConnectedComponent() {
  releaseMemory();
}


//! Operator overload.
ConnectedComponent& ConnectedComponent::operator=(
  const ConnectedComponent& i_cc) {

  if (&i_cc == this) {
    return *this;
  }

  releaseMemory();

  new (this) ConnectedComponent(i_cc);
  return *this;
}


//! Extract connected components of the current pixel.
void ConnectedComponent::extract(
  const Image& i_im,
  Image& io_mask,
  const int p_i,
  const int p_j,
  const bool p_useDiag) {

  //! Clear the connected component
  m_elements.clear();

  //! Update the level
  m_level = i_im.getPtr(0, p_i)[p_j];

  //! Initialize the queue
  std::queue<std::pair<int, int> > q;
  q.push(make_pair(p_i, p_j));

  //! Update the mask at the current position
  io_mask.getPtr(0, p_i)[p_j] = 1;

  //! Look for connected components
  while(!q.empty()) {

    //! Get the position of the pixel
    const std::pair<int, int> s = q.front();
    const int i = s.first;
    const int j = s.second;

    //! Remove this position from the queue
    q.pop();

    //! Add this position to the connected components
    m_elements.push_back(s);

    //! For convenience
    const float* iIT = i_im.getPtr(0, i - 1);
    const float* iIC = i_im.getPtr(0, i);
    const float* iIB = i_im.getPtr(0, i + 1);
    float* oMT = io_mask.getPtr(0, i - 1);
    float* oMC = io_mask.getPtr(0, i);
    float* oMB = io_mask.getPtr(0, i + 1);

    //! Look for potential connected neighbor non already processed.
    //! Left neighbor
    if (m_level == iIC[j - 1] && oMC[j - 1] == 0) {
      oMC[j - 1] = 1;
      q.push(make_pair(i, j - 1));
    }

    //! Right neighbor
    if (m_level == iIC[j + 1] && oMC[j + 1] == 0) {
      oMC[j + 1] = 1;
      q.push(make_pair(i, j + 1));
    }

    //! Top neighbor
    if (m_level == iIT[j] && oMT[j] == 0) {
      oMT[j] = 1;
      q.push(make_pair(i - 1, j));
    }

    //! Bottom neighbor
    if (m_level == iIB[j] && oMB[j] == 0) {
      oMB[j] = 1;
      q.push(make_pair(i + 1, j));
    }

    //! Diagonal elements
    if (p_useDiag) {

      //! Top left
      if (m_level == iIT[j - 1] && oMT[j - 1] == 0) {
        oMT[j - 1] = 1;
        q.push(make_pair(i - 1, j - 1));
      }

      //! Top right
      if (m_level == iIT[j + 1] && oMT[j + 1] == 0) {
        oMT[j + 1] = 1;
        q.push(make_pair(i - 1, j + 1));
      }

      //! Bottom left
      if (m_level == iIB[j - 1] && oMB[j - 1] == 0) {
        oMB[j - 1] = 1;
        q.push(make_pair(i + 1, j - 1));
      }

      //! Bottom left
      if (m_level == iIB[j + 1] && oMB[j + 1] == 0) {
        oMB[j + 1] = 1;
        q.push(make_pair(i + 1, j + 1));
      }
    }
  }

  //! Update the flag
  m_toBeProcess = false;
}


//! Process the regularization of the image according to the current set of
//! connected pixels.
void ConnectedComponent::process(
  const Image& i_im,
  Image &o_im,
  const bool p_useDiag) {

  //! Initialization
  const size_t width  = i_im.width ();
  const size_t height = i_im.height();
  const float epsilon = 1.f / (255.f * float(width * height));

  //! The image is_processed is used during connected component extraction
  Image isProcessed(width, height, 1, i_im.border());

  //! All points that have a lower neighbor are set to distance 0 and we put
  //! them into the queue
  std::queue<std::pair<int, int> > q[2];
  for (size_t l = 0; l < m_elements.size(); l++) {

    //! Get the current position
    const std::pair<int, int> s = m_elements[l];
    const int i = s.first;
    const int j = s.second;

    //! Check if it has lower neighbor
    if (i_im.hasLowerNeighbor(i, j, 0, p_useDiag)) {
      isProcessed.getPtr(0, i)[j] = 1;
      q[0].push(s);
    }
  }

  //! q[0] contains all pixel that touches the border
  //! They are at distance 0
  int indexCurr = 0;
  int indexNext = 1;
  int dist      = 0;
  bool shallContinue = true;

  //! While it exists connected components
  while (shallContinue) {
    while (!q[indexCurr].empty()) {

      //! Get the current position
      const std::pair<int, int> s = q[indexCurr].front();
      const int i = s.first;
      const int j = s.second;

      //! For convenience
      const float* iIT = i_im.getPtr(0, i - 1);
      const float* iIC = i_im.getPtr(0, i);
      const float* iIB = i_im.getPtr(0, i + 1);
      float* oI = o_im.getPtr(0, i);
      float* iPT = isProcessed.getPtr(0, i - 1);
      float* iPC = isProcessed.getPtr(0, i);
      float* iPB = isProcessed.getPtr(0, i + 1);

      //! Remove the current pixel from the queue
      q[indexCurr].pop();

      //! Set new distance
      oI[j] = iIC[j] + epsilon * dist;

      //! Left neighbor
      if (iPC[j - 1] == 0 && iIC[j - 1] == m_level) {
        iPC[j - 1] = 1;
        q[indexNext].push(std::make_pair(i, j - 1));
      }

      //! Right neighbor
      if (iPC[j + 1] == 0 && iIC[j + 1] == m_level) {
        iPC[j + 1] = 1;
        q[indexNext].push(std::make_pair(i, j + 1));
      }

      //! Top neighbor
      if (iPT[j] == 0 && iIT[j] == m_level) {
        iPT[j] = 1;
        q[indexNext].push(std::make_pair(i - 1, j));
      }

      //! Bottom neighbor
      if (iPB[j] == 0 && iIB[j] == m_level) {
        iPB[j] = 1;
        q[indexNext].push(std::make_pair(i + 1, j));
      }

      //! Diagonal elements
      if (p_useDiag) {

        //! Top left
        if (iPT[j - 1] == 0 && iIT[j - 1] == m_level) {
          iPT[j - 1] = 1;
          q[indexNext].push(std::make_pair(i - 1, j - 1));
        }

        //! Top right
        if (iPT[j + 1] == 0 && iIT[j + 1] == m_level) {
          iPT[j + 1] = 1;
          q[indexNext].push(std::make_pair(i - 1, j + 1));
        }

        //! Bottom left
        if (iPB[j - 1] == 0 && iIB[j - 1] == m_level) {
          iPB[j - 1] = 1;
          q[indexNext].push(std::make_pair(i + 1, j - 1));
        }

        //! Bottom right
        if (iPB[j + 1] == 0 && iIB[j + 1] == m_level) {
          iPB[j + 1] = 1;
          q[indexNext].push(std::make_pair(i + 1, j + 1));
        }
      }
    }

    //! Check if it exists at least one pixel to process next
    if (q[indexNext].empty()) {
      shallContinue = false;
    }

    //! Change the queue index
    std::swap(indexNext, indexCurr);

    //! Update the current distance
    dist += 1;
  }
}


//! Release memory
void ConnectedComponent::releaseMemory() {
  m_elements.clear();
}












