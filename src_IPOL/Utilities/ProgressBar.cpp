/**
 * @file ProgressBar.cpp
 *
 * @brief Small class to display the progression of the a loop.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <iostream>
#include <sstream>


//! Local includes
#include "ProgressBar.h"


using namespace std;


//! Default constructor
ProgressBar::ProgressBar() :
  m_current  (0),
  m_total    (1),
  m_nbDisplay(70) {

}


//! Initialization constructor
ProgressBar::ProgressBar(
  const size_t p_nbTotal,
  const size_t p_nbDisplay) :
  m_current  (0),
  m_total    (p_nbTotal),
  m_nbDisplay(p_nbDisplay) {

}


//! Default destructor
ProgressBar::~ProgressBar(
  ) {
}


//! Initialization
void ProgressBar::init(
  const size_t p_nbTotal,
  const size_t p_current) {

  m_total   = p_nbTotal;
  m_current = p_current;
}


//! Update and show.
void ProgressBar::update(
  const size_t p_current) {

  //! Check that the end has been reached
  m_current = std::min(int(p_current), std::max(int(m_total) - 1, 0));
  this->display();

  //! If we are done
  if (p_current >= m_total - 1) {
    cout << endl;
  }
}


//!
void ProgressBar::display(
  ) {

  //! Get the percentage
  const float perc = 100.f * float(m_current + 1) / std::max(float(m_total), 1.f);
  ostringstream ossPerc;
  ossPerc << (int) perc;
  const string sPerc = (perc < 10 ? "00" + ossPerc.str() :
                        perc < 100 ? "0" + ossPerc.str() : ossPerc.str());

  //! Get the number of '-' and '='
  const size_t nbDone = (size_t) (perc * float(m_nbDisplay) / 100.f);
  const size_t nbToDo = m_nbDisplay - nbDone;
  const string sBar = std::string(nbDone, '=') + std::string(nbToDo, '-');

  //! Display the current state
  cout << "\r" << "[" << sBar << "] - " << sPerc << "%" << std::flush;
}

















