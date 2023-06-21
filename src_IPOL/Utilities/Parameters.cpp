/**
 * @file Parameters.cpp
 *
 * @brief Parameters class: will initialize and contains the parameters used
 *        for the main function.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <cstdio>


//! Local includes
#include "Parameters.h"


using namespace std;


//! Default constructor
Parameters::Parameters(
  ) :
  //! Miscellaneous
  m_size         (43 + 2 * 15),

  //! Parameters
  m_inpName      (NULL),       // Path of the input image            -i  [%s]
  m_outName      (NULL),       // Path of the ouput image            -o  [%s]
  m_verbose      (false),      // Activate the verbose mode          -v
  m_timeStep     (1.f),        // Time step                          -t  [%f]
  m_rain         (1.f  * 1e-7),// Rain                               -r  [%f]
  m_rainMatrix   (NULL),       // Rain matrix                        -rm [%s]
  m_expM         (0.5f),       // Exp M                              -m  [%f]
  m_expN         (1.0f),       // Exp N                              -n  [%f]
  m_e            (1.f),        // Erosion e                          -e  [%f]
  m_creep        (1.f  * 1e-3),// Erosion creep                      -c  [%f]
  m_s            (0.f  * 1e-6),// Erosion s                          -s  [%f]
  m_nbIter       (1000),       // Maximal number of iterations       -ni [%d]
  m_nbWater      (40000),      // Number of iterations for water     -nw [%d]
  m_oceanLevel   (0.f),        // Ocean level                        -lo [%f]
  m_percLand     (15.f * 1e-2),// Percentage of landscape to remove  -p  [%f]
  m_initWaterLvl (0.f  * 1e-3),// Initial water level                -lw [%f]
  m_alpha        (1.f),        // Water divergence efficiency "alpha"-ca [%f]
#ifdef _OPENMP
  m_nbThreads    (4),          // Number of threads to use           -d  [%d]
#else
  m_nbThreads    (1),
#endif // _OPEN
  m_useMultiScale(false),      // Activate the multiscale water init -wm
  m_boostWater   (0.f),        // Activate the boost for water init  -bw [%s]
  m_boostIter    (0.f),        // Activate the boost for main iter   -bi [%s]
  m_visualizeMesh(false),      // Save 3D meshes of the landscape    -f
  m_visualizeTopo(NULL),       // Visualize topography with colors   -g  [%s]
  m_zoom         (1.f),        // Apply a zoom in the input image    -z  [%f]
  m_computeTWI   (false),      // Compute the TWI                    -j
  m_blur         (0.f),        // Apply a Gaussian blur              -k  [%f]
  m_waterColor   (1.f),        // Threshold to see the water         -wc [%f]
  m_noise        (0.f),        // Noise to add                       -q  [%f]
  m_autoParams   (false),      // Activate auto parameters selection -a
  m_stabCurv     (false),      // Activate the stability curve       -cs
  m_autoUplift   (false),      // Activate the automatic uplift      -au
  m_autoE        (0),          // Activate the auto estimation of E  -ae [%d]
  m_addNoise     (0.f),        // Add noise every iterNoise iter     -qa [%f]
  m_iterNoise    (1),          // Nb of iterations between addNoise  -qi [%d]
  m_iterWater    (0),          // Nb of iterations of water only     -wi [%d]
  m_rainWater    (m_rain),     // Rain for the water only            -rw [%f]
  m_borderWater  (0),          // Border conditions for the water    -cw [%d]
  m_borderLand   (0),          // Border conditions for the landscape-cl [%d]
  m_deltaM       (0.1f),       // Increment for expM selection       -dm [%f]
  m_deltaN       (0.1f),       // Increment for expN selection       -dn [%f]
  m_deltaS       (0.1f * 1e-6),// Increment for sediment selection   -ds [%f]
  m_noZoomOut    (false)       // Don't apply the zoom out if zoom   -nz
  {

}


//! Default destructor
Parameters::~Parameters(
  ) {
}


//! Copy constructor.
Parameters::Parameters(
  const Parameters& i_params) :
  //! Miscalleneous
  m_size(i_params.m_size),

  //! Parameters
  m_inpName      (i_params.inpName()),
  m_outName      (i_params.outName()),
  m_verbose      (i_params.verbose()),
  m_timeStep     (i_params.timeStep()),
  m_rain         (i_params.rain()),
  m_rainMatrix   (i_params.rainMatrix()),
  m_expM         (i_params.expM()),
  m_expN         (i_params.expN()),
  m_e            (i_params.e()),
  m_creep        (i_params.creep()),
  m_s            (i_params.s()),
  m_nbIter       (i_params.nbIter()),
  m_nbWater      (i_params.nbWater()),
  m_oceanLevel   (i_params.oceanLevel()),
  m_percLand     (i_params.percLand()),
  m_initWaterLvl (i_params.initWaterLvl()),
  m_alpha        (i_params.alpha()),
  m_nbThreads    (i_params.nbThreads()),
  m_useMultiScale(i_params.useMultiScale()),
  m_boostWater   (i_params.boostWater()),
  m_boostIter    (i_params.boostIter()),
  m_visualizeMesh(i_params.visualizeMesh()),
  m_visualizeTopo(i_params.visualizeTopo()),
  m_zoom         (i_params.zoom()),
  m_computeTWI   (i_params.computeTWI()),
  m_blur         (i_params.blur()),
  m_waterColor   (i_params.waterColor()),
  m_noise        (i_params.noise()),
  m_autoParams   (i_params.autoParams()),
  m_stabCurv     (i_params.stabCurv()),
  m_autoUplift   (i_params.autoUplift()),
  m_autoE        (i_params.autoE()),
  m_addNoise     (i_params.addNoise()),
  m_iterNoise    (i_params.iterNoise()),
  m_iterWater    (i_params.iterWater()),
  m_rainWater    (i_params.rainWater()),
  m_borderWater  (i_params.borderWater()),
  m_borderLand   (i_params.borderLand()),
  m_deltaM       (i_params.deltaM()),
  m_deltaN       (i_params.deltaN()),
  m_deltaS       (i_params.deltaS()),
  m_noZoomOut    (i_params.noZoomOut()) {

}

//! Operator overload.
Parameters& Parameters::operator=(
  const Parameters& i_params) {

  if (&i_params == this) {
    return *this;
  }

  new (this) Parameters(i_params);
  return *this;
}


//! Read the input arguments.
int Parameters::checkArgs(
  const int i_argc,
  char** i_argv) {

  //! If no arguments entered
  if (i_argc <= 1) {
    cout << "No arguments detected. Type -h or -help for some help." << endl;
    return EXIT_FAILURE;
  }

  //! Detect the parameter entered
  for (int n = 1; n < i_argc; n++) {

    //! Get the current argument
    const string sarg = i_argv[n];

    //! If help is required
    if (sarg.find("-h") == 0 || sarg.find("-help") == 0) {

      //! Copyright
      this->printCopyright();

      //! Synopsis
      this->printSynopsis();

      //! Description
      this->printDescription();

      //! Input parameters
      this->printInput();

      //! Output parameters
      this->printOutput();

      //! Other options
      this->printOptional();

      return EXIT_FAILURE;
    }

    //! Input path
    if (sarg == "-i") {
      if (n + 1 < i_argc) {
        m_inpName = i_argv[++n];
      }
    }

    //! Output path
    if (sarg == "-o") {
      if (n + 1 < i_argc) {
        m_outName = i_argv[++n];
      }
    }

    //! Time step
    if (sarg == "-t") {
      if (n + 1 < i_argc) {
        m_timeStep = atof(i_argv[++n]);
      }
    }

    //! Rain
    if (sarg == "-r") {
      if (n + 1 < i_argc) {
        m_rain = atof(i_argv[++n]);
      }
    }

    //! Rain matrix
    m_rainName = NULL; // CORRECCION!![Miguel]
    if (sarg == "-rm") {
      if (n + 1 < i_argc) {
        m_rainName = i_argv[++n];
        m_rainMatrix = NULL;
      }
    }

    //! ExpM
    if (sarg == "-m") {
      if (n + 1 < i_argc) {
        m_expM = atof(i_argv[++n]);
      }
    }

    //! ExpN
    if (sarg == "-n") {
      if (n + 1 < i_argc) {
        m_expN = atof(i_argv[++n]);
      }
    }

    //! Erosion E
    if (sarg == "-e") {
      if (n + 1 < i_argc) {
        m_e = atof(i_argv[++n]);
      }
    }

    //! Creep
    if (sarg == "-c") {
      if (n + 1 < i_argc) {
        m_creep = atof(i_argv[++n]);
      }
    }

    //! Erosion S
    if (sarg == "-s") {
      if (n + 1 < i_argc) {
        m_s = atof(i_argv[++n]);
      }
    }

    //! Number of iterations
    if (sarg == "-ni") {
      if (n + 1 < i_argc) {
        m_nbIter = atoi(i_argv[++n]);
      }
    }

    //! Number of iteration for the water initialization
    if (sarg == "-nw") {
      if (n + 1 < i_argc) {
        m_nbWater = atoi(i_argv[++n]);
      }
    }

    //! Ocean level
    if (sarg == "-lo") {
      if (n + 1 < i_argc) {
        m_oceanLevel = atof(i_argv[++n]);
      }
    }

    //! Percentage of landscape to remove
    if (sarg == "-p") {
      if (n + 1 < i_argc) {
        m_percLand = atof(i_argv[++n]) * 1e-2;
      }
    }

    //! Initialization of the water level
    if (sarg == "-lw") {
      if (n + 1 < i_argc) {
        m_initWaterLvl = atof(i_argv[++n]);
      }
    }

    //! Water divergence efficiency "alpha"
    if (sarg == "-ca") {
      if (n + 1 < i_argc) {
        m_alpha = atof(i_argv[++n]);
        printf("Using coefficient alpha =  %f\n", m_alpha);
      }
    }

    //! Number of threads
    if (sarg == "-d") {
      if (n + 1 < i_argc) {
        m_nbThreads = atoi(i_argv[++n]);
      }
    }

    //! Multiscale water initialization
    if (sarg == "-wm") {
      m_useMultiScale = true;
    }

    //! Boost for water init
    if (sarg == "-bw") {
      if (n + 1 < i_argc) {
        m_boostWater = atof(i_argv[++n]);
      }
    }

    //! Boost for main iteration
    if (sarg == "-bi") {
      if (n + 1 < i_argc) {
        m_boostIter = atof(i_argv[++n]);
      }
    }

    //! Visualize mesh
    if (sarg == "-f") {
      m_visualizeMesh = true;
    }

    //! Visualize topo
    if (sarg == "-g") {
      if (n + 1 < i_argc) {
        m_visualizeTopo = i_argv[++n];
      }
    }

    //! Zoom
    if (sarg == "-z") {
      if (n + 1 < i_argc) {
        m_zoom = atof(i_argv[++n]);
      }
    }

    //! Topographic Wetness Index
    if (sarg == "-j") {
      m_computeTWI = true;
    }

    //! Automatic selection of the parameters
    if (sarg == "-a") {
      m_autoParams = true;
    }

    //! Gaussian blur
    if (sarg == "-k") {
      if (n + 1 < i_argc) {
        m_blur = atof(i_argv[++n]);
      }
    }

    //! Water color threshold
    if (sarg == "-wc") {
      if (n + 1 < i_argc) {
        m_waterColor = atof(i_argv[++n]);
      }
    }

    //! Noise
    if (sarg == "-q") {
      if (n + 1 < i_argc) {
        m_noise = atof(i_argv[++n]) * 1e-2;
      }
    }

    //! Stability curve
    if (sarg == "-cs") {
      m_stabCurv = true;
    }

    //! Automatic uplift
    if (sarg == "-au") {
      if (n + 1 < i_argc) {
        m_autoUplift = atoi(i_argv[++n]);
      }
    }

    //! Auto estimation of E
    if (sarg == "-ae") {
      m_autoE = true;
    }

    //! Add noise
    if (sarg == "-qa") {
      if (n + 1 < i_argc) {
        m_addNoise = atof(i_argv[++n]) * 1e-2;
      }
    }

    //! Iter noise
    if (sarg == "-qi") {
      if (n + 1 < i_argc) {
        m_iterNoise = atoi(i_argv[++n]);
      }
    }

    //! Iter water only
    if (sarg == "-wi") {
      if (n + 1 < i_argc) {
        m_iterWater = atoi(i_argv[++n]);
      }
    }

    //! Rain for the water only
    if (sarg == "-rw") {
      if (n + 1 < i_argc) {
        m_rainWater = atof(i_argv[++n]);
      }
    }

    //! Border conditions for the water
    if (sarg == "-cw") {
      if (n + 1 < i_argc) {
        m_borderWater = atoi(i_argv[++n]);
      }
    }

    //! Border conditions for the landscape
    if (sarg == "-cl") {
      if (n + 1 < i_argc) {
        m_borderLand = atoi(i_argv[++n]);
      }
    }

    //! Increment for expM selection
    if (sarg == "-dm") {
      if (n + 1 < i_argc) {
        m_deltaM = atof(i_argv[++n]);
      }
    }

    //! Increment for expN selection
    if (sarg == "-dn") {
      if (n + 1 < i_argc) {
        m_deltaN = atof(i_argv[++n]);
      }
    }

    //! Increment for sediment selection
    if (sarg == "-ds") {
      if (n + 1 < i_argc) {
        m_deltaS = atof(i_argv[++n]) * 1e-6;
      }
    }

    //! Don't apply the zoom out if zoom
    if (sarg == "-nz") {
      m_noZoomOut = true;
    }

    //! Verbose option
    if (sarg == "-v") {
      m_verbose = true;
    }
  }

  //! Check that the mandatory arguments have been given
  if (m_inpName == NULL) {
    cout << "The input image path must be precised. Use \n";
    cout << "    -i path/image.ext" << endl;
    return EXIT_FAILURE;
  }
  if (m_outName == NULL) {
    cout << "The output image path must be precised. Use \n";
    cout << "    -o path/image.ext" << endl;
    return EXIT_FAILURE;
  }

  //! To avoid numerical issue that slow down a lot
  m_rain = std::max(m_rain, 1e-16f);

  //! For convenience
  if (m_autoE) {
    m_boostIter  = 0.f;
    m_autoUplift = false;
  }

  return EXIT_SUCCESS;
}


//! Contain the copyright.
void Parameters::printCopyright(
  ) const {

  const string s0 = "* -------------------------------- *";
  const string s1 = "| Copyright (c) 2015-2016 by IPOL. |";
  const string s2 = "|       All Rights Reserved.       |";
  const string st = std::string((m_size - s0.size()) / 2, ' ');

  cout << endl;
  cout << st << s0 << st << endl;
  cout << st << s1 << st << endl;
  cout << st << s2 << st << endl;
  cout << st << s0 << st << endl;
  cout << endl;
}


//! Print the synopsis of this program.
void Parameters::printSynopsis(
  ) const {

  //! Header
  const string header = "SYNOPSIS:";

  //! Name of the program
  const string name = "./LEM";

  //! Parameters
  string sentence = name;
  sentence += " [-i input image path]";
  sentence += " [-o output image path]";
  sentence += " [-t time step]";
  sentence += " [-r rain]";
  sentence += " [-rm rainMatrix]";
  sentence += " [-m expM]";
  sentence += " [-n expM]";
  sentence += " [-e erosion E]";
  sentence += " [-c creep]";
  sentence += " [-s erosion S]";
  sentence += " [-ni number of iterations]";
  sentence += " [-nw number of water iterations]";
  sentence += " [-lo ocean level]";
  sentence += " [-p landscape percentage]";
  sentence += " [-lw water level init]";
  sentence += " [-ca Water divergence efficiency 'alpha']";
  sentence += " [-a eta]";
  sentence += " [-d number of threads]";
  sentence += " [-wm activate multiscale]";
  sentence += " [-bw activate water boost]";
  sentence += " [-bi activate init boost]";
  sentence += " [-f visualize mesh]";
  sentence += " [-g visualize topo]";
  sentence += " [-z zoom]";
  sentence += " [-j TWI]";
  sentence += " [-k blur]";
  sentence += " [-a auto params]";
  sentence += " [-wc water color threshold]";
  sentence += " [-q  noise]";
  sentence += " [-cs stability curve]";
  sentence += " [-au auto uplift]";
  sentence += " [-ae auto E]";
  sentence += " [-qa add Noise]";
  sentence += " [-qi iter Noise]";
  sentence += " [-qi iter Water only]";
  sentence += " [-rw rain for Water only]";
  sentence += " [-cw border conditions for water]";
  sentence += " [-cl border conditions for land]";
  sentence += " [-dm expM increment]";
  sentence += " [-dn expN increment]";
  sentence += " [-ds sediment increment]";
  sentence += " [-nz deactivate the zoom out]";
  sentence += " [-v verbose]";

  //! Print the synopsis
  cout << header << endl << endl;

  //! Print the content
  this->printLine(sentence, "  ");

  //! Finish the synopsis
  cout << endl << endl;
}


//! Print the description of this program.
void Parameters::printDescription(
  ) const {

  //! Header
  const string header = "DESCRIPTION:";

  //! Content
  const string s1 = "This algorithm implements a landscape evolution model.";

  //! Print the description
  cout << header << endl << endl;

  //! Print the content
  this->printLine(s1, "  ");

  //! Finish the description
  cout << endl << endl;
}


//! Print the input parameters list.
void Parameters::printInput(
  ) const {

  //! Header
  const string header = "INPUT OPTIONS:";

  //! List of options
  string sentence;
  sentence += "-i [%s] : path to the input image";

  //! Print the header
  cout << header << endl << endl;

  //! Print the options
  this->printLine(sentence, "    ");

  //! Finish the list
  cout << endl << endl;
}


//! Print the output parameters list.
void Parameters::printOutput(
  ) const {

  //! Header
  const string header = "OUTPUT OPTIONS (results):";

  //! List of options
  string sentence;
  sentence += "-o [%s] : path to the output image";

  //! Print the header
  cout << header << endl << endl;

  //! Print the options
  this->printLine(sentence, "    ");

  //! Finish the list
  cout << endl << endl;
}


//! Print the optional parameters list.
void Parameters::printOptional(
  ) const {

  //! For convenience
  const string s4 = std::string(4, ' ');
  const string s7 = std::string(7, ' ');

  //! Print header
  this->printWord("OTHER OPTIONS:", "Default");
  cout << endl;

  //! Time Step
  this->printWord("-t  (optional)", "1.0", s4);
  this->printLine("Time step", s7);

  //! Rain
  this->printWord("-r  (optional)", "1.0", s4);
  this->printLine("Rain parameter for the hole filling and water "\
    "initialization", s7);

  //! Exp M
  this->printWord("-m  (optional)", "0.5", s4);
  this->printLine("Exponent M", s7);

  //! Exp N
  this->printWord("-n  (optional)", "1.0", s4);
  this->printLine("Exponent N", s7);

  //! Erosion E
  this->printWord("-e  (optional)", "1.0", s4);
  this->printLine("Erosion E", s7);

  //! Erosion creep
  this->printWord("-c  (optional)", "1.0", s4);
  this->printLine("Erosion creep", s7);

  //! Erosion S
  this->printWord("-s  (optional)", "0.0", s4);
  this->printLine("Erosion S", s7);

  //! Maximal number of iterations
  this->printWord("-ni (optional)", "1000", s4);
  this->printLine("Maximal number of iterations", s7);

  //! Number of iterations for the water
  this->printWord("-nw (optional)", "40000", s4);
  this->printLine("Number of iterations for the water", s7);

  //! Ocean level
  this->printWord("-lo (optional)", "0.0", s4);
  this->printLine("Ocean level", s7);

  //! Landscape percentage
  this->printWord("-p  (optional)", "15", s4);
  this->printLine("Percentage of landscape to remove (x1e-2)", s7);

  //! Initial water of level
  this->printWord("-lw (optional)", "0.0", s4);
  this->printLine("Initial water of level", s7);
  
  //! Water divergence efficiency "alpha"
  this->printWord("-ca (optional)", "1.0", s4);
  this->printLine("Water divergence efficiency 'alpha'", s7);

  //! Eta
  this->printWord("-a  (optional)", "0.01", s4);
  this->printLine("Eta parameter for the hole filling and water initialization",
                  s7);

  //! Number of threads
  this->printWord("-d  (optional)", "4", s4);
  this->printLine("If OMP is used, determine the number of threads to use.",
                  s7);

  //! Activate the multiscale water initialization
  this->printWord("-wm (optional)", "False", s4);
  this->printLine("Activate the multiscale water initialization. Three scales \
                  will be used.", s7);

  //! Boost for water init
  this->printWord("-bw (optional)", "0.f", s4);
  this->printLine("Activate the boost for water init if bw > 0. The initial \
                  time step will be equal to bw.", s7);

  //! Boost for main iterations
  this->printWord("-bi (optional)", "0.f", s4);
  this->printLine("Activate the boost for main iterations if bi > 0. The \
                  initial time step will be equal to bi.", s7);

  //! Meshes
  this->printWord("-f  (optional)", "False", s4);
  this->printLine("Save 3D meshes for visualization.", s7);

  //! Topography
  this->printWord("-g  (optional)", "NULL", s4);
  this->printLine("If a path to the folder containing both water and ground \
                  topological palette is provided, visualize topography of the \
                  landscape with false color.", s7);

  //! Zoom
  this->printWord("-z  (optional)", "1.0", s4);
  this->printLine("If > 1 a zoom-in of this factor will first be applied on \
                   the input image.", s7);

  //! Topographic Wetness Index
  this->printWord("-j  (optional)", "False", s4);
  this->printLine("Activate the computation of the Topographic Wetness Index."
                  , s7);

  //! Automatic selection of the parameters
  this->printWord("-a  (optional)", "False", s4);
  this->printLine("Activate the automatic selection of the parameters.", s7);

  //! Gaussian blur
  this->printWord("-k  (optional)", "0.0", s4);
  this->printLine("If k > 0, apply a Gaussian blur with sigma = k before the\
  landscape evolution.", s7);

  //! Water color threshold
  this->printWord("-wc (optional)", "1.0", s4);
  this->printLine("Threshold (wc * mean(water)) for topo colorization below \
  which the ground appears as ocean.", s7);

  //! Noise
  this->printWord("-q  (optional)", "0.0", s4);
  this->printLine("If q > 0, apply a Gaussian noise of standard deviation (q / \
  100) over the normalized landscape.", s7);

  //! Stability curve
  this->printWord("-cs (optional)", "False", s4);
  this->printLine("Activate the stability curve computation and saving." ,s7);

  //! Automatic uplift
  this->printWord("-au (optional)", "0", s4);
  this->printLine("Activate the automatic uplift: 0 no uplift, 1 uniform uplift\
  , 2 local uplift.", s7);

  //! Auto estimation of E
  this->printWord("-ae (optional)", "False", s4);
  this->printLine("Activate the automatic estimation of E for a given number\
  of iteration and percentage to reach.", s7);

  //! Add Noise
  this->printWord("-qa (optional)", "0.0", s4);
  this->printLine("If qa > 0, apply a Gaussian noise of standard deviation (q /\
   100) over the normalized landscape every qi iterations.", s7);

  //! Iter Noise
  this->printWord("-qi (optional)", "1", s4);
  this->printLine("If qa > 0, will apply the noise every qi iterations.", s7);

  //! Iter Water only
  this->printWord("-wi (optional)", "0", s4);
  this->printLine("Number of iterations of water only equation to apply for \
  every iterations of the main iterations.", s7);

  //! Rain for water only
  this->printWord("-rw (optional)", "1.0", s4);
  this->printLine("Rain for the water only option.", s7);

  //! Border conditions for water
  this->printWord("-cw (optional)", "0", s4);
  this->printLine("Border conditions for the water. 0: 44 456 66, 1: 65 456 54\
  , 2: 23 456 78, 3: 00 456 00.", s7);

  //! Border conditions for landscape
  this->printWord("-cl (optional)", "0", s4);
  this->printLine("Border conditions for the landscape. 0: 44 456 66, 1: 65 456\
   54, 2: 23 456 78, 3: 00 456 00.", s7);

  //! Increment for expM selection
  this->printWord("-dm (optional)", "0.1", s4);
  this->printLine("Increment for expM values to test during the automatic \
  selection of parameter values.", s7);

  //! Increment for expM selection
  this->printWord("-dn (optional)", "0.1", s4);
  this->printLine("Increment for expN values to test during the automatic \
  selection of parameter values.", s7);

  //! Increment for expM selection
  this->printWord("-ds (optional)", "0.1 (x1e-6)", s4);
  this->printLine("Increment for sediment values to test during the automatic \
  selection of parameter values.", s7);

  //! Deactivate the zoom-out
  this->printWord("-nz (optional)", "False", s4);
  this->printLine("Deactivate the final zoom-out when the zoom option is chosen\
  .", s7);

  //! Verbose
  this->printWord("-v  (optional)", "False", s4);
  this->printLine("Activate the verbose mode.", s7);
  cout << endl;

  //! Help
  this->printWord("-h  (optional)", "False", s4);
  this->printLine("Print the help.", s7);
  cout << endl << endl;
}

//! Print a line of a specific m_size maximum number of characters
void Parameters::printLine(
  const std::string &i_sentence,
  const std::string &i_pad) const {

  //! Initializations
  istringstream iss(i_sentence);

  //! Get all the words contained in this sentence
  vector<string> line;
  do {
    string word;
    iss >> word;
    line.push_back(word);
  } while (iss);

  //! Print the line
  size_t nb = i_pad.size();
  cout << i_pad;
  for (size_t n = 0; n < line.size(); n++) {

    //! Get the number of characters contained in the current word
    const size_t nc = line[n].size();

    //! Print the current word
    if (nb + nc <= m_size) {
      cout << line[n];
      nb += nc;
    }
    else {
      nb = nc + i_pad.size();
      cout << endl;
      cout << i_pad << line[n];
    }

    //! Print the space character
    if (nb < m_size) {
      cout << " ";
      nb++;
    }
    else {
      cout << endl << i_pad;
      nb = i_pad.size();
    }
  }

  //! Finish the line
  cout << endl;
}


//! Print a line of a specific m_size maximum number of characters,
//! and add a word in the last position.
void Parameters::printWord(
  const std::string &i_line,
  const std::string &i_word,
  const std::string &i_pad) const {

  //! Check the number of characters
  const size_t nbL = i_line.size();
  const size_t nbW = i_word.size();
  const size_t nbP = i_pad.size();

  if (nbL + nbW + nbP < m_size) {
    const string s = std::string(m_size - nbL - nbW - nbP, ' ');
    cout << i_pad << i_line << s << i_word << endl;
  }
  else {
    const string s = std::string(m_size - nbP - nbW, ' ');
    cout << i_pad << i_line << endl;
    cout << i_pad << s << i_word << endl;
  }
}


//! Print the list of parameters
void Parameters::print() const {

  const string s4 = std::string(4, ' ');
  /*cout << "Path to the input image: " << endl;
  cout << "  " << m_inpName << endl;
  cout << "Path to the output folder: " << endl;
  cout << "  " << m_outName << endl;*/
  this->printWord("List of parameters", "current value");

  //! Time step
  ostringstream oss1;
  oss1 << m_timeStep;
  this->printWord(" - Time step:", oss1.str(), s4);

  //! Rain
  ostringstream oss2;
  oss2 << m_rain;
  this->printWord(" - Rain:", oss2.str(), s4);

  //! ExpM
  ostringstream oss3;
  oss3 << m_expM;
  this->printWord(" - ExpM:", oss3.str(), s4);

  //! ExpM
  ostringstream oss14;
  oss14 << m_expN;
  this->printWord(" - ExpN:", oss14.str(), s4);

  //! Erosion e
  ostringstream oss4;
  oss4 << m_e;
  this->printWord(" - Erosion e:", oss4.str(), s4);

  //! Creep
  ostringstream oss5;
  oss5 << m_creep;
  this->printWord(" - Creep:", oss5.str(), s4);

  //! Erosion s
  ostringstream oss6;
  oss6 << m_s;
  this->printWord(" - Erosion s:", oss6.str(), s4);

  //! Nb Iter
  ostringstream oss7;
  oss7 << m_nbIter;
  this->printWord(" - NbIter:", oss7.str(), s4);

  //! Nb Water
  ostringstream oss8;
  oss8 << m_nbWater;
  this->printWord(" - NbWater:", oss8.str(), s4);

  //! Ocean level
  ostringstream oss9;
  oss9 << m_oceanLevel;
  this->printWord(" - Ocean level:", oss9.str(), s4);

  //! Percentage
  ostringstream oss11;
  oss11 << m_percLand;
  this->printWord(" - Percentage to remove:", oss11.str(), s4);

  //! Init water level
  ostringstream oss12;
  oss12 << m_initWaterLvl;
  this->printWord(" - Init water level:", oss12.str(), s4);

  //! Water divergence efficiency "alpha"
  ostringstream oss13;
  oss13 << m_alpha;
  this->printWord(" - alpha:", oss13.str(), s4);

  //! Number of threads
  ostringstream oss15;
  oss15 << m_nbThreads;
  this->printWord(" - NbThreads:", oss15.str(), s4);
}




