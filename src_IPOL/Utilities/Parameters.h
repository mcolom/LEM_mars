#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED


//! Global includes
#include <iostream>
#include <string.h>


//! Local includes


/**
 * @brief Class for the parameters of the main function.
 **/
class Parameters {

  //! Methods
  public:

    /**
     * @brief Default constructor.
     **/
    Parameters();

    /**
     * @brief Default destructor.
     **/
    ~ Parameters();

    /**
     * @brief Copy constructor.
     **/
    Parameters(
      const Parameters& i_params);

    /**
     * @brief Operator overload.
     **/
    Parameters& operator=(const Parameters& i_params);

    /**
     * @brief Getters.
     **/
    char*  inpName      () const {return m_inpName      ;}
    char*  outName      () const {return m_outName      ;}
    bool   verbose      () const {return m_verbose      ;}
    float  timeStep     () const {return m_timeStep     ;}
    float  rain         () const {return m_rain         ;}
    float*  rainMatrix  () const {return m_rainMatrix   ;}
    const char*  rainName  () const {return m_rainName  ;}
    float  expM         () const {return m_expM         ;}
    float  expN         () const {return m_expN         ;}
    float  e            () const {return m_e            ;}
    float  creep        () const {return m_creep        ;}
    float  s            () const {return m_s            ;}
    int    nbIter       () const {return m_nbIter       ;}
    int    nbWater      () const {return m_nbWater      ;}
    float  oceanLevel   () const {return m_oceanLevel   ;}
    float  percLand     () const {return m_percLand     ;}
    float  initWaterLvl () const {return m_initWaterLvl ;}
    float  alpha        () const {return m_alpha        ;}
    int    nbThreads    () const {return m_nbThreads    ;}
    bool   useMultiScale() const {return m_useMultiScale;}
    float  boostWater   () const {return m_boostWater   ;}
    float  boostIter    () const {return m_boostIter    ;}
    bool   visualizeMesh() const {return m_visualizeMesh;}
    char*  visualizeTopo() const {return m_visualizeTopo;}
    float  zoom         () const {return m_zoom         ;}
    bool   computeTWI   () const {return m_computeTWI   ;}
    bool   autoParams   () const {return m_autoParams   ;}
    float  blur         () const {return m_blur         ;}
    float  noise        () const {return m_noise        ;}
    float  waterColor   () const {return m_waterColor   ;}
    bool   stabCurv     () const {return m_stabCurv     ;}
    int    autoUplift   () const {return m_autoUplift   ;}
    bool   autoE        () const {return m_autoE        ;}
    float  addNoise     () const {return m_addNoise     ;}
    int    iterNoise    () const {return m_iterNoise    ;}
    int    iterWater    () const {return m_iterWater    ;}
    float  rainWater    () const {return m_rainWater    ;}
    int    borderWater  () const {return m_borderWater  ;}
    int    borderLand   () const {return m_borderLand   ;}
    bool   noZoomOut    () const {return m_noZoomOut    ;}
    float  deltaM       () const {return m_deltaM       ;}
    float  deltaN       () const {return m_deltaN       ;}
    float  deltaS       () const {return m_deltaS       ;}

    /**
     * @brief Setters.
     **/
    void setE(const float i_e)  {m_e         = i_e;}
    void setS(const float i_s)  {m_s         = i_s;}
    void setC(const float i_c)  {m_creep     = i_c;}
    void setP(const float i_p)  {m_percLand  = i_p;}
    void setM(const float i_m)  {m_expM      = i_m;}
    void setN(const float i_n)  {m_expN      = i_n;}
    void setI(const float i_i)  {m_nbIter    = i_i;}
    void setB(const float i_b)  {m_boostIter = i_b;}
    void setV(const bool  i_v)  {m_verbose   = i_v;}
    void setR(const float i_r)  {m_rain      = i_r;}
    void setRMatrix(float* i_r) {m_rainMatrix  = (float*)i_r;}
    void setW(const bool  i_w)  {m_useMultiScale = i_w;}

    /**
     * @brief Setters.
     **/
    void setTimeStep(const float i_timeStep) {m_timeStep = i_timeStep;}

    /**
     * @brief Read the input arguments. Detect the "-h" option, and print the
     *        informations. Otherwise, get the input arguments, if valid.
     **/
    int checkArgs(
      const int i_argc,
      char** i_argv);

    /**
     * @brief Print the list of parameters.
     **/
    void print() const;


  //! Miscellaneous functions
  private:

    /**
     * @brief Contain the copyright.
     **/
    void printCopyright() const;

    /**
     * @brief Print the synopsis of this program.
     **/
    void printSynopsis() const;

    /**
     * @brief Print the description of this program.
     **/
    void printDescription() const;

    /**
     * @brief Print the input parameters list.
     **/
    void printInput() const;

    /**
     * @brief Print the output parameters list.
     **/
    void printOutput() const;

    /**
     * @brief Print the optional parameters list.
     **/
    void printOptional() const;

    /**
     * @brief Print a line of a specific m_size maximum number of characters,
     *        with a pad start.
     **/
    void printLine(
      const std::string &i_sentence,
      const std::string &i_pad) const;

    /**
     * @brief Print a line of a specific m_size maximum number of characters,
     *        and add a word in the last position.
     **/
    void printWord(
      const std::string &i_line,
      const std::string &i_word,
      const std::string &i_pad = "") const;

  //! Data members
  private:

    //! Misellaneaous
    size_t m_size;          // Maximum number of characters to print by line

    //! Parameters
    char*  m_inpName;       // Path of the input image            -i  [%s]
    char*  m_outName;       // Path of the ouput image            -o  [%s]
    bool   m_verbose;       // Activate the verbose mode          -v
    float  m_timeStep;      // Time step                          -t  [%f]
    float  m_rain;          // Rain                               -r  [%f]
    float* m_rainMatrix;    // Rain matrix
    const char* m_rainName; // Path of the rain image             -rm [%s]
    float  m_expM;          // Exp M                              -m  [%f]
    float  m_expN;          // Exp N                              -n  [%f]
    float  m_e;             // Erosion e                          -e  [%f]
    float  m_creep;         // Erosion creep                      -c  [%f]
    float  m_s;             // Erosion s                          -s  [%f]
    int    m_nbIter;        // Maximal number of iterations       -ni [%d]
    int    m_nbWater;       // Number of iterations for water     -nw [%d]
    float  m_oceanLevel;    // Ocean level                        -lo [%f]
    float  m_percLand;      // Percentage of landscape to remove  -p  [%f]
    float  m_initWaterLvl;  // Initial water level                -lw [%f]
    float  m_alpha;         // Water divergence efficiency "alpha"-ca [%f]
    int    m_nbThreads;     // Number of threads to use           -d  [%d]
    bool   m_useMultiScale; // Activate the multiscale water init -wm
    float  m_boostWater;    // Activate the boost for water init  -bw [%s]
    float  m_boostIter;     // Activate the boost for main iter   -bi [%s]
    bool   m_visualizeMesh; // Save 3D mesh of the landscape      -f
    char*  m_visualizeTopo; // Visualize topography with colors   -g  [%s]
    float  m_zoom;          // Apply a zoom in the input          -z  [%f]
    bool   m_computeTWI;    // Compute the TWI index              -j
    float  m_blur;          // Apply a Gaussian blur              -k  [%f]
    float  m_waterColor;    // Threshold to see the water         -wc [%f]
    float  m_noise;         // Noise to add over the landscape    -q  [%f]
    bool   m_autoParams;    // Activate auto parameters selection -a
    bool   m_stabCurv;      // Activate the stability curve       -cs
    int    m_autoUplift;    // Activate the automatic uplift      -au [%d]
    bool   m_autoE;         // Activate the auto estimation of e  -ae
    float  m_addNoise;      // Add noise every iterNoise iteration-qa [%f]
    int    m_iterNoise;     // Nb of iterations between add noise -qi [%d]
    int    m_iterWater;     // Nb of iterations of water only     -wi [%d]
    float  m_rainWater;     // Rain for the water only            -rw [%f]
    int    m_borderWater;   // Border conditions for water        -cw [%d]
    int    m_borderLand;    // Border conditions for landscape    -cl [%d]
    float  m_deltaM;        // Increment for expM selection       -dm [%f]
    float  m_deltaN;        // Increment for expN selection       -dn [%f]
    float  m_deltaS;        // Increment for sediment selection   -ds [%f]
    bool   m_noZoomOut;     // Don't apply the zoom out if zoom   -nz
};
#else
class Parameters;

#endif // PARAMETERS_H_INCLUDED
