#ifndef SOHNVAD_H_INC
#define SOHNVAD_H_INC

#include <vad.h>
#include <armadillo>

namespace t2r {

/*!
	* \brief Class implements the VAD algorithm.
  *
  * The algorithm used is documented in the paper:
	*
	* Jongseo Sohn, Nam Soo Kim, Wonyong Sung,
	* A Statistical Model-Based Voice Activity Detection,
	* IEEE Signal Processing Letters, vol. 6., no. 1, January 1999
  *
  * For noise background estimation the algorithm for fast noise tracking proposed by Q. Zhang
  * in the paper is used:
  * 
  * Qiquan Zhang, Mingijang Wang, Yun Lu, Muhammad Idrees, Lu Zhang,
  * Fast Nonstationary Noise Tracking Based on Log-Spectral Power
  * MMSE Estimator and Temporal Recursive Averaging,
  * IEEE Access, vol. 7., pp. 80985-80999, 2019
  *
	*/
class Vad {

  public:

    Vad(int sampling_rate);

    /*!
     * \brief Computes the speech presence probability
     *
     * \param[in] signal_frame signal frame for the VAD decision
     *
     * \return the speech presence probability
     *
     */
    double operator()(arma::Col<float> signal_frame);

  private:

    arma::Col<float> window_;

    int frame_length_;

    FastNoisePSDTracking fast_noise_tracking_;

    Vad_ vad_;

    bool start_ = true; 
    
    int sampling_rate_;
};

}

#endif 
