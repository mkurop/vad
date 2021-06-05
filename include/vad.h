#ifndef VAD_H_INC
#define VAD_H_INC

#include <armadillo>
#include <cmath>

#include <log-probability.h>
#include <fast-noise-psd-tracking.h>
#include <string>

/*!
	* \brief Class implements the VAD algorithm from the paper:
	*
	* Jongseo Sohn, Nam Soo Kim, Wonyong Sung,
	* A Statistical Model-Based Voice Activity Detection,
	* IEEE Signal Processing Letters, vol. 6., no. 1, January 1999
	*/


class Vad {
  bool start;
  int mag_len;
  arma::Row<float> X_prev;
  double ksi_min, a01, a10, a00, a11, aa;

 public:
  double logsp, lognp, spprb;

  /*!
   * \brief Class implements the VAD algorithm
   *
   * \param samplingRate sampling rate in [Hz], defaults to 16000Hz
   * \param frame frame length before zero adding, defaults to 400 samples
   * \param frame2 frame length after zero adding, defaults to 512 samples
   * \param a01 log natural base probability of noise -> speech transition
   * \param a10 log natural base probability of speech -> noise transition
   */
  Vad(int samplingRate = 16000, int frame = 400, int frame2 = 512, double a01 = -.3, double a10 = -.3); 

  /*!
   * \brief Method for speech presence probability
   *
   * \param X the noisy speech PSD
   * \param lambdaN the estimated noise PSD
   *
   * \return speech presence probability
   */
  double vad(const arma::Row<float> &X, const arma::Row<float> &lambdaN); 
    
};

#endif
