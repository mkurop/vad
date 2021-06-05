
#ifndef SPECTRAL_SUBTRACTION_H_
#define SPECTRAL_SUBTRACTION_H_

#include <armadillo>

//! Function for the simplest possible spectral subtraction algorithm
/*!
 * \param[in] noisy_psd noisy speech PSD of \f$K/2+1\f$ length where \f$K\f$ is the frame length
 * \param[in] noise_psd noise estimated PSD
 * \param[in] subtraction_floor the minimum of the enhanced speech magnitude spectrum
 * \param[in] oversubtraction oversubtraction parameter
 */
arma::Col<float> SpectralSubtraction(arma::Col<float> &noisy_psd, arma::Col<float> &noise_psd, float subtraction_floor, float oversubtraction);


#endif
