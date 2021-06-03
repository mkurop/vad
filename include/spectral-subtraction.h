
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
arma::Col<float> spectralSubtraction(arma::Col<float> &noisy_psd, arma::Col<float> &noise_psd, float subtraction_floor, float oversubtraction){

  arma::Col<float> enhanced_speech_magnitude_spectrum = \
      arma::sqrt(noisy_psd) - oversubtraction*arma::sqrt(noise_psd)

  enhanced_speech_magnitude_spectrum.clamp(subtraction_floor,arma::max(enhanced_speech_magnitude_spectrum))

  return enhanced_speech_magnitude_spectrum 
}


#endif
