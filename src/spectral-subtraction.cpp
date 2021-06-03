#include <spectral-subtraction.h>

arma::Col<float> spectralSubtraction(arma::Col<float> &noisy_psd, arma::Col<float> &noisy_speech_angular_spectrum, \
                                     arma::Col<float> &noise_psd, float subtraction_floor, float oversubtraction){

  arma::Col<float> enhanced_speech_magnitude_spectrum = \
      arma::sqrt(noisy_psd) - oversubtraction*arma::sqrt(noise_psd)

  enhanced_speech_magnitude_spectrum.clamp(subtraction_floor,arma::max(enhanced_speech_magnitude_spectrum))

  return enhanced_speech_magnitude_spectrum 
}
