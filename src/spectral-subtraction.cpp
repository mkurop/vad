#include <spectral-subtraction.h>
#include <fast-noise-psd-tracking.h>

arma::Col<float> SpectralSubtraction(arma::Col<float> &noisy_psd, arma::Col<float> &noise_psd, float subtraction_floor, float oversubtraction){

  arma::Col<float> enhanced_speech_magnitude_spectrum = \
      arma::sqrt(noisy_psd) - oversubtraction*arma::sqrt(noise_psd);

  enhanced_speech_magnitude_spectrum = amax(enhanced_speech_magnitude_spectrum,subtraction_floor);
   
  return enhanced_speech_magnitude_spectrum % enhanced_speech_magnitude_spectrum; // return power spectral density
}
