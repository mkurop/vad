#include <sohnvad.h>
#include <vad.h>
#include <hamming.h>
#include <fast-noise-psd-tracking.h>
#include <spectral-subtraction.h>

t2r::Vad::Vad(int sampling_rate) : sampling_rate_(sampling_rate), fast_noise_tracking_(sampling_rate){
}

double t2r::Vad::operator()(arma::Col<float> signal_frame){

   if (start_){
    frame_length_ = signal_frame.n_rows;
    Vad_ vad_(sampling_rate_, frame_length_, frame_length_);
    window_ = hamming(frame_length_);
   } 
  
  // Dither
  arma::arma_rng::set_seed(1000);
  arma::Col<float> dither = arma::randn<arma::Col<float>>(signal_frame.size());
  dither *= 1.e-7;
  signal_frame += dither;

  // Window frame
  signal_frame %= window_;


  // Power spectrum
  arma::Col<std::complex<float>> signal_frame_spectrum = arma::fft(signal_frame);
  arma::Col<float> power_density_spectrum_of_signal_frame = (arma::real(signal_frame_spectrum % arma::conj(signal_frame_spectrum)));

  power_density_spectrum_of_signal_frame = power_density_spectrum_of_signal_frame.head_rows(frame_length_/2+1);

  // Noise PSD
  arma::Col<float> power_density_spectrum_of_noise(frame_length_/2+1);

  
  // Previous frame clean speech PSD
  arma::Col<float> previous_frame_clean_speech_power_density_spectrum(frame_length_/2+1);

  if(start_){
    fast_noise_tracking_.setPrevFrameSpeechPSDEstimate(power_density_spectrum_of_signal_frame);
    start_ = false;
  }else{
    fast_noise_tracking_.setPrevFrameSpeechPSDEstimate(previous_frame_clean_speech_power_density_spectrum);
  }  

  fast_noise_tracking_.noisePowRunning(power_density_spectrum_of_signal_frame);

  power_density_spectrum_of_noise = fast_noise_tracking_.getNoisePsd();

  previous_frame_clean_speech_power_density_spectrum = SpectralSubtraction(power_density_spectrum_of_signal_frame,power_density_spectrum_of_noise,0.01,2.0);

  double speech_presence_probability = vad_.vad(power_density_spectrum_of_signal_frame.t(), power_density_spectrum_of_noise.t());

  return speech_presence_probability;
}
