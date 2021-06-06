#include "spectral-subtraction.h"
#include <algorithm>
#include <fast-noise-psd-tracking.h>
#include <sohnvad.h>
#include <iostream>
#include <ostream>
#include <vad.h>
#include <hamming.h>
#include <AudioFile.h> 
#include <armadillo>

const int kFrame = 400;
const int kStep = 200;

arma::Col<float> GetFrame(AudioFile<float> &audioFile, int start, int frame){

  arma::Col<float> signal_frame(frame);

  const int kChannel = 0;

  for(int i = start; i < start+frame; i++){
    signal_frame[i-start] = audioFile.samples[kChannel][i];
  }

  return signal_frame;

}

int main(int argc, char *argv[]){

  arma::Col<float> window = hamming(kFrame);

  AudioFile<float> audioFile;

  audioFile.load("../../data/SI2290.wav");

  int start = 0;

  FastNoisePSDTracking noiseTracking(audioFile.getSampleRate());

  Vad_ vad(audioFile.getSampleRate(), kFrame, kFrame);

  t2r::Vad sohnvad(audioFile.getSampleRate());


  // VAD decision file
  std::ofstream vad_decision_text_file;

  vad_decision_text_file.open("../../data/vad.dat");

  std::vector<float>  vad_decisions{};

  // Vector of signal_frame powers
  std::vector<float> powers_of_the_signal_frame{};

  while(start + kFrame < audioFile.getNumSamplesPerChannel()){

    arma::Col<float> signal_frame_raw = GetFrame(audioFile,start,kFrame);

    // Dither
    arma::arma_rng::set_seed(1000);
    arma::Col<float> dither = arma::randn<arma::Col<float>>(signal_frame_raw.size());
    dither *= 1.e-7;
    arma::Col<float> signal_frame = signal_frame_raw + dither;

    // Window frame
    signal_frame %= window;

    // Power spectrum
    arma::Col<std::complex<float>> signal_frame_spectrum = arma::fft(signal_frame);
    arma::Col<float> power_density_spectrum_of_signal_frame = (arma::real(signal_frame_spectrum % arma::conj(signal_frame_spectrum)));

    power_density_spectrum_of_signal_frame = power_density_spectrum_of_signal_frame.head_rows(kFrame/2+1);

    // Signal frame power
    double power_of_the_signal_frame = arma::sum(signal_frame%signal_frame);
    powers_of_the_signal_frame.push_back(power_of_the_signal_frame);
    
    // Noise PSD
    arma::Col<float> power_density_spectrum_of_noise(kFrame/2+1);

    // Previous frame clean speech PSD
    arma::Col<float> previous_frame_clean_speech_power_density_spectrum(kFrame/2+1);
    
    if(start == 0){
      noiseTracking.setPrevFrameSpeechPSDEstimate(power_density_spectrum_of_signal_frame);
    }else{
      noiseTracking.setPrevFrameSpeechPSDEstimate(previous_frame_clean_speech_power_density_spectrum);
    }  

    noiseTracking.noisePowRunning(power_density_spectrum_of_signal_frame);
    power_density_spectrum_of_noise = noiseTracking.getNoisePsd();
    previous_frame_clean_speech_power_density_spectrum = SpectralSubtraction(power_density_spectrum_of_signal_frame,power_density_spectrum_of_noise,0.01,2.0);

    double spprb = vad.vad(power_density_spectrum_of_signal_frame.t(), power_density_spectrum_of_noise.t());

    vad_decisions.push_back(spprb);

    std::cout << "speech probability: " << spprb << ", " << sohnvad(signal_frame_raw) << std::endl;


    start += kFrame;
  }

  float max_power = *std::max_element(powers_of_the_signal_frame.begin(),powers_of_the_signal_frame.end());
  // for(int i = 0; i < vad_decisions.size(); i++){
  //   vad_decision_text_file << vad_decisions[i] << ", "<<  powers_of_the_signal_frame[i]/max_power*1000. << std::endl;
  // }
  vad_decision_text_file.close();


  return 0;
}


