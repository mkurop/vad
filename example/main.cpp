#include "spectral-subtraction.h"
#include <fast-noise-psd-tracking.h>
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

  Vad vad(audioFile.getSampleRate(), kFrame, kFrame);
  

  // VAD decision file
  std::ofstream vad_decision_text_file;

  vad_decision_text_file.open("../../data/vad.dat");

  while(start + kFrame < audioFile.getNumSamplesPerChannel()){

    arma::Col<float> signal_frame = GetFrame(audioFile,start,kFrame);

    // Dither
    arma::Col<float> dither = arma::randn<arma::Col<float>>(signal_frame.size());
    dither *= 1.e-7;
    signal_frame += dither;

    // Window frame
    signal_frame %= window;

    // Power spectrum
    arma::Col<std::complex<float>> signal_frame_spectrum = arma::fft(signal_frame);
    arma::Col<float> power_density_spectrum_of_signal_frame = (arma::real(signal_frame_spectrum % arma::conj(signal_frame_spectrum)));

    power_density_spectrum_of_signal_frame = power_density_spectrum_of_signal_frame.head_rows(kFrame/2+1);
    
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

    std::cout << spprb << std::endl;

    vad_decision_text_file << (spprb > 0.5 ? 1 : 0) << std::endl;

    start += kFrame;
  }

  vad_decision_text_file.close();


  return 0;
}


