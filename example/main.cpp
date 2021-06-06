#include <sohnvad.h>
#include <AudioFile.h> 
#include <armadillo>

const int kFrame = 400;
const int kStep = kFrame;

/*!
 * \brief Form the frame from the input speech file
 *
 * \param[in] audioFile AudioFile object
 * \param[in] start beginning sample of the frame
 * \param[in] frame the frame length in samples
 * \return signal_frame Armadillo column vector containing the samples in the signal frame
 *
 */
arma::Col<float> GetFrame(const AudioFile<float> &audioFile, const int start, const int frame){

  arma::Col<float> signal_frame(frame);

  const int kChannel = 0;

  for(int i = start; i < start+frame; i++){
    signal_frame[i-start] = audioFile.samples[kChannel][i];
  }

  return signal_frame;

}

int main(int argc, char *argv[]){
  
  //!< Declare the audioFile object
  AudioFile<float> audioFile;

  // !< Load test file 
  audioFile.load("../../data/SI2290.wav");

  int start = 0;

  t2r::Vad sohnvad(audioFile.getSampleRate());

  while(start + kFrame < audioFile.getNumSamplesPerChannel()){

    arma::Col<float> signal_frame_raw = GetFrame(audioFile,start,kFrame);

    double speech_signal_presence_probability = sohnvad(signal_frame_raw);

    std::cout << "speech probability: " << speech_signal_presence_probability << std::endl;

    start += kStep;
  }

  return 0;

}


