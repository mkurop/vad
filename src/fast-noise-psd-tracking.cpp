#include <fast-noise-psd-tracking.h>

static float db2linpower(const float db) {

  float aux;

  aux = pow(10, db / 10.);

  return aux;
}


FastNoisePSDTracking::FastNoisePSDTracking(float sample_rate) : sample_rate(sample_rate) {
  state.frmCntr = 0;
  queue = NDeque<arma::Col<float>>(3);
}
void FastNoisePSDTracking::setPrevFrameSpeechPSDEstimate(const arma::Col<float> &speech_psd) {
  speech_psd = true;
  this->speech_psd = speech_psd;
}
void FastNoisePSDTracking::noisePowRunning(const arma::Col<float> &noisy_per) {
  if (state.frm_cntr == 0) {
    state.frm_cntr++;
    state.noise_pow = noisy_per;
    K = noisy_per.n_rows * 2 - 2;
    setPsi();
    setP();
    return;
  }
  if (!speech_psd) {
    throw std::logic_error(
        "Speech spectrum not set, program exits ... run first the method "
        "setPrevFrameSpeechPSDEstimate ...");
  }
  //< compute a posteriori SNR
  arma::Col<float> gamma = noisy_per / state.noise_pow;
  //< compute a priori SNR
  arma::Col<float> ksi = amax(alpha_ns * speech_psd / state.noise_pow +
                                  (1 - alpha_ns) * (gamma - 1),
                              db2linpower(-15.));
  //< compute noise estimate
  arma::Col<float> aux1 = gamma / (ksi % (ksi + 1));
  //< compute gamma
  arma::Col<float> aux2 = arma::zeros<arma::Col<float>>(K / 2 + 1);
  for (size_t k = 0; k < (K / 2 + 1); k++) {
    //      double aux = boost::math::tgamma(1e-15,aux1[k]);
    //      double v = aux1[k];
    aux2[k] = std::exp(
        std::min(boost::math::tgamma(1e-15, aux1[k]), std::log(1e5)));
  }
  arma::Col<float> aux3 = 1 / (1 + ksi);
  arma::Col<float> G = aux3 % aux3 % aux2;
  arma::Col<float> N = G % noisy_per;
  //< smoothing of gamma
  queue.push_back(gamma);
  arma::Col<float> gamma_smoothed = smoothGamma();

  //< the I variable (speech presence)
  arma::Col<float> I =
      arma::conv_to<arma::Col<float>>::from(gamma_smoothed > psi);

  //< compute speech presence probability
  p = alpha_p * p + (1 - alpha_p) * I;

  //< compute noise estimate
  state.noise_pow = p % state.noise_pow +
                   (1 - p) % (alpha_n * state.noise_pow + (1 - alpha_n) * N);

  speech_psd = false;
}

const arma::Col<float> FastNoisePSDTracking::getNoisePsd() { return state.noisePow; }

const FastNsePowState &FastNoisePSDTracking::getState() { return state; }
