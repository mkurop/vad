#include <fast-noise-psd-tracking.h>

static float db2linpower(const float db) {

  float aux;

  aux = pow(10, db / 10.);

  return aux;
}

arma::Col<float> amax(const arma::Col<float>& V, float val)
{
	return clamp(V, val, std::max(V.max(), val));
}

arma::Col<float> amin(const arma::Col<float> &V, float val) {
  return clamp(V, std::min(V.min(), val), val);
}

arma::Row<float> amaxR(const arma::Row<float> &V, float val) {
  return clamp(V, val, std::max(V.max(), val));
}

arma::Row<float> aminR(const arma::Row<float> &V, float val) {
  return clamp(V, std::min(V.min(), val), val);
}

FastNoisePSDTracking::FastNoisePSDTracking(float sample_rate) : sample_rate(sample_rate) {
  state.frm_cntr = 0;
  queue = NDeque<arma::Col<float>>(3);
}

void FastNoisePSDTracking::setPrevFrameSpeechPSDEstimate(const arma::Col<float> &speech_psd) {
  speech_psd_set_flag = true;
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
  if (!speech_psd_set_flag) {
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

  speech_psd_set_flag = false;
}

arma::Col<float> FastNoisePSDTracking::smoothGamma(){
  arma::Col<float> aux = arma::zeros<arma::Col<float>>(K/2+1);
  arma::Col<float> out = arma::zeros<arma::Col<float>>(K/2+1);
  for(auto it = queue.begin(); it != queue.end(); ++it){
    aux += 1./queue.size() * *it;
  }
  for(int i = 0; i < K/2+1; i++){
    int l = ((i-deltak) < 0) ? 0 : i-deltak;
    int h = ((i+deltak) > K/2) ? K/2 : i+deltak;
    out[i] = arma::mean( aux.rows(l,h) );
  }
  return out;
}

void FastNoisePSDTracking::setPsi() {
  psi = arma::zeros<arma::Col<float>>(K/2+1);
  for(size_t i = 0; i < (K/2+1); i++){
    float f = i/(float)K*sample_rate;
    if(f < 1000.){
      psi[i] = 5.;
      continue;
    }
    if((f >= 1000.) && (f < 3000.)){
      psi[i] = 6.5;
      continue;
    }
    if(f > 3000.){
      psi[i] = 8.;
      continue;
    }
  }
}

void FastNoisePSDTracking::setP() {
  p = arma::zeros<arma::Col<float>>(K/2+1);
}
