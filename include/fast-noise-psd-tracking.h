//
// Created by marcin on 31.01.2020.
//

#ifndef ZOOMER_NOISEPSDTRACKING_FASTNOISEPSDTRACKING_H_
#define ZOOMER_NOISEPSDTRACKING_FASTNOISEPSDTRACKING_H_

#include <ndeque.h>
#include <algorithm>
#include <armadillo>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <cstdint>

class FastNsePowState {

public:

  int frm_cntr; 
  arma::Col<float> noise_pow, ph1mean;
  
};

static float db2linpower(const float db); 

arma::Col<float> amax(const arma::Col<float>& V, float val);

arma::Col<float> amin(const arma::Col<float> &V, float val); 

arma::Row<float> amaxR(const arma::Row<float> &V, float val); 

arma::Row<float> aminR(const arma::Row<float> &V, float val); 

//! Class for tracking the power spectral density of the noise contaminatig the speech signal
/*!
 * The implementation is based on the article:
 * Qiquan Zhang et al. "Fast Nonstationary Noise Tracking Based on Log-Spectral \
 * Power MMSE Estimator and Temporal Recursive Averaging", Open Access Article
 *
 * \author Marcin Kuropatwiński
 * \date 02.06.2021
 * \copyright MIT licence
 */

class FastNoisePSDTracking {

private:

  const float alpha_n = 0.8;
  const float alpha_ns = 0.98;
  const float alpha_p = 0.2;
  const int32_t deltak = 1;
  const uint32_t deltal = 2;
  uint32_t K{0}; //< length of the the spectrum
  arma::Col<float> psi;
  arma::Col<float> p;
  arma::Col<float> speech_psd;
  bool speech_psd_set_flag = false;
  NDeque<arma::Col<float>> queue;
  float sample_rate;
  FastNsePowState state;

public:
  //! A constructor.
  /*!
   * The constructor for the object performing the noise tracking.
   *
   * \param[in] sample_rate sampling rate in Herz
   */
  FastNoisePSDTracking(float sample_rate);

  //! Method for setting the previous speech power density spectrum (PSD)   
  /*!
   * The argument to this method is PSD of the speech esimated in previous frame
   *
   *  \param[in] armadillo vector containing the speech PSD 
   */
  void setPrevFrameSpeechPSDEstimate(const arma::Col<float> &speech_psd); 

  //! Method running the main noise tracing procedure
  /*!
   * Takes the noisy speech PSD and stores the estimated noise PSD in the member variable.
   *
   * \param[in] the noisy speech PSD of lenght K/2+1 where K is the frame length
   */

  void noisePowRunning(const arma::Col<float> &noisy_per); 

  //! Getter for the estimated noise PSD
  /*!
   * Returns the estimate of the noise
   *
   * \return the noise PSD in the current frame
   */

  const arma::Col<float> getNoisePsd() { return state.noise_pow; }


private:
  
  arma::Col<float> smoothGamma();
    
  const FastNsePowState &getState() { return state; }

  void setPsi(); 

  void setP();
};

#endif // ZOOMER_NOISEPSDTRACKING_FASTNOISEPSDTRACKING_H_
