#ifndef VAD_H_INC
#define VAD_H_INC

#include <armadillo>
#include <cmath>

#include <log-probability.h>
#include <fast-noise-psd-tracking.h>
#include <string>

/*!
	* \brief Class implements the VAD algorithm from the paper:
	*
	* Jongseo Sohn, Nam Soo Kim, Wonyong Sung,
	* A Statistical Model-Based Voice Activity Detection,
	* IEEE Signal Processing Letters, vol. 6., no. 1, January 1999
	*/


class Vad {
  bool start;
  int mag_len;
  arma::Row<float> X_prev;
  double ksi_min, a01, a10, a00, a11, aa;
  NsePow np;
  FilterBankEnergies fbe;

 public:
  double logsp, lognp, spprb;

  /*!
   * \brief Class implements the VAD algorithm
   *
   * \param samplingRate sampling rate in [Hz], defaults to 16000Hz
   * \param frame frame length before zero adding, defaults to 400 samples
   * \param frame2 frame length after zero adding, defaults to 512 samples
   * \param a01 log natural base probability of noise -> speech transition
   * \param a10 log natural base probability of speech -> noise transition
   */
  Vad(int samplingRate = 16000, int frame = 400, int frame2 = 512, double a01 = -.3, double a10 = -.3) : fbe(
      samplingRate,
      frame,
      frame2), start(true), np(frame2) {
    this->mag_len = frame2 / 2 + 1;
    this->ksi_min = pow(10.0, -250.0 / 10.0);
    this->a01 = a01;
    this->a10 = a10;
    this->a00 = (LogProbability(1.) - LogProbability::fromLog(a01)).getLog(); // logAddExp1minus1(0, a01);
    this->a11 = (LogProbability(1.) - LogProbability::fromLog(a10)).getLog(); // logAddExp1minus1(0, a10);
    this->aa = 0.9;
    this->X_prev = arma::zeros<arma::Row<float>>(1, mag_len);
  }

  /*!
   * \brief Method for speech presence probability
   *
   * \param X the noisy speech PSD
   * \param lambdaN the estimated noise PSD
   *
   * \return speech presence probability
   */
  double vad(const arma::Row<float> &X, const arma::Row<float> &lambdaN) {
    arma::Row<float> _lambdaN = amaxR(lambdaN, 1e-3);

    arma::Row<float> gammak = aminR(X / _lambdaN, 1000.0);

    arma::Row<float> ksi;
    if (start) {
      ksi = aa + (1 - aa) * amaxR(gammak - 1, 0.0);
    } else {
      ksi = aa * aminR(X_prev / _lambdaN, 100.0) + (1 - aa) * amaxR(gammak - 1, 0.0);
      ksi = amaxR(ksi, ksi_min);
    }

    arma::Row<float> l = gammak % ksi / (1 + ksi) - log(1 + ksi);

    X_prev = X;

    _lambdaN(0) *= .5;
    _lambdaN(_lambdaN.n_cols - 1) *= .5;
    gammak(0) *= .5;
    gammak(gammak.n_cols - 1) *= .5;
    X_prev(0) *= .5;
    X_prev(X_prev.n_cols - 1) *= .5;

    l(0) *= .5;
    l(l.n_cols - 1) *= .5;

    double lpXH0 = (2. * sum(-log(M_PI * _lambdaN) - gammak));
    double lpXH1dpXH0 = (2. * sum(l));
    double lpXH1 = lpXH1dpXH0 + lpXH0;

    if (start) {
      logsp = lpXH1;
      lognp = lpXH0;
      start = false;
    } else {
      double oldLogSP = logsp;
      logsp = logAddExp(logsp + a11, lognp + a01) + lpXH1;
      lognp = logAddExp(oldLogSP + a10, lognp + a00) + lpXH0;
    }

    spprb = (LogProbability::fromLog(logsp) / (LogProbability::fromLog(logsp) + LogProbability::fromLog(lognp)))
        .getLinear();

    if (isnan(spprb)) {
      std::string what = "Speech presence probability is NaN."
      throw std::runtime_error(what);
    }

    return spprb;
  }
};

#endif
