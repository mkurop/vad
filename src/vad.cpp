#include <vad.h>

Vad_::Vad_(int samplingRate, int frame, int frame2, double a01, double a10) : start(true) {
  this->mag_len = frame2 / 2 + 1;
  this->ksi_min = pow(10.0, -250.0 / 10.0);
  this->a01 = a01;
  this->a10 = a10;
  this->a00 = (LogProbability(1.) - LogProbability::fromLog(a01)).getLog(); // logAddExp1minus1(0, a01);
  this->a11 = (LogProbability(1.) - LogProbability::fromLog(a10)).getLog(); // logAddExp1minus1(0, a10);
  this->aa = 0.9;
  this->X_prev = arma::zeros<arma::Row<float>>(1, mag_len);
}

double Vad_::vad(const arma::Row<float> &X, const arma::Row<float> &lambdaN) {
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
    std::string what = "Speech presence probability is NaN.";
    throw std::runtime_error(what);
  }

  return spprb;
}
