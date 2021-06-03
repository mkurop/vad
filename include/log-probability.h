#ifndef LogProbability_INCLUDED
#define LogProbability_INCLUDED

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <random>


inline double logAddExp(double a, double b) {
  double max_exp = (a > b ? a : b);
  double sum = exp(a - max_exp) + exp(b - max_exp);
  return log(sum) + max_exp;
}

static double logProbAdd(double logProb1, double logProb2) {
  if (std::isnan(logProb1) || std::isnan(logProb2))
    throw std::runtime_error("logProbAdd - nan");

  if (logProb1 == std::numeric_limits<double>::lowest())
    return logProb2;
  else
    return logAddExp(logProb2, logProb1);
}
static double logProbSub(double logProb1, double logProb2) {
  if (std::isnan(logProb1) || std::isnan(logProb2))
    throw std::runtime_error("logProbSub - nan");

  if (logProb1 > logProb2) {
    return logProb1 + std::log1p(-std::exp(logProb2 - logProb1));
  } else if (logProb1 == logProb2) {
    return std::numeric_limits<double>::lowest();
  } else {
    LOG(error) << "Probability < 0.";
    throw (std::underflow_error("logProb1ability cannot store values below 0"));
  }
}
static double logProbMul(double logProb1, double logProb2) {
  if (std::isnan(logProb1) || std::isnan(logProb2)) {
    std::cout << "logProb1 = " << logProb1 << std::endl;
    std::cout << "logProb2 = " << logProb2 << std::endl;
    throw std::runtime_error("logProbMul - nan");
  }

  if (logProb2 == std::numeric_limits<double>::lowest())
    return std::numeric_limits<double>::lowest();
  return logProb1 + logProb2;
}
static double logProbDiv(double logProb1, double logProb2) {
  if (std::isnan(logProb1) || std::isnan(logProb2))
    throw std::runtime_error("logProbDiv - nan");

  if (logProb2 == std::numeric_limits<double>::lowest()) {
    std::string what{POSITION + ": divide by linear 0. detected"};
    throw std::runtime_error(what);
  }
  return logProb1 - logProb2;
}
static double logProbPow(double logProb, double r) {
  if (std::isnan(logProb))
    throw std::runtime_error("logProbPow - nan");

  return logProb * r;
}

/**
 * Allows storing values in their logarithmic states but using them like linear values
 */
class LogProbability {
 private:
  double logProb = std::numeric_limits<double>::lowest();

 public:
  LogProbability() = default;
  LogProbability(double probability) {
    if (probability > 0.) logProb = std::log(probability);
    else if (probability < 0.) {
      std::string what{POSITION + ": LogProbability cannot store values below 0, like " + std::to_string(probability)};
      throw std::runtime_error(what);
    }
  }

  bool sanityCheck() { return std::isnan(logProb) || std::isinf(logProb); }

  explicit operator float() const { return (float) std::exp(logProb); }
  operator double() const { return std::exp(logProb); }

  // operators
  LogProbability &operator+=(LogProbability r) { return logProb = logProbAdd(logProb, r.logProb), *this; }
  LogProbability &operator-=(LogProbability r) { return logProb = logProbSub(logProb, r.logProb), *this; }
  LogProbability &operator*=(LogProbability r) { return logProb = logProbMul(logProb, r.logProb), *this; }
  LogProbability &operator/=(LogProbability r) { return logProb = logProbDiv(logProb, r.logProb), *this; }
  LogProbability &operator^=(double r) { return logProb = logProbPow(logProb, r), *this; }

  friend LogProbability operator+(LogProbability lhs, LogProbability rhs) { return LogProbability(lhs) += rhs; }
  friend LogProbability operator+(LogProbability lhs, double rhs) { return lhs += LogProbability(rhs); }
  friend LogProbability operator+(LogProbability lhs, int rhs) { return lhs += LogProbability(rhs); }
  friend LogProbability operator+(double lhs, LogProbability rhs) { return LogProbability(lhs) += rhs; }
  friend LogProbability operator+(int lhs, LogProbability rhs) { return LogProbability(lhs) += rhs; }

  friend LogProbability operator-(LogProbability lhs, LogProbability rhs) { return LogProbability(lhs) -= rhs; }
  friend LogProbability operator-(LogProbability lhs, double rhs) { return lhs -= LogProbability(rhs); }
  friend LogProbability operator-(LogProbability lhs, int rhs) { return lhs -= LogProbability(rhs); }
  friend LogProbability operator-(double lhs, LogProbability rhs) { return LogProbability(lhs) -= rhs; }
  friend LogProbability operator-(int lhs, LogProbability rhs) { return LogProbability(lhs) -= rhs; }

  friend LogProbability operator*(LogProbability lhs, LogProbability rhs) { return LogProbability(lhs) *= rhs; }
  friend LogProbability operator*(LogProbability lhs, double rhs) { return lhs *= LogProbability(rhs); }
  friend LogProbability operator*(LogProbability lhs, int rhs) { return lhs *= LogProbability(rhs); }
  friend LogProbability operator*(double lhs, LogProbability rhs) { return LogProbability(lhs) *= rhs; }
  friend LogProbability operator*(int lhs, LogProbability rhs) { return LogProbability(lhs) *= rhs; }

  friend LogProbability operator/(LogProbability lhs, LogProbability rhs) { return LogProbability(lhs) /= rhs; }
  friend LogProbability operator/(LogProbability lhs, double rhs) { return lhs /= LogProbability(rhs); }
  friend LogProbability operator/(LogProbability lhs, int rhs) { return lhs /= LogProbability(rhs); }
  friend LogProbability operator/(double lhs, LogProbability rhs) { return LogProbability(lhs) /= rhs; }
  friend LogProbability operator/(int lhs, LogProbability rhs) { return LogProbability(lhs) /= rhs; }

  friend LogProbability operator^(LogProbability lhs, double rhs) { return lhs ^= rhs; }

  bool close(LogProbability o, double precision = 1.e-4) {
    auto maxXY = std::max(std::fabs(o.logProb), std::fabs(logProb));
    return std::fabs(o.logProb - logProb) / maxXY <= precision;
  }

  friend bool operator==(LogProbability lhs, LogProbability rhs) { return lhs.logProb == rhs.logProb; }
  friend bool operator==(LogProbability lhs, double rhs) { return lhs == LogProbability(rhs); }
  friend bool operator==(LogProbability lhs, int rhs) { return lhs == LogProbability(rhs); }
  friend bool operator==(double lhs, LogProbability rhs) { return LogProbability(lhs) == rhs; }
  friend bool operator==(int lhs, LogProbability rhs) { return LogProbability(lhs) == rhs; }

  friend bool operator!=(LogProbability lhs, LogProbability rhs) { return lhs.logProb != rhs.logProb; }
  friend bool operator!=(LogProbability lhs, double rhs) { return lhs != LogProbability(rhs); }
  friend bool operator!=(LogProbability lhs, int rhs) { return lhs != LogProbability(rhs); }
  friend bool operator!=(double lhs, LogProbability rhs) { return LogProbability(lhs) != rhs; }
  friend bool operator!=(int lhs, LogProbability rhs) { return LogProbability(lhs) != rhs; }

  friend bool operator>(LogProbability lhs, LogProbability rhs) { return lhs.logProb > rhs.logProb; }
  friend bool operator>(LogProbability lhs, double rhs) { return lhs > LogProbability(rhs); }
  friend bool operator>(LogProbability lhs, int rhs) { return lhs > LogProbability(rhs); }
  friend bool operator>(double lhs, LogProbability rhs) { return LogProbability(lhs) > rhs; }
  friend bool operator>(int lhs, LogProbability rhs) { return LogProbability(lhs) > rhs; }

  friend bool operator>=(LogProbability lhs, LogProbability rhs) { return lhs.logProb >= rhs.logProb; }
  friend bool operator>=(LogProbability lhs, double rhs) { return lhs >= LogProbability(rhs); }
  friend bool operator>=(LogProbability lhs, int rhs) { return lhs >= LogProbability(rhs); }
  friend bool operator>=(double lhs, LogProbability rhs) { return LogProbability(lhs) >= rhs; }
  friend bool operator>=(int lhs, LogProbability rhs) { return LogProbability(lhs) >= rhs; }

  friend bool operator<(LogProbability lhs, LogProbability rhs) { return lhs.logProb < rhs.logProb; }
  friend bool operator<(LogProbability lhs, double rhs) { return lhs < LogProbability(rhs); }
  friend bool operator<(LogProbability lhs, int rhs) { return lhs < LogProbability(rhs); }
  friend bool operator<(double lhs, LogProbability rhs) { return LogProbability(lhs) < rhs; }
  friend bool operator<(int lhs, LogProbability rhs) { return LogProbability(lhs) < rhs; }

  friend bool operator<=(LogProbability lhs, LogProbability rhs) { return lhs.logProb <= rhs.logProb; }
  friend bool operator<=(LogProbability lhs, double rhs) { return lhs <= LogProbability(rhs); }
  friend bool operator<=(LogProbability lhs, int rhs) { return lhs <= LogProbability(rhs); }
  friend bool operator<=(double lhs, LogProbability rhs) { return LogProbability(lhs) <= rhs; }
  friend bool operator<=(int lhs, LogProbability rhs) { return LogProbability(lhs) <= rhs; }

  // utils
  double get(bool logarithmic = true) const {
    if (logarithmic) return logProb;
    return std::exp(logProb);
  };

  double getLog() const { return logProb; }
  double getLinear() const { return std::exp(logProb); }

  void setLog(double logPrb) {
    logProb = logPrb;
  }

  void randDensity() {
    rand(0., 10.);
  }

  void randProbability() {
    rand(0., 1.);
  }

  void rand(double begin, double end) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(begin, end);
    auto aux = dist(gen);
    operator=(aux);
  }

  std::string strLog() const { return std::to_string(logProb); }
  std::string strLinear() const { return std::to_string(std::exp(logProb)); }

  std::ostream &save(std::ostream &os) const { return oswrite(os, logProb); }
  std::istream &load(std::istream &is) { return isread(is, logProb); }

  friend std::ostream &operator<<(std::ostream &os, const LogProbability &v) { return v.save(os); }
  friend std::istream &operator>>(std::istream &is, LogProbability &v) { return v.load(is); }

  // factory
  static LogProbability fromLinear(double srcLinear) {
    return LogProbability(srcLinear);
  }

  static LogProbability fromLog(double srcLog) {
    LogProbability lp;
    lp.logProb = srcLog;
    return lp;
  }
};

#endif
