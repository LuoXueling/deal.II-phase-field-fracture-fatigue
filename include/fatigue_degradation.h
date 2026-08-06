//
// Created by xlluo on 24-8-5.
//

#ifndef CRACKS_FATIGUE_DEGRADATION_H
#define CRACKS_FATIGUE_DEGRADATION_H

#include "controller.h"
#include "dealii_includes.h"
#include <complex>
using namespace dealii;

template<int dim>
class FatigueAccumulation {
public:
  FatigueAccumulation(Controller<dim> &ctl) {
  };

  void step(const std::shared_ptr<PointHistory> &lqph_q, double phasefield,
            double degrade, double degrade_derivative,
            double degrade_second_derivative, Controller<dim> &ctl) {
    double increm = increment(lqph_q, phasefield, degrade, degrade_derivative,
                              degrade_second_derivative, ctl);
    lqph_q->update("Fatigue history", std::max(increm, 0.0), "accumulate");
    record(lqph_q, phasefield, degrade, degrade_derivative,
           degrade_second_derivative, ctl);
  }

  virtual double increment(const std::shared_ptr<PointHistory> &lqph,
                           double phasefield, double degrade,
                           double degrade_derivative,
                           double degrade_second_derivative,
                           Controller<dim> &ctl) {
    AssertThrow(false, ExcNotImplemented());
  };

  virtual void record(const std::shared_ptr<PointHistory> &lqph,
                      double phasefield, double degrade,
                      double degrade_derivative,
                      double degrade_second_derivative, Controller<dim> &ctl) {
    // If one is using the latest fatigue history to record something, it's
    // necessary to record it in this function, because if solving the phase
    // field first and then the displacement, and the multistep staggered is not
    // triggered, the record in increment() is still using the old fatigue
    // variable updated by old elastic energy. The issue happens in VHCF.
  };
};

/**
 * The default fatigue threshold, Gc/(12*l_phi).
 *
 * This is the value Carrara Eq. 47 reduces to once the extra 0.5 needed to
 * reproduce the published results is folded in:
 *   0.25 * E * eps^2  with  eps^2 = Gc/(3*l_phi*E)   ==>   Gc/(12*l_phi)
 * It is shared by the Carrara/Kristensen asymptotic degradations (alpha_t) and
 * by CarraraMeanEffectAccumulation (alpha_n), which historically each spelled
 * the same formula out separately.
 */
template<int dim>
inline double default_fatigue_alpha_t(Controller<dim> &ctl) {
  double epsilon_at2 =
      std::sqrt(ctl.params.Gc / (3 * ctl.params.l_phi * ctl.params.E));
  return 0.5 * 0.5 * epsilon_at2 * ctl.params.E * epsilon_at2;
}

/**
 * Resolves the fatigue threshold for one consumer, in priority order:
 *   1. the global "Fatigue alpha_t" entry, if set
 *   2. the scheme's own parameter string, if it supplies a value
 *   3. `fallback` -- the consumer's historical hardcoded formulation
 *
 * The global deliberately outranks the per-scheme strings. Its whole purpose is
 * to hold every consumer (degradation alpha_t, CarraraMeanEffect alpha_n, and
 * JonasCycleJump's slot) at one common threshold; if a per-scheme value could
 * override it, setting it would silently leave those consumers disagreeing,
 * which is the inconsistency it exists to prevent. Leaving it empty reproduces
 * the old behaviour exactly, whatever each consumer computed for itself.
 */
template<int dim>
inline double resolve_fatigue_alpha_t(Controller<dim> &ctl, double fallback,
                                      const std::string &own_parameters = "") {
  if (ctl.params.fatigue_alpha_t != "") {
    std::istringstream iss(ctl.params.fatigue_alpha_t);
    double value;
    AssertThrow(static_cast<bool>(iss >> value),
                ExcInternalError("'Fatigue alpha_t' is not a number: " +
                                 ctl.params.fatigue_alpha_t));
    return value;
  }
  if (own_parameters != "") {
    // Mirrors the original `iss >> alpha_t`: a non-empty string is taken as
    // authoritative, and the branch is entered on emptiness alone, not on
    // whether the parse succeeds. Since C++11 a failed extraction zeroes the
    // target, so an unparseable value gives 0 -- surfaced here as an error
    // rather than silently degrading with alpha_t = 0.
    std::istringstream iss(own_parameters);
    double value;
    AssertThrow(static_cast<bool>(iss >> value),
                ExcInternalError("Leading entry of '" + own_parameters +
                                 "' is not a number."));
    return value;
  }
  return fallback;
}

/**
 * Per-cycle increment laws.
 *
 * These are the elementary "how much fatigue history does one resolved cycle
 * add" rules, factored out so that an acceleration algorithm can be told which
 * one to use for its resolved-cycle branch via the "Fatigue increment"
 * parameter. The three differ in how they scale with stress amplitude, which is
 * what sets the Paris exponent:
 *   CarraraNoMeanEffect : dpsi*degrade            (~sigma^2 dsigma)
 *   Kristensen          : dpsi, undegraded        (~sigma^2 dsigma)
 *   CarraraMeanEffect   : dpsi*psi/alpha_n        (~sigma^4 dsigma)
 * Only positive increments contribute, matching the original inline code.
 */
inline double carrara_no_mean_effect_increment(
  const std::shared_ptr<PointHistory> &lqph, double degrade) {
  double dpsi = lqph->get_increment_latest("Positive elastic energy", 0.0) *
                degrade;
  return (dpsi > 0 ? 1.0 : 0.0) * dpsi;
}

inline double kristensen_increment(
  const std::shared_ptr<PointHistory> &lqph) {
  // Note: deliberately undegraded -- this is what distinguishes it from
  // CarraraNoMeanEffect.
  double dpsi = lqph->get_increment_latest("Positive elastic energy", 0.0);
  return (dpsi > 0 ? 1.0 : 0.0) * dpsi;
}

inline double carrara_mean_effect_increment(
  const std::shared_ptr<PointHistory> &lqph, double degrade, double alpha_n) {
  double dpsi = lqph->get_increment_latest("Positive elastic energy", 0.0) *
                degrade;
  double psi = lqph->get_latest("Positive elastic energy", 0.0) * degrade;
  return (dpsi > 0 ? 1.0 : 0.0) * dpsi * psi / alpha_n;
}

/**
 * alpha_n for CarraraMeanEffect-style increments selected via the "Fatigue
 * increment" parameter. This path cannot read "Fatigue accumulation
 * parameters" -- the host acceleration algorithm already owns that string --
 * so alpha_n comes from "Fatigue increment parameters" when given, else the
 * global "Fatigue alpha_t", else the shared default.
 *
 * Note the precedence differs from resolve_fatigue_alpha_t's usual ordering:
 * here the scheme's own string wins over the global, because "Fatigue
 * increment parameters" exists only to set alpha_n independently of the
 * degradation threshold. Leaving it empty preserves the previous behaviour.
 */
template<int dim>
inline double carrara_default_alpha_n(Controller<dim> &ctl) {
  // increment() calls this per quadrature point per cycle, so the string parse
  // is resolved exactly once and every later call returns the cached value.
  // Everything it depends on -- the parameter strings and the material
  // constants behind the default -- is fixed after construction, so a single
  // resolution is correct. The static initialiser is thread-safe (C++11): a
  // concurrent caller blocks until it completes, then reads the same value.
  static const double alpha_n = [&ctl]() {
    double value;
    std::string source;
    if (ctl.params.fatigue_increment_parameters != "") {
      std::istringstream iss(ctl.params.fatigue_increment_parameters);
      AssertThrow(static_cast<bool>(iss >> value),
                  ExcInternalError(
                    "Leading entry of 'Fatigue increment parameters' (alpha_n) "
                    "is not a number: " +
                    ctl.params.fatigue_increment_parameters));
      source = "Fatigue increment parameters";
    } else {
      value = resolve_fatigue_alpha_t(ctl, default_fatigue_alpha_t(ctl));
      source = (ctl.params.fatigue_alpha_t != "") ? "Fatigue alpha_t"
                                                  : "internal default";
    }
    // Printed once per run, from inside the one-time initialiser -- this is
    // also the marker that the cache resolved exactly once.
    ctl.dcout << "Using alpha_n: " << value << " from " << source << std::endl;
    return value;
  }();
  return alpha_n;
}

/**
 * Resolves the "Fatigue increment" parameter for an acceleration algorithm.
 *
 * `own` is the law the algorithm uses natively; it is returned for "Auto" so
 * that existing parameter files keep their exact present behaviour.
 */
template<int dim>
inline std::string resolve_fatigue_increment(const std::string &own,
                                             Controller<dim> &ctl) {
  const std::string &requested = ctl.params.fatigue_increment;
  return (requested == "" || requested == "Auto") ? own : requested;
}

/**
 * Evaluates the selected increment law. Shared by every acceleration algorithm
 * so they all interpret "Fatigue increment" identically.
 */
template<int dim>
inline double evaluate_fatigue_increment(
  const std::string &law, const std::shared_ptr<PointHistory> &lqph,
  double degrade, Controller<dim> &ctl) {
  if (law == "Kristensen")
    return kristensen_increment(lqph);
  if (law == "CarraraMeanEffect")
    return carrara_mean_effect_increment(lqph, degrade,
                                         carrara_default_alpha_n(ctl));
  if (law == "CarraraNoMeanEffect")
    return carrara_no_mean_effect_increment(lqph, degrade);
  AssertThrow(false,
              ExcInternalError("Unknown 'Fatigue increment' value: " + law));
  return 0.0;
}

/**
 * Guards the accumulations that are definitionally tied to one increment law.
 * CarraraNoMeanEffect/Kristensen/CarraraMeanEffect each accept only their own,
 * so a mismatched "Fatigue increment" is a hard error rather than a silently
 * ignored setting.
 */
template<int dim>
inline void assert_fixed_fatigue_increment(const std::string &own,
                                           Controller<dim> &ctl) {
  const std::string &requested = ctl.params.fatigue_increment;
  AssertThrow(requested == "" || requested == "Auto" || requested == own,
              ExcInternalError(
                own + "Accumulation only accepts the " + own +
                " increment, but 'Fatigue increment' is set to '" + requested +
                "'. Use an acceleration-algorithm accumulation (Cojocaru, Li, "
                "Jonas, Yang, Jaccon, ...) to select a different increment."));
}

template<int dim>
class CarraraNoMeanEffectAccumulation : public FatigueAccumulation<dim> {
public:
  CarraraNoMeanEffectAccumulation(Controller<dim> &ctl)
    : FatigueAccumulation<dim>(ctl) {
    assert_fixed_fatigue_increment("CarraraNoMeanEffect", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    return carrara_no_mean_effect_increment(lqph, degrade);
  };
};

template<int dim>
class KristensenAccumulation : public FatigueAccumulation<dim> {
public:
  KristensenAccumulation(Controller<dim> &ctl)
    : FatigueAccumulation<dim>(ctl) {
    assert_fixed_fatigue_increment("Kristensen", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    return kristensen_increment(lqph);
  };
};

template<int dim>
class KristensenCLAAccumulation : public FatigueAccumulation<dim> {
public:
  KristensenCLAAccumulation(Controller<dim> &ctl)
    : FatigueAccumulation<dim>(ctl) {
    AssertThrow(ctl.params.adaptive_timestep == "KristensenCLA",
                ExcInternalError("KristensenCLATimeStep must be used "
                  "with KristensenCLAAccumulation."));
    AssertThrow(
      ctl.params.fatigue_accumulation_parameters != "",
      ExcInternalError(
        "Parameters of KristensenCLAAccumulation is not assigned."));
    std::istringstream iss(ctl.params.fatigue_accumulation_parameters);
    iss >> R;
    AssertThrow(
      R >= 0 || (R < 0 && ctl.params.degradation == "hybridnotension"),
      ExcInternalError("Cannot use KristensenCLAAccumulation when "
        "R<0 while hybridnotension split is not used"));
    increment_law = resolve_fatigue_increment("CarraraNoMeanEffect", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    if (ctl.current_timestep != ctl.params.timestep_size_2) {
      return evaluate_fatigue_increment(increment_law, lqph, degrade, ctl);
    } else {
      double psi = lqph->get_latest("Positive elastic energy", 0.0) * degrade;
      double n_jump = ctl.get_info("N jump", 1);
      double increm = psi * (1 - R * R * (R >= 0 ? 1 : 0)) * n_jump;
      return increm;
    }
  };
  double R;
  std::string increment_law;
};

template<int dim>
class CojocaruAccumulation : public FatigueAccumulation<dim> {
public:
  CojocaruAccumulation(Controller<dim> &ctl) : FatigueAccumulation<dim>(ctl) {
    AssertThrow(ctl.params.fatigue_accumulation_parameters != "",
                ExcInternalError(
                  "Parameters of CojocaruCLAAccumulation is not assigned."));
    std::istringstream iss(ctl.params.fatigue_accumulation_parameters);
    iss >> R >> q_jump;
    increment_law = resolve_fatigue_increment("CarraraNoMeanEffect", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    int n_jumps = static_cast<int>(ctl.get_info("N jump", 0.0));
    double increm;
    if (n_jumps == 0 || ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      increm = evaluate_fatigue_increment(increment_law, lqph, degrade, ctl);
    } else {
      double s12 = lqph->get_initial("s12", 0.0);
      double s23 = lqph->get_initial("s23", 0.0);
      increm = s12 * n_jumps + (s12 - s23) * std::pow(n_jumps, 2) / 2.0;
    }
    return increm;
  };

  void record(const std::shared_ptr<PointHistory> &lqph, double phasefield,
              double degrade, double degrade_derivative,
              double degrade_second_derivative, Controller<dim> &ctl) {
    // Determine the number of jumps
    double subcycle = ctl.get_info("Subcycle", 0.0);
    if (std::abs(subcycle - 1) < 1e-8) {
      lqph->update("y3", lqph->get_latest("Fatigue history", 0.0));
    } else if (std::abs(subcycle - 2) < 1e-8) {
      double y3 = lqph->get_latest("y3", 0.0);
      double y2 = lqph->get_latest("Fatigue history", 0.0);
      lqph->update("y2", y2);
      lqph->update("s23", y2 - y3);
    } else if (std::abs(subcycle - 3) < 1e-8) {
      double s23 = lqph->get_latest("s23", 0.0);
      double y1 = lqph->get_latest("Fatigue history", 0.0);
      double y2 = lqph->get_latest("y2", 0.0);
      double s12 = y1 - y2;
      lqph->update("s12", s12);
      double max_jump = ctl.get_info("Maximum jump", 1.0e8);
      double n_jump_local = max_jump;
      if (phasefield < 0.95 && phasefield > 0.001 &&
          std::abs(s12 - s23) / std::abs(s12) > 1e-5) {
        n_jump_local = q_jump * s12 / std::abs(s12 - s23);
      }
      lqph->update("n_jump_local", n_jump_local);
    }
  }

  double R, q_jump;
  std::string increment_law;
};

template<int dim>
class CojocaruCLAAccumulation : public CojocaruAccumulation<dim> {
public:
  CojocaruCLAAccumulation(Controller<dim> &ctl)
    : CojocaruAccumulation<dim>(ctl) {
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    int n_jumps = static_cast<int>(ctl.get_info("N jump", 0.0));
    double increm;
    if (n_jumps == 0 || ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      if (ctl.current_timestep != ctl.params.timestep_size_2) {
        increm = evaluate_fatigue_increment(this->increment_law, lqph, degrade,
                                            ctl);
      } else {
        double psi = lqph->get_latest("Positive elastic energy", 0.0) * degrade;
        increm = psi * (1 - this->R * this->R * (this->R >= 0 ? 1 : 0));
      }
    } else {
      double s12 = lqph->get_initial("s12", 0.0);
      double s23 = lqph->get_initial("s23", 0.0);
      increm = s12 * n_jumps + (s12 - s23) * std::pow(n_jumps, 2) / 2.0;
    }
    return increm;
  };
};

template<int dim>
class LiAccumulation : public FatigueAccumulation<dim> {
public:
  LiAccumulation(Controller<dim> &ctl) : FatigueAccumulation<dim>(ctl) {
    AssertThrow(ctl.params.fatigue_accumulation_parameters != "",
                ExcInternalError(
                  "Parameters of LiAccumulation is not assigned."));
    std::istringstream iss(ctl.params.fatigue_accumulation_parameters);
    iss >> R >> chi_cr;
    increment_law = resolve_fatigue_increment("CarraraNoMeanEffect", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    int n_jumps = static_cast<int>(ctl.get_info("N jump", 0.0));
    double increm;
    if (n_jumps == 0 || ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      increm = evaluate_fatigue_increment(increment_law, lqph, degrade, ctl);
    } else {
      double y3 = lqph->get_initial("y3", 0.0);
      double y2 = lqph->get_initial("y2", 0.0);
      double y1 = lqph->get_initial("y1", 0.0);
      double y0 = lqph->get_initial("y0", 0.0);
      // Closed form of y_{k+1} = 11/4 y_k - 11/4 y_{k-1} + 5/4 y_{k-2} - 1/4
      // y_{k-3}, replacing an O(n_jumps) loop per quadrature point. Its roots
      // are a defective 1 (giving A + B n) and a pair of modulus 1/2 (the
      // transient), so y_n = A + B n + 2 Re(C lambda^n). Do not evaluate this by
      // matrix powers: the defective unit root amplifies its splitting error
      // linearly in n (~3% at n = 1e5).
      const double A_sec = y0 + 0.25 * y1 - 0.5 * y2 + 0.25 * y3;
      (void)A_sec; // cancels in y_n - y_0.
      const double B_sec = 2.0 * y0 - 3.5 * y1 + 2.0 * y2 - 0.5 * y3;
      increm = B_sec * static_cast<double>(n_jumps);
      // The transient decays as 2^-n; past ~60 jumps dropping it costs <1e-4
      // relative, so only small jumps pay for the complex arithmetic.
      if (n_jumps < 60) {
        const std::complex<double> lambda(0.375, 0.330718913883073824);
        const std::complex<double> C =
            std::complex<double>(0.0, 0.377964473009227227) * y0 +
            std::complex<double>(-0.125, -0.897665623396914665) * y1 +
            std::complex<double>(0.25, 0.661437827766147648) * y2 +
            std::complex<double>(-0.125, -0.14173667737846021) * y3;
        increm += 2.0 * (C * std::pow(lambda, n_jumps)).real() - 2.0 * C.real();
      }
    }
    return increm;
  };

  void record(const std::shared_ptr<PointHistory> &lqph, double phasefield,
              double degrade, double degrade_derivative,
              double degrade_second_derivative, Controller<dim> &ctl) {
    // Determine the number of jumps
    double subcycle = ctl.get_info("Subcycle", 0.0);
    double max_jump = ctl.get_info("Maximum jump", 1.0e8);
    if (std::abs(subcycle - 1) < 1e-8) {
      lqph->update("y3", lqph->get_latest("Fatigue history", 0.0));
    } else if (std::abs(subcycle - 2) < 1e-8) {
      lqph->update("y2", lqph->get_latest("Fatigue history", 0.0));
      lqph->update("phi2", lqph->get_latest("Phase field", 0.0));
    } else if (std::abs(subcycle - 3) < 1e-8) {
      lqph->update("y1", lqph->get_latest("Fatigue history", 0.0));
      lqph->update("phi1", lqph->get_latest("Phase field", 0.0));
    } else if (std::abs(subcycle - 4) < 1e-8) {
      lqph->update("y0", lqph->get_latest("Fatigue history", 0.0));
      double phi2 = std::max(lqph->get_initial("phi2", 0.0), 1e-10);
      double phi1 = std::max(lqph->get_initial("phi1", 0.0), 1e-10);
      double phi0 = std::max(lqph->get_latest("Phase field", 0.0), 1e-10);
      if (phi0 < 1e-1 && (phi0 <= phi1 * (1 + 1e-6) || phi1 <= phi2 * (1 + 1e-6))) {
        // The point is subject to minor numerical error, or they are not updated.
        lqph->update("n_jump_local", max_jump);
      } else {
        double n_jump_local = 1;
        if (phi0 - 2 * phi1 + phi2 < 0) {
          lqph->update("n_jump_local", max_jump);
        } else {
          // f(n) = (n+1) d1 + (n+1)^2 d2/2 is strictly increasing (d1 > 0 and
          // d2 >= 0 are guaranteed above), so the original O(n) scan for the
          // smallest n >= 1 with f(n) > rhs reduces to bisection.
          const double d1 = phi0 - phi1;
          const double d2 = phi0 - 2 * phi1 + phi2;
          const double rhs = d1 / (d1 + d2) * chi_cr * phi0;
          auto f = [d1, d2](double n) {
            return (n + 1) * d1 + (n + 1) * (n + 1) * 0.5 * d2;
          };
          if (f(1.0) > rhs) {
            n_jump_local = 1.0;
          } else {
            double lo = 1.0, hi = 2.0;
            bool capped = false;
            while (f(hi) <= rhs) {
              lo = hi;
              hi *= 2.0;
              if (hi >= max_jump) {
                n_jump_local = max_jump;
                capped = true;
                break;
              }
            }
            if (!capped) {
              // Invariant: f(lo) <= rhs < f(hi); converge on the smallest such n.
              while (hi - lo > 1.0) {
                double mid = std::floor((lo + hi) / 2.0);
                if (f(mid) <= rhs) lo = mid; else hi = mid;
              }
              n_jump_local = std::min(hi, max_jump);
            }
          }
          lqph->update("n_jump_local", n_jump_local);
        }
      }
    }
  }

  double R, chi_cr;
  std::string increment_law;
};

template<int dim>
class LiCLAAccumulation : public LiAccumulation<dim> {
public:
  LiCLAAccumulation(Controller<dim> &ctl)
    : LiAccumulation<dim>(ctl) {
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    int n_jumps = static_cast<int>(ctl.get_info("N jump", 0.0));
    double increm;
    if (n_jumps == 0 || ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      if (ctl.current_timestep != ctl.params.timestep_size_2) {
        increm = evaluate_fatigue_increment(this->increment_law, lqph, degrade,
                                            ctl);
      } else {
        double psi = lqph->get_latest("Positive elastic energy", 0.0) * degrade;
        increm = psi * (1 - this->R * this->R * (this->R >= 0 ? 1 : 0));
      }
    } else {
      double y3 = lqph->get_initial("y3", 0.0);
      double y2 = lqph->get_initial("y2", 0.0);
      double y1 = lqph->get_initial("y1", 0.0);
      double y0 = lqph->get_initial("y0", 0.0);
      // Closed form of y_{k+1} = 11/4 y_k - 11/4 y_{k-1} + 5/4 y_{k-2} - 1/4
      // y_{k-3}, replacing an O(n_jumps) loop per quadrature point. Its roots
      // are a defective 1 (giving A + B n) and a pair of modulus 1/2 (the
      // transient), so y_n = A + B n + 2 Re(C lambda^n). Do not evaluate this by
      // matrix powers: the defective unit root amplifies its splitting error
      // linearly in n (~3% at n = 1e5).
      const double A_sec = y0 + 0.25 * y1 - 0.5 * y2 + 0.25 * y3;
      (void)A_sec; // cancels in y_n - y_0.
      const double B_sec = 2.0 * y0 - 3.5 * y1 + 2.0 * y2 - 0.5 * y3;
      increm = B_sec * static_cast<double>(n_jumps);
      // The transient decays as 2^-n; past ~60 jumps dropping it costs <1e-4
      // relative, so only small jumps pay for the complex arithmetic.
      if (n_jumps < 60) {
        const std::complex<double> lambda(0.375, 0.330718913883073824);
        const std::complex<double> C =
            std::complex<double>(0.0, 0.377964473009227227) * y0 +
            std::complex<double>(-0.125, -0.897665623396914665) * y1 +
            std::complex<double>(0.25, 0.661437827766147648) * y2 +
            std::complex<double>(-0.125, -0.14173667737846021) * y3;
        increm += 2.0 * (C * std::pow(lambda, n_jumps)).real() - 2.0 * C.real();
      }
    }
    return increm;
  };
};

template<int dim>
class JonasAccumulation : public FatigueAccumulation<dim> {
public:
  JonasAccumulation(Controller<dim> &ctl) : FatigueAccumulation<dim>(ctl) {
    increment_law = resolve_fatigue_increment("CarraraNoMeanEffect", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    int n_jumps = static_cast<int>(ctl.get_info("N jump", 0.0));
    double increm = 0;
    double trial_cycle = ctl.get_info("Trial cycle", 0.0);

    if ((n_jumps == 0 && std::abs(trial_cycle) < 1e-8) ||
        ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      increm = evaluate_fatigue_increment(increment_law, lqph, degrade, ctl);
    } else if (n_jumps > 0) {
      double y1 = lqph->get_initial("y1", 0.0);
      double y2 = lqph->get_initial("y2", 0.0);
      double y3 = lqph->get_initial("y3", 0.0);
      double y4 = lqph->get_initial("y4", 0.0);
      increm =
          1.0 / 6.0 * (-2.0 * y1 + 9.0 * y2 - 18.0 * y3 + 11.0 * y4) * n_jumps +
          0.5 * (-y1 + 4.0 * y2 - 5.0 * y3 + 2.0 * y4) * std::pow(n_jumps, 2.0);
    }
    return increm;
  };

  void record(const std::shared_ptr<PointHistory> &lqph, double phasefield,
              double degrade, double degrade_derivative,
              double degrade_second_derivative, Controller<dim> &ctl) {
    double subcycle = ctl.get_info("Subcycle", 0.0);
    if (std::fmod(subcycle, 1) < 1e-8 && subcycle < 5 && subcycle > 1e-8) {
      // y1, y2, y3, and y4
      lqph->update(
        "y" + std::to_string(static_cast<unsigned int>(std::round(subcycle))),
        lqph->get_latest("Fatigue history", 0.0));
    }
  }

  std::string increment_law;
};

template<int dim>
class JonasCLAAccumulation : public JonasAccumulation<dim> {
public:
  JonasCLAAccumulation(Controller<dim> &ctl) : JonasAccumulation<dim>(ctl) {
    AssertThrow(ctl.params.fatigue_accumulation_parameters != "",
                ExcInternalError(
                  "Parameters of JonasCLAAccumulation is not assigned."));
    std::istringstream iss(ctl.params.fatigue_accumulation_parameters);
    iss >> R;
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    int n_jumps = static_cast<int>(ctl.get_info("N jump", 0.0));
    double increm = 0;
    double trial_cycle = ctl.get_info("Trial cycle", 0.0);

    if ((n_jumps == 0 && std::abs(trial_cycle) < 1e-8) ||
        ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      if (ctl.current_timestep != ctl.params.timestep_size_2) {
        increm = evaluate_fatigue_increment(this->increment_law, lqph, degrade,
                                            ctl);
      } else {
        double psi = lqph->get_latest("Positive elastic energy", 0.0) * degrade;
        increm = psi * (1 - R * R * (R >= 0 ? 1 : 0));
      }
    } else if (n_jumps > 0) {
      double y1 = lqph->get_initial("y1", 0.0);
      double y2 = lqph->get_initial("y2", 0.0);
      double y3 = lqph->get_initial("y3", 0.0);
      double y4 = lqph->get_initial("y4", 0.0);
      increm =
          1.0 / 6.0 * (-2.0 * y1 + 9.0 * y2 - 18.0 * y3 + 11.0 * y4) * n_jumps +
          0.5 * (-y1 + 4.0 * y2 - 5.0 * y3 + 2.0 * y4) * std::pow(n_jumps, 2.0);
    }
    return increm;
  };
  double R;
};

template<int dim>
class JonasNodegradeAccumulation : public JonasAccumulation<dim> {
public:
  JonasNodegradeAccumulation(Controller<dim> &ctl)
    : JonasAccumulation<dim>(ctl) {
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    return JonasAccumulation<dim>::increment(
      lqph, phasefield, ctl.params.constant_k + 1, -2.0, 0.0, ctl);
  };
};

template<int dim>
class YangAccumulation : public FatigueAccumulation<dim> {
public:
  YangAccumulation(Controller<dim> &ctl) : FatigueAccumulation<dim>(ctl) {
    increment_law = resolve_fatigue_increment("CarraraNoMeanEffect", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    double n_jumps = ctl.get_info("N jump", 0.0);
    double increm = 0;
    double subcycle = ctl.get_info("Subcycle", 0.0);

    if (std::abs(subcycle - 0) > 1e-8 || std::abs(n_jumps) < 1e-8 ||
        ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      increm = evaluate_fatigue_increment(increment_law, lqph, degrade, ctl);
    } else {
      double last_jump = ctl.get_info("Last jump", 0.0);
      double new_increment = lqph->get_initial("Fast increment", 0.0);
      double diff = lqph->get_initial("Fast increment diff", 0.0);
      double extra_increm = new_increment + n_jumps / last_jump * diff;
      increm = (n_jumps - 1) * (extra_increm + new_increment) / 2;
    }
    return increm;
  };

  void record(const std::shared_ptr<PointHistory> &lqph, double phasefield,
              double degrade, double degrade_derivative,
              double degrade_second_derivative, Controller<dim> &ctl) {
    double subcycle = ctl.get_info("Subcycle", 0.0);

    if (std::fmod(subcycle, 1) < 1e-8) {
      lqph->update("Initial history", lqph->get_latest("Fatigue history", 0.0));
      if (std::abs(subcycle - 1) < 1e-8) {
        double new_increment = lqph->get_latest("Fatigue history", 0.0) -
                               lqph->get_initial("Initial history", 0.0);
        double old_increment =
            lqph->get_initial("Fast increment", new_increment);
        lqph->update("Fast increment", new_increment);
        lqph->update("Fast increment diff", new_increment - old_increment);
      }
    }
  }

  std::string increment_law;
};

template<int dim>
class JacconAccumulation : public FatigueAccumulation<dim> {
public:
  JacconAccumulation(Controller<dim> &ctl) : FatigueAccumulation<dim>(ctl) {
    increment_law = resolve_fatigue_increment("CarraraNoMeanEffect", ctl);
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    double n_jumps = ctl.get_info("N jump", 0.0);
    double increm = 0;
    double subcycle = ctl.get_info("Subcycle", 0.0);

    if (std::abs(subcycle - 0) > 1e-8 ||
        ctl.current_timestep != ctl.params.timestep_size_2) {
      // Regular accumulation
      increm = evaluate_fatigue_increment(increment_law, lqph, degrade, ctl);
    } else {
      double n_trials = ctl.get_info("N trials", 0);
      if (std::abs(n_trials) < 1e-8) {
        increm = n_jumps * lqph->get_initial("Trial increment");
      } else {
        double last_increm = lqph->get_independent_initial("increm", 0.0);
        double residual = lqph->get_independent_latest(
          "Residual", 0.0); // pointhistory is not finalized when it has not
        // converged, so we need the latest one.
        increm = last_increm - (-1) * residual;
      }
      lqph->update_independent("increm", increm);
      increm -= lqph->get_initial("Trial increment"); // we are at cycle N+1
    }
    return increm;
  };

  void record(const std::shared_ptr<PointHistory> &lqph, double phasefield,
              double degrade, double degrade_derivative,
              double degrade_second_derivative, Controller<dim> &ctl) {
    double n_jumps = ctl.get_info("N jump", 0.0);
    double subcycle = ctl.get_info("Subcycle", 0.0);
    if (std::fmod(subcycle, 1) < 1e-8) {
      double alpha_n1 = lqph->get_latest("Fatigue history", 0.0);
      lqph->update("alpha_n1", alpha_n1);
      if (std::abs(subcycle - 1) < 1e-8) {
        double alpha_n0 = lqph->get_initial("alpha_n1", 0.0);
        double new_increment = alpha_n1 - alpha_n0;
        lqph->update("Trial increment", new_increment);
        lqph->update("alpha_0", alpha_n0);
        lqph->update("alpha_1", alpha_n1);
      }
    }
    if (std::abs(subcycle - 1) < 1e-8) {
      double alpha_n0 = lqph->get_initial("alpha_n1", 0.0);
      double alpha_n1 = lqph->get_latest("alpha_n1", 0.0);
      double alpha_0 = lqph->get_initial("alpha_0", 0.0);
      double alpha_1 = lqph->get_initial("alpha_1", 0.0);
      double residual = alpha_0 * (1 - n_jumps / 2) +
                        (alpha_n1 + alpha_1) * (n_jumps / 2) -
                        alpha_n0 * (1 + n_jumps / 2);
      lqph->update_independent("Residual", residual);
    }
  }

  std::string increment_law;
};

template<int dim>
class JacconNodegradeAccumulation : public JacconAccumulation<dim> {
public:
  JacconNodegradeAccumulation(Controller<dim> &ctl)
    : JacconAccumulation<dim>(ctl) {
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    return JacconAccumulation<dim>::increment(
      lqph, phasefield, ctl.params.constant_k + 1, -2.0, 0.0, ctl);
  };
};

template<int dim>
class CarraraMeanEffectAccumulation : public FatigueAccumulation<dim> {
public:
  CarraraMeanEffectAccumulation(Controller<dim> &ctl)
    : FatigueAccumulation<dim>(ctl) {
    assert_fixed_fatigue_increment("CarraraMeanEffect", ctl);
    // Eq. 47 does not match any of alpha_t claimed in the result section
    // We have to multiply another 0.5 to reproduce the results.
    alpha_n = resolve_fatigue_alpha_t(ctl, default_fatigue_alpha_t(ctl),
                                      ctl.params.fatigue_accumulation_parameters);
    if (ctl.params.fatigue_alpha_t != "")
      ctl.dcout << "Using alpha_n: " << alpha_n << " from Fatigue alpha_t"
          << std::endl;
    else if (ctl.params.fatigue_accumulation_parameters != "")
      ctl.dcout << "Using alpha_n: " << alpha_n << "from configuration"
          << std::endl;
    else
      ctl.dcout << "Using alpha_n: " << alpha_n << std::endl;
  };

  double increment(const std::shared_ptr<PointHistory> &lqph, double phasefield,
                   double degrade, double degrade_derivative,
                   double degrade_second_derivative,
                   Controller<dim> &ctl) override {
    return carrara_mean_effect_increment(lqph, degrade, alpha_n);
  };

private:
  double alpha_n;
};

template<int dim>
std::unique_ptr<FatigueAccumulation<dim> >
select_fatigue_accumulation(std::string method, Controller<dim> &ctl) {
  if (method == "CarraraNoMeanEffect")
    return std::make_unique<CarraraNoMeanEffectAccumulation<dim> >(ctl);
  else if (method == "CarraraMeanEffect")
    return std::make_unique<CarraraMeanEffectAccumulation<dim> >(ctl);
  else if (method == "Kristensen")
    return std::make_unique<KristensenAccumulation<dim> >(ctl);
  else if (method == "KristensenCLA")
    return std::make_unique<KristensenCLAAccumulation<dim> >(ctl);
  else if (method == "Cojocaru")
    return std::make_unique<CojocaruAccumulation<dim> >(ctl);
  else if (method == "CojocaruCLA")
    return std::make_unique<CojocaruCLAAccumulation<dim> >(ctl);
  else if (method == "Li")
    return std::make_unique<LiAccumulation<dim> >(ctl);
  else if (method == "LiCLA")
    return std::make_unique<LiCLAAccumulation<dim> >(ctl);
  else if (method == "Jonas")
    return std::make_unique<JonasAccumulation<dim> >(ctl);
  else if (method == "JonasCLA")
    return std::make_unique<JonasCLAAccumulation<dim> >(ctl);
  else if (method == "JonasNodegrade")
    return std::make_unique<JonasNodegradeAccumulation<dim> >(ctl);
  else if (method == "Yang")
    return std::make_unique<YangAccumulation<dim> >(ctl);
  else if (method == "Jaccon")
    return std::make_unique<JacconAccumulation<dim> >(ctl);
  else if (method == "JacconNodegrade")
    return std::make_unique<JacconNodegradeAccumulation<dim> >(ctl);
  else
    AssertThrow(false, ExcNotImplemented());
}

template<int dim>
class FatigueDegradation {
public:
  FatigueDegradation(Controller<dim> &ctl) {
  };

  virtual double degradation_value(const std::shared_ptr<PointHistory> &lqph,
                                   double phasefield, double degrade,
                                   Controller<dim> &ctl) {
    AssertThrow(false, ExcNotImplemented());
  };
};

// https://www.sciencedirect.com/science/article/pii/S0045782519306218
template<int dim>
class CarraraAsymptoticFatigueDegradation : public FatigueDegradation<dim> {
public:
  CarraraAsymptoticFatigueDegradation(Controller<dim> &ctl)
    : FatigueDegradation<dim>(ctl) {
    // Eq. 47 does not match any of alpha_t claimed in the result section
    // We have to multiply another 0.5 to reproduce the results.
    alpha_t = resolve_fatigue_alpha_t(ctl, default_fatigue_alpha_t(ctl),
                                      ctl.params.fatigue_degradation_parameters);
    if (ctl.params.fatigue_alpha_t != "")
      ctl.dcout << "Using alpha_t: " << alpha_t << " from Fatigue alpha_t"
          << std::endl;
    else if (ctl.params.fatigue_degradation_parameters != "")
      ctl.dcout << "Using alpha_t: " << alpha_t << " from configuration"
          << std::endl;
    else
      ctl.dcout << "Using alpha_t: " << alpha_t << std::endl;
  };

  double degradation_value(const std::shared_ptr<PointHistory> &lqph,
                           double phasefield, double phasefield_degrade,
                           Controller<dim> &ctl) override {
    double degrade;
    double alpha = lqph->get_latest("Fatigue history", 0.0);
    if (alpha <= alpha_t) {
      degrade = 1;
    } else {
      degrade = std::pow(2 * alpha_t / (alpha + alpha_t), 2);
    }
    return degrade;
  };

  double alpha_t;
};

// https://www.sciencedirect.com/science/article/pii/S0045782519306218
template<int dim>
class CarraraLogarithmicFatigueDegradation : public FatigueDegradation<dim> {
public:
  CarraraLogarithmicFatigueDegradation(Controller<dim> &ctl)
    : FatigueDegradation<dim>(ctl) {
    // kappa has no default, so the parameter string stays mandatory even when
    // the global "Fatigue alpha_t" is set -- the global only overrides the
    // leading alpha_t slot, kappa is always read positionally from here.
    AssertThrow(
      ctl.params.fatigue_degradation_parameters != "",
      ExcInternalError("Parameters of CarraraLogarithmicFatigueDegradation "
        "is not assigned."));
    // Only kappa is taken from this parse; the leading alpha_t slot is read by
    // resolve_fatigue_alpha_t below, the same way every other scheme does it.
    std::istringstream iss(ctl.params.fatigue_degradation_parameters);
    double parsed_alpha_t;
    AssertThrow(static_cast<bool>(iss >> parsed_alpha_t >> kappa),
                ExcInternalError(
                  "Parameters of CarraraLogarithmicFatigueDegradation must be "
                  "'<alpha_t> <kappa>', but got: " +
                  ctl.params.fatigue_degradation_parameters));
    alpha_t = resolve_fatigue_alpha_t(ctl, default_fatigue_alpha_t(ctl),
                                      ctl.params.fatigue_degradation_parameters);
    if (ctl.params.fatigue_alpha_t != "")
      ctl.dcout << "Using alpha_t: " << alpha_t << " from Fatigue alpha_t"
          << " and kappa: " << kappa << std::endl;
    else
      ctl.dcout << "Using alpha_t: " << alpha_t << " and kappa: " << kappa
          << std::endl;
  };

  double degradation_value(const std::shared_ptr<PointHistory> &lqph,
                           double phasefield, double phasefield_degrade,
                           Controller<dim> &ctl) override {
    double degrade;
    double alpha = lqph->get_latest("Fatigue history", 0.0);
    if (alpha <= alpha_t) {
      degrade = 1;
    } else if (alpha <= (alpha_t * std::pow(10, 1 / kappa))) {
      degrade = std::pow(1 - kappa * std::log10(alpha / alpha_t), 2);
    } else {
      degrade = 0;
    }
    return degrade;
  };

  double alpha_t, kappa;
};

template<int dim>
class KristensenAsymptoticFatigueDegradation
    : public CarraraAsymptoticFatigueDegradation<dim> {
public:
  KristensenAsymptoticFatigueDegradation(Controller<dim> &ctl)
    : CarraraAsymptoticFatigueDegradation<dim>(ctl) {
    // Deliberately ignores "Fatigue degradation parameters" -- this scheme
    // defines alpha_t from Gc and l_phi. The global "Fatigue alpha_t" does
    // override it, so the threshold can be set independently of Gc.
    this->alpha_t = resolve_fatigue_alpha_t(
      ctl, ctl.params.Gc / (12 * ctl.params.l_phi));
    if (ctl.params.fatigue_alpha_t != "")
      ctl.dcout << "Using alpha_t: " << this->alpha_t << " from Fatigue alpha_t"
          << std::endl;
  };
};

template<int dim>
class CojocaruAsymptoticFatigueDegradation : public FatigueDegradation<dim> {
public:
  CojocaruAsymptoticFatigueDegradation(Controller<dim> &ctl)
    : FatigueDegradation<dim>(ctl) {
    // alpha_t is the only parameter here, so the global "Fatigue alpha_t" can
    // supply it outright; the parameter string is required only without it.
    AssertThrow(ctl.params.fatigue_degradation_parameters != "" ||
                  ctl.params.fatigue_alpha_t != "",
                ExcInternalError(
                  "Parameters of CojocaruCLAAccumulation is not assigned."));
    alpha_t = resolve_fatigue_alpha_t(
      ctl, 0.0, ctl.params.fatigue_degradation_parameters);
    if (ctl.params.fatigue_alpha_t != "")
      ctl.dcout << "Using alpha_t: " << alpha_t << " from Fatigue alpha_t"
          << std::endl;
    else
      ctl.dcout << "Using alpha_t: " << alpha_t << " from configuration"
          << std::endl;
  };

  double degradation_value(const std::shared_ptr<PointHistory> &lqph,
                           double phasefield, double phasefield_degrade,
                           Controller<dim> &ctl) override {
    double alpha = lqph->get_latest("Fatigue history", 0.0);
    double degrade = std::pow(alpha_t / (alpha + alpha_t), 2);
    return degrade;
  };
  double alpha_t;
};

template<int dim>
std::unique_ptr<FatigueDegradation<dim> >
select_fatigue_degradation(std::string method, Controller<dim> &ctl) {
  if (method == "CarraraAsymptotic")
    return std::make_unique<CarraraAsymptoticFatigueDegradation<dim> >(ctl);
  else if (method == "CarraraLogarithmic")
    return std::make_unique<CarraraLogarithmicFatigueDegradation<dim> >(ctl);
  else if (method == "KristensenAsymptotic")
    return std::make_unique<KristensenAsymptoticFatigueDegradation<dim> >(ctl);
  else if (method == "CojocaruAsymptotic")
    return std::make_unique<CojocaruAsymptoticFatigueDegradation<dim> >(ctl);
  else
    AssertThrow(false, ExcNotImplemented());
}

#endif // CRACKS_FATIGUE_DEGRADATION_H
