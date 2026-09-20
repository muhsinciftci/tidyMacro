#ifndef RNG_TIDYMACRO_H
#define RNG_TIDYMACRO_H

// Self-contained pseudo-random number generator for the Bayesian sign- and
// narrative-restriction routines.
//
// Rationale
//   * R's RNG (and arma::randn) is NOT thread safe and must never be touched
//     from inside an OpenMP region.  Every worker therefore owns one RNG
//     object seeded deterministically from (base_seed, draw_index), so results
//     are reproducible and independent of the thread count.
//   * The generator and every distribution are implemented here rather than
//     taken from <random>, because the standard library's distribution
//     implementations differ between libstdc++ and libc++; hand-rolling them
//     keeps draws bit-identical across platforms.
//
// Generator: xoshiro256++ (Blackman & Vigna), seeded through splitmix64.

#include <cstdint>
#include <cmath>

namespace tidymacro {

class RNG {
public:
    explicit RNG(std::uint64_t seed) : has_spare_(false), spare_(0.0) {
        std::uint64_t z = seed + 0x9E3779B97F4A7C15ULL;
        for (int i = 0; i < 4; ++i) {
            s_[i] = splitmix64(z);
        }
        // Discard a short prefix so nearby seeds decorrelate immediately.
        for (int i = 0; i < 16; ++i) (void)next_u64();
    }

    // Uniform on (0, 1): 53 significant bits, never exactly 0 or 1.
    inline double unif() {
        const double u = static_cast<double>(next_u64() >> 11) * (1.0 / 9007199254740992.0);
        return (u <= 0.0) ? 1e-300 : ((u >= 1.0) ? 0.99999999999999989 : u);
    }

    // Standard normal via the Marsaglia polar method (one cached spare).
    inline double norm() {
        if (has_spare_) { has_spare_ = false; return spare_; }
        double u, v, s;
        do {
            u = 2.0 * unif() - 1.0;
            v = 2.0 * unif() - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double f = std::sqrt(-2.0 * std::log(s) / s);
        spare_     = v * f;
        has_spare_ = true;
        return u * f;
    }

    // Gamma(shape, scale = 1) via Marsaglia & Tsang (2000), with the
    // shape < 1 boosting step.
    inline double gamma(double shape) {
        if (shape < 1.0) {
            const double u = unif();
            return gamma(shape + 1.0) * std::pow(u, 1.0 / shape);
        }
        const double d = shape - 1.0 / 3.0;
        const double c = 1.0 / std::sqrt(9.0 * d);
        for (;;) {
            double x, v;
            do {
                x = norm();
                v = 1.0 + c * x;
            } while (v <= 0.0);
            v = v * v * v;
            const double u  = unif();
            const double x2 = x * x;
            if (u < 1.0 - 0.0331 * x2 * x2) return d * v;
            if (std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v))) return d * v;
        }
    }

    // Chi-square with df degrees of freedom (df > 0, need not be an integer).
    inline double chisq(double df) { return 2.0 * gamma(0.5 * df); }

private:
    static inline std::uint64_t splitmix64(std::uint64_t& z) {
        z += 0x9E3779B97F4A7C15ULL;
        std::uint64_t r = z;
        r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9ULL;
        r = (r ^ (r >> 27)) * 0x94D049BB133111EBULL;
        return r ^ (r >> 31);
    }

    static inline std::uint64_t rotl(std::uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    inline std::uint64_t next_u64() {
        const std::uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    std::uint64_t s_[4];
    bool   has_spare_;
    double spare_;
};

} // namespace tidymacro

#endif // RNG_TIDYMACRO_H
