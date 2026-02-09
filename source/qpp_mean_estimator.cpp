#include "nozzle_qns/qpp_mean_estimator.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace nozzle_qns {

namespace {

    using namespace qpp::literals;

    inline double clamp01(double x) {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }

    idx ceil_log2(idx x) {
        if (x <= 1) return 0;
        idx p = 0;
        idx v = 1;
        while (v < x) {
            v <<= 1;
            ++p;
        }
        return p;
    }

    std::vector<double> pad_to_pow2(const std::vector<double>& g, idx N) {
        std::vector<double> out(N, 0.0);
        const idx K = static_cast<idx>(g.size());
        for (idx i = 0; i < std::min<idx>(K, N); ++i)
            out[i] = clamp01(g[i]);
        return out; // remaining padded with 0
    }

    std::vector<idx> make_dims(idx nq) {
        return std::vector<idx>(static_cast<std::size_t>(nq), 2);
    }

    inline bool get_qubit_bit(idx basis, idx q, idx nq) {
        const idx pos = (nq - 1 - q);
        return ((basis >> pos) & 1ULL) != 0ULL;
    }

    // Map a global basis index -> integer y from the phase register qubits.
    idx phase_outcome_from_basis(idx basis, const std::vector<idx>& phase_regs, idx nq_total) {
        idx y = 0;
        const idx t = static_cast<idx>(phase_regs.size());
        for (idx r = 0; r < t; ++r) {
            const bool bit = get_qubit_bit(basis, phase_regs[r], nq_total);
            y = (y << 1) | (bit ? 1ULL : 0ULL); // phase_regs[0] is MSB
        }
        return y;
    }

    // Sample measurement outcome y on phase register (by marginalizing probabilities)
    idx sample_phase_y(const qpp::ket& st,
                    const std::vector<idx>& phase_regs,
                    idx nq_total,
                    std::mt19937_64& rng) {
        const idx t = static_cast<idx>(phase_regs.size());
        const idx M = 1ULL << t;

        std::vector<double> probs(static_cast<std::size_t>(M), 0.0);
        const idx dim_total = 1ULL << nq_total;

        for (idx basis = 0; basis < dim_total; ++basis) {
            const idx y = phase_outcome_from_basis(basis, phase_regs, nq_total);
            const double p = std::norm(st(static_cast<qpp::idx>(basis)));
            probs[static_cast<std::size_t>(y)] += p;
        }

        std::discrete_distribution<idx> dist(probs.begin(), probs.end());
        return dist(rng);
    }

    // Decode y to mu estimate using standard AE decoding:
    // φ = y / 2^t, θ ≈ π * min(φ, 1-φ), mu = sin^2 θ
    double decode_mu_from_y(idx y, idx t) {
        const double M = static_cast<double>(1ULL << t);
        const double phi = static_cast<double>(y) / M;
        const double theta = static_cast<double>(pi) * std::min(phi, 1.0 - phi);
        const double s = std::sin(theta);
        return clamp01(s * s);
    }

    // -----Gate-based building blocks on a ket (qpp)-----

    // Apply the oracle O on (index_regs + ancilla):
    // For each j, apply multi-controlled Ry(theta_j) on ancilla,
    // controlled on index_regs == |j>, with theta_j = 2 asin sqrt(g[j]).
    void apply_oracle_O(qpp::ket& state,
                        const std::vector<idx>& index_regs, // MSB..LSB ordering
                        idx ancilla_reg,
                        const std::vector<double>& gpad) {
        const idx n = static_cast<idx>(index_regs.size());
        const idx N = 1ULL << n;
        ensure(static_cast<idx>(gpad.size()) == N, "apply_oracle_O: gpad size mismatch");

        for (idx j = 0; j < N; ++j) {
            const double gj = clamp01(gpad[j]);
            const double theta = -2.0 * std::acos(std::sqrt(gj)); //2.0 * std::asin(std::sqrt(gj));
            if (std::abs(theta) < 1e-12) continue;

            std::vector<idx> flipped;
            flipped.reserve(static_cast<std::size_t>(n));

            for (idx bit = 0; bit < n; ++bit) {
                const bool is_one = ((j >> bit) & 1ULL) != 0ULL;
                const idx q = index_regs[n - 1 - bit];
                if (!is_one) {
                    state = qpp::apply(state, qpp::gt.X, {q});
                    flipped.push_back(q);
                }
            }

            // Multi-controlled RY(theta) on ancilla, controlled by all index_regs
            state = qpp::applyCTRL(state, qpp::gt.RY(theta), index_regs, {ancilla_reg});

            // Uncompute X flips
            for (idx q : flipped)
                state = qpp::apply(state, qpp::gt.X, {q});
        }
    }

    void apply_oracle_O_dag(qpp::ket& state,
                            const std::vector<idx>& index_regs, // MSB..LSB
                            idx ancilla_reg,
                            const std::vector<double>& gpad) {
        // Same but with -theta
        const idx n = static_cast<idx>(index_regs.size());
        const idx N = 1ULL << n;
        ensure(static_cast<idx>(gpad.size()) == N, "apply_oracle_O_dag: gpad size mismatch");

        for (idx j = N; j-- > 0;) {
            const double gj = clamp01(gpad[j]);
            const double theta = 2.0 * std::acos(std::sqrt(gj));//-2.0 * std::asin(std::sqrt(gj));
            if (std::abs(theta) < 1e-12) continue;

            std::vector<idx> flipped;
            flipped.reserve(static_cast<std::size_t>(n));

            for (idx bit = 0; bit < n; ++bit) {
                const bool is_one = ((j >> bit) & 1ULL) != 0ULL;
                const idx q = index_regs[n - 1 - bit];
                if (!is_one) {
                    state = qpp::apply(state, qpp::gt.X, {q});
                    flipped.push_back(q);
                }
            }

            state = qpp::applyCTRL(state, qpp::gt.RY(theta), index_regs, {ancilla_reg});

            for (idx q : flipped)
                state = qpp::apply(state, qpp::gt.X, {q});
        }
    }

    // A = O (QFT ⊗ I) acting on (index_regs + ancilla)
    void apply_A(qpp::ket& state,
                const std::vector<idx>& index_regs,
                idx ancilla_reg,
                const std::vector<double>& gpad) {
        // QFT on index register
        state = qpp::applyQFT(state, index_regs);
        // Oracle
        apply_oracle_O(state, index_regs, ancilla_reg, gpad);
    }

    // A^{-1} = A^\dagger = (QFT)^\dagger O^\dagger
    void apply_A_dag(qpp::ket& state,
                    const std::vector<idx>& index_regs,
                    idx ancilla_reg,
                    const std::vector<double>& gpad) {
        apply_oracle_O_dag(state, index_regs, ancilla_reg, gpad);
        state = qpp::applyTFQ(state, index_regs);
    }

    // S_chi: flip phase of “good” states (ancilla |1>) -> Z on ancilla
    void apply_Schi(qpp::ket& state, idx ancilla_reg) {
        state = qpp::apply(state, qpp::gt.Z, {ancilla_reg});
    }

    // S0: Reflection about |0...0>_index |1>_anc
    void apply_S0(qpp::ket& state, 
                const std::vector<qpp::idx>& index_regs, 
                qpp::idx ancilla_reg) {
        
        // 1. Flip only the index registers (0 -> 1)
        for (qpp::idx q : index_regs)
            state = qpp::apply(state, qpp::gt.X, {q});

        // 2. Reflection: Multi-Controlled Z
        state = qpp::applyCTRL(state, qpp::gt.Z, index_regs, {ancilla_reg});

        // 3. Uncompute (Flip index registers back 1 -> 0)
        for (qpp::idx q : index_regs)
            state = qpp::apply(state, qpp::gt.X, {q});
    }

    // Apply Grover iterate Q on (index+anc) (phase register untouched):
    // Q = -A S0 A^{-1} S_chi
    void apply_Q(qpp::ket& state,
                const std::vector<idx>& index_regs,
                idx ancilla_reg,
                const std::vector<double>& gpad) {
        apply_Schi(state, ancilla_reg);
        apply_A_dag(state, index_regs, ancilla_reg, gpad);
        apply_S0(state, index_regs, ancilla_reg);
        apply_A(state, index_regs, ancilla_reg, gpad);
        state = qpp::apply(state, qpp::gt.RZ(2 * M_PI), {ancilla_reg});
    }

    // Build a dense system-unitary for Q on (index+anc) only (dimension 2N) by simulating apply_Q on basis vectors of the system.

    qpp::cmat build_Q_system_unitary(idx n, const std::vector<double>& gpad) {
        const idx N = 1ULL << n;
        const idx nq_sys = n + 1;
        const idx dim_sys = 1ULL << nq_sys;
        ensure(static_cast<idx>(gpad.size()) == N, "build_Q_system_unitary: gpad size mismatch");

        // Local indexing for system simulation:
        // index regs: 0..n-1 (MSB..LSB ordering), anc = n
        std::vector<idx> index_regs(n);
        for (idx i = 0; i < n; ++i) index_regs[i] = i;
        const idx anc = n;

        const std::vector<idx> dims_sys = make_dims(nq_sys);

        qpp::cmat Qsys = qpp::cmat::Zero(static_cast<qpp::idx>(dim_sys),
                                        static_cast<qpp::idx>(dim_sys));

        // Build columns: Q |e_k>
        for (idx k = 0; k < dim_sys; ++k) {
            qpp::ket basis = qpp::ket::Zero(static_cast<qpp::idx>(dim_sys));
            basis(static_cast<qpp::idx>(k)) = 1.0;

            apply_Q(basis, index_regs, anc, gpad);

            // set as k-th column
            Qsys.col(static_cast<qpp::idx>(k)) = basis;
        }

        return Qsys;
    }

    // One QAEA run
    double qaea_once_mu(const std::vector<double>& g,
                        const QppMeanEstimator::Options& opt,
                        double eps,
                        std::mt19937_64& rng) {
        
        // 1. Setup and Dimensions
        const idx K = static_cast<idx>(g.size());
        if (K == 0) throw std::runtime_error("QppMeanEstimator: empty g");

        const idx n = (opt.index_qubits > 0) ? opt.index_qubits : ceil_log2(K);
        const idx N = 1ULL << n;
        auto gpad = pad_to_pow2(g, N);

        // Determine phase register size t
        idx t = n; 
        if (eps > 0.0) {
            const idx t_eps = ceil_log2(static_cast<idx>(std::ceil((2.0 * M_PI) / eps)));
            t = std::max(t, t_eps);
        }

        const idx nq_total = t + n + 1;
        const auto dims_total = make_dims(nq_total);

        // 2. Register Mapping
        // Phase: 0..t-1  | Index: t..t+n-1 | Ancilla: t+n
        std::vector<idx> phase_regs(t);
        std::iota(phase_regs.begin(), phase_regs.end(), 0);

        std::vector<idx> index_regs(n);
        std::iota(index_regs.begin(), index_regs.end(), t);

        const idx anc = t + n;
        
        // System registers for Q application (Index + Ancilla)
        std::vector<idx> sys_regs = index_regs;
        sys_regs.push_back(anc);

        // 3. Build Base Operator Q
        qpp::cmat Q_power = build_Q_system_unitary(n, gpad);

        // 4. State Initialization
        // |0...0>_phase |0...0>_index |1>_anc
        const idx dim_total = 1ULL << nq_total;
        qpp::ket st = qpp::ket::Zero(static_cast<qpp::idx>(dim_total));
        st(0) = 1.0; 
        st = qpp::apply(st, qpp::gt.X, {anc}); 

        // Superposition on Phase Register
        st = qpp::applyQFT(st, phase_regs);

        // Prepare eigenstate |psi> on System (A |0...01>)
        apply_A(st, index_regs, anc, gpad);

        // 5. Phase Estimation Loop (The Optimization)
        // We iterate t times. In each step r:
        //   1. Apply Controlled-Q^(2^r)
        //   2. Square Q to get Q^(2^(r+1)) for the next step
        
        for (idx r = 0; r < t; ++r) {
            const idx control = phase_regs[t - 1 - r];

            // Apply Controlled-CurrentPower
            st = qpp::applyCTRL(st, Q_power, {control}, sys_regs, dims_total);

            if (r < t - 1) {
                Q_power = Q_power * Q_power;
            }
        }

        // 6. Inverse QFT
        st = qpp::applyTFQ(st, phase_regs);

        // Measure
        const idx y = sample_phase_y(st, phase_regs, nq_total, rng);

        return decode_mu_from_y(y, t);
    }

    // Repeat O(log(1/delta)) times and take median
    double qaea_mu(const std::vector<double>& g,
                const QppMeanEstimator::Options& opt,
                double eps,
                double delta,
                std::mt19937_64& rng) {
        idx R = 1;
        if (delta > 0.0 && delta < 1.0) {
            R = static_cast<idx>(std::ceil(std::log(1.0 / delta)));
            R = std::max<idx>(1, R);
        }
        std::vector<double> samples;
        samples.reserve(static_cast<std::size_t>(R));
        for (idx i = 0; i < R; ++i) {
            samples.push_back(qaea_once_mu(g, opt, eps, rng));
        }
        std::nth_element(samples.begin(), samples.begin() + samples.size() / 2, samples.end());
        return clamp01(samples[samples.size() / 2]);
    }

    } // unnamed namespace

    // Public API
    MeanEstimate QppMeanEstimator::estimate_mean(const std::vector<double>& g1,
                                                const std::vector<double>& g2,
                                                const std::vector<double>& g3,
                                                double eps_target,
                                                double delta_target) {
        const double eps = (eps_target > 0.0) ? eps_target : opt_.default_eps;
        const double delta = (delta_target > 0.0) ? delta_target : opt_.default_delta;

        const double G1 = qaea_mu(g1, opt_, eps, delta, rng_);
        const double G2 = qaea_mu(g2, opt_, eps, delta, rng_);
        const double G3 = qaea_mu(g3, opt_, eps, delta, rng_);

        return MeanEstimate{G1, G2, G3, eps, 1.0 - delta};
    }

} // namespace nozzle_qns
