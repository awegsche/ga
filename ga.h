#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <iterator>
#include <ostream>
#include <random>
#include <vector>

using std::cin, std::cout, std::cerr, std::endl;
using std::ostream;
using std::vector;

// -------------------------------------------------------------------------------------------------
// ---- genetic algorithm --------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
namespace ga
{

    // CRTP
    template <typename Derived>
    class Nucl
    {
    public:
        template <typename R>
        static auto random(R& rng) -> Derived
        {
            return Derived::random(rng);
        }

        static void crossover_inplace(Derived const& a, Derived const& b, Derived* c) {
            return Derived::crossover_inplace(a, b, c);
        }

        Derived const& derived() const
        {
            return *static_cast<Derived const*>(this);
        }

        template <typename D>
        friend auto operator<<(ostream& os, Nucl<D> const& nucl) -> ostream&;

    private:
        Derived& derived()
        {
            return *static_cast<Derived*>(this);
        }
    };

    /// The Genom class
    ///
    /// Provides genetic operators such as crossover and mutation.
    /// Provides insight in its performance after the simulator has calculated a scorer.
    template <typename N, typename S>
    class Genom
    {
    public:
        // -----------------------------------------------------------------------------------------
        // ---- Init -------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------

        /**
         * @brief Creates a random Genom
         * 
         * @tparam R (Psuedo) random number generator to be used
         * @param n size of the genom (= count of nucleotides)
         * @param rng random number generator
         * @return Genom<N, S> 
         */
        template <typename R>
        static auto random(size_t n, R& rng) -> Genom<N, S>
        {
            Genom g;
            g.nucleotides.reserve(n);

            for (auto i = 0; i < n; i++)
                g.nucleotides.push_back(Nucl<N>::random(rng));

            return g;
        }

        // -----------------------------------------------------------------------------------------
        // ---- Genetic operators ------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------

        /**
         * @brief Crossover (a [genetic operator](https://en.wikipedia.org/wiki/Genetic_operator))
         * Inplace version, avoids allocation
         * 
         *
         * @param a the first parent
         * @param b the second parent
         * @return the child genom after the crossover process
         */
        static auto crossover_inplace(Genom const& a, Genom const& b, Genom* c, size_t n)
        {
            c->nucleotides.resize(a.nucleotides.size());
            for (auto i = 0; i < n; i++)
            {
                Nucl<N>::crossover_inplace(a.nucleotides[i], b.nucleotides[i], &c->nucleotides[i]);
            }
            for (auto i = n; i < a.nucleotides.size(); i++)
            {
                Nucl<N>::crossover_inplace(b.nucleotides[i], a.nucleotides[i], &c->nucleotides[i]);
            }
        }

        /**
         * @brief Mutation (a [genetic operator](https://en.wikipedia.org/wiki/Genetic_operator))
         * 
         * @tparam R 
         * @param n 
         * @param rng 
         * @return auto 
         */
        template <typename R>
        auto mutate(size_t n, R& rng)
        {
            nucleotides[n].mutate(rng);
        }

        /**
         * @brief Shifts all the nucleotides left by `n` places.
         * For real time simulations where the starting point continously shifts this can be used
         * to keep the genoms synched with the simulation
         * 
         * @param n 
         * @return auto 
         */
        auto shift(size_t n = 1)
        {
            std::rotate(nucleotides.begin(), nucleotides.begin() + n, nucleotides.end());
        }

        // -----------------------------------------------------------------------------------------
        // ---- Accessors --------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------

        /**
         * @brief returns the number of nucleotides
         * 
         * @return size_t 
         */
        auto size() const -> size_t
        {
            return nucleotides.size();
        }

        /**
         * @brief gives const access to the nucleotides
         * 
         * @return vector<N> const& 
         */
        auto get_nucleotides() const -> vector<N> const&
        {
            return nucleotides;
        }

        /**
         * @brief Sets the scorer. This function is likely to be called by the simulator after
         * simulating the effect of this genom.
         * 
         * @param s 
         * @return auto 
         */
        auto set_scorer(S&& s)
        {
            scorer = s;
        }

        /**
         * @brief Retrieves the scorer.
         *
         * @attention The simulator has to calculate and assign a scorer first.
         * 
         * @return S const& 
         */
        auto get_scorer() const -> S const& { return scorer; }

        auto get_scorer_mut() -> S& { return scorer; }

        /**
         * @brief Convenience function to retrieve the score of this genom
         * 
         * @return float 
         */
        auto score() const -> float { return scorer.score; }

        template <typename A, typename B>
        friend auto operator<<(ostream& os, Genom<A, B> const& g) -> ostream&;

        auto operator<(Genom const& other) const -> bool
        {
            return scorer.score < other.scorer.score;
        }

    private:
        vector<N> nucleotides;
        S scorer;
    };

    template <typename N, typename S, typename R>
    class GenePool
    {
    public:
        GenePool(size_t genoms_count,
                 size_t genom_len,
                 float take_value)
            : genoms(), n_genoms(genoms_count), rng(std::default_random_engine{}()), distr(0, genom_len), generation(0), take(genoms_count * take_value), retain(genoms_count * (0.5 - take_value))
        {
            genoms.reserve(n_genoms);
            for (int i = 0; i < n_genoms; i++)
                genoms.push_back(Genom<N, S>::random(genom_len, rng));

            shuffle_indices.reserve(n_genoms);
            for (size_t i = 0; i < n_genoms; i++)
                shuffle_indices.push_back(i);

            next_generation.reserve(n_genoms);
        }

        template <typename Simul>
        auto simulate(Simul const& simulator)
        {
            for (auto& genom : genoms)
            {
                simulator.simulate(genom);
            }
            std::sort(genoms.rbegin(), genoms.rend());
            generation++;
        }

        auto shift(size_t n = 1)
        {
            for (auto& genom : genoms)
                genom.shift(n);
        }

        auto reset()
        {
            generation = 0;
        }

        auto select()
        {
            std::shuffle(shuffle_indices.begin() + take, shuffle_indices.end(), rng);
            auto half = n_genoms / 2;
            for (size_t i = 0; i < half; i++)
                shuffle_indices[i + half] = shuffle_indices[i];

            next_generation.resize(n_genoms);

            std::shuffle(shuffle_indices.begin() + half, shuffle_indices.end(), rng);
            for (size_t i = 0; i < half; i++)
                Genom<N, S>::crossover_inplace(genoms[shuffle_indices[i]],
                        genoms[shuffle_indices[i + half]],
                        &next_generation[i],
                        distr(rng));
                //next_generation.push_back(Genom<N, S>::crossover(genoms[shuffle_indices[i]],
                //                                                 genoms[shuffle_indices[i + half]],
                //                                                 distr(rng)));

            std::shuffle(shuffle_indices.begin() + half, shuffle_indices.end(), rng);
            for (size_t i = 0; i < half; i++)
                Genom<N, S>::crossover_inplace(genoms[shuffle_indices[i]],
                        genoms[shuffle_indices[i + half]],
                        &next_generation[i+half],
                        distr(rng));
                //next_generation.push_back(Genom<N, S>::crossover(genoms[shuffle_indices[i]],
                //                                                 genoms[shuffle_indices[i + half]],
                //                                                 distr(rng)));

            for (auto& genom : next_generation)
                genom.mutate(distr(rng), rng);
            genoms.swap(next_generation);
        }

        auto generations() const -> size_t { return generation; }

        auto best() const -> Genom<N, S> const& { return genoms[0]; }

        auto begin() const -> auto
        {
            return genoms.begin();
        }
        auto end() const -> auto
        {
            return genoms.end();
        }

        template <typename N_, typename S_, typename R_>
        friend auto operator<<(ostream& os, GenePool<N_, S_, R_> const& pool) -> ostream&;

    private:
        size_t n_genoms;
        size_t take;
        size_t retain;
        vector<Genom<N, S>> genoms;
        R rng;
        std::uniform_int_distribution<size_t> distr;
        size_t generation;

        // state (to avoid allocations)
        vector<size_t> shuffle_indices;
        vector<Genom<N, S>> next_generation;
    };

    //template<typename Derived, typename N, typename S>
    //class Simu {
    //public:
    //    auto simulate(Genom<N, S>& genom) {
    //        derived().simulate(genom);
    //    }
    //private:
    //    Derived& derived() {
    //        return *static_cast<Derived*>(this);
    //    }
    //};

    // ---------------------------------------------------------------------------------------------
    // ---- ostream operators ----------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    template <typename D>
    auto operator<<(ostream& os, Nucl<D> const& nucl) -> ostream&
    {
        return os << nucl.derived();
    }

    template <typename N, typename S, typename R>
    auto operator<<(ostream& os, GenePool<N, S, R> const& pool) -> ostream&
    {
        os << "GenePool [" << pool.n_genoms << " genoms]\n{\n"
           << "  generations: " << pool.generation << "\n"
           << "  genoms:\n";
        if (pool.genoms.size() > 20)
        {
            for (int i = 0; i < 5; i++)
                os << "    " << pool.genoms[i] << "\n";
            os << "    ...\n";
            for (int i = pool.genoms.size() - 5; i < pool.genoms.size(); i++)
                os << "    " << pool.genoms[i] << "\n";
        }
        else
        {
            for (auto& genom : pool.genoms)
                os << "    " << genom << "\n";
        }
        os << "\n}";
        return os;
    }

    template <typename N, typename S>
    auto operator<<(ostream& os, Genom<N, S> const& g) -> ostream&
    {
        os << "[ ";
        if (g.nucleotides.size() < 30)
        {
            for (auto const& nucl : g.nucleotides)
                os << nucl << ", ";
        }
        else
        {
            for (auto i = 0; i < 5; i++)
                os << g.nucleotides[i] << ", ";
            os << "... ";
            for (auto i = g.nucleotides.size() - 5; i < g.nucleotides.size(); i++)
                os << g.nucleotides[i] << ", ";
        }
        os << "]";
        return os;
    }
}  // namespace ga
