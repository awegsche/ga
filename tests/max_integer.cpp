#include <gtest/gtest.h>
#include <numeric>
#include <vector>
#include <random>

#include "../ga.h"

std::uniform_int_distribution<int32_t> dist(-5, 5);

class IntNucl : public ga::Nucl<IntNucl> {
    public:

        // ---- Initialisation --------------------------------------------------------------------
        //
        IntNucl() : value_(0) {}
        IntNucl(int32_t value) : value_(value) {}

        /// All nucleotides need a static random generation initialisation
        template<typename R>
        static auto random(R& rng) -> IntNucl {
            return IntNucl(dist(rng));
        }

        // ---- Genetic Operators -----------------------------------------------------------------
        //
        template<typename R>
        auto mutate(R& rng) {
            value_ += dist(rng);
        }

        static void crossover_inplace(const IntNucl& a, const IntNucl& b, IntNucl* c) {
            c->value_ = (a.value_ + b.value_)/2;
        }

        // ---- Accessors -------------------------------------------------------------------------
        //
        auto value() const -> int32_t {
            return value_;
        }

        // ---- Printing (needed if you want to print the whole pool) -----------------------------
        //
        friend ostream& operator<<(ostream& os, const IntNucl& nucl) {
            return os << nucl.value();
        }

    private:
        int32_t value_;
};

class Scorer {
    public:
        // public score field is needed for the sorting of genoms
        float score;

    public: 
        Scorer() : score(0) {}
        Scorer(float score) : score(score) {}
};

class Simulator {
    public:
        Simulator(std::vector<int> integers) : integers_(integers) {
            max_ = *std::max_element(integers_.begin(), integers_.end());
        }


        // simulation function, required
        void simulate(ga::Genom<IntNucl, Scorer> &genom) const {
            int m = get_value(genom);
            genom.set_scorer(m == max_ ? Scorer(10.0f) : Scorer(1.0f / std::abs(m - max_)));
        };

        auto get_value(ga::Genom<IntNucl, Scorer> const &genom) const -> int32_t {
            return std::accumulate(genom.get_nucleotides().begin(), genom.get_nucleotides().end(), 0,
                    [](int sum, const IntNucl& nucl) {
                    return sum + nucl.value();
                    });
        }

    private:
        std::vector<int> integers_;
        int max_ = 0;

};

TEST(GetMaxInteger, BasicTest) {

    std::cout << "test" << std::endl;

    ga::GenePool<IntNucl, Scorer, std::default_random_engine> pool{10, 7, 0.3};
    Simulator sim{{1, 2, 3, 4, 5, 6}};

    pool.simulate(sim);
    std::cout << pool << std::endl;

    for(int i = 0; i < 100; ++i){
        pool.select();
        pool.simulate(sim);

        if (pool.best().get_scorer().score == 10.0f)
            break;

        std::cout << pool << std::endl;
    }

    ASSERT_EQ(sim.get_value(pool.best()) , 6);
    ASSERT_EQ(pool.best().get_scorer().score , 10.0f);
}
