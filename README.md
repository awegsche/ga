# ga

C++ genetic algorithm template library.

For a general introduction into GAs, see [wikipedia](https://en.wikipedia.org/wiki/Genetic_algorithm).

## Example

To illustrate the usage, let us try to find the maximum element in a list of integers, using a
genetic algorithm. (_Note_: this example is implemented in the [max_integer](/tests/max_integer.cpp)
test.)

The genes will be a list of integers and the result will be the sum of these integers.

We need to define the following classes:

- A nucleotide (this is the atomic element of the gene, i.e. an integer in our case)

```cpp

// dist represents the range of change of each mutation
std::uniform_int_distribution<int32_t> dist(-5, 5);

class IntNucl : public ga::Nucl<IntNucl> {
    public:

        // ---- Initialisation --------------------------------------------------------------------
        //
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

        static auto crossover(const IntNucl& a, const IntNucl& b) -> IntNucl {
            return IntNucl((a.value_ + b.value_)/2);
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

```

- A scorer class (this holds the score of each genom and ideally also the result)

```cpp

class Scorer {
    public:
        // public score field is needed for the sorting of genoms
        float score;

    public: 
        Scorer() : score(0) {}
        Scorer(float score) : score(score) {}
};

```

- The simulator, this class performs the heavy lifting: calculating the propagation of the simulation.
```cpp

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
```
