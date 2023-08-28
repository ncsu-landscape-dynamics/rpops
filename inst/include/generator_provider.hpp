/*
 * Simple random number generator provider class.
 *
 * Copyright (C) 2023 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
 * This file is part of PoPS.

 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef POPS_SIMPLE_GENERATOR_HPP
#define POPS_SIMPLE_GENERATOR_HPP

#include <random>
#include <map>
#include <string>
#include <exception>

#include "config.hpp"

namespace pops {

/**
 * Interface for generator providers.
 *
 * Now used only internally to switch between providers single
 * and multiple seeds.
 */
template<typename Generator>
class RandomNumberGeneratorProviderInterface
{
public:
    virtual void seed(unsigned seed) = 0;
    virtual void seed(const std::map<std::string, unsigned>& seeds) = 0;
    virtual void seed(const Config& config) = 0;
    virtual Generator& disperser_generation() = 0;
    virtual Generator& natural_dispersal() = 0;
    virtual Generator& anthropogenic_dispersal() = 0;
    virtual Generator& establishment() = 0;
    virtual Generator& weather() = 0;
    virtual Generator& movement() = 0;
    virtual Generator& overpopulation() = 0;
    virtual Generator& survival_rate() = 0;
    virtual Generator& soil() = 0;
    virtual ~RandomNumberGeneratorProviderInterface() = default;
};

/**
 * Provider which always supplies only one generator based on one seed.
 *
 * Satisfies UniformRandomBitGenerator, so it can be used in place of
 * standard generators. This makes testing of other components easier
 * as this generator object can be used directly or, more importantly,
 * standard generator can be used in its place.
 */
template<typename Generator>
class SingleGeneratorProvider : public RandomNumberGeneratorProviderInterface<Generator>
{
public:
    /**
     * @brief Seeds the underlying generator
     * @param seed for the underlying generator
     */
    SingleGeneratorProvider(unsigned seed)
    {
        this->seed(seed);
    }

    /* Re-seed the generator */
    void seed(unsigned seed)
    {
        general_generator_.seed(seed);
    }

    /* This overload always throws std::invalid_argument because only one seed is
     * supported. */
    void seed(const std::map<std::string, unsigned>& seeds)
    {
        UNUSED(seeds);
        throw std::invalid_argument(
            "Multiple seeds are not supported by SimpleGeneratorProvider (only one seed is supported)");
    }

    /* Re-seeds the generator
     *
     * Throws std::invalid_argument if configuration contains multiple seeds.
     */
    void seed(const Config& config)
    {
        if (config.multiple_random_seeds) {
            throw std::invalid_argument(
                "Config cannot have multiple_random_seeds set for SimpleGeneratorProvider (only random_seed is supported)");
        }
        if (!config.random_seeds.empty()) {
            throw std::invalid_argument(
                "Config cannot have random_seeds set for SimpleGeneratorProvider (only random_seed is supported)");
        }
        general_generator_.seed(config.random_seed);
    }

    Generator& general()
    {
        return general_generator_;
    }

    Generator& disperser_generation()
    {
        return general();
    }

    Generator& natural_dispersal()
    {
        return general();
    }

    Generator& anthropogenic_dispersal()
    {
        return general();
    }

    Generator& establishment()
    {
        return general();
    }

    Generator& weather()
    {
        return general();
    }

    Generator& movement()
    {
        return general();
    }

    Generator& overpopulation()
    {
        return general();
    }

    Generator& survival_rate()
    {
        return general();
    }

    Generator& soil()
    {
        return general();
    }

    // API to behave like the underlying generator.

    using result_type = typename Generator::result_type;

    static result_type min()
    {
        return Generator::min();
    }

    static result_type max()
    {
        return Generator::max();
    }

    result_type operator()()
    {
        return general_generator_();
    }

    void discard(unsigned long long n)
    {
        general_generator_.discard(n);
    }

private:
    Generator general_generator_;
};

/** Default generator provider to be used for tests and development */
using DefaultSingleGeneratorProvider =
    SingleGeneratorProvider<std::default_random_engine>;

/** Generator provider providing multiple isolated generators
 *
 * All ways of seeding the provider result in multiple independent generators.
 */
template<typename Generator>
class MultiRandomNumberGeneratorProvider
    : public RandomNumberGeneratorProviderInterface<Generator>
{
public:
    /**
     * Seeds first generator with the seed and then each subsequent generator with
     * seed += 1.
     */
    MultiRandomNumberGeneratorProvider(unsigned seed)
    {
        this->seed(seed);
    }

    /**
     * Seeds generators by name.
     */
    MultiRandomNumberGeneratorProvider(const std::map<std::string, unsigned>& seeds)
    {
        this->seed(seeds);
    }

    /**
     * Seeds generators according to the configuration. If named seeds are available,
     * generators are initialized by name, otherwise single seed is incremented for
     * each generator.
     */
    MultiRandomNumberGeneratorProvider(const Config& config)
    {
        this->seed(config);
    }

    /** Re-seed with single value incremented for each generator. */
    void seed(unsigned seed)
    {
        disperser_generation_generator_.seed(seed++);
        natural_dispersal_generator_.seed(seed++);
        anthropogenic_dispersal_generator_.seed(seed++);
        establishment_generator_.seed(seed++);
        weather_generator_.seed(seed++);
        movement_generator_.seed(seed++);
        overpopulation_generator_.seed(seed++);
        survival_rate_generator_.seed(seed++);
        soil_generator_.seed(seed);
    }

    /** Re-seed generators by name. */
    void seed(const std::map<std::string, unsigned>& seeds)
    {
        this->set_seed_by_name(
            seeds, "disperser_generation", disperser_generation_generator_);
        this->set_seed_by_name(
            seeds, "natural_dispersal", natural_dispersal_generator_);
        this->set_seed_by_name(
            seeds, "anthropogenic_dispersal", anthropogenic_dispersal_generator_);
        this->set_seed_by_name(seeds, "establishment", establishment_generator_);
        this->set_seed_by_name(seeds, "weather", weather_generator_);
        this->set_seed_by_name(seeds, "movement", movement_generator_);
        this->set_seed_by_name(seeds, "overpopulation", overpopulation_generator_);
        this->set_seed_by_name(seeds, "survival_rate", survival_rate_generator_);
        this->set_seed_by_name(seeds, "soil", soil_generator_);
    }

    /** Re-seed using named seeds, otherwise increment single seed */
    void seed(const Config& config)
    {
        if (!config.random_seeds.empty()) {
            this->seed(config.random_seeds);
        }
        else {
            this->seed(config.random_seed);
        }
    }

    Generator& disperser_generation()
    {
        return disperser_generation_generator_;
    }

    Generator& natural_dispersal()
    {
        return natural_dispersal_generator_;
    }

    Generator& anthropogenic_dispersal()
    {
        return anthropogenic_dispersal_generator_;
    }

    Generator& establishment()
    {
        return establishment_generator_;
    }

    Generator& weather()
    {
        return weather_generator_;
    }

    Generator& movement()
    {
        return movement_generator_;
    }

    Generator& overpopulation()
    {
        return overpopulation_generator_;
    }

    Generator& survival_rate()
    {
        return survival_rate_generator_;
    }

    Generator& soil()
    {
        return soil_generator_;
    }

private:
    /** Seed a given generator by value associated with the key */
    void set_seed_by_name(
        const std::map<std::string, unsigned>& seeds,
        const char* key,
        Generator& generator)
    {
        try {
            generator.seed(seeds.at(key));
        }
        catch (const std::out_of_range&) {
            throw std::invalid_argument(
                std::string("Seed '") + key + "' missing from the seeds configuration");
        }
    }

    Generator disperser_generation_generator_;
    Generator natural_dispersal_generator_;
    Generator anthropogenic_dispersal_generator_;
    Generator establishment_generator_;
    Generator weather_generator_;
    Generator movement_generator_;
    Generator overpopulation_generator_;
    Generator survival_rate_generator_;
    Generator soil_generator_;
};

/**
 * Generator provider which can switch between using single seed and
 * multiple seeds.
 *
 * Depending on how the object is seeded, result is single generator or
 * multiple independent generators.
 *
 * Satisfies UniformRandomBitGenerator, so it can be used in place of
 * standard generators. This makes testing of other components easier
 * as this generator object can be used directly or, more importantly,
 * standard generator can be used in its place.
 *
 * However, unlike the simple generator for single seed, this will throw
 * an exception if used directly as UniformRandomBitGenerator, but the
 * object was seeded with multiple seeds.
 */
template<typename Generator>
class RandomNumberGeneratorProvider
{
public:
    /**
     * Seeds first generator with the seed and then each subsequent generator with
     * seed += 1. *multi* decides if single generator is created or if multiple
     * isolated generators are created and seed is incremented as needed.
     */
    RandomNumberGeneratorProvider(unsigned seed, bool multi = false) : mutli_(multi)
    {
        if (multi) {
            impl.reset(new MultiRandomNumberGeneratorProvider<Generator>(seed));
        }
        else {
            impl.reset(new SingleGeneratorProvider<Generator>(seed));
        }
    }

    /** Creates multiple isolated generators based on named seeds */
    RandomNumberGeneratorProvider(const std::map<std::string, unsigned>& seeds)
        : impl(new MultiRandomNumberGeneratorProvider<Generator>(seeds)), mutli_(true)
    {}

    /**
     * Result can be multiple independent generators or a single generator
     * based on the configuration.
     */
    RandomNumberGeneratorProvider(const Config& config) : impl(nullptr)
    {
        if (config.multiple_random_seeds) {
            impl.reset(new MultiRandomNumberGeneratorProvider<Generator>(config));
            mutli_ = true;
        }
        else {
            impl.reset(new SingleGeneratorProvider<Generator>(config.random_seed));
            mutli_ = false;
        }
    }

    Generator& disperser_generation()
    {
        return impl->disperser_generation();
    }

    Generator& natural_dispersal()
    {
        return impl->natural_dispersal();
    }

    Generator& anthropogenic_dispersal()
    {
        return impl->anthropogenic_dispersal();
    }

    Generator& establishment()
    {
        return impl->establishment();
    }

    Generator& weather()
    {
        return impl->weather();
    }

    Generator& movement()
    {
        return impl->movement();
    }

    Generator& overpopulation()
    {
        return impl->overpopulation();
    }

    Generator& survival_rate()
    {
        return impl->survival_rate();
    }

    Generator& soil()
    {
        return impl->soil();
    }

    // API to behave like the underlying generator.

    using result_type = typename Generator::result_type;

    static result_type min()
    {
        return Generator::min();
    }

    static result_type max()
    {
        return Generator::max();
    }

    /* Throws std::runtime_error if using multiple isolated generators */
    result_type operator()()
    {
        if (mutli_) {
            std::runtime_error(
                "RandomNumberGeneratorProvider used as a single generator "
                "but it is set to provide multiple isolated generators");
        }
        return impl->disperser_generation().operator()();
    }

    /* Throws std::runtime_error if using multiple isolated generators */
    void discard(unsigned long long n)
    {
        if (mutli_) {
            std::runtime_error(
                "RandomNumberGeneratorProvider used as a single generator "
                "but it is set to provide multiple isolated generators");
        }
        impl->disperser_generation().discard(n);
    }

private:
    std::unique_ptr<RandomNumberGeneratorProviderInterface<Generator>> impl;
    bool mutli_ = false;
};

/**
 * Try constructing a provider
 *
 * Throws the underlying exception. Does nothing on success.
 */
void validate_random_number_generator_provider_config(const Config& config)
{
    MultiRandomNumberGeneratorProvider<std::default_random_engine> provider(config);
    UNUSED(provider);
}

/**
 * Try constructing a provider
 *
 * Throws the underlying exception. Does nothing on success.
 */
void validate_random_number_generator_provider_seeds(
    const std::map<std::string, unsigned>& seeds)
{
    MultiRandomNumberGeneratorProvider<std::default_random_engine> provider(seeds);
    UNUSED(provider);
}

}  // namespace pops

#endif  // POPS_SIMPLE_GENERATOR_HPP
