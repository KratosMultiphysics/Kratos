#ifndef VIENNACL_GENERATOR_AUTOTUNE_HPP
#define VIENNACL_GENERATOR_AUTOTUNE_HPP


/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/** @file viennacl/generator/autotune.hpp
 *
 * @brief User interface for the autotuning procedure
*/

#include <ctime>
#include <iomanip>
#include <cmath>
#include <iterator>

#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/infos.hpp"

#include "viennacl/scheduler/forwards.h"
#include "viennacl/generator/generate.hpp"

#include "viennacl/tools/timer.hpp"

namespace viennacl{

  namespace generator{

    namespace autotune{

      /** @brief class for a tuning parameter */
      class tuning_param{
        public:

          /** @brief The constructor
           *
           *  @param values    The set of values which this particular tuning parameter can take
           */
          tuning_param(std::vector<int> const & values) : values_(values){ reset();  }

          /** @brief Returns true if the parameter has reached its maximum value */
          bool is_max() const { return current_ ==  (values_.size()-1); }

          /** @brief Increments the parameter */
          bool inc(){
            ++current_ ;
            if(current_ < values_.size() )
              return false;
            reset();
            return true;
          }

          /** @brief Returns the current value of the parameter */
          int current() const{ return values_[current_]; }

          /** @brief Resets the parameter to its minimum value */
          void reset() { current_ = 0; }

        private:
          std::vector<int> values_;
          unsigned int current_;
      };

      /** @brief Tuning configuration
       *
       *  ConfigType must have a profile_type typedef
       *  ConfigType must implement is_invalid that returns whether or not a given parameter is invalid
       *  ConfigType must implement create_profile that creates a profile_type given a set of parameters
       *
       *  Parameters are stored in a std::map<std::string, viennacl::generator::autotune::tuning_param>
       */
      template<class ConfigType>
      class tuning_config{
        private:
          /** @brief Storage type of the parameters */
          typedef std::map<std::string, viennacl::generator::autotune::tuning_param> params_t;

        public:
          typedef ConfigType config_type;

          /** @brief Accessor for profile_type */
          typedef typename config_type::profile_type profile_type;

          /** @brief Add a tuning parameter to the config */
          void add_tuning_param(std::string const & name, std::vector<int> const & values){
              params_.insert(std::make_pair(name,values));
          }

          /** @brief Returns true if the tuning config has still not explored all its possibilities */
          bool has_next() const{
              bool res = false;
              for(typename params_t::const_iterator it = params_.begin() ; it != params_.end() ; ++it)
                  res = res || !it->second.is_max();
              return res;
          }

          /** @brief Update the parameters of the config */
          void update(){
              for(typename params_t::iterator it = params_.begin() ; it != params_.end() ; ++it)
                  if(it->second.inc()==false)
                      break;
          }

          /** @brief Returns true if the compilation/execution of the underlying profile has an undefined behavior */
          bool is_invalid(viennacl::ocl::device const & dev) const{
              return config_type::is_invalid(dev,params_);
          }

          /** @brief Returns the current profile */
          typename config_type::profile_type get_current(){
              return config_type::create_profile(params_);
          }

          /** @brief Reset the config */
          void reset(){
              for(params_t::iterator it = params_.begin() ; it != params_.end() ; ++it){
                  it->second.reset();
              }
          }

        private:
          params_t params_;
      };

      /** @brief Add the timing value for a given profile and an statement */
      template<class ProfileT>
      double benchmark_impl(viennacl::scheduler::statement const & statement, code_generator::forced_profile_key_type key, ProfileT const & prof, unsigned int n_runs){

        tools::timer t;
        std::list<viennacl::ocl::kernel *> kernels;
        viennacl::generator::code_generator gen;
        gen.force_profile(key, prof);
        gen.add(statement, statement.array()[0]);
        viennacl::generator::get_configured_program(gen, kernels, true);
        viennacl::generator::enqueue(gen);
        viennacl::backend::finish();
        t.start();

        for(unsigned int i = 0 ; i < n_runs ; ++i)
            viennacl::generator::enqueue(gen);
        viennacl::backend::finish();
        return (double)t.get()/n_runs;
      }


      /** @brief Fills a timing map for a given statement and a benchmark configuration
       *
       * @tparam OpT         type of the statement
       * @tparam ConfigType  type of the benchmark configuration
       * @param timings      the timings to fill
       * @param op           the given statement
       * @param key          a key for forcing a particular kernel profile (i.e. to pick profile A for a device which would usually use profile B)
       * @param config       the given configuration
       * @param n_runs       Number of runs for the benchmark
       * @param out          Pointer to output file stream for writing to file (if not NULL)
       */
      template<class ConfigType>
      void benchmark(std::map<double, typename ConfigType::profile_type> * timings, scheduler::statement const & op, code_generator::forced_profile_key_type const & key, tuning_config<ConfigType> & config, unsigned int n_runs, std::ofstream * out){
        viennacl::ocl::device const & dev = viennacl::ocl::current_device();
        unsigned int n_conf = 0;
        while(config.has_next()){
          config.update();
          typename ConfigType::profile_type const & profile = config.get_current();
          if(config.is_invalid(dev) || profile.is_slow(dev))
              continue;
          ++n_conf;
        }
        config.reset();

        unsigned int n = 0;
        while(config.has_next()){
          config.update();
          typename ConfigType::profile_type const & profile = config.get_current();
          if(config.is_invalid(dev) || profile.is_slow(dev))
              continue;
          double percent = (double)n++*100/n_conf;
          double exec_time = benchmark_impl(op,key,profile,n_runs);
          timings->insert(std::make_pair(exec_time, profile));
          std::cout << '\r' << "Autotuning..." << "[" << std::setprecision(2) << std::setfill (' ') << std::setw(6) << std::fixed  << percent << "%" << "]"
                    << " | Best : " << timings->begin()->second << " => " << std::scientific << std::right << std::setprecision(2) << timings->begin()->first << std::flush;
          if(out)
            *out << std::setprecision(3) << std::scientific << exec_time << "," << profile.csv_representation() << std::endl ;
        }
        std::cout << '\r' << "Autotuning..." << "[100.00%]" << std::endl;
      }

    }

  }

}
#endif // AUTOTUNE_HPP
