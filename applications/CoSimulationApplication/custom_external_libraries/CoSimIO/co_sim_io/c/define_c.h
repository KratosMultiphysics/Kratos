/*   ______     _____ _           ________
    / ____/___ / ___/(_)___ ___  /  _/ __ |
   / /   / __ \\__ \/ / __ `__ \ / // / / /
  / /___/ /_/ /__/ / / / / / / // // /_/ /
  \____/\____/____/_/_/ /_/ /_/___/\____/
  Kratos CoSimulationApplication

  License:         BSD License, see license.txt

  Main authors:    Philipp Bucher (https://github.com/philbucher)
*/

#ifndef CO_SIM_IO_C_DEFINE_INCLUDED
#define CO_SIM_IO_C_DEFINE_INCLUDED

/*nodiscard is part of C23:
https://en.cppreference.com/w/c/language/attributes/nodiscard
hence using custom solution
*/
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define CO_SIM_IO_NODISCARD __attribute__((warn_unused_result))
#else
#define CO_SIM_IO_NODISCARD
#endif

#endif /* CO_SIM_IO_C_DEFINE_INCLUDED */
