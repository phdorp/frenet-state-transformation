#ifndef CONSTEXPR_TYPES_H
#define CONSTEXPR_TYPES_H

#include <array>

namespace FrenetTransform {
namespace Internal {
/**
 * @brief Struct containing an integral type for the definition of test
 * parameters.
 *
 * @tparam Value integral value.
 */
template <int Value> struct Integral { static constexpr int s_val{Value}; };

/**
 * @brief Struct containing an array for for the definition of test parameters.
 *
 * @tparam ValType array element type.
 * @tparam NumVals number of array elements.
 * @tparam Vals array elements.
 */
template <typename ValType, int NumVals, std::array<ValType, NumVals> Vals>
struct ConstVals {
  static constexpr std::array<ValType, NumVals> s_vals{Vals};
};
}; // namespace Internal
}; // namespace FrenetTransform

#endif