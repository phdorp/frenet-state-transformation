#ifndef CONSTEXPR_TYPES_H
#define CONSTEXPR_TYPES_H

#include <array>

namespace FrenetTransform
{
    namespace Testing
    {
        template <int Value>
        struct Integral { static constexpr int s_val { Value }; };

        template <typename ValType, int NumVals, std::array<ValType, NumVals> Vals>
        struct ConstVals { static constexpr std::array<ValType, NumVals> s_vals { Vals }; };
    };
};

#endif