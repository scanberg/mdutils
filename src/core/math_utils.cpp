#include <core/math_utils.h>

namespace math {

void generate_halton_sequence(float* dst, int count, int base) {
    for (int i = 0; i < count; i++) {
        dst[i] = halton(i + 1, base);
    }
}

}  // namespace math
