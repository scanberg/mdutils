#include <core/math_utils.h>

namespace math {

void generate_halton_sequence(float* dst, int count, int base) {
    for (int i = 0; i < count; i++) {
        dst[i] = halton(i + 1, base);
    }
}

void generate_halton_sequence(vec2* dst, int count, int base_x, int base_y) {
    for (int i = 0; i < count; i++) {
        dst[i].x = halton(i + 1, base_x);
        dst[i].y = halton(i + 1, base_y);
    }
}

}  // namespace math
