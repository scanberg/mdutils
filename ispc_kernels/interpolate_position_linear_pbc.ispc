export void interpolate_position_linear(uniform float<4> dst_pos[], uniform float<4> prev_pos[], uniform float<4> next_pos[], uniform int count, uniform float t, uniform ) {
	uniform float one_minus_t = 1.0f - t;
    foreach (i = 0 ... count) {
        dst_pos[i] = prev_pos[i] * one_minus_t + next_pos[i] * t;
    }
}
