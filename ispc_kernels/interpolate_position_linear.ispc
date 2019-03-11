// Interpolates positions linearly

export void interpolate_position_linear(uniform float out_x[], uniform float out_y[], uniform float out_z[],
										uniform float in_x0[], uniform float in_y0[], uniform float in_z0[],
										uniform float in_x1[], uniform float in_y1[], uniform float in_z1[],
										uniform int count, uniform float t) {
	uniform float one_minus_t = 1.0f - t;
    foreach (i = 0 ... count) {
        out_x[i] = in_x0[i] * one_minus_t + in_x1[i] * t;
		out_y[i] = in_y0[i] * one_minus_t + in_y1[i] * t;
		out_z[i] = in_z0[i] * one_minus_t + in_z1[i] * t;
    }
}
