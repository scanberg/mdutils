

export void interpolate_position_linear(uniform float out_x[], uniform float out_y[], uniform float out_z[],
										uniform float in_x0[], uniform float in_y0[], uniform float in_z0[],
										uniform float in_x1[], uniform float in_y1[], uniform float in_z1[],
										uniform int count, uniform float t,
										uniform float ext_x, uniform float ext_y, uniform float ext_z) {
	uniform float one_minus_t = 1.0f - t;
    foreach (i = 0 ... count) {
		float dx = in_x1[i] - in_x0[i];
		float dy = in_y1[i] - in_y0[i];
		float dz = in_z1[i] - in_z0[i];

		float abs_dx = abs(dx);
		float abs_dy = abs(dy);
		float abs_dz = abs(dz);

    }
}
