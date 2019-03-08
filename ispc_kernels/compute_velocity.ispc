export void compute_velocity(uniform float vel_out[], uniform float curr_pos[], uniform float prev_pos[], uniform float dt, uniform int count) {
    foreach (i = 0 ... count) {
        float cp = curr_pos[i];
        float pp = prev_pos[i];
        vel_out[i] = (cp - pp) * dt;
    }
}
