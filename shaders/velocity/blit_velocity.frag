#version 150 core

uniform mat4 u_curr_clip_to_prev_clip_mat;

in  vec2 uv;
out vec4 out_ss_vel;

void main() {
    vec2 p_uv = uv;
    vec4 p_cs = vec4(uv * 2.0 - 1.0, 1.0, 1.0);
    
	vec4 q_cs = u_curr_clip_to_prev_clip_mat * p_cs;
    vec2 q_uv = (q_cs.xy / q_cs.w) * 0.5 + 0.5;

    out_ss_vel = vec4(p_uv - q_uv, 0, 0);
}