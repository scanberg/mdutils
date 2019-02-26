#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;

uniform float u_radius_scale = 1.0;
uniform vec4 u_color;

layout (location = 0) in vec3  in_position;
layout (location = 1) in float in_radius;
layout (location = 2) in uint  in_mask;

out VS_GS {
    flat vec2 axis_a;
    flat vec2 axis_b;
    flat vec2 center;
    flat float z;
    flat float radius;
} out_geom;

// From Inigo Quilez!
void proj_sphere(in vec4 sphere, 
                 in float fle,
                 out vec2 axis_a,
                 out vec2 axis_b,
                 out vec2 center) {
    vec3  o = sphere.xyz;
    float r2 = sphere.w*sphere.w;
    float z2 = o.z*o.z; 
    float l2 = dot(o,o);
    
    // axis
    axis_a = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-z2)))*vec2( o.x,o.y);
    axis_b = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-l2)))*vec2(-o.y,o.x);
    center = -fle*o.z*o.xy/(z2-r2);
}

void main() {
    vec3 pos = in_position;
    float rad = in_radius * u_radius_scale;
    uint mask = in_mask;
    if (mask == 0U) {
    	out_geom.radius = 0.0;
	} else {
		vec4 view_coord = u_view_mat * vec4(pos, 1.0);
		vec4 view_sphere = vec4(view_coord.xyz, rad);

		float fle = u_proj_mat[1][1]; // Focal length
		float inv_ar = u_proj_mat[0][0] / u_proj_mat[1][1]; // 1.0 / aspect_ratio
		float z = -u_proj_mat[2][2] - u_proj_mat[3][2] / (view_coord.z + rad); // Compute view depth (with bias of radius) 

		vec2 axis_a;
		vec2 axis_b;
		vec2 center;
		proj_sphere(view_sphere, fle, axis_a, axis_b, center);

		vec2 scl = vec2(inv_ar, 1.0);

		out_geom.axis_a = axis_a * scl;
		out_geom.axis_b = axis_b * scl;
		out_geom.center = center * scl;
		out_geom.z = z;
		out_geom.radius = rad;

		gl_Position = view_coord;
	}
}