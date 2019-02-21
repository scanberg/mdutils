#version 150 core
#pragma optionNV(unroll all)

// From http://blog.tuxedolabs.com/2018/05/04/bokeh-depth-of-field-in-single-pass.html

uniform sampler2D uHalfRes; // Half res color (rgb) and coc (a)
uniform sampler2D uColor;   // Image to be processed 
uniform sampler2D uDepth;   // Linear depth, where 1.0 == far plane 

uniform vec2 uPixelSize;    // The size of a pixel: vec2(1.0/width, 1.0/height) 
uniform float uFocusPoint; 
uniform float uFocusScale;

const float GOLDEN_ANGLE = 2.39996323; 
const float MAX_BLUR_SIZE = 15.0; 
const float RAD_SCALE = 1.0; // Smaller = nicer blur, larger = faster

#define APPROX

float getBlurSize(float depth, float focusPoint, float focusScale)
{
	float coc = clamp((1.0 / focusPoint - 1.0 / depth)*focusScale, -1.0, 1.0);
	return abs(coc) * MAX_BLUR_SIZE;
}

vec3 depthOfField(vec2 tex_coord, float focus_point, float focus_scale)
{
	float center_depth  = texture(uDepth, tex_coord).r;
	vec3  center_color  = texture(uColor, tex_coord).rgb;
	float center_coc    = getBlurSize(center_depth, focus_point, focus_scale);
	vec4  color_coc_sum = vec4(center_color, center_coc);

	float contrib_sum   = 1.0;
	float radius        = RAD_SCALE;
	int	  offset_idx    = 0;
	float ang           = 0.0;


	for (; radius < MAX_BLUR_SIZE; ang += GOLDEN_ANGLE)
	{
		vec2 tc = tex_coord + vec2(cos(ang), sin(ang)) * uPixelSize * radius;
		float sample_depth = texture(uDepth, tc).r;
		vec3  sample_color = texture(uColor, tc).rgb;
		float sample_coc   = getBlurSize(sample_depth, focus_point, focus_scale);
		if (sample_depth > center_depth)
			sample_coc = min(sample_coc, center_coc*2.0);

		vec4 sample_color_coc = vec4(sample_color, sample_coc);

		color_coc_sum     += mix(color_coc_sum / contrib_sum, sample_color_coc, smoothstep(radius-0.5, radius+0.5, sample_color_coc.a));
		contrib_sum       += 1.0;
		radius            += RAD_SCALE/radius;

#ifdef APPROX
		if (color_coc_sum.a / contrib_sum > 10.0) break;
#endif
	}

#ifdef APPROX
const float HALF_RES_RAD_SCALE = 2.0;
	for (; radius < MAX_BLUR_SIZE; ang += GOLDEN_ANGLE) {
		vec2 tc = tex_coord + vec2(cos(ang), sin(ang)) * uPixelSize * radius;
		vec4  sample_color_coc = texture(uHalfRes, tc) * vec4(1, 1, 1, MAX_BLUR_SIZE);
		
		float sample_depth = texture(uDepth, tc).r;
		if (sample_depth > center_depth) sample_color_coc.a = clamp(sample_color_coc.a, 0.0, center_coc*2.0);

		color_coc_sum     += mix(color_coc_sum / contrib_sum, sample_color_coc, smoothstep(radius-0.5, radius+0.5, sample_color_coc.a));
		contrib_sum       += 1.0;
		radius            += HALF_RES_RAD_SCALE/radius;
	}
#endif
	return color_coc_sum.rgb /= contrib_sum;
}

in vec2 tc;
out vec4 out_frag;

void main() {
	vec3 dof = depthOfField(tc, uFocusPoint, uFocusScale);
	out_frag = vec4(dof, 1);
}