/*-----------------------------------------------------------------------
  Copyright (c) 2014, NVIDIA. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Neither the name of its contributors may be used to endorse
     or promote products derived from this software without specific
     prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-----------------------------------------------------------------------*/

// Shaders for HBAO are taken from nVidias examples and are copyright protected as stated above

#include "postprocessing_utils.h"
#include <core/types.h>
#include <core/common.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <gfx/gl_utils.h>
#include <stdio.h>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace postprocessing {

// @TODO: Use half-res render targets for SSAO
// @TODO: Use shared textures for all postprocessing operations
// @TODO: Use some kind of unified pipeline for all post processing operations

static struct {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint v_shader_fs_quad = 0;
    GLuint tex_width = 0;
    GLuint tex_height = 0;

	struct {
		GLuint fbo = 0;
		GLuint tex_color[2] = { 0,0 };
		GLuint tex_temporal_buffer[2] = { 0,0 };	// These are dedicated and cannot be use as intermediate buffers by other shaders
	} targets;

	struct {
		GLuint fbo = 0;
		GLuint tex_tilemax = 0;
		GLuint tex_neighbormax = 0;
		int32 tex_width = 0;
		int32 tex_height = 0;
	} velocity;

    struct {
        GLuint fbo = 0;
        GLuint program = 0;
        struct {
            GLuint color_coc = 0;
        } tex;
        struct {
            GLint tex_depth = -1;
            GLint tex_color = -1;
            GLint focus_point = -1;
            GLint focus_scale = -1;
        } uniform_loc;
    } half_res;

    struct {
        GLuint fbo = 0;
        GLuint texture = 0;
        GLuint program_persp = 0;
        GLuint program_ortho = 0;
        struct {
            GLint clip_info = -1;
            GLint tex_depth = -1;
        } uniform_loc;
    } linear_depth;

    struct {
        GLuint tex_random = 0;
        GLuint ubo_hbao_data = 0;

        struct {
            GLuint fbo = 0;
            GLuint texture = 0;
            GLuint program = 0;

            struct {
                GLint control_buffer = -1;
                GLint tex_linear_depth = -1;
                GLint tex_normal = -1;
                GLint tex_random = -1;
            } uniform_loc;
        } hbao;

        struct {
            GLuint fbo = 0;
            GLuint texture = 0;
            GLuint program_first = 0;
            GLuint program_second = 0;
            struct {
                GLint sharpness = -1;
                GLint inv_res_dir = -1;
                GLint texture = -1;
            } uniform_loc;
        } blur;
    } ssao;

    struct {
        GLuint program = 0;
        struct {
            GLint tex_half_res = -1;
            GLint tex_color = -1;
            GLint tex_depth = -1;
            GLint pixel_size = -1;
            GLint focus_point = -1;
            GLint focus_scale = -1;
        } uniform_loc;
    } bokeh_dof;

    struct {
        GLuint program = 0;
    } bloom;

    struct {
        GLuint program = 0;
        struct {
            GLint mode = -1;
            GLint tex_color = -1;
        } uniform_loc;
    } tonemapping;

	struct {
		GLuint program = 0;
		struct {
			GLint tex_linear_depth = -1;
			GLint tex_main = -1;
			GLint tex_prev = -1;
			GLint tex_vel = -1;
			GLint tex_vel_neighbormax = -1;

			GLint texel_size = -1;

			GLint sin_time = -1;
			GLint feedback_min = -1;
			GLint feedback_max = -1;
			GLint motion_scale = -1;

			GLint jitter_uv = -1;
		} uniform_loc;
	} temporal;

} gl;

static const char* v_shader_src_fs_quad = R"(
#version 150 core

out vec2 tc;

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	tc = gl_Position.xy * 0.5 + 0.5;
}
)";

static const char* f_shader_src_linearize_depth = R"(
#ifndef PERSPECTIVE
#define PERSPECTIVE 1
#endif

// z_n * z_f,  z_n - z_f,  z_f, *not used*
uniform vec4 u_clip_info;
uniform sampler2D u_tex_depth;

float ReconstructCSZ(float d, vec4 clip_info) {
#ifdef PERSPECTIVE
    return (clip_info[0] / (d*clip_info[1] + clip_info[2]));
#else
    return (clip_info[1] + clip_info[2] - d*clip_info[1]);
#endif
}

out vec4 out_frag;

void main() {
  float d = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
  out_frag = vec4(ReconstructCSZ(d, u_clip_info));
}
)";

static bool setup_program(GLuint* program, const char* name, const char* f_shader_src, const char* defines = nullptr) {
    ASSERT(program);
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    auto f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    if (defines) {
        const char* sources[2] = {defines, f_shader_src};
        glShaderSource(f_shader, 2, sources, 0);
    } else {
        glShaderSource(f_shader, 1, &f_shader_src, 0);
    }

    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Error while compiling %s shader:\n%s", name, buffer);
        return false;
    }

    if (!*program) {
        *program = glCreateProgram();
    }

    glAttachShader(*program, gl.v_shader_fs_quad);
    glAttachShader(*program, f_shader);
    glLinkProgram(*program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, *program)) {
        LOG_ERROR("Error while linking %s program:\n%s", name, buffer);
        return false;
    }

    glDetachShader(*program, gl.v_shader_fs_quad);
    glDetachShader(*program, f_shader);
    glDeleteShader(f_shader);

    return true;
}

static bool is_orthographic_proj_matrix(const mat4& proj_mat) { return math::length2(vec2(proj_mat[3])) > 0.f; }

namespace ssao {
#ifndef AO_RANDOM_TEX_SIZE
#define AO_RANDOM_TEX_SIZE 4
#endif

#ifndef AO_MAX_SAMPLES
#define AO_MAX_SAMPLES 1
#endif

#ifndef AO_DIRS
#define AO_DIRS 8
#endif

#ifndef AO_SAMPLES
#define AO_SAMPLES 4
#endif

#ifndef AO_BLUR
#define AO_BLUR 1
#endif

static GLuint fbo_hbao = 0;
static GLuint fbo_blur = 0;

static GLuint tex_random = 0;
static GLuint tex_hbao = 0;
static GLuint tex_blur = 0;

static GLuint prog_hbao = 0;
static GLuint prog_blur_first = 0;
static GLuint prog_blur_second = 0;

static GLint uniform_block_index_hbao_control_buffer = -1;
static GLint uniform_loc_hbao_tex_linear_depth = -1;
static GLint uniform_loc_hbao_tex_normal = -1;
static GLint uniform_loc_hbao_tex_random = -1;

static GLint uniform_loc_blur_sharpness = -1;
static GLint uniform_loc_blur_inv_res_dir = -1;
static GLint uniform_loc_blur_tex_ao = -1;
static GLint uniform_loc_blur_tex_linear_depth = -1;

static GLuint ubo_hbao_data = 0;

static GLuint tex_width;
static GLuint tex_height;

struct HBAOData {
    float radius_to_screen;
    float r2;
    float neg_inv_r2;
    float n_dot_v_bias;

    vec2 inv_full_res;
    vec2 inv_quarter_res;

    float ao_multiplier;
    float pow_exponent;
    vec2 _pad0;

    vec4 proj_info;

    vec2 proj_scale;
    int proj_ortho;
    float rnd;
};

static const char* f_shader_src_hbao = R"(
#pragma optionNV(unroll all)

#ifndef AO_BLUR
#define AO_BLUR 0
#endif

#define M_PI 3.14159265

#ifndef AO_STEPS
#define AO_STEPS 4
#endif

#ifndef AO_DIRS
#define AO_DIRS 8
#endif

#ifndef AO_USE_NORMAL
#define AO_USE_NORMAL 1
#endif

#ifndef AO_RANDOM_TEX_SIZE
#define AO_RANDOM_TEX_SIZE 4
#endif

struct HBAOData {
  float   radius_to_screen;
  float   r2;
  float   neg_inv_r2;
  float   n_dot_v_bias;
 
  vec2    inv_full_res;
  vec2    inv_quarter_res;
  
  float   ao_multiplier;
  float   pow_exponent;
  vec2    _pad0;
  
  vec4    proj_info;

  vec2    proj_scale;
  int     proj_ortho;
  float   rnd;
  
  //vec4    offsets[AO_RANDOM_TEX_SIZE*AO_RANDOM_TEX_SIZE];
  //vec4    jitters[AO_RANDOM_TEX_SIZE*AO_RANDOM_TEX_SIZE];
};

// tweakables
const float NUM_STEPS = AO_STEPS;
const float NUM_DIRECTIONS = AO_DIRS; // tex_random/jitter initialization depends on this

layout(std140) uniform u_control_buffer {
  HBAOData control;
};

uniform sampler2D u_tex_linear_depth;
uniform sampler2D u_tex_normal;
uniform sampler2D u_tex_random;

in vec2 tc;
out vec4 out_frag;

void OutputColor(vec4 color) {
  out_frag = color;
}

//----------------------------------------------------------------------------------

vec3 UVToView(vec2 uv, float eye_z) {
  return vec3((uv * control.proj_info.xy + control.proj_info.zw) * (control.proj_ortho != 0 ? 1. : eye_z), eye_z);
}

vec3 FetchViewPos(vec2 uv, float lod) {
  float ViewDepth = textureLod(u_tex_linear_depth, uv, lod).x;
  return UVToView(uv, ViewDepth);
}

vec3 MinDiff(vec3 P, vec3 Pr, vec3 Pl) {
  vec3 V1 = Pr - P;
  vec3 V2 = P - Pl;
  return (dot(V1,V1) < dot(V2,V2)) ? V1 : V2;
}

vec3 DecodeNormal(vec2 enc) {
    vec2 fenc = enc*4-2;
    float f = dot(fenc,fenc);
    float g = sqrt(1-f/4.0);
    vec3 n;
    n.xy = fenc*g;
    n.z = 1-f/2.0;
    return n;
}

vec3 FetchViewNormal(vec2 uv) {
	vec2 enc = textureLod(u_tex_normal, uv, 0).xy;
	vec3 n = DecodeNormal(enc);
	return n * vec3(1,1,-1);
}

//----------------------------------------------------------------------------------
float Falloff(float DistanceSquare) {
  // 1 scalar mad instruction
  return DistanceSquare * control.neg_inv_r2 + 1.0;
}

//----------------------------------------------------------------------------------
// P = view-space position at the kernel center
// N = view-space normal at the kernel center
// S = view-space position of the current sample
//----------------------------------------------------------------------------------
float ComputeAO(vec3 P, vec3 N, vec3 S) {
  vec3 V = S - P;
  float VdotV = dot(V, V);
  float NdotV = dot(N, V) * 1.0/sqrt(VdotV);

  // Use saturate(x) instead of max(x,0.f) because that is faster on Kepler
  return clamp(NdotV - control.n_dot_v_bias,0,1) * clamp(Falloff(VdotV),0,1);
}

//----------------------------------------------------------------------------------
vec2 RotateDirection(vec2 Dir, vec2 CosSin) {
  return vec2(Dir.x*CosSin.x - Dir.y*CosSin.y,
              Dir.x*CosSin.y + Dir.y*CosSin.x);
}

//----------------------------------------------------------------------------------
vec4 GetJitter() {
  // (cos(Alpha),sin(Alpha),rand1,rand2)
  return textureLod(u_tex_random, (gl_FragCoord.xy / AO_RANDOM_TEX_SIZE), 0);
}

//----------------------------------------------------------------------------------
float ComputeCoarseAO(vec2 FullResUV, float RadiusPixels, vec4 Rand, vec3 ViewPosition, vec3 ViewNormal) {
  // Divide by NUM_STEPS+1 so that the farthest samples are not fully attenuated
  float StepSizePixels = RadiusPixels / (NUM_STEPS + 1);

  float AngleOffset = control.rnd * 2.0 * M_PI;
  const float Alpha = 2.0 * M_PI / NUM_DIRECTIONS;
  float AO = 0;

  for (float DirectionIndex = 0; DirectionIndex < NUM_DIRECTIONS; ++DirectionIndex)
  {
    float Angle = AngleOffset + Alpha * DirectionIndex;

    // Compute normalized 2D direction
    vec2 Direction = RotateDirection(vec2(cos(Angle), sin(Angle)), Rand.xy);

    // Jitter starting sample within the first step
    float RayPixels = (Rand.z * StepSizePixels + 1.0);

    for (float StepIndex = 0; StepIndex < NUM_STEPS; ++StepIndex)
    {
      vec2 SnappedUV = round(RayPixels * Direction) * control.inv_full_res + FullResUV;
      vec3 S = FetchViewPos(SnappedUV, StepIndex);
	  vec3 N = ViewNormal;

      RayPixels += StepSizePixels;

      AO += ComputeAO(ViewPosition, N, S);
    }
  }

  AO *= control.ao_multiplier / (NUM_DIRECTIONS * NUM_STEPS);

  return clamp(1.0 - AO, 0, 1);
}

//----------------------------------------------------------------------------------
void main() {
  vec2 uv = tc;
  vec3 ViewPosition = FetchViewPos(uv, 0);

  // Reconstruct view-space normal from nearest neighbors
#if AO_USE_NORMAL
  vec3 ViewNormal = FetchViewNormal(uv);
#else
  vec3 ViewNormal = vec3(0,0,-1);
#endif

  // Compute projection of disk of radius control.R into screen space
  float RadiusPixels = control.radius_to_screen / (control.proj_ortho != 0 ? 1.0 : ViewPosition.z);

  // Get jitter vector for the current full-res pixel
  vec4 Rand = GetJitter();

  float AO = ComputeCoarseAO(uv, RadiusPixels, Rand, ViewPosition, ViewNormal);

#if AO_BLUR
  OutputColor(vec4(pow(AO, control.pow_exponent), ViewPosition.z, 0, 1));
#else
  OutputColor(vec4(vec3(pow(AO, control.pow_exponent)), 1));
#endif 
}
)";

static const char* f_shader_src_hbao_blur = R"(
#pragma optionNV(unroll all)

const float KERNEL_RADIUS = 5;
  
uniform float u_sharpness;
uniform vec2  u_inv_res_dir; // either set x to 1/width or y to 1/height
uniform sampler2D u_tex_ao;
uniform sampler2D u_tex_linear_depth;

in vec2 tc;
out vec4 out_frag;

#ifndef AO_BLUR_PRESENT
#define AO_BLUR_PRESENT 1
#endif

//-------------------------------------------------------------------------

float BlurFunction(vec2 uv, float r, float center_c, float center_d, inout float w_total)
{
  float c = texture(u_tex_ao, uv).x;
  float d = texture(u_tex_linear_depth, uv).x;
  
  const float BlurSigma = float(KERNEL_RADIUS) * 0.5;
  const float BlurFalloff = 1.0 / (2.0*BlurSigma*BlurSigma);
  
  float ddiff = (d - center_d) * u_sharpness;
  float w = exp2(-r*r*BlurFalloff - ddiff*ddiff);
  w_total += w;

  return c*w;
}

void main()
{
  float center_c = texture(u_tex_ao, tc).x;
  float center_d = texture(u_tex_linear_depth, tc).x;
  
  float c_total = center_c;
  float w_total = 1.0;
  
  for (float r = 1; r <= KERNEL_RADIUS; ++r)
  {
    vec2 uv = tc + u_inv_res_dir * r;
    c_total += BlurFunction(uv, r, center_c, center_d, w_total);  
  }
  
  for (float r = 1; r <= KERNEL_RADIUS; ++r)
  {
    vec2 uv = tc - u_inv_res_dir * r;
    c_total += BlurFunction(uv, r, center_c, center_d, w_total);  
  }
  
#if AO_BLUR_PRESENT
  out_frag = vec4(vec3(c_total/w_total), 1);
#else
  out_frag = vec4(c_total/w_total, center_d, 0, 1);
#endif
}
)";

void setup_ubo_hbao_data(GLuint ubo, int width, int height, const mat4& proj_mat, float intensity, float radius, float bias) {
    ASSERT(ubo);
    constexpr float METERS_TO_VIEWSPACE = 1.f;
    const float* proj_data = &proj_mat[0][0];

    static float rnd_data[16] = {0};
    static int idx = 0;
    if (rnd_data[0] == 0.f) {
        math::generate_halton_sequence(rnd_data, 16, 2);
    }
    idx = (idx + 1) % 16;
    float rnd = rnd_data[idx];

    bool is_ortho = is_orthographic_proj_matrix(proj_mat);

    vec4 proj_info;
    float proj_scl;
    if (!is_ortho) {
        proj_info = vec4(2.0f / (proj_data[4 * 0 + 0]),                          // (x) * (R - L)/N
                         2.0f / (proj_data[4 * 1 + 1]),                          // (y) * (T - B)/N
                         -(1.0f - proj_data[4 * 2 + 0]) / proj_data[4 * 0 + 0],  // L/N
                         -(1.0f + proj_data[4 * 2 + 1]) / proj_data[4 * 1 + 1]   // B/N
        );

        // proj_scl = float(height) / (math::tan(fovy * 0.5f) * 2.0f);
        proj_scl = float(height) * proj_data[4 * 1 + 1] * 0.5f;
        proj_scl = float(height) / (math::tan(math::PI / 8.f) * 2.0f);
    } else {
        proj_info = vec4(2.0f / (proj_data[4 * 0 + 0]),                          // ((x) * R - L)
                         2.0f / (proj_data[4 * 1 + 1]),                          // ((y) * T - B)
                         -(1.0f + proj_data[4 * 3 + 0]) / proj_data[4 * 0 + 0],  // L
                         -(1.0f - proj_data[4 * 3 + 1]) / proj_data[4 * 1 + 1]   // B
        );
        proj_scl = float(height) / proj_info[1];
    }

    float r = radius * METERS_TO_VIEWSPACE;

    HBAOData data;
    data.radius_to_screen = r * 0.5f * proj_scl;
    data.r2 = r * r;
    data.neg_inv_r2 = -1.f / (r * r);
    data.n_dot_v_bias = math::clamp(bias, 0.f, 1.f - math::EPSILON);

    data.inv_full_res = vec2(1.f / float(width), 1.f / float(height));
    data.inv_quarter_res = vec2(1.f / float((width + 3) / 4), 1.f / float((height + 3) / 4));

    data.ao_multiplier = 1.f / (1.f - data.n_dot_v_bias);
    data.pow_exponent = math::max(intensity, 0.f);

    data.proj_info = proj_info;
    // @NOTE: Is one needed?
    data.proj_scale = vec2(proj_scl);
    data.proj_ortho = is_ortho ? 1 : 0;
    data.rnd = rnd;

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(HBAOData), &data, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void initialize_rnd_tex(GLuint rnd_tex, int num_direction) {
    ASSERT(AO_MAX_SAMPLES == 1);
    constexpr int buffer_size = AO_RANDOM_TEX_SIZE * AO_RANDOM_TEX_SIZE * AO_MAX_SAMPLES;
    signed short buffer[buffer_size * 4];

    vec2 rnd_vals[buffer_size];
    math::generate_halton_sequence(rnd_vals, buffer_size, 2, 3);
    /*
for (int i = 0; i < buffer_size; i++) {
    rnd_vals[i] = vec2(math::rnd(), math::rnd());
}
    */

    for (int i = 0; i < buffer_size; i++) {
#define SCALE ((1 << 15))
        float rand1 = rnd_vals[i].x;
        float rand2 = rnd_vals[i].y;
        float angle = 2.f * math::PI * rand1 / (float)num_direction;

        buffer[i * 4 + 0] = (signed short)(SCALE * math::cos(angle));
        buffer[i * 4 + 1] = (signed short)(SCALE * math::sin(angle));
        buffer[i * 4 + 2] = (signed short)(SCALE * rand2);
        buffer[i * 4 + 3] = (signed short)(SCALE * 0);
#undef SCALE
    }

    // @TODO: If MSAA and AO_MAX_SAMPLES > 1, then this probably has to go into a texture array
    glBindTexture(GL_TEXTURE_2D, rnd_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16_SNORM, AO_RANDOM_TEX_SIZE, AO_RANDOM_TEX_SIZE, 0, GL_RGBA, GL_SHORT, buffer);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

float compute_sharpness(float radius) { return 30.f / math::sqrt(radius); }

void initialize(int width, int height) {
    // @TODO: dynamically generate this
    const char* defines_hbao = R"(
        #version 150 core
        #define AO_RANDOM_TEX_SIZE 4
        #define PERSPECTIVE 1
        #define AO_BLUR 1
        #define AO_STEPS 4
        #define AO_DIRS 8
        #define AO_USE_NORMAL 1
    )";

    const char* defines_blur_first = R"(
        #version 150 core
        #define AO_BLUR_PRESENT 0
    )";

    const char* defines_blur_second = R"(
        #version 150 core
        #define AO_BLUR_PRESENT 1
    )";

    setup_program(&prog_hbao, "hbao", f_shader_src_hbao, defines_hbao);
    setup_program(&prog_blur_first, "hbao first blur", f_shader_src_hbao_blur, defines_blur_first);
    setup_program(&prog_blur_second, "hbao second blur", f_shader_src_hbao_blur, defines_blur_second);

    uniform_block_index_hbao_control_buffer = glGetUniformBlockIndex(prog_hbao, "u_control_buffer");
    uniform_loc_hbao_tex_linear_depth = glGetUniformLocation(prog_hbao, "u_tex_linear_depth");
    uniform_loc_hbao_tex_normal = glGetUniformLocation(prog_hbao, "u_tex_normal");
    uniform_loc_hbao_tex_random = glGetUniformLocation(prog_hbao, "u_tex_random");

    ASSERT(uniform_block_index_hbao_control_buffer != -1);
    ASSERT(uniform_loc_hbao_tex_linear_depth != -1);
    ASSERT(uniform_loc_hbao_tex_normal != -1);
    ASSERT(uniform_loc_hbao_tex_random != -1);

    uniform_loc_blur_sharpness = glGetUniformLocation(prog_blur_second, "u_sharpness");
    uniform_loc_blur_inv_res_dir = glGetUniformLocation(prog_blur_second, "u_inv_res_dir");
    uniform_loc_blur_tex_ao = glGetUniformLocation(prog_blur_second, "u_tex_ao");
	uniform_loc_blur_tex_linear_depth = glGetUniformLocation(prog_blur_second, "u_tex_linear_depth");

    ASSERT(uniform_loc_blur_sharpness != -1);
    ASSERT(uniform_loc_blur_inv_res_dir != -1);
    ASSERT(uniform_loc_blur_tex_ao != -1);
	ASSERT(uniform_loc_blur_tex_linear_depth != -1);

    if (!fbo_hbao) glGenFramebuffers(1, &fbo_hbao);
    if (!fbo_blur) glGenFramebuffers(1, &fbo_blur);

    if (!tex_random) glGenTextures(1, &tex_random);
    if (!tex_hbao) glGenTextures(1, &tex_hbao);
    if (!tex_blur) glGenTextures(1, &tex_blur);

    initialize_rnd_tex(tex_random, AO_DIRS);

    glBindTexture(GL_TEXTURE_2D, tex_hbao);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, tex_blur);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, fbo_hbao);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_hbao, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, fbo_blur);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_blur, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    if (!ubo_hbao_data) glGenBuffers(1, &ubo_hbao_data);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_hbao_data);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(HBAOData), nullptr, GL_DYNAMIC_DRAW);

    tex_width = width;
    tex_height = height;
}

void shutdown() {
    if (fbo_hbao) glDeleteFramebuffers(1, &fbo_hbao);
    if (fbo_blur) glDeleteFramebuffers(1, &fbo_blur);

    if (tex_random) glDeleteTextures(1, &tex_random);
    if (tex_blur) glDeleteTextures(1, &tex_blur);

    if (ubo_hbao_data) glDeleteBuffers(1, &ubo_hbao_data);

    if (prog_hbao) glDeleteProgram(prog_hbao);
    if (prog_blur_second) glDeleteProgram(prog_blur_second);
    if (prog_blur_first) glDeleteProgram(prog_blur_first);
}

}  // namespace ssao

namespace highlight {

static struct {
    GLuint program = 0;
    GLuint selection_texture = 0;
    struct {
        GLint texture_atom_idx = -1;
        GLint buffer_selection = -1;
    } uniform_loc;
} highlight;

void initialize() {
    String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/highlight.frag");
    defer { FREE(f_shader_src); };
    setup_program(&highlight.program, "highlight", f_shader_src);
    if (!highlight.selection_texture) glGenTextures(1, &highlight.selection_texture);
    highlight.uniform_loc.texture_atom_idx = glGetUniformLocation(highlight.program, "u_texture_atom_idx");
    highlight.uniform_loc.buffer_selection = glGetUniformLocation(highlight.program, "u_buffer_selection");
}

void shutdown() {
    if (highlight.program) glDeleteProgram(highlight.program);
}
}  // namespace highlight

namespace deferred {

static struct {
    GLuint program = 0;
    struct {
        GLint texture_depth = -1;
        GLint texture_color = -1;
        GLint texture_normal = -1;
        GLint inv_proj_mat = -1;
    } uniform_loc;
} deferred;

void initialize() {
    String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/deferred_shading.frag");
    defer { FREE(f_shader_src); };
    setup_program(&deferred.program, "deferred", f_shader_src);
    deferred.uniform_loc.texture_depth = glGetUniformLocation(deferred.program, "u_texture_depth");
    deferred.uniform_loc.texture_color = glGetUniformLocation(deferred.program, "u_texture_color");
    deferred.uniform_loc.texture_normal = glGetUniformLocation(deferred.program, "u_texture_normal");
    deferred.uniform_loc.inv_proj_mat = glGetUniformLocation(deferred.program, "u_inv_proj_mat");
}

void shutdown() {
    if (deferred.program) glDeleteProgram(deferred.program);
}
}  // namespace deferred

namespace tonemapping {

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
    } uniform_loc;
} passthrough;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} exposure_gamma;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} hejl_dawsson;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} filmic;

void initialize() {
    {
        // PASSTHROUGH
        String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/tonemap/passthrough.frag");
        defer { FREE(f_shader_src); };

        setup_program(&passthrough.program, "tonemap_passthrough", f_shader_src);
        passthrough.uniform_loc.texture = glGetUniformLocation(passthrough.program, "u_texture");
    }
    {
        // EXPOSURE GAMMA
        String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/tonemap/exposure_gamma.frag");
        defer { FREE(f_shader_src); };

        setup_program(&exposure_gamma.program, "tonemap_exposure_gamma", f_shader_src);
        exposure_gamma.uniform_loc.texture = glGetUniformLocation(exposure_gamma.program, "u_texture");
        exposure_gamma.uniform_loc.exposure = glGetUniformLocation(exposure_gamma.program, "u_exposure");
        exposure_gamma.uniform_loc.gamma = glGetUniformLocation(exposure_gamma.program, "u_gamma");
    }
    {
        // UNCHARTED
        String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/tonemap/uncharted.frag");
        defer { FREE(f_shader_src); };

        setup_program(&filmic.program, "tonemap_filmic", f_shader_src);
        filmic.uniform_loc.texture = glGetUniformLocation(filmic.program, "u_texture");
        filmic.uniform_loc.exposure = glGetUniformLocation(filmic.program, "u_exposure");
        filmic.uniform_loc.gamma = glGetUniformLocation(filmic.program, "u_gamma");
    }
}

void shutdown() {
    if (passthrough.program) glDeleteProgram(passthrough.program);
    if (exposure_gamma.program) glDeleteProgram(exposure_gamma.program);
    if (filmic.program) glDeleteProgram(filmic.program);
}

}  // namespace tonemapping

namespace blit {
static GLuint program = 0;
static GLint uniform_loc_texture = -1;
static const char* f_shader_src = R"(
#version 150 core

uniform sampler2D u_texture;

in vec2 tc;
out vec4 out_frag;

void main() {
	out_frag = texture(u_texture, tc);
}
)";

void initialize() {
    if (!program) setup_program(&program, "render_texture", f_shader_src);
    uniform_loc_texture = glGetUniformLocation(program, "u_texture");
}

void shutdown() {
    if (program) glDeleteProgram(program);
}
}  // namespace blit

namespace velocity {
#define VEL_TILE_SIZE 20

	struct {
		GLuint program = 0;
		struct {
			GLint curr_clip_to_prev_clip_mat = -1;
		} uniform_loc;
	} blit_velocity;

	struct {
		GLuint program = 0;
		struct {
			GLint tex_vel;
			GLint tex_vel_texel_size;
		} uniform_loc;
	} blit_tilemax;

	struct {
		GLuint program = 0;
		struct {
			GLint tex_vel;
			GLint tex_vel_texel_size;
		} uniform_loc;
	} blit_neighbormax;

	void initialize() {
		{
			String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/velocity/blit_velocity.frag");
			defer{ free_string(&f_shader_src); };
			setup_program(&blit_velocity.program, "screen-space velocity", f_shader_src);
			blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat = glGetUniformLocation(blit_velocity.program, "u_curr_clip_to_prev_clip_mat");
		}
		{
			//const char* defines = { "#version 150 core\n#define TILE_SIZE " TOSTRING(VEL_TILE_SIZE) };
			String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/velocity/blit_tilemax.frag");
			defer{ free_string(&f_shader_src); };
			setup_program(&blit_tilemax.program, "tilemax", f_shader_src);
			blit_tilemax.uniform_loc.tex_vel = glGetUniformLocation(blit_tilemax.program, "u_tex_vel");
			blit_tilemax.uniform_loc.tex_vel_texel_size = glGetUniformLocation(blit_tilemax.program, "u_tex_vel_texel_size");
		}
		{
			String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/velocity/blit_neighbormax.frag");
			defer{ free_string(&f_shader_src); };
			setup_program(&blit_neighbormax.program, "neighbormax", f_shader_src);
			blit_neighbormax.uniform_loc.tex_vel = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel");
			blit_neighbormax.uniform_loc.tex_vel_texel_size = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel_texel_size");
		}
	}

	void shutdown() {
		if (blit_velocity.program) glDeleteProgram(blit_velocity.program);
	}
}

namespace temporal {
void initialize() {
	{
		String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/temporal.frag");
		defer{ free_string(&f_shader_src); };
		setup_program(&gl.temporal.program, "temporal aa + motion-blur", f_shader_src);

		gl.temporal.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.temporal.program, "u_tex_linear_depth");
		gl.temporal.uniform_loc.tex_main = glGetUniformLocation(gl.temporal.program, "u_tex_main");
		gl.temporal.uniform_loc.tex_prev = glGetUniformLocation(gl.temporal.program, "u_tex_prev");
		gl.temporal.uniform_loc.tex_vel = glGetUniformLocation(gl.temporal.program, "u_tex_vel");
		gl.temporal.uniform_loc.tex_vel_neighbormax = glGetUniformLocation(gl.temporal.program, "u_tex_vel_neighbormax");
		gl.temporal.uniform_loc.texel_size = glGetUniformLocation(gl.temporal.program, "u_texel_size");
		gl.temporal.uniform_loc.jitter_uv = glGetUniformLocation(gl.temporal.program, "u_jitter_uv");
		gl.temporal.uniform_loc.sin_time = glGetUniformLocation(gl.temporal.program, "u_sin_time");
		gl.temporal.uniform_loc.feedback_min = glGetUniformLocation(gl.temporal.program, "u_feedback_min");
		gl.temporal.uniform_loc.feedback_max = glGetUniformLocation(gl.temporal.program, "u_feedback_max");
		gl.temporal.uniform_loc.motion_scale = glGetUniformLocation(gl.temporal.program, "u_motion_scale");
	}
}

void shutdown() {

}
};

void initialize(int width, int height) {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    if (!gl.vao) glGenVertexArrays(1, &gl.vao);
    if (!gl.vbo) glGenBuffers(1, &gl.vbo);

    glBindVertexArray(gl.vao);
    glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid*)0);
    glBindVertexArray(0);

    if (!gl.v_shader_fs_quad) {
        gl.v_shader_fs_quad = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(gl.v_shader_fs_quad, 1, &v_shader_src_fs_quad, 0);
        glCompileShader(gl.v_shader_fs_quad);
        if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, gl.v_shader_fs_quad)) {
            LOG_ERROR("Error while compiling postprocessing fs-quad vertex shader:\n%s", buffer);
        }
    }

    // LINEARIZE DEPTH
    if (!gl.linear_depth.program_persp) setup_program(&gl.linear_depth.program_persp, "linearize depth persp", f_shader_src_linearize_depth, "#version 150 core\n#define PERSPECTIVE 1");
    if (!gl.linear_depth.program_ortho) setup_program(&gl.linear_depth.program_ortho, "linearize depth ortho", f_shader_src_linearize_depth, "#version 150 core\n#define PERSPECTIVE 0");

    if (!gl.linear_depth.texture) glGenTextures(1, &gl.linear_depth.texture);
    glBindTexture(GL_TEXTURE_2D, gl.linear_depth.texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.linear_depth.fbo) {
        glGenFramebuffers(1, &gl.linear_depth.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 0);
        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            LOG_ERROR("Something went wrong");
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    gl.linear_depth.uniform_loc.clip_info = glGetUniformLocation(gl.linear_depth.program_persp, "u_clip_info");
    gl.linear_depth.uniform_loc.tex_depth = glGetUniformLocation(gl.linear_depth.program_persp, "u_tex_depth");

	// COLOR
	if (!gl.targets.tex_color[0]) glGenTextures(2, gl.targets.tex_color);
	glBindTexture(GL_TEXTURE_2D, gl.targets.tex_color[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, gl.targets.tex_color[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);

	if (!gl.targets.tex_temporal_buffer[0]) glGenTextures(2, gl.targets.tex_temporal_buffer);
	glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);

	if (!gl.targets.fbo) {
		glGenFramebuffers(1, &gl.targets.fbo);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.targets.tex_color[0], 0);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.targets.tex_color[1], 0);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[0], 0);
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[1], 0);

		GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
		if (status != GL_FRAMEBUFFER_COMPLETE) {
			LOG_ERROR("Something went wrong");
		}
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	}

    // HALF RES
	{
        String src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/dof_half_res_prepass.frag");
        defer { FREE(src); };
        setup_program(&gl.half_res.program, "dof pre-pass", src);
        if (gl.half_res.program) {
            gl.half_res.uniform_loc.tex_depth = glGetUniformLocation(gl.half_res.program, "u_tex_depth");
            gl.half_res.uniform_loc.tex_color = glGetUniformLocation(gl.half_res.program, "u_tex_color");
            gl.half_res.uniform_loc.focus_point = glGetUniformLocation(gl.half_res.program, "u_focus_point");
            gl.half_res.uniform_loc.focus_scale = glGetUniformLocation(gl.half_res.program, "u_focus_scale");
        }
	}

    if (!gl.half_res.tex.color_coc) {
        glGenTextures(1, &gl.half_res.tex.color_coc);
    }
    glBindTexture(GL_TEXTURE_2D, gl.half_res.tex.color_coc);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width / 2, height / 2, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.half_res.fbo) {
        glGenFramebuffers(1, &gl.half_res.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.half_res.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.half_res.tex.color_coc, 0);
        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            LOG_ERROR("Something went wrong");
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    // DOF
    {
        String src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/dof.frag");
        defer { FREE(src); };
        if (setup_program(&gl.bokeh_dof.program, "bokeh dof", src)) {
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "uHalfRes");
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "uColor");
            gl.bokeh_dof.uniform_loc.tex_depth = glGetUniformLocation(gl.bokeh_dof.program, "uDepth");
            gl.bokeh_dof.uniform_loc.pixel_size = glGetUniformLocation(gl.bokeh_dof.program, "uPixelSize");
            gl.bokeh_dof.uniform_loc.focus_point = glGetUniformLocation(gl.bokeh_dof.program, "uFocusPoint");
            gl.bokeh_dof.uniform_loc.focus_scale = glGetUniformLocation(gl.bokeh_dof.program, "uFocusScale");
        }
    }

	// Velocity
	{
		if (!gl.velocity.tex_tilemax) {
			glGenTextures(1, &gl.velocity.tex_tilemax);
		}

		if (!gl.velocity.tex_neighbormax) {
			glGenTextures(1, &gl.velocity.tex_neighbormax);
		}

		gl.velocity.tex_width = width / VEL_TILE_SIZE;
		gl.velocity.tex_height = height / VEL_TILE_SIZE;

		glBindTexture(GL_TEXTURE_2D, gl.velocity.tex_tilemax);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, gl.velocity.tex_width, gl.velocity.tex_height, 0, GL_RG, GL_FLOAT, nullptr);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_2D, 0);

		glBindTexture(GL_TEXTURE_2D, gl.velocity.tex_neighbormax);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, gl.velocity.tex_width, gl.velocity.tex_height, 0, GL_RG, GL_FLOAT, nullptr);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_2D, 0);

		if (!gl.velocity.fbo) {
			glGenFramebuffers(1, &gl.velocity.fbo);
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.velocity.fbo);
			glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.velocity.tex_tilemax, 0);
			glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.velocity.tex_neighbormax, 0);
			GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
			if (status != GL_FRAMEBUFFER_COMPLETE) {
				LOG_ERROR("Something went wrong");
			}
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		}


	}

    gl.tex_width = width;
    gl.tex_height = height;

    ssao::initialize(width, height);
    deferred::initialize();
    highlight::initialize();
    tonemapping::initialize();
	velocity::initialize();
	temporal::initialize();
    blit::initialize();
}

void shutdown() {
    ssao::shutdown();
    deferred::shutdown();
    highlight::shutdown();
    tonemapping::shutdown();
	velocity::shutdown();
	temporal::shutdown();
    blit::shutdown();

    if (gl.vao) glDeleteVertexArrays(1, &gl.vao);
    if (gl.vbo) glDeleteBuffers(1, &gl.vbo);
    if (gl.v_shader_fs_quad) glDeleteShader(gl.v_shader_fs_quad);
}

void compute_linear_depth(GLuint depth_tex, float near_plane, float far_plane, bool orthographic = false) {
    const vec4 clip_info(near_plane * far_plane, near_plane - far_plane, far_plane, 0);

	glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);

    if (orthographic)
        glUseProgram(gl.linear_depth.program_ortho);
    else
        glUseProgram(gl.linear_depth.program_persp);
    glUniform1i(gl.linear_depth.uniform_loc.tex_depth, 0);
    glUniform4fv(gl.linear_depth.uniform_loc.clip_info, 1, &clip_info[0]);

    // ASSUME THAT THE APPROPRIATE FS_QUAD VAO IS BOUND
	glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
	glBindVertexArray(0);
}

void apply_ssao(GLuint linear_depth_tex, GLuint normal_tex, const mat4& proj_matrix, float intensity, float radius, float bias) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(normal_tex));

	const float sharpness = ssao::compute_sharpness(radius);
	const vec2 inv_res = vec2(1.f / ssao::tex_width, 1.f / ssao::tex_height);

	GLint last_fbo;
	GLint last_viewport[4];
	glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
	glGetIntegerv(GL_VIEWPORT, last_viewport);

	glBindVertexArray(gl.vao);

    ssao::setup_ubo_hbao_data(ssao::ubo_hbao_data, ssao::tex_width, ssao::tex_height, proj_matrix, intensity, radius, bias);

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ssao::fbo_hbao);
	glViewport(0, 0, ssao::tex_width, ssao::tex_height);

    // RENDER HBAO
    PUSH_GPU_SECTION("HBAO")
    glUseProgram(ssao::prog_hbao);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, ssao::tex_random);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ssao::ubo_hbao_data);
    glUniformBlockBinding(ssao::prog_hbao, ssao::uniform_block_index_hbao_control_buffer, 0);
    glUniform1i(ssao::uniform_loc_hbao_tex_linear_depth, 0);
    glUniform1i(ssao::uniform_loc_hbao_tex_normal, 1);
    glUniform1i(ssao::uniform_loc_hbao_tex_random, 2);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    POP_GPU_SECTION()

    // BLUR FIRST
    PUSH_GPU_SECTION("BLUR 1st")
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ssao::fbo_blur);

    glUseProgram(ssao::prog_blur_first);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, ssao::tex_hbao);

    glUniform1i(ssao::uniform_loc_blur_tex_linear_depth, 0);
	glUniform1i(ssao::uniform_loc_blur_tex_ao, 1);
    glUniform1f(ssao::uniform_loc_blur_sharpness, sharpness);
    glUniform2f(ssao::uniform_loc_blur_inv_res_dir, inv_res.x, 0);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    POP_GPU_SECTION()

    // BLUR SECOND
    PUSH_GPU_SECTION("BLUR 2nd")
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
	glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glUseProgram(ssao::prog_blur_second);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, ssao::tex_blur);

	glUniform1i(ssao::uniform_loc_blur_tex_linear_depth, 0);
	glUniform1i(ssao::uniform_loc_blur_tex_ao, 1);
    glUniform1f(ssao::uniform_loc_blur_sharpness, sharpness);
    glUniform2f(ssao::uniform_loc_blur_inv_res_dir, 0, inv_res.y);

	glEnable(GL_BLEND);
	glBlendFunc(GL_ZERO, GL_SRC_COLOR);
	glColorMask(1, 1, 1, 0);

    glDrawArrays(GL_TRIANGLES, 0, 3);

	glDisable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);
	glColorMask(1, 1, 1, 1);

    glBindVertexArray(0);
    POP_GPU_SECTION()
}

void shade_deferred(GLuint depth_tex, GLuint color_tex, GLuint normal_tex, const mat4& inv_proj_matrix) {
    ASSERT(glIsTexture(depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(normal_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, normal_tex);

    glUseProgram(deferred::deferred.program);
    glUniform1i(deferred::deferred.uniform_loc.texture_depth, 0);
    glUniform1i(deferred::deferred.uniform_loc.texture_color, 1);
    glUniform1i(deferred::deferred.uniform_loc.texture_normal, 2);
    glUniformMatrix4fv(deferred::deferred.uniform_loc.inv_proj_mat, 1, GL_FALSE, &inv_proj_matrix[0][0]);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void highlight_selection(GLuint atom_idx_tex, GLuint selection_buffer) {
    ASSERT(glIsTexture(atom_idx_tex));
    ASSERT(glIsBuffer(selection_buffer));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, atom_idx_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, highlight::highlight.selection_texture);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8, selection_buffer);

    glUseProgram(highlight::highlight.program);
    glUniform1i(highlight::highlight.uniform_loc.texture_atom_idx, 0);
    glUniform1i(highlight::highlight.uniform_loc.buffer_selection, 1);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void half_res_color_coc(GLuint linear_depth_tex, GLuint color_tex, float focus_point, float focus_scale) {
    PUSH_GPU_SECTION("DOF Prepass");
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    glViewport(0, 0, gl.tex_width / 2, gl.tex_height / 2);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.half_res.fbo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(gl.half_res.program);

    glUniform1i(gl.half_res.uniform_loc.tex_depth, 0);
    glUniform1i(gl.half_res.uniform_loc.tex_color, 1);
    glUniform1f(gl.half_res.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.half_res.uniform_loc.focus_scale, focus_scale);

    // ASSUME THAT THE APPROPRIATE FS_QUAD VAO IS BOUND
    glDrawArrays(GL_TRIANGLES, 0, 3);

    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    POP_GPU_SECTION();
}

void apply_dof(GLuint linear_depth_tex, GLuint color_tex, const mat4& proj_matrix, float focus_point, float focus_scale) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));

    const float n = proj_matrix[3][2] / (proj_matrix[2][2] - 1.f);
    const float f = (proj_matrix[2][2] - 1.f) * n / (proj_matrix[2][2] + 1.f);
    const bool ortho = is_orthographic_proj_matrix(proj_matrix);
    const vec2 pixel_size = vec2(1.f / gl.tex_width, 1.f / gl.tex_height);

    GLint last_fbo;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);

    glBindVertexArray(gl.vao);

    half_res_color_coc(linear_depth_tex, color_tex, focus_point, focus_scale);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.half_res.tex.color_coc);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);

    glUseProgram(gl.bokeh_dof.program);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_half_res, 0);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_depth, 1);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_color, 2);
    glUniform2f(gl.bokeh_dof.uniform_loc.pixel_size, pixel_size.x, pixel_size.y);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_scale, focus_scale);

    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
    glActiveTexture(GL_TEXTURE0);
}

void apply_tonemapping(GLuint color_tex, Tonemapping tonemapping, float exposure, float gamma) {
    ASSERT(glIsTexture(color_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    switch (tonemapping) {
        case Tonemapping_ExposureGamma:
            glUseProgram(tonemapping::exposure_gamma.program);
            glUniform1i(tonemapping::exposure_gamma.uniform_loc.texture, 0);
            glUniform1f(tonemapping::exposure_gamma.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::exposure_gamma.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_Filmic:
            glUseProgram(tonemapping::filmic.program);
            glUniform1i(tonemapping::filmic.uniform_loc.texture, 0);
            glUniform1f(tonemapping::filmic.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::filmic.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_Passthrough:
        default:
            glUseProgram(tonemapping::passthrough.program);
            glUniform1i(tonemapping::passthrough.uniform_loc.texture, 0);
            break;
    }

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_velocity(const ViewParam& view_param) {
	mat4 curr_clip_to_prev_clip_mat = view_param.previous.matrix.view_proj * view_param.matrix.inverse.view_proj;

	glUseProgram(velocity::blit_velocity.program);
	glUniformMatrix4fv(velocity::blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat, 1, GL_FALSE, &curr_clip_to_prev_clip_mat[0][0]);
	glBindVertexArray(gl.vao);
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glBindVertexArray(0);
	glUseProgram(0);
}

void blit_tilemax(GLuint velocity_tex, int tex_width, int tex_height) {
	ASSERT(glIsTexture(velocity_tex));
	const vec2 texel_size = { 1.f / tex_width, 1.f / tex_height };

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, velocity_tex);

	glUseProgram(velocity::blit_tilemax.program);
	glUniform1i(velocity::blit_tilemax.uniform_loc.tex_vel, 0);
	glUniform2fv(velocity::blit_tilemax.uniform_loc.tex_vel_texel_size, 1, &texel_size[0]);
	glBindVertexArray(gl.vao);
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glBindVertexArray(0);
	glUseProgram(0);
}

void blit_neighbormax(GLuint velocity_tex, int tex_width, int tex_height) {
	ASSERT(glIsTexture(velocity_tex));
	const vec2 texel_size = { 1.f / tex_width, 1.f / tex_height };

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, velocity_tex);

	glUseProgram(velocity::blit_neighbormax.program);
	glUniform1i(velocity::blit_neighbormax.uniform_loc.tex_vel, 0);
	glUniform2fv(velocity::blit_neighbormax.uniform_loc.tex_vel_texel_size, 1, &texel_size[0]);
	glBindVertexArray(gl.vao);
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glBindVertexArray(0);
	glUseProgram(0);
}

void apply_temporal_aa(GLuint linear_depth_tex, GLuint color_tex, GLuint velocity_tex, GLuint velocity_neighbormax_tex, const vec2& curr_jitter, const vec2& prev_jitter, float feedback_min, float feedback_max, float motion_scale) {
	ASSERT(glIsTexture(linear_depth_tex));
	ASSERT(glIsTexture(color_tex));
	ASSERT(glIsTexture(velocity_tex));
	ASSERT(glIsTexture(velocity_neighbormax_tex));

	static int target = 0;
	target = (target + 1) % 2;

	static float time = 0.f;
	time += 0.016f;

	if (time > 100.f) {
		time -= 100.f;
	}

	int dst_buf = target;
	int src_buf = (target + 1) % 2;

	const vec2 res = { gl.tex_width, gl.tex_height };
	const vec4 texel_size = vec4(1.f / res, res);
	const vec4 jitter_uv = vec4(curr_jitter / res, prev_jitter / res);
	const float sin_time = time;

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, color_tex);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[src_buf]);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, velocity_tex);

	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, velocity_neighbormax_tex);

	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

	GLenum draw_buffers[2];
	draw_buffers[0] = GL_COLOR_ATTACHMENT2 + dst_buf; // tex_temporal_buffer[0 or 1]
	draw_buffers[1] = GL_COLOR_ATTACHMENT1; // tex_final

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
	glViewport(0, 0, gl.tex_width, gl.tex_height);
	glDrawBuffers(2, draw_buffers);

	glUseProgram(gl.temporal.program);

	// Uniforms
	glUniform1i(gl.temporal.uniform_loc.tex_main, 0);
	glUniform1i(gl.temporal.uniform_loc.tex_prev, 1);
	glUniform1i(gl.temporal.uniform_loc.tex_vel, 2);
	glUniform1i(gl.temporal.uniform_loc.tex_vel_neighbormax, 3);
	glUniform1i(gl.temporal.uniform_loc.tex_linear_depth, 4);

	glUniform4fv(gl.temporal.uniform_loc.texel_size, 1, &texel_size[0]);
	glUniform4fv(gl.temporal.uniform_loc.jitter_uv, 1, &jitter_uv[0]);

	glUniform1f(gl.temporal.uniform_loc.sin_time, sin_time);
	glUniform1f(gl.temporal.uniform_loc.feedback_min, feedback_min);
	glUniform1f(gl.temporal.uniform_loc.feedback_max, feedback_max);
	glUniform1f(gl.temporal.uniform_loc.motion_scale, motion_scale);

	glBindVertexArray(gl.vao);
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glBindVertexArray(0);
	glUseProgram(0);
}

void blit_texture(GLuint tex) {
    ASSERT(glIsTexture(tex));
    glUseProgram(blit::program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUniform1i(blit::uniform_loc_texture, 0);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void apply_postprocessing(const PostProcessingDesc& desc, const ViewParam& view_param, GLuint depth_tex, GLuint color_tex, GLuint normal_tex, GLuint velocity_tex) {
	ASSERT(glIsTexture(depth_tex));
	ASSERT(glIsTexture(color_tex));
	ASSERT(glIsTexture(normal_tex));
	ASSERT(glIsTexture(velocity_tex));

	const auto near_dist = view_param.matrix.proj[3][2] / (view_param.matrix.proj[2][2] - 1.f);
	const auto far_dist = (view_param.matrix.proj[2][2] - 1.f) * near_dist / (view_param.matrix.proj[2][2] + 1.f);
	const auto ortho = is_orthographic_proj_matrix(view_param.matrix.proj);

	GLint last_fbo;
	GLint last_viewport[4];
	GLint last_draw_buffer;
	glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
	glGetIntegerv(GL_VIEWPORT, last_viewport);
	glGetIntegerv(GL_DRAW_BUFFER, &last_draw_buffer);

	glViewport(0, 0, gl.tex_width, gl.tex_height);
	glBindVertexArray(gl.vao);

	PUSH_GPU_SECTION("Linearize Depth") {
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo);
		compute_linear_depth(depth_tex, near_dist, far_dist, ortho);
	}
	POP_GPU_SECTION()

	if (desc.temporal_reprojection.enabled) {
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.velocity.fbo);
		glViewport(0, 0, gl.velocity.tex_width, gl.velocity.tex_height);

		PUSH_GPU_SECTION("Velocity: Tilemax") {
			glDrawBuffer(GL_COLOR_ATTACHMENT0);
			blit_tilemax(velocity_tex, gl.tex_width, gl.tex_height);
		}
		POP_GPU_SECTION()

		PUSH_GPU_SECTION("Velocity: Neighbormax") {
			glDrawBuffer(GL_COLOR_ATTACHMENT1);
			blit_neighbormax(gl.velocity.tex_tilemax, gl.velocity.tex_width, gl.velocity.tex_height);
		}
		POP_GPU_SECTION()
	}

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
	glDrawBuffer(GL_COLOR_ATTACHMENT0);
	glViewport(0, 0, gl.tex_width, gl.tex_height);

	if (desc.tonemapping.enabled) {
		PUSH_GPU_SECTION("Tonemapping")
		apply_tonemapping(color_tex, desc.tonemapping.mode, desc.tonemapping.exposure, desc.tonemapping.gamma);
		POP_GPU_SECTION()
	}

	if (desc.ambient_occlusion.enabled) {
		PUSH_GPU_SECTION("SSAO")
		apply_ssao(gl.linear_depth.texture, normal_tex, view_param.matrix.proj, desc.ambient_occlusion.intensity, desc.ambient_occlusion.radius, desc.ambient_occlusion.bias);
		POP_GPU_SECTION()
	}
	
	if (desc.depth_of_field.enabled) {
		glDrawBuffer(GL_COLOR_ATTACHMENT1);
		PUSH_GPU_SECTION("DOF")
		apply_dof(gl.linear_depth.texture, gl.targets.tex_color[0], view_param.matrix.proj, desc.depth_of_field.focus_depth, desc.depth_of_field.focus_scale);
		POP_GPU_SECTION()
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		blit_texture(gl.targets.tex_color[1]);
	}

	if (desc.temporal_reprojection.enabled) {
		PUSH_GPU_SECTION("Temporal AA + Motion-Blur")
		apply_temporal_aa(gl.linear_depth.texture, gl.targets.tex_color[0], velocity_tex, gl.velocity.tex_neighbormax, view_param.jitter, view_param.previous.jitter, desc.temporal_reprojection.feedback_min, desc.temporal_reprojection.feedback_max, desc.temporal_reprojection.motion_blur.motion_scale);
		POP_GPU_SECTION()
	}

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
	glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
	glDrawBuffer(last_draw_buffer);
	blit_texture(gl.targets.tex_color[1]);
}

}  // namespace postprocessing
