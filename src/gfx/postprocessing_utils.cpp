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
        GLuint tex_color[2] = {0, 0};
        GLuint tex_temporal_buffer[2] = {0, 0};  // These are dedicated and cannot be use as intermediate buffers by other shaders
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
            GLuint program_persp = 0;
            GLuint program_ortho = 0;

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
            GLuint program = 0;
            struct {
                GLint sharpness = -1;
                GLint inv_res_dir = -1;
                GLint tex_ao = -1;
                GLint tex_linear_depth = -1;
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
            GLint time = -1;
        } uniform_loc;

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
        struct {
            GLuint program = 0;
            struct {
                GLint tex_linear_depth = -1;
                GLint tex_main = -1;
                GLint tex_prev = -1;
                GLint tex_vel = -1;
                GLint tex_vel_neighbormax = -1;
                GLint texel_size = -1;
                GLint time = -1;
                GLint feedback_min = -1;
                GLint feedback_max = -1;
                GLint motion_scale = -1;
                GLint jitter_uv = -1;
            } uniform_loc;
        } with_motion_blur;
        struct {
            GLuint program = 0;
            struct {
                GLint tex_linear_depth = -1;
                GLint tex_main = -1;
                GLint tex_prev = -1;
                GLint tex_vel = -1;
                GLint texel_size = -1;
                GLint time = -1;
                GLint feedback_min = -1;
                GLint feedback_max = -1;
                GLint motion_scale = -1;
                GLint jitter_uv = -1;
            } uniform_loc;
        } no_motion_blur;
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

static const char* f_shader_src_mip_map_min_depth = R"(
uniform sampler2D u_tex_depth;

out vec4 out_frag;

void main() {
	float d00 = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy) + ivec2(0,0), 0).x;
	float d01 = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy) + ivec2(0,1), 0).x;
	float d10 = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy) + ivec2(1,0), 0).x;
	float d11 = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy) + ivec2(1,1), 0).x;

	float dmin0 = min(d00, d01);
	float dmin1 = min(d10, d11);

	out_frag = vec4(min(dmin0, dmin1));
}
)";

static bool setup_program(GLuint* program, CString name, CString f_shader_src, CString defines = {}) {
    ASSERT(program);
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    auto f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    if (defines) {
        StringBuffer<64> version_str;
        if (compare_n(f_shader_src, "#version", 8)) {
            version_str = extract_line(f_shader_src);
        }
        const char* sources[5] = {version_str, "\n", defines, "\n", f_shader_src};
        glShaderSource(f_shader, 5, sources, 0);
    } else {
        const char* sources[1] = {f_shader_src};
        glShaderSource(f_shader, 1, sources, 0);
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

struct HBAOData {
    float radius_to_screen;
    float neg_inv_r2;
    float n_dot_v_bias;
    float time;

    float ao_multiplier;
    float pow_exponent;
    vec2 inv_full_res;

    vec4 proj_info;

    vec4 sample_pattern[32];
};

void setup_ubo_hbao_data(GLuint ubo, int width, int height, const mat4& proj_mat, float intensity, float radius, float bias, float time) {
    ASSERT(ubo);

    // From intel ASSAO
    static constexpr float SAMPLE_PATTERN[] = {
        0.78488064,  0.56661671,  1.500000, -0.126083, 0.26022232,  -0.29575172, 1.500000, -1.064030, 0.10459357,  0.08372527,  1.110000, -2.730563, -0.68286800, 0.04963045,  1.090000, -0.498827,
        -0.13570161, -0.64190155, 1.250000, -0.532765, -0.26193795, -0.08205118, 0.670000, -1.783245, -0.61177456, 0.66664219,  0.710000, -0.044234, 0.43675563,  0.25119025,  0.610000, -1.167283,
        0.07884444,  0.86618668,  0.640000, -0.459002, -0.12790935, -0.29869005, 0.600000, -1.729424, -0.04031125, 0.02413622,  0.600000, -4.792042, 0.16201244,  -0.52851415, 0.790000, -1.067055,
        -0.70991218, 0.47301072,  0.640000, -0.335236, 0.03277707,  -0.22349690, 0.600000, -1.982384, 0.68921727,  0.36800742,  0.630000, -0.266718, 0.29251814,  0.37775412,  0.610000, -1.422520,
        -0.12224089, 0.96582592,  0.600000, -0.426142, 0.11071457,  -0.16131058, 0.600000, -2.165947, 0.46562141,  -0.59747696, 0.600000, -0.189760, -0.51548797, 0.11804193,  0.600000, -1.246800,
        0.89141309,  -0.42090443, 0.600000, 0.028192,  -0.32402530, -0.01591529, 0.600000, -1.543018, 0.60771245,  0.41635221,  0.600000, -0.605411, 0.02379565,  -0.08239821, 0.600000, -3.809046,
        0.48951152,  -0.23657045, 0.600000, -1.189011, -0.17611565, -0.81696892, 0.600000, -0.513724, -0.33930185, -0.20732205, 0.600000, -1.698047, -0.91974425, 0.05403209,  0.600000, 0.062246,
        -0.15064627, -0.14949332, 0.600000, -1.896062, 0.53180975,  -0.35210401, 0.600000, -0.758838, 0.41487166,  0.81442589,  0.600000, -0.505648, -0.24106961, -0.32721516, 0.600000, -1.665244};
    constexpr float METERS_TO_VIEWSPACE = 1.f;

    vec4 proj_info;
    float proj_scl;

    const float* proj_data = &proj_mat[0][0];
    if (!is_orthographic_proj_matrix(proj_mat)) {
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
    data.neg_inv_r2 = -1.f / (r * r);
    data.n_dot_v_bias = math::clamp(bias, 0.f, 1.f - math::EPSILON);
    data.time = time;
    data.ao_multiplier = 1.f / (1.f - data.n_dot_v_bias);
    data.pow_exponent = math::max(intensity, 0.f);
    data.inv_full_res = vec2(1.f / float(width), 1.f / float(height));
    data.proj_info = proj_info;
    memcpy(&data.sample_pattern, SAMPLE_PATTERN, sizeof(SAMPLE_PATTERN));

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(HBAOData), &data, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void initialize_rnd_tex(GLuint rnd_tex) {
    constexpr int buffer_size = AO_RANDOM_TEX_SIZE * AO_RANDOM_TEX_SIZE;
    signed short buffer[buffer_size * 4];

    vec2 rnd_vals[buffer_size];
    math::generate_halton_sequence(rnd_vals, buffer_size, 2, 3);

    for (int i = 0; i < buffer_size; i++) {
#define SCALE ((1 << 15))
        float rand1 = rnd_vals[i].x;
        float rand2 = rnd_vals[i].y;
        float angle = 2.f * math::PI * rand1;

        buffer[i * 4 + 0] = (signed short)(SCALE * math::cos(angle));
        buffer[i * 4 + 1] = (signed short)(SCALE * math::sin(angle));
        buffer[i * 4 + 2] = (signed short)(SCALE * rand2);
        buffer[i * 4 + 3] = (signed short)(SCALE * 0);
#undef SCALE
    }

    glBindTexture(GL_TEXTURE_2D, rnd_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16_SNORM, AO_RANDOM_TEX_SIZE, AO_RANDOM_TEX_SIZE, 0, GL_RGBA, GL_SHORT, buffer);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

float compute_sharpness(float radius) { return 30.f / math::sqrt(radius); }

void initialize(int width, int height) {

    String f_shader_src_ssao = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/ssao/ssao.frag");
    String f_shader_src_blur = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/ssao/blur.frag");
    defer {
        free_string(&f_shader_src_ssao);
        free_string(&f_shader_src_blur);
    };

    setup_program(&gl.ssao.hbao.program_persp, "ssao perspective", f_shader_src_ssao, "#define AO_PERSPECTIVE 1");
    setup_program(&gl.ssao.hbao.program_ortho, "ssao orthographic", f_shader_src_ssao, "#define AO_PERSPECTIVE 0");
    setup_program(&gl.ssao.blur.program, "ssao blur", f_shader_src_blur);

    gl.ssao.hbao.uniform_loc.control_buffer = glGetUniformBlockIndex(gl.ssao.hbao.program_persp, "u_control_buffer");
    gl.ssao.hbao.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.ssao.hbao.program_persp, "u_tex_linear_depth");
    gl.ssao.hbao.uniform_loc.tex_normal = glGetUniformLocation(gl.ssao.hbao.program_persp, "u_tex_normal");
    gl.ssao.hbao.uniform_loc.tex_random = glGetUniformLocation(gl.ssao.hbao.program_persp, "u_tex_random");

    gl.ssao.blur.uniform_loc.sharpness = glGetUniformLocation(gl.ssao.blur.program, "u_sharpness");
    gl.ssao.blur.uniform_loc.inv_res_dir = glGetUniformLocation(gl.ssao.blur.program, "u_inv_res_dir");
    gl.ssao.blur.uniform_loc.tex_ao = glGetUniformLocation(gl.ssao.blur.program, "u_tex_ao");
    gl.ssao.blur.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.ssao.blur.program, "u_tex_linear_depth");

    if (!gl.ssao.hbao.fbo) glGenFramebuffers(1, &gl.ssao.hbao.fbo);
    if (!gl.ssao.blur.fbo) glGenFramebuffers(1, &gl.ssao.blur.fbo);

    if (!gl.ssao.tex_random) glGenTextures(1, &gl.ssao.tex_random);
    if (!gl.ssao.hbao.texture) glGenTextures(1, &gl.ssao.hbao.texture);
    if (!gl.ssao.blur.texture) glGenTextures(1, &gl.ssao.blur.texture);

    if (!gl.ssao.ubo_hbao_data) glGenBuffers(1, &gl.ssao.ubo_hbao_data);

    initialize_rnd_tex(gl.ssao.tex_random);

    glBindTexture(GL_TEXTURE_2D, gl.ssao.hbao.texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gl.ssao.blur.texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, gl.ssao.hbao.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.ssao.hbao.texture, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, gl.ssao.blur.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.ssao.blur.texture, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    glBindBuffer(GL_UNIFORM_BUFFER, gl.ssao.ubo_hbao_data);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(HBAOData), nullptr, GL_DYNAMIC_DRAW);
}

void shutdown() {
    if (gl.ssao.hbao.fbo) glDeleteFramebuffers(1, &gl.ssao.hbao.fbo);
    if (gl.ssao.blur.fbo) glDeleteFramebuffers(1, &gl.ssao.blur.fbo);

    if (gl.ssao.tex_random) glDeleteTextures(1, &gl.ssao.tex_random);
    if (gl.ssao.hbao.texture) glDeleteTextures(1, &gl.ssao.hbao.texture);
    if (gl.ssao.blur.texture) glDeleteTextures(1, &gl.ssao.blur.texture);

    if (gl.ssao.ubo_hbao_data) glDeleteBuffers(1, &gl.ssao.ubo_hbao_data);

    if (gl.ssao.hbao.program_persp) glDeleteProgram(gl.ssao.hbao.program_persp);
    if (gl.ssao.hbao.program_ortho) glDeleteProgram(gl.ssao.hbao.program_ortho);
    if (gl.ssao.blur.program) glDeleteProgram(gl.ssao.blur.program);
}

}  // namespace ssao

namespace highlight {

static struct {
    GLuint program = 0;
    GLuint selection_texture = 0;
    struct {
        GLint texture_atom_idx = -1;
        GLint buffer_selection = -1;
        GLint highlight = -1;
        GLint selection = -1;
        GLint outline = -1;
    } uniform_loc;
} highlight;

void initialize() {
    String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/highlight.frag");
    defer { FREE(f_shader_src); };
    setup_program(&highlight.program, "highlight", f_shader_src);
    if (!highlight.selection_texture) glGenTextures(1, &highlight.selection_texture);
    highlight.uniform_loc.texture_atom_idx = glGetUniformLocation(highlight.program, "u_texture_atom_idx");
    highlight.uniform_loc.buffer_selection = glGetUniformLocation(highlight.program, "u_buffer_selection");
    highlight.uniform_loc.highlight = glGetUniformLocation(highlight.program, "u_highlight");
    highlight.uniform_loc.selection = glGetUniformLocation(highlight.program, "u_selection");
    highlight.uniform_loc.outline = glGetUniformLocation(highlight.program, "u_outline");
}

void shutdown() {
    if (highlight.program) glDeleteProgram(highlight.program);
}
}  // namespace highlight

namespace hsv {

static struct {
    GLuint program = 0;
    struct {
        GLint texture_color = -1;
        GLint hsv_scale = -1;
    } uniform_loc;
} gl;

void initialize() {
    String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/scale_hsv.frag");
    defer { FREE(f_shader_src); };
    setup_program(&gl.program, "scale hsv", f_shader_src);
    gl.uniform_loc.texture_color = glGetUniformLocation(gl.program, "u_texture_atom_color");
    gl.uniform_loc.hsv_scale = glGetUniformLocation(gl.program, "u_hsv_scale");
}

void shutdown() {
    if (gl.program) glDeleteProgram(gl.program);
}
}  // namespace hsv

namespace deferred {

static struct {
    GLuint program = 0;
    struct {
        GLint texture_depth = -1;
        GLint texture_color = -1;
        GLint texture_normal = -1;
        GLint inv_proj_mat = -1;
        GLint time = -1;
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
    deferred.uniform_loc.time = glGetUniformLocation(deferred.program, "u_time");
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

namespace dof {
void initialize(int32 width, int32 height) {
    {
        String src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/dof/dof_half_res_prepass.frag");
        defer { FREE(src); };
        setup_program(&gl.bokeh_dof.half_res.program, "dof pre-pass", src);
        if (gl.bokeh_dof.half_res.program) {
            gl.bokeh_dof.half_res.uniform_loc.tex_depth = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_tex_depth");
            gl.bokeh_dof.half_res.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_tex_color");
            gl.bokeh_dof.half_res.uniform_loc.focus_point = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_focus_point");
            gl.bokeh_dof.half_res.uniform_loc.focus_scale = glGetUniformLocation(gl.bokeh_dof.half_res.program, "u_focus_scale");
        }
    }

    if (!gl.bokeh_dof.half_res.tex.color_coc) {
        glGenTextures(1, &gl.bokeh_dof.half_res.tex.color_coc);
    }
    glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.half_res.tex.color_coc);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width / 2, height / 2, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.bokeh_dof.half_res.fbo) {
        glGenFramebuffers(1, &gl.bokeh_dof.half_res.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.bokeh_dof.half_res.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.half_res.tex.color_coc, 0);
        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            LOG_ERROR("Something went wrong");
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    // DOF
    {
        String src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/dof/dof.frag");
        defer { FREE(src); };
        if (setup_program(&gl.bokeh_dof.program, "bokeh dof", src)) {
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "u_half_res");
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "u_tex_color");
            gl.bokeh_dof.uniform_loc.tex_depth = glGetUniformLocation(gl.bokeh_dof.program, "u_tex_depth");
            gl.bokeh_dof.uniform_loc.pixel_size = glGetUniformLocation(gl.bokeh_dof.program, "u_texel_size");
            gl.bokeh_dof.uniform_loc.focus_point = glGetUniformLocation(gl.bokeh_dof.program, "u_focus_depth");
            gl.bokeh_dof.uniform_loc.focus_scale = glGetUniformLocation(gl.bokeh_dof.program, "u_focus_scale");
            gl.bokeh_dof.uniform_loc.time = glGetUniformLocation(gl.bokeh_dof.program, "u_time");
        }
    }
}

void shutdown() {}
}  // namespace dof

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
		GLint tex_depth = -1;
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

void initialize(int32 width, int32 height) {
    {
        String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/velocity/blit_velocity.frag");
        defer { free_string(&f_shader_src); };
        setup_program(&blit_velocity.program, "screen-space velocity", f_shader_src);
		blit_velocity.uniform_loc.tex_depth = glGetUniformLocation(blit_velocity.program, "u_tex_depth");
        blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat = glGetUniformLocation(blit_velocity.program, "u_curr_clip_to_prev_clip_mat");
    }
    {
        const char* defines = {"#define TILE_SIZE " TOSTRING(VEL_TILE_SIZE)};
        String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/velocity/blit_tilemax.frag");
        defer { free_string(&f_shader_src); };
        setup_program(&blit_tilemax.program, "tilemax", f_shader_src, defines);
        blit_tilemax.uniform_loc.tex_vel = glGetUniformLocation(blit_tilemax.program, "u_tex_vel");
        blit_tilemax.uniform_loc.tex_vel_texel_size = glGetUniformLocation(blit_tilemax.program, "u_tex_vel_texel_size");
    }
    {
        String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/velocity/blit_neighbormax.frag");
        defer { free_string(&f_shader_src); };
        setup_program(&blit_neighbormax.program, "neighbormax", f_shader_src);
        blit_neighbormax.uniform_loc.tex_vel = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel");
        blit_neighbormax.uniform_loc.tex_vel_texel_size = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel_texel_size");
    }

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

void shutdown() {
    if (blit_velocity.program) glDeleteProgram(blit_velocity.program);
    if (gl.velocity.tex_tilemax) glDeleteTextures(1, &gl.velocity.tex_tilemax);
    if (gl.velocity.tex_neighbormax) glDeleteTextures(1, &gl.velocity.tex_neighbormax);
    if (gl.velocity.fbo) glDeleteFramebuffers(1, &gl.velocity.fbo);
}
}  // namespace velocity

namespace temporal {
void initialize() {
    {
        String f_shader_src = allocate_and_read_textfile(MDUTILS_SHADER_DIR "/temporal.frag");
        defer { free_string(&f_shader_src); };
        setup_program(&gl.temporal.with_motion_blur.program, "temporal aa + motion-blur", f_shader_src);
        setup_program(&gl.temporal.no_motion_blur.program, "temporal aa", f_shader_src, "#define USE_MOTION_BLUR 0\n");

        gl.temporal.with_motion_blur.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_linear_depth");
        gl.temporal.with_motion_blur.uniform_loc.tex_main = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_main");
        gl.temporal.with_motion_blur.uniform_loc.tex_prev = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_prev");
        gl.temporal.with_motion_blur.uniform_loc.tex_vel = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_vel");
        gl.temporal.with_motion_blur.uniform_loc.tex_vel_neighbormax = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_vel_neighbormax");
        gl.temporal.with_motion_blur.uniform_loc.texel_size = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_texel_size");
        gl.temporal.with_motion_blur.uniform_loc.jitter_uv = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_jitter_uv");
        gl.temporal.with_motion_blur.uniform_loc.time = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_time");
        gl.temporal.with_motion_blur.uniform_loc.feedback_min = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_feedback_min");
        gl.temporal.with_motion_blur.uniform_loc.feedback_max = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_feedback_max");
        gl.temporal.with_motion_blur.uniform_loc.motion_scale = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_motion_scale");

        gl.temporal.no_motion_blur.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_linear_depth");
        gl.temporal.no_motion_blur.uniform_loc.tex_main = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_main");
        gl.temporal.no_motion_blur.uniform_loc.tex_prev = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_prev");
        gl.temporal.no_motion_blur.uniform_loc.tex_vel = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_vel");
        gl.temporal.no_motion_blur.uniform_loc.texel_size = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_texel_size");
        gl.temporal.no_motion_blur.uniform_loc.jitter_uv = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_jitter_uv");
        gl.temporal.no_motion_blur.uniform_loc.time = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_time");
        gl.temporal.no_motion_blur.uniform_loc.feedback_min = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_feedback_min");
        gl.temporal.no_motion_blur.uniform_loc.feedback_max = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_feedback_max");
        gl.temporal.no_motion_blur.uniform_loc.motion_scale = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_motion_scale");
    }
}

void shutdown() {}
};  // namespace temporal

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
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
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
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_color[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.targets.tex_temporal_buffer[0]) glGenTextures(2, gl.targets.tex_temporal_buffer);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
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

        GLenum buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3};
        glDrawBuffers(4, buffers);
        glClear(GL_COLOR_BUFFER_BIT);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    gl.tex_width = width;
    gl.tex_height = height;

    ssao::initialize(width, height);
    dof::initialize(width, height);
    velocity::initialize(width, height);
    deferred::initialize();
    highlight::initialize();
    hsv::initialize();
    tonemapping::initialize();
    temporal::initialize();
    blit::initialize();
}

void shutdown() {
    ssao::shutdown();
    dof::shutdown();
    velocity::shutdown();
    deferred::shutdown();
    highlight::shutdown();
    hsv::shutdown();
    tonemapping::shutdown();
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

void apply_ssao(GLuint linear_depth_tex, GLuint normal_tex, const mat4& proj_matrix, float intensity, float radius, float bias, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(normal_tex));

    const bool ortho = is_orthographic_proj_matrix(proj_matrix);
    const float sharpness = ssao::compute_sharpness(radius);
    const vec2 inv_res = vec2(1.f / gl.tex_width, 1.f / gl.tex_height);

    GLint last_fbo;
    GLint last_viewport[4];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);

    glBindVertexArray(gl.vao);

    ssao::setup_ubo_hbao_data(gl.ssao.ubo_hbao_data, gl.tex_width, gl.tex_height, proj_matrix, intensity, radius, bias, time);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.ssao.hbao.fbo);
    glViewport(0, 0, gl.tex_width, gl.tex_height);

    // RENDER HBAO
    GLuint program = ortho ? gl.ssao.hbao.program_ortho : gl.ssao.hbao.program_persp;

    PUSH_GPU_SECTION("HBAO")
    glUseProgram(program);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, gl.ssao.tex_random);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ssao.ubo_hbao_data);
    glUniformBlockBinding(program, gl.ssao.hbao.uniform_loc.control_buffer, 0);
    glUniform1i(gl.ssao.hbao.uniform_loc.tex_linear_depth, 0);
    glUniform1i(gl.ssao.hbao.uniform_loc.tex_normal, 1);
    glUniform1i(gl.ssao.hbao.uniform_loc.tex_random, 2);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    POP_GPU_SECTION()

    glUseProgram(gl.ssao.blur.program);

    glUniform1i(gl.ssao.blur.uniform_loc.tex_linear_depth, 0);
    glUniform1i(gl.ssao.blur.uniform_loc.tex_ao, 1);
    glUniform1f(gl.ssao.blur.uniform_loc.sharpness, sharpness);
    glUniform2f(gl.ssao.blur.uniform_loc.inv_res_dir, inv_res.x, 0);

    glActiveTexture(GL_TEXTURE1);

    // BLUR FIRST
    PUSH_GPU_SECTION("BLUR 1st")
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.ssao.blur.fbo);
    glBindTexture(GL_TEXTURE_2D, gl.ssao.hbao.texture);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    POP_GPU_SECTION()

    // BLUR SECOND
    PUSH_GPU_SECTION("BLUR 2nd")
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, gl.ssao.blur.texture);
    glUniform2f(gl.ssao.blur.uniform_loc.inv_res_dir, 0, inv_res.y);

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

void shade_deferred(GLuint depth_tex, GLuint color_tex, GLuint normal_tex, const mat4& inv_proj_matrix, float time) {
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
    glUniform1f(deferred::deferred.uniform_loc.time, time);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void highlight_selection(GLuint atom_idx_tex, GLuint selection_buffer, const vec3& highlight, const vec3& selection, const vec3& outline) {
    ASSERT(glIsTexture(atom_idx_tex));
    ASSERT(glIsBuffer(selection_buffer));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, atom_idx_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, highlight::highlight.selection_texture);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, selection_buffer);

    glUseProgram(highlight::highlight.program);
    glUniform1i(highlight::highlight.uniform_loc.texture_atom_idx, 0);
    glUniform1i(highlight::highlight.uniform_loc.buffer_selection, 1);
    glUniform3fv(highlight::highlight.uniform_loc.highlight, 1, &highlight[0]);
    glUniform3fv(highlight::highlight.uniform_loc.selection, 1, &selection[0]);
    glUniform3fv(highlight::highlight.uniform_loc.outline, 1, &outline[0]);
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

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.bokeh_dof.half_res.fbo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(gl.bokeh_dof.half_res.program);

    glUniform1i(gl.bokeh_dof.half_res.uniform_loc.tex_depth, 0);
    glUniform1i(gl.bokeh_dof.half_res.uniform_loc.tex_color, 1);
    glUniform1f(gl.bokeh_dof.half_res.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.bokeh_dof.half_res.uniform_loc.focus_scale, focus_scale);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    POP_GPU_SECTION();
}

void apply_dof(GLuint linear_depth_tex, GLuint color_tex, float focus_point, float focus_scale, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));

    const vec2 pixel_size = vec2(1.f / gl.tex_width, 1.f / gl.tex_height);

    GLint last_fbo;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);

    half_res_color_coc(linear_depth_tex, color_tex, focus_point, focus_scale);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.half_res.tex.color_coc);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(gl.bokeh_dof.program);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_half_res, 0);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_depth, 1);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_color, 2);
    glUniform2f(gl.bokeh_dof.uniform_loc.pixel_size, pixel_size.x, pixel_size.y);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_scale, focus_scale);
    glUniform1f(gl.bokeh_dof.uniform_loc.time, time);

    glBindVertexArray(gl.vao);
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

void blit_static_velocity(GLuint depth_tex, const ViewParam& view_param) {
    mat4 curr_clip_to_prev_clip_mat = view_param.previous.matrix.view_proj * view_param.matrix.inverse.view_proj;

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, depth_tex);

    glUseProgram(velocity::blit_velocity.program);
	glUniform1i(velocity::blit_velocity.uniform_loc.tex_depth, 0);
    glUniformMatrix4fv(velocity::blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat, 1, GL_FALSE, &curr_clip_to_prev_clip_mat[0][0]);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_tilemax(GLuint velocity_tex, int tex_width, int tex_height) {
    ASSERT(glIsTexture(velocity_tex));
    const vec2 texel_size = {1.f / tex_width, 1.f / tex_height};

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
    const vec2 texel_size = {1.f / tex_width, 1.f / tex_height};

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

void apply_temporal_aa(GLuint linear_depth_tex, GLuint color_tex, GLuint velocity_tex, GLuint velocity_neighbormax_tex, const vec2& curr_jitter, const vec2& prev_jitter, float feedback_min,
                       float feedback_max, float motion_scale, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(velocity_neighbormax_tex));

    static int target = 0;
    target = (target + 1) % 2;

    const int dst_buf = target;
    const int src_buf = (target + 1) % 2;

    const vec2 res = {gl.tex_width, gl.tex_height};
    const vec4 texel_size = vec4(1.f / res, res);
    const vec4 jitter_uv = vec4(curr_jitter / res, prev_jitter / res);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, gl.targets.tex_temporal_buffer[src_buf]);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_2D, velocity_neighbormax_tex);

    GLint bound_buffer;
    glGetIntegerv(GL_DRAW_BUFFER, &bound_buffer);

    GLenum draw_buffers[2];
    draw_buffers[0] = GL_COLOR_ATTACHMENT2 + dst_buf;  // tex_temporal_buffer[0 or 1]
    draw_buffers[1] = bound_buffer;                    // assume that this is part of the same fbo

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
    glViewport(0, 0, gl.tex_width, gl.tex_height);
    glDrawBuffers(2, draw_buffers);

    if (motion_scale != 0.f) {
        glUseProgram(gl.temporal.with_motion_blur.program);

        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_linear_depth, 0);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_main, 1);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_prev, 2);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_vel, 3);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_vel_neighbormax, 4);

        glUniform4fv(gl.temporal.with_motion_blur.uniform_loc.texel_size, 1, &texel_size[0]);
        glUniform4fv(gl.temporal.with_motion_blur.uniform_loc.jitter_uv, 1, &jitter_uv[0]);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.time, time);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.feedback_min, feedback_min);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.feedback_max, feedback_max);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.motion_scale, motion_scale);
    } else {
        glUseProgram(gl.temporal.no_motion_blur.program);

        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_linear_depth, 0);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_main, 1);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_prev, 2);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_vel, 3);

        glUniform4fv(gl.temporal.no_motion_blur.uniform_loc.texel_size, 1, &texel_size[0]);
        glUniform4fv(gl.temporal.no_motion_blur.uniform_loc.jitter_uv, 1, &jitter_uv[0]);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.time, time);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.feedback_min, feedback_min);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.feedback_max, feedback_max);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.motion_scale, motion_scale);
    }

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void scale_hsv(GLuint color_tex, vec3 hsv_scale) {
    GLint last_fbo;
    GLint last_viewport[4];
    GLint last_draw_buffer;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    glGetIntegerv(GL_DRAW_BUFFER, &last_draw_buffer);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, gl.tex_width, gl.tex_height);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(hsv::gl.program);
    glUniform1i(hsv::gl.uniform_loc.texture_color, 0);
    glUniform3fv(hsv::gl.uniform_loc.hsv_scale, 1, &hsv_scale[0]);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glDrawBuffer(last_draw_buffer);

    blit_texture(gl.targets.tex_color[0]);
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

void shade_and_postprocess(const Descriptor& desc, const ViewParam& view_param) {
    ASSERT(glIsTexture(desc.input_textures.depth));
    ASSERT(glIsTexture(desc.input_textures.color));
    ASSERT(glIsTexture(desc.input_textures.normal));
    ASSERT(glIsTexture(desc.input_textures.velocity));

    // For seeding noise
    static float time = 0.f;
    time = time + 0.016f;
    if (time > 100.f) time -= 100.f;

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
        compute_linear_depth(desc.input_textures.depth, near_dist, far_dist, ortho);
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Generate Linear Depth Mipmaps") {
        glBindTexture(GL_TEXTURE_2D, gl.linear_depth.texture);
        glGenerateMipmap(GL_TEXTURE_2D);
    }
    POP_GPU_SECTION()

    if (desc.temporal_reprojection.enabled) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.velocity.fbo);
        glViewport(0, 0, gl.velocity.tex_width, gl.velocity.tex_height);

        PUSH_GPU_SECTION("Velocity: Tilemax") {
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            blit_tilemax(desc.input_textures.velocity, gl.tex_width, gl.tex_height);
        }
        POP_GPU_SECTION()

        PUSH_GPU_SECTION("Velocity: Neighbormax") {
            glDrawBuffer(GL_COLOR_ATTACHMENT1);
            blit_neighbormax(gl.velocity.tex_tilemax, gl.velocity.tex_width, gl.velocity.tex_height);
        }
        POP_GPU_SECTION()
    }

    GLenum dst_buffer = GL_COLOR_ATTACHMENT1;
    GLuint src_texture = gl.targets.tex_color[0];

    auto swap_target = [&dst_buffer, &src_texture]() {
        dst_buffer = dst_buffer == GL_COLOR_ATTACHMENT0 ? GL_COLOR_ATTACHMENT1 : GL_COLOR_ATTACHMENT0;
        src_texture = src_texture == gl.targets.tex_color[0] ? gl.targets.tex_color[1] : gl.targets.tex_color[0];
    };

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.targets.fbo);
    glDrawBuffer(dst_buffer);
    glViewport(0, 0, gl.tex_width, gl.tex_height);

    PUSH_GPU_SECTION("Clear HDR")
    glClearColor(desc.background.intensity.r, desc.background.intensity.g, desc.background.intensity.b, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Shade")
    shade_deferred(desc.input_textures.depth, desc.input_textures.color, desc.input_textures.normal, view_param.matrix.inverse.proj, time);
    POP_GPU_SECTION()

    if (desc.ambient_occlusion.enabled) {
        PUSH_GPU_SECTION("SSAO")
        apply_ssao(gl.linear_depth.texture, desc.input_textures.normal, view_param.matrix.proj, desc.ambient_occlusion.intensity, desc.ambient_occlusion.radius, desc.ambient_occlusion.bias, time);
        POP_GPU_SECTION()
    }

    if (desc.input_textures.emissive) {
        PUSH_GPU_SECTION("Add Emissive")
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
        blit_texture(desc.input_textures.emissive);
        glDisable(GL_BLEND);
        POP_GPU_SECTION()
    }

    if (desc.depth_of_field.enabled) {
        swap_target();
        glDrawBuffer(dst_buffer);
        PUSH_GPU_SECTION("DOF")
        apply_dof(gl.linear_depth.texture, src_texture, desc.depth_of_field.focus_depth, desc.depth_of_field.focus_scale, time);
        POP_GPU_SECTION()
    }

    PUSH_GPU_SECTION("Tonemapping") {
        swap_target();
        glDrawBuffer(dst_buffer);
        const auto tonemapper = desc.tonemapping.enabled ? desc.tonemapping.mode : Tonemapping_Passthrough;
        apply_tonemapping(src_texture, tonemapper, desc.tonemapping.exposure, desc.tonemapping.gamma);
    }
    POP_GPU_SECTION()

    if (desc.input_textures.post_tonemap) {
        PUSH_GPU_SECTION("Add Post Tonemap")
        glEnable(GL_BLEND);
        glColorMask(1, 1, 1, 1);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        blit_texture(desc.input_textures.post_tonemap);
        glDisable(GL_BLEND);
        POP_GPU_SECTION()
    }

    if (desc.temporal_reprojection.enabled) {
        swap_target();
        glDrawBuffer(dst_buffer);
        const float feedback_min = desc.temporal_reprojection.feedback_min;
        const float feedback_max = desc.temporal_reprojection.feedback_max;
        const float motion_scale = desc.temporal_reprojection.motion_blur.enabled ? desc.temporal_reprojection.motion_blur.motion_scale : 0.f;
        if (motion_scale != 0.f)
            PUSH_GPU_SECTION("Temporal AA + Motion Blur")
        else
            PUSH_GPU_SECTION("Temporal AA")
        apply_temporal_aa(gl.linear_depth.texture, src_texture, desc.input_textures.velocity, gl.velocity.tex_neighbormax, view_param.jitter, view_param.previous.jitter, feedback_min, feedback_max,
                          motion_scale, time);
        POP_GPU_SECTION()
    }

    // Activate backbuffer or whatever was bound before
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glDrawBuffer(last_draw_buffer);

    swap_target();
    blit_texture(src_texture);

    glColorMask(1, 1, 1, 1);
}

}  // namespace postprocessing
