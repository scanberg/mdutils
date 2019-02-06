//#define NOMINMAX

#include "molecule_draw.h"
#include <core/common.h>
#include <core/log.h>
#include <core/string_utils.h>
#include <gfx/gl_utils.h>
#include <mol/molecule_utils.h>

namespace draw {
static GLuint vao = 0;

namespace vdw {
static GLuint program = 0;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_inv_proj_mat = -1;
static GLint uniform_loc_curr_view_to_prev_clip_mat = -1;
static GLint uniform_loc_radius_scale = -1;

static void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/vdw.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/vdw.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/vdw.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_inv_proj_mat = glGetUniformLocation(program, "u_inv_proj_mat");
    uniform_loc_radius_scale = glGetUniformLocation(program, "u_radius_scale");
	uniform_loc_curr_view_to_prev_clip_mat = glGetUniformLocation(program, "u_curr_view_to_prev_clip_mat");
}

static void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace vdw

namespace licorice {
static GLuint program = 0;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_radius = -1;

static void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/licorice.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/licorice.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/licorice.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_radius = glGetUniformLocation(program, "u_radius");
}

static void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace licorice

namespace ribbon {
static GLuint program = 0;
static GLuint atom_color_tex = 0;

static GLint uniform_loc_normal_mat = 0;
static GLint uniform_loc_view_proj_mat = 0;
static GLint uniform_loc_scale = 0;

void intitialize() {
    GLuint v_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/ribbons.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/ribbons.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/ribbons.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);

    uniform_loc_normal_mat = glGetUniformLocation(program, "u_normal_mat");
    uniform_loc_view_proj_mat = glGetUniformLocation(program, "u_view_proj_mat");
    uniform_loc_scale = glGetUniformLocation(program, "u_scale");

    if (!atom_color_tex) {
        glGenTextures(1, &atom_color_tex);
    }
}

void shutdown() {
    if (program) glDeleteProgram(program);
    if (atom_color_tex) glDeleteTextures(1, &atom_color_tex);
}

}  // namespace ribbon

namespace cartoon {
static GLuint program = 0;
static GLuint atom_color_tex = 0;

static GLint uniform_loc_normal_mat = 0;
static GLint uniform_loc_view_proj_mat = 0;
static GLint uniform_loc_scale = 0;
static GLint uniform_loc_atom_color_tex = 0;

void intitialize() {
    GLuint v_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/cartoon.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/cartoon.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/cartoon.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);

    uniform_loc_normal_mat = glGetUniformLocation(program, "u_normal_mat");
    uniform_loc_view_proj_mat = glGetUniformLocation(program, "u_view_proj_mat");
    uniform_loc_scale = glGetUniformLocation(program, "u_scale");
    uniform_loc_atom_color_tex = glGetUniformLocation(program, "u_atom_color_tex");

    if (!atom_color_tex) {
        glGenTextures(1, &atom_color_tex);
    }
}

void shutdown() {
    if (program) glDeleteProgram(program);
    if (atom_color_tex) glDeleteTextures(1, &atom_color_tex);
}

}  // namespace cartoon

namespace backbone_spline {
static GLuint extract_control_points_program = 0;
static GLuint compute_spline_program = 0;
static GLuint draw_spline_program = 0;

static GLint uniform_loc_ramachandran_tex = 0;

static GLint uniform_loc_num_subdivisions = 0;
static GLint uniform_loc_tension = 0;

static GLint uniform_loc_view_proj_mat = 0;
static GLint uniform_loc_s_color = 0;
static GLint uniform_loc_v_color = 0;
static GLint uniform_loc_t_color = 0;

void initialize() {
    {
        GLuint v_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/backbone_control_points.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/backbone_control_points.geom", GL_GEOMETRY_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader};
        const GLchar* feedback_varyings[] = {"out_control_point", "out_support_vector_xy", "out_support_vector_z_tangent_vector_x", "out_tangent_vector_yz", "out_classification", "out_atom_index"};
        if (!extract_control_points_program) extract_control_points_program = glCreateProgram();
        gl::attach_link_detach_with_transform_feedback(extract_control_points_program, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);

        uniform_loc_ramachandran_tex = glGetUniformLocation(extract_control_points_program, "u_ramachandran_tex");
    }

    {
        GLuint v_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/backbone_compute_spline.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/backbone_compute_spline.geom", GL_GEOMETRY_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader};
        const GLchar* feedback_varyings[] = {"out_control_point", "out_support_vector_xy", "out_support_vector_z_tangent_vector_x", "out_tangent_vector_yz", "out_classification", "out_atom_index"};
        if (!compute_spline_program) compute_spline_program = glCreateProgram();
        gl::attach_link_detach_with_transform_feedback(compute_spline_program, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);

        uniform_loc_num_subdivisions = glGetUniformLocation(compute_spline_program, "u_num_subdivisions");
        uniform_loc_tension = glGetUniformLocation(compute_spline_program, "u_tension");
    }

    {
        GLuint v_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/backbone_draw_spline.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/backbone_draw_spline.geom", GL_GEOMETRY_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(MDUTILS_SHADER_DIR "/backbone_draw_spline.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
            glDeleteShader(f_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader, f_shader};
        if (!draw_spline_program) draw_spline_program = glCreateProgram();
        gl::attach_link_detach(draw_spline_program, shaders);

        uniform_loc_view_proj_mat = glGetUniformLocation(draw_spline_program, "u_view_proj_mat");
        uniform_loc_s_color = glGetUniformLocation(draw_spline_program, "u_s_color");
        uniform_loc_v_color = glGetUniformLocation(draw_spline_program, "u_v_color");
        uniform_loc_t_color = glGetUniformLocation(draw_spline_program, "u_t_color");
    }
}

void shutdown() {
    if (extract_control_points_program) glDeleteProgram(extract_control_points_program);
    if (compute_spline_program) glDeleteProgram(compute_spline_program);
    if (draw_spline_program) glDeleteProgram(draw_spline_program);
}
}  // namespace backbone_spline

void initialize() {
    if (!vao) glGenVertexArrays(1, &vao);
    vdw::initialize();
    licorice::initialize();
    ribbon::intitialize();
    cartoon::intitialize();
    backbone_spline::initialize();
}

void shutdown() {
    if (vao) glDeleteVertexArrays(1, &vao);
    vdw::shutdown();
    licorice::shutdown();
    ribbon::shutdown();
    cartoon::shutdown();
    backbone_spline::shutdown();
}

void draw_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, int32 atom_count, const ViewParam& view_param, float radius_scale) {
	ASSERT(glIsBuffer(atom_position_buffer));
	ASSERT(glIsBuffer(atom_radius_buffer));
	ASSERT(glIsBuffer(atom_color_buffer));
	ASSERT(glIsBuffer(atom_velocity_buffer));

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);
    glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, atom_radius_buffer);
	glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(AtomRadius), (const GLvoid*)0);
	glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, atom_color_buffer);
    glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(AtomColor), (const GLvoid*)0);
    glEnableVertexAttribArray(2);

	glBindBuffer(GL_ARRAY_BUFFER, atom_velocity_buffer);
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(AtomVelocity), (const GLvoid*)0);
	glEnableVertexAttribArray(3);

    glUseProgram(vdw::program);

	mat4 curr_view_to_prev_clip_mat = view_param.previous.matrix.view_proj * view_param.matrix.inverse.view;

    // Uniforms
    glUniform1f(vdw::uniform_loc_radius_scale, radius_scale);
    glUniformMatrix4fv(vdw::uniform_loc_view_mat, 1, GL_FALSE, &view_param.matrix.view[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_proj_mat, 1, GL_FALSE, &view_param.matrix.proj[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_inv_proj_mat, 1, GL_FALSE, &view_param.matrix.inverse.proj[0][0]);
	glUniformMatrix4fv(vdw::uniform_loc_curr_view_to_prev_clip_mat, 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);

    glDrawArrays(GL_POINTS, 0, atom_count);

    glBindVertexArray(0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint bond_buffer, int32 bond_count, const ViewParam& view_param, float radius_scale) {
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glBindBuffer(GL_ARRAY_BUFFER, atom_color_buffer);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(AtomColor), (const GLvoid*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bond_buffer);

    glUseProgram(licorice::program);
    glUniformMatrix4fv(licorice::uniform_loc_view_mat, 1, GL_FALSE, &view_param.matrix.view[0][0]);
    glUniformMatrix4fv(licorice::uniform_loc_proj_mat, 1, GL_FALSE, &view_param.matrix.proj[0][0]);
    glUniform1f(licorice::uniform_loc_radius, 0.25f * radius_scale);

    glDrawElements(GL_LINES, bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void compute_backbone_control_points(GLuint dst_buffer, GLuint atom_position_buffer, GLuint backbone_index_buffer, int num_backbone_indices, GLuint ramachandran_tex) {
    glEnable(GL_RASTERIZER_DISCARD);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, backbone_index_buffer);

    glUseProgram(backbone_spline::extract_control_points_program);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ramachandran_tex);
    glUniform1i(backbone_spline::uniform_loc_ramachandran_tex, 0);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, dst_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawElements(GL_TRIANGLES_ADJACENCY, num_backbone_indices, GL_UNSIGNED_INT, 0);
    glEndTransformFeedback();

    glUseProgram(0);

    glDisable(GL_RASTERIZER_DISCARD);
}

void compute_backbone_spline(GLuint dst_buffer, GLuint control_point_buffer, GLuint control_point_index_buffer, int num_control_point_indices, int num_subdivisions, float tension) {
    glEnable(GL_RASTERIZER_DISCARD);
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, control_point_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, classification));
    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, atom_index));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, control_point_index_buffer);

    glUseProgram(backbone_spline::compute_spline_program);

    glUniform1i(backbone_spline::uniform_loc_num_subdivisions, num_subdivisions);
    glUniform1f(backbone_spline::uniform_loc_tension, tension);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, dst_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawElements(GL_LINE_STRIP_ADJACENCY, num_control_point_indices, GL_UNSIGNED_INT, 0);
    glEndTransformFeedback();

    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
    glDisable(GL_RASTERIZER_DISCARD);
}

void draw_spline(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const ViewParam& view_param, uint32 s_color, uint32 v_color, uint32 t_color) {
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, spline_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, tangent_vector));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spline_index_buffer);

    glUseProgram(backbone_spline::draw_spline_program);
    glUniformMatrix4fv(backbone_spline::uniform_loc_view_proj_mat, 1, GL_FALSE, &(view_param.matrix.view_proj)[0][0]);
    glUniform4fv(backbone_spline::uniform_loc_s_color, 1, &math::convert_color(s_color)[0]);
    glUniform4fv(backbone_spline::uniform_loc_v_color, 1, &math::convert_color(v_color)[0]);
    glUniform4fv(backbone_spline::uniform_loc_t_color, 1, &math::convert_color(t_color)[0]);
    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, int32 num_spline_indices, const ViewParam& view_param) {
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, spline_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, tangent_vector));
    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, atom_index));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spline_index_buffer);

    glBindTexture(GL_TEXTURE_BUFFER, ribbon::atom_color_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    glUseProgram(ribbon::program);
    glUniformMatrix4fv(ribbon::uniform_loc_normal_mat, 1, GL_FALSE, &view_param.matrix.norm[0][0]);
    glUniformMatrix4fv(ribbon::uniform_loc_view_proj_mat, 1, GL_FALSE, &view_param.matrix.view_proj[0][0]);
    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

void draw_cartoon(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, int32 num_spline_indices, const ViewParam& view_param) {
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, spline_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, tangent_vector));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, classification));
    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, atom_index));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spline_index_buffer);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, cartoon::atom_color_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    glUseProgram(cartoon::program);
    glUniform1i(cartoon::uniform_loc_atom_color_tex, 0);
    glUniformMatrix4fv(cartoon::uniform_loc_normal_mat, 1, GL_FALSE, &view_param.matrix.norm[0][0]);
    glUniformMatrix4fv(cartoon::uniform_loc_view_proj_mat, 1, GL_FALSE, &view_param.matrix.view_proj[0][0]);
    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

}  // namespace draw
