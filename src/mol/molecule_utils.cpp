#include "molecule_utils.h"

#include <core/common.h>
#include <core/simd.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/spatial_hash.h>

#include <mol/element.h>
#include <mol/element_utils.h>
#include <mol/aminoacid.h>
#include <mol/aminoacid_utils.h>
#include <mol/trajectory_utils.h>

#include <svd3/svd3.h>
#include <ctype.h>

inline SIMD_TYPE_F apply_pbc(const SIMD_TYPE_F x, const SIMD_TYPE_F box_ext) {
    const SIMD_TYPE_F add = simd::bit_and(simd::cmp_lt(x, SIMD_ZERO_F), box_ext);
    const SIMD_TYPE_F sub = simd::bit_and(simd::cmp_gt(x, box_ext), box_ext);
    const SIMD_TYPE_F res = simd::bit_and(x, simd::sub(add, sub));
    return res;
}

template <typename T>
inline T de_periodize(T pos, T ref_pos, T box_ext) {
    const T delta = pos - ref_pos;
    const T signed_mask = math::sign(delta) * math::step(box_ext * T(0.5), math::abs(delta));
    const T res = pos - box_ext * signed_mask;
    return res;
}

inline SIMD_TYPE_F de_periodize(const SIMD_TYPE_F pos, const SIMD_TYPE_F ref_pos, const SIMD_TYPE_F box_ext) {
    const SIMD_TYPE_F half_box = simd::mul(box_ext, SIMD_SET_F(0.5f));
    const SIMD_TYPE_F delta = simd::sub(pos, ref_pos);
    const SIMD_TYPE_F signed_mask = simd::mul(simd::sign(delta), simd::step(half_box, simd::abs(delta)));
    const SIMD_TYPE_F res = simd::sub(pos, simd::mul(box_ext, signed_mask));
    return res;
}

void translate_ref(soa_vec3 in_out, i64 count, const vec3& translation) {
    for (i64 i = 0; i < count; i++) {
        in_out.x[i] += translation.x;
        in_out.y[i] += translation.y;
        in_out.z[i] += translation.z;
    }
}

void translate(soa_vec3 in_out, i64 count, const vec3& translation) {
    i64 i = 0;

    const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    if (simd_count > 0) {
        SIMD_TYPE_F t_x = SIMD_SET_F(translation.x);
        SIMD_TYPE_F t_y = SIMD_SET_F(translation.y);
        SIMD_TYPE_F t_z = SIMD_SET_F(translation.z);

        for (; i < simd_count; i += SIMD_WIDTH) {
            SIMD_TYPE_F p_x = SIMD_LOAD_F(in_out.x + i);
            SIMD_TYPE_F p_y = SIMD_LOAD_F(in_out.y + i);
            SIMD_TYPE_F p_z = SIMD_LOAD_F(in_out.z + i);

            p_x = simd::add(p_x, t_x);
            p_y = simd::add(p_y, t_y);
            p_z = simd::add(p_z, t_z);

            SIMD_STORE(in_out.x + i, p_x);
            SIMD_STORE(in_out.y + i, p_y);
            SIMD_STORE(in_out.z + i, p_z);
        }
    }

    for (; i < count; i++) {
        in_out.x[i] += translation.x;
        in_out.y[i] += translation.y;
        in_out.z[i] += translation.z;
    }
}

void translate(soa_vec3 out, const soa_vec3 in, i64 count, const vec3& translation) {
    i64 i = 0;

    const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    if (simd_count > 0) {
        SIMD_TYPE_F t_x = SIMD_SET_F(translation.x);
        SIMD_TYPE_F t_y = SIMD_SET_F(translation.y);
        SIMD_TYPE_F t_z = SIMD_SET_F(translation.z);

        for (; i < simd_count; i += SIMD_WIDTH) {
            SIMD_TYPE_F p_x = SIMD_LOAD_F(in.x + i);
            SIMD_TYPE_F p_y = SIMD_LOAD_F(in.y + i);
            SIMD_TYPE_F p_z = SIMD_LOAD_F(in.z + i);

            p_x = simd::add(p_x, t_x);
            p_y = simd::add(p_y, t_y);
            p_z = simd::add(p_z, t_z);

            SIMD_STORE(out.x + i, p_x);
            SIMD_STORE(out.y + i, p_y);
            SIMD_STORE(out.z + i, p_z);
        }
    }

    for (; i < count; i++) {
        out.x[i] += translation.x;
        out.y[i] += translation.y;
        out.z[i] += translation.z;
    }
}

void transform_ref(soa_vec3 in_out, i64 count, const mat4& transformation, float w_comp) {
    for (i64 i = 0; i < count; i++) {
        vec4 v = {in_out.x[i], in_out.y[i], in_out.z[i], w_comp};
        v = transformation * v;
        in_out.x[i] = v.x;
        in_out.y[i] = v.y;
        in_out.z[i] = v.z;
    }
}

void transform(soa_vec3 in_out, i64 count, const mat4& M, float w_comp) {
    const SIMD_TYPE_F m11 = SIMD_SET_F(M[0][0]);
    const SIMD_TYPE_F m12 = SIMD_SET_F(M[0][1]);
    const SIMD_TYPE_F m13 = SIMD_SET_F(M[0][2]);

    const SIMD_TYPE_F m21 = SIMD_SET_F(M[1][0]);
    const SIMD_TYPE_F m22 = SIMD_SET_F(M[1][1]);
    const SIMD_TYPE_F m23 = SIMD_SET_F(M[1][2]);

    const SIMD_TYPE_F m31 = SIMD_SET_F(M[2][0]);
    const SIMD_TYPE_F m32 = SIMD_SET_F(M[2][1]);
    const SIMD_TYPE_F m33 = SIMD_SET_F(M[2][2]);

    const SIMD_TYPE_F m41 = SIMD_SET_F(M[3][0]);
    const SIMD_TYPE_F m42 = SIMD_SET_F(M[3][1]);
    const SIMD_TYPE_F m43 = SIMD_SET_F(M[3][2]);

    const SIMD_TYPE_F w = SIMD_SET_F(w_comp);

    i64 i = 0;
    const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    for (; i < simd_count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x = SIMD_LOAD_F(in_out.x + i);
        const SIMD_TYPE_F y = SIMD_LOAD_F(in_out.y + i);
        const SIMD_TYPE_F z = SIMD_LOAD_F(in_out.z + i);

        const SIMD_TYPE_F m11x = simd::mul(m11, x);
        const SIMD_TYPE_F m21y = simd::mul(m21, y);
        const SIMD_TYPE_F m31z = simd::mul(m31, z);
        const SIMD_TYPE_F m41w = simd::mul(m41, w);

        const SIMD_TYPE_F m12x = simd::mul(m12, x);
        const SIMD_TYPE_F m22y = simd::mul(m22, y);
        const SIMD_TYPE_F m32z = simd::mul(m32, z);
        const SIMD_TYPE_F m42w = simd::mul(m42, w);

        const SIMD_TYPE_F m13x = simd::mul(m13, x);
        const SIMD_TYPE_F m23y = simd::mul(m23, y);
        const SIMD_TYPE_F m33z = simd::mul(m33, z);
        const SIMD_TYPE_F m43w = simd::mul(m43, w);

        const SIMD_TYPE_F res_x = simd::add(simd::add(m11x, m21y), simd::add(m31z, m41w));
        const SIMD_TYPE_F res_y = simd::add(simd::add(m12x, m22y), simd::add(m32z, m42w));
        const SIMD_TYPE_F res_z = simd::add(simd::add(m13x, m23y), simd::add(m33z, m43w));

        SIMD_STORE(in_out.x + i, res_x);
        SIMD_STORE(in_out.y + i, res_y);
        SIMD_STORE(in_out.z + i, res_z);
    }

    for (; i < count; i++) {
        const float x = in_out.x[i];
        const float y = in_out.y[i];
        const float z = in_out.z[i];

        in_out.x[i] = x * M[0][0] + y * M[1][0] + z * M[2][0] + w_comp * M[3][0];
        in_out.y[i] = x * M[0][1] + y * M[1][1] + z * M[2][1] + w_comp * M[3][1];
        in_out.z[i] = x * M[0][2] + y * M[1][2] + z * M[2][2] + w_comp * M[3][2];
    }
}

void transform(soa_vec3 out, const soa_vec3 in, i64 count, const mat4& M, float w_comp) {
    const SIMD_TYPE_F m11 = SIMD_SET_F(M[0][0]);
    const SIMD_TYPE_F m12 = SIMD_SET_F(M[0][1]);
    const SIMD_TYPE_F m13 = SIMD_SET_F(M[0][2]);

    const SIMD_TYPE_F m21 = SIMD_SET_F(M[1][0]);
    const SIMD_TYPE_F m22 = SIMD_SET_F(M[1][1]);
    const SIMD_TYPE_F m23 = SIMD_SET_F(M[1][2]);

    const SIMD_TYPE_F m31 = SIMD_SET_F(M[2][0]);
    const SIMD_TYPE_F m32 = SIMD_SET_F(M[2][1]);
    const SIMD_TYPE_F m33 = SIMD_SET_F(M[2][2]);

    const SIMD_TYPE_F m41 = SIMD_SET_F(M[3][0]);
    const SIMD_TYPE_F m42 = SIMD_SET_F(M[3][1]);
    const SIMD_TYPE_F m43 = SIMD_SET_F(M[3][2]);

    const SIMD_TYPE_F w = SIMD_SET_F(w_comp);

    i64 i = 0;
    const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    for (; i < simd_count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x = SIMD_LOAD_F(in.x + i);
        const SIMD_TYPE_F y = SIMD_LOAD_F(in.y + i);
        const SIMD_TYPE_F z = SIMD_LOAD_F(in.z + i);

        const SIMD_TYPE_F m11x = simd::mul(m11, x);
        const SIMD_TYPE_F m21y = simd::mul(m21, y);
        const SIMD_TYPE_F m31z = simd::mul(m31, z);
        const SIMD_TYPE_F m41w = simd::mul(m41, w);

        const SIMD_TYPE_F m12x = simd::mul(m12, x);
        const SIMD_TYPE_F m22y = simd::mul(m22, y);
        const SIMD_TYPE_F m32z = simd::mul(m32, z);
        const SIMD_TYPE_F m42w = simd::mul(m42, w);

        const SIMD_TYPE_F m13x = simd::mul(m13, x);
        const SIMD_TYPE_F m23y = simd::mul(m23, y);
        const SIMD_TYPE_F m33z = simd::mul(m33, z);
        const SIMD_TYPE_F m43w = simd::mul(m43, w);

        const SIMD_TYPE_F res_x = simd::add(simd::add(m11x, m21y), simd::add(m31z, m41w));
        const SIMD_TYPE_F res_y = simd::add(simd::add(m12x, m22y), simd::add(m32z, m42w));
        const SIMD_TYPE_F res_z = simd::add(simd::add(m13x, m23y), simd::add(m33z, m43w));

        SIMD_STORE(out.x + i, res_x);
        SIMD_STORE(out.y + i, res_y);
        SIMD_STORE(out.z + i, res_z);
    }

    for (; i < count; i++) {
        const float x = in.x[i];
        const float y = in.y[i];
        const float z = in.z[i];

        out.x[i] = x * M[0][0] + y * M[1][0] + z * M[2][0] + w_comp * M[3][0];
        out.y[i] = x * M[0][1] + y * M[1][1] + z * M[2][1] + w_comp * M[3][1];
        out.z[i] = x * M[0][2] + y * M[1][2] + z * M[2][2] + w_comp * M[3][2];
    }
}

void homogeneous_transform(soa_vec3 in_out, i64 count, const mat4& M, float w_comp) {
    for (i64 i = 0; i < count; i++) {
        float& x = in_out.x[i];
        float& y = in_out.y[i];
        float& z = in_out.z[i];

        const float w_scl = 1.0f / (x * M[0][3] + y * M[1][3] + z * M[2][3] + w_comp * M[3][3]);
        in_out.z[i] = (x * M[0][0] + y * M[1][0] + z * M[2][0] + w_comp * M[3][0]) * w_scl;
        in_out.y[i] = (x * M[0][1] + y * M[1][1] + z * M[2][1] + w_comp * M[3][1]) * w_scl;
        in_out.z[i] = (x * M[0][2] + y * M[1][2] + z * M[2][2] + w_comp * M[3][2]) * w_scl;
    }
}

void homogeneous_transform(soa_vec3 out, const soa_vec3 in, i64 count, const mat4& M, float w_comp) {
    for (i64 i = 0; i < count; i++) {
        const float x = in.x[i];
        const float y = in.y[i];
        const float z = in.z[i];

        const float w_scl = 1.0f / (x * M[0][3] + y * M[1][3] + z * M[2][3] + w_comp * M[3][3]);
        out.z[i] = (x * M[0][0] + y * M[1][0] + z * M[2][0] + w_comp * M[3][0]) * w_scl;
        out.y[i] = (x * M[0][1] + y * M[1][1] + z * M[2][1] + w_comp * M[3][1]) * w_scl;
        out.z[i] = (x * M[0][2] + y * M[1][2] + z * M[2][2] + w_comp * M[3][2]) * w_scl;
    }
}

AABB compute_aabb(const soa_vec3 in_pos, i64 count) {
    if (count == 0) {
        return {};
    }

    AABB aabb;

    i64 i = 0;
    if (count > SIMD_WIDTH) {  // @NOTE: There is probably some number where this makes most sense
        SIMD_TYPE_F min_x = SIMD_LOAD_F(in_pos.x);
        SIMD_TYPE_F min_y = SIMD_LOAD_F(in_pos.y);
        SIMD_TYPE_F min_z = SIMD_LOAD_F(in_pos.z);

        SIMD_TYPE_F max_x = min_x;
        SIMD_TYPE_F max_y = min_y;
        SIMD_TYPE_F max_z = min_z;

        i += SIMD_WIDTH;
        const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            const SIMD_TYPE_F x = SIMD_LOAD_F(in_pos.x + i);
            const SIMD_TYPE_F y = SIMD_LOAD_F(in_pos.y + i);
            const SIMD_TYPE_F z = SIMD_LOAD_F(in_pos.z + i);

            min_x = simd::min(min_x, x);
            min_y = simd::min(min_y, y);
            min_z = simd::min(min_z, z);

            max_x = simd::max(max_x, x);
            max_y = simd::max(max_y, y);
            max_z = simd::max(max_z, z);
        }

        aabb.min = {simd::horizontal_min(min_x), simd::horizontal_min(min_y), simd::horizontal_min(min_z)};
        aabb.max = {simd::horizontal_max(max_x), simd::horizontal_max(max_y), simd::horizontal_max(max_z)};
    }

    for (; i < count; i++) {
        const vec3 p = {in_pos.x[i], in_pos.y[i], in_pos.z[i]};
        aabb.min = math::min(aabb.min, p);
        aabb.max = math::max(aabb.max, p);
    }

    return aabb;
}

AABB compute_aabb(const soa_vec3 in_pos, const float in_r[], i64 count) {
    if (count == 0) {
        return {};
    }

    AABB aabb;

    i64 i = 0;
    if (count > SIMD_WIDTH) {  // @NOTE: There is probably some number where this makes most sense
        SIMD_TYPE_F x = SIMD_LOAD_F(in_pos.x);
        SIMD_TYPE_F y = SIMD_LOAD_F(in_pos.y);
        SIMD_TYPE_F z = SIMD_LOAD_F(in_pos.z);
        SIMD_TYPE_F r = SIMD_LOAD_F(in_r);

        SIMD_TYPE_F min_x = simd::sub(x, r);
        SIMD_TYPE_F min_y = simd::sub(y, r);
        SIMD_TYPE_F min_z = simd::sub(z, r);

        SIMD_TYPE_F max_x = simd::add(x, r);
        SIMD_TYPE_F max_y = simd::add(y, r);
        SIMD_TYPE_F max_z = simd::add(z, r);

        i += SIMD_WIDTH;
        const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            x = SIMD_LOAD_F(in_pos.x + i);
            y = SIMD_LOAD_F(in_pos.y + i);
            z = SIMD_LOAD_F(in_pos.z + i);
            r = SIMD_LOAD_F(in_r + i);

            min_x = simd::min(min_x, simd::sub(x, r));
            min_y = simd::min(min_y, simd::sub(y, r));
            min_z = simd::min(min_z, simd::sub(z, r));

            max_x = simd::max(max_x, simd::add(x, r));
            max_y = simd::max(max_y, simd::add(y, r));
            max_z = simd::max(max_z, simd::add(z, r));
        }

        aabb.min = {simd::horizontal_min(min_x), simd::horizontal_min(min_y), simd::horizontal_min(min_z)};
        aabb.max = {simd::horizontal_max(max_x), simd::horizontal_max(max_y), simd::horizontal_max(max_z)};
    }

    for (; i < count; i++) {
        const vec3 p = {in_pos.x[i], in_pos.y[i], in_pos.z[i]};
        aabb.min = math::min(aabb.min, p);
        aabb.max = math::max(aabb.max, p);
    }

    return aabb;
}

vec3 compute_com(const float in_x[], const float in_y[], const float in_z[], i64 count) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_x[0], in_y[0], in_z[0]};

    vec3 sum{0};
    i64 i = 0;

    if (count > SIMD_WIDTH) {
        const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
        SIMD_TYPE_F x = SIMD_LOAD_F(in_x);
        SIMD_TYPE_F y = SIMD_LOAD_F(in_y);
        SIMD_TYPE_F z = SIMD_LOAD_F(in_z);
        i += SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            x = simd::add(x, SIMD_LOAD_F(in_x + i));
            y = simd::add(y, SIMD_LOAD_F(in_y + i));
            z = simd::add(z, SIMD_LOAD_F(in_z + i));
        }
        sum.x = simd::horizontal_add(x);
        sum.y = simd::horizontal_add(y);
        sum.z = simd::horizontal_add(z);
    }

    for (; i < count; i++) {
        sum.x += in_x[i];
        sum.y += in_y[i];
        sum.z += in_z[i];
    }
    return sum / (float)count;
}

vec3 compute_com(const float in_x[], const float in_y[], const float in_z[], const float in_m[], i64 count) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_x[0], in_y[0], in_z[0]};

    vec3 vec_sum{0, 0, 0};
    float mass_sum = 0.0f;
    i64 i = 0;

    const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    if (simd_count > SIMD_WIDTH) {
        SIMD_TYPE_F m = SIMD_LOAD_F(in_m);
        SIMD_TYPE_F x = simd::mul(SIMD_LOAD_F(in_x), m);
        SIMD_TYPE_F y = simd::mul(SIMD_LOAD_F(in_y), m);
        SIMD_TYPE_F z = simd::mul(SIMD_LOAD_F(in_z), m);
        i += SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            const SIMD_TYPE_F mass = SIMD_LOAD_F(in_m + i);
            x = simd::add(x, simd::mul(SIMD_LOAD_F(in_x + i), mass));
            y = simd::add(y, simd::mul(SIMD_LOAD_F(in_y + i), mass));
            z = simd::add(z, simd::mul(SIMD_LOAD_F(in_z + i), mass));
            m = simd::add(m, mass);
        }
        vec_sum.x = simd::horizontal_add(x);
        vec_sum.y = simd::horizontal_add(y);
        vec_sum.z = simd::horizontal_add(z);
        mass_sum = simd::horizontal_add(m);
    }

    for (; i < count; i++) {
        const float mass = in_m[i];
        vec_sum.x += in_x[i] * mass;
        vec_sum.y += in_y[i] * mass;
        vec_sum.z += in_z[i] * mass;
        mass_sum += mass;
    }
    return vec_sum / mass_sum;
}

vec3 compute_com_periodic_ref(const soa_vec3 in_pos, const float in_m[], i64 count, const mat3& box) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_pos.x[0], in_pos.y[0], in_pos.z[0]};

    const vec3 box_ext = box * vec3(1.0f);

    float mass_sum = in_m[0];
    vec3 vec_sum = vec3(in_pos.x[0], in_pos.y[0], in_pos.z[0]) * in_m[0];

    for (i64 i = 1; i < count; i++) {
        const vec3 pos = {in_pos.x[i], in_pos.y[i], in_pos.z[i]};
        const float mass = in_m[i];
        const vec3 com = vec_sum / mass_sum;
        vec_sum += de_periodize(pos, com, box_ext) * mass;
        mass_sum += mass;
    }
    return vec_sum / mass_sum;
}

vec3 compute_com_periodic(const soa_vec3 in_pos, const float in_mass[], i64 count, const mat3& box) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_pos.x[0], in_pos.y[0], in_pos.z[0]};

    const vec3 box_ext = box * vec3(1.0f);

    vec3 vec_sum = {0, 0, 0};
    float mass_sum = 0;
    i64 i = 0;

    if (count > SIMD_WIDTH) {
        const i64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
        const SIMD_TYPE_F box_ext_x = SIMD_SET_F(box[0][0]);
        const SIMD_TYPE_F box_ext_y = SIMD_SET_F(box[1][1]);
        const SIMD_TYPE_F box_ext_z = SIMD_SET_F(box[2][2]);

        SIMD_TYPE_F m_sum = SIMD_LOAD_F(in_mass);
        SIMD_TYPE_F x_sum = simd::mul(SIMD_LOAD_F(in_pos.x), m_sum);
        SIMD_TYPE_F y_sum = simd::mul(SIMD_LOAD_F(in_pos.y), m_sum);
        SIMD_TYPE_F z_sum = simd::mul(SIMD_LOAD_F(in_pos.z), m_sum);

        i += SIMD_WIDTH;

        for (; i < simd_count; i += SIMD_WIDTH) {
            const SIMD_TYPE_F m = SIMD_LOAD_F(in_mass + i);
            const SIMD_TYPE_F x = SIMD_LOAD_F(in_pos.x + i);
            const SIMD_TYPE_F y = SIMD_LOAD_F(in_pos.y + i);
            const SIMD_TYPE_F z = SIMD_LOAD_F(in_pos.z + i);

            const SIMD_TYPE_F x_com = simd::div(x_sum, m_sum);
            const SIMD_TYPE_F y_com = simd::div(y_sum, m_sum);
            const SIMD_TYPE_F z_com = simd::div(z_sum, m_sum);

            const SIMD_TYPE_F x_dp = de_periodize(x, x_com, box_ext_x);
            const SIMD_TYPE_F y_dp = de_periodize(y, y_com, box_ext_y);
            const SIMD_TYPE_F z_dp = de_periodize(z, z_com, box_ext_z);

            x_sum = simd::add(x_sum, simd::mul(x_dp, m));
            y_sum = simd::add(y_sum, simd::mul(y_dp, m));
            z_sum = simd::add(z_sum, simd::mul(z_dp, m));
            m_sum = simd::add(m_sum, m);
        }

        float x[SIMD_WIDTH];
        float y[SIMD_WIDTH];
        float z[SIMD_WIDTH];
        float m[SIMD_WIDTH];

        SIMD_STORE(x, x_sum);
        SIMD_STORE(y, y_sum);
        SIMD_STORE(z, z_sum);
        SIMD_STORE(m, m_sum);

        vec_sum = {x[0], y[0], z[0]};
        mass_sum = m[0];

        for (int j = 1; j < SIMD_WIDTH; j++) {
            const vec3 sum_p = {x[j], y[j], z[j]};
            const float sum_m = m[j];
            const vec3 com = vec_sum / mass_sum;
            vec_sum.x += de_periodize(sum_p.x / sum_m, com.x, box_ext.x) * sum_m;
            vec_sum.y += de_periodize(sum_p.y / sum_m, com.y, box_ext.x) * sum_m;
            vec_sum.z += de_periodize(sum_p.z / sum_m, com.z, box_ext.x) * sum_m;
            mass_sum += sum_m;
        }
    }

    for (; i < count; i++) {
        const float mass = in_mass[i];
        const vec3 pos = {in_pos.x[i], in_pos.y[i], in_pos.z[i]};
        const vec3 com = vec_sum / mass_sum;
        vec_sum.x += de_periodize(pos.x, com.x, box_ext.x) * mass;
        vec_sum.y += de_periodize(pos.y, com.y, box_ext.y) * mass;
        vec_sum.z += de_periodize(pos.z, com.z, box_ext.z) * mass;
        mass_sum += mass;
    }

    return vec_sum / mass_sum;
}

mat3 compute_covariance_matrix(const soa_vec3 in_pos, const float mass[], i64 count, const vec3& com) {
    mat3 A{0};
    float mass_sum = 0.0f;
    for (i64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = in_pos.x[i] - com.x;
        const float qy = in_pos.y[i] - com.y;
        const float qz = in_pos.z[i] - com.z;
        const float m = mass[i];
        mass_sum += m;

        A[0][0] += m * qx * qx;
        A[0][1] += m * qy * qx;
        A[0][2] += m * qz * qx;
        A[1][0] += m * qx * qy;
        A[1][1] += m * qy * qy;
        A[1][2] += m * qz * qy;
        A[2][0] += m * qx * qz;
        A[2][1] += m * qy * qz;
        A[2][2] += m * qz * qz;
    }

    return A / mass_sum;
}

#define ARGS(M) M[0][0], M[1][0], M[2][0], M[0][1], M[1][1], M[2][1], M[0][2], M[1][2], M[2][2]
EigenFrame compute_eigen_frame(const soa_vec3 in_pos, const float in_mass[], i64 count) {
    const vec3 com = compute_com(in_pos, count);
    const mat3 M = compute_covariance_matrix(in_pos, in_mass, count, com);
    mat3 U, S, V;
    svd(ARGS(M), ARGS(U), ARGS(S), ARGS(V));
    const mat3 Ut = glm::transpose(U);

    const float e_val[] = {math::sqrt(S[0][0]), math::sqrt(S[1][1]), math::sqrt(S[2][2])};
    const vec3 e_vec[] = {Ut[0], Ut[1], Ut[2]};
    int l[3] = {0, 1, 2};

    const auto swap = [](int& x, int& y) {
        int tmp = x;
        x = y;
        y = tmp;
    };

    if (e_val[l[0]] < e_val[l[1]]) swap(l[0], l[1]);
    if (e_val[l[1]] < e_val[l[2]]) swap(l[1], l[2]);
    if (e_val[l[0]] < e_val[l[1]]) swap(l[0], l[1]);

    EigenFrame ef;
    ef.values[0] = e_val[l[0]];
    ef.values[1] = e_val[l[1]];
    ef.values[2] = e_val[l[2]];

    ef.vectors[0] = e_vec[l[0]];
    ef.vectors[1] = e_vec[l[1]];
    ef.vectors[2] = e_vec[l[2]];

    return ef;
}
#undef ARGS

void recenter_trajectory(MoleculeDynamic* dynamic, AtomRange range) {
    ASSERT(dynamic);
    if (!dynamic->operator bool()) {
        LOG_ERROR("Dynamic is not valid.");
        return;
    }

    const auto& mol = dynamic->molecule;
    auto& traj = dynamic->trajectory;
    const i64 count = range.ext();
    const float* mass = mol.atom.mass + range.beg;

    for (auto& frame : traj.frame_buffer) {
        const soa_vec3 range_pos = frame.atom_position + range.beg;
        const vec3 com = compute_com_periodic(range_pos, mass, count, frame.box);
        const vec3 translation = frame.box * vec3(0.5f) - com;
        translate(frame.atom_position, mol.atom.count, translation);
        apply_pbc(frame.atom_position, mol.atom.mass, mol.residue.atom_range, mol.residue.count, frame.box);
    }
}

void linear_interpolation_ref(soa_vec3 out, const soa_vec3 in[2], i64 count, float t) {
    for (i64 i = 0; i < count; i++) {
        out.x[i] = in[0].x[i] * (1.0f - t) + in[1].x[i] * t;
        out.y[i] = in[0].y[i] * (1.0f - t) + in[1].y[i] * t;
        out.z[i] = in[0].z[i] * (1.0f - t) + in[1].z[i] * t;
    }
}

void linear_interpolation(soa_vec3 out, const soa_vec3 in[2], i64 count, float t) {
    for (i64 i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in[0].x + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in[0].y + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in[0].z + i);

        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in[1].x + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in[1].y + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in[1].z + i);

        const SIMD_TYPE_F x = simd::lerp(x0, x1, t);
        const SIMD_TYPE_F y = simd::lerp(y0, y1, t);
        const SIMD_TYPE_F z = simd::lerp(z0, z1, t);

        simd::store(out.x + i, x);
        simd::store(out.y + i, y);
        simd::store(out.z + i, z);
    }
}

void linear_interpolation_pbc_ref(soa_vec3 out, const soa_vec3 in[2], i64 count, float t, const mat3& sim_box) {
    const float box_ext_x = sim_box[0][0];
    const float box_ext_y = sim_box[1][1];
    const float box_ext_z = sim_box[2][2];

    for (i64 i = 0; i < count; i++) {
        float x0 = in[0].x[i];
        float y0 = in[0].y[i];
        float z0 = in[0].z[i];

        float x1 = in[1].x[i];
        float y1 = in[1].x[i];
        float z1 = in[1].x[i];

        x1 = de_periodize(x1, x0, box_ext_x);
        y1 = de_periodize(y1, y0, box_ext_y);
        z1 = de_periodize(z1, z0, box_ext_z);

        const float x = x0 * (1.0f - t) + x1 * t;
        const float y = y0 * (1.0f - t) + y1 * t;
        const float z = z0 * (1.0f - t) + z1 * t;

        out.x[i] = x;
        out.y[i] = y;
        out.z[i] = z;
    }
}

void linear_interpolation_pbc(soa_vec3 out, const soa_vec3 in[2], i64 count, float t, const mat3& sim_box) {
    const SIMD_TYPE_F box_ext_x = SIMD_SET_F(sim_box[0][0]);
    const SIMD_TYPE_F box_ext_y = SIMD_SET_F(sim_box[1][1]);
    const SIMD_TYPE_F box_ext_z = SIMD_SET_F(sim_box[2][2]);

    for (i64 i = 0; i < count; i += SIMD_WIDTH) {
        SIMD_TYPE_F x0 = SIMD_LOAD_F(in[0].x + i);
        SIMD_TYPE_F y0 = SIMD_LOAD_F(in[0].y + i);
        SIMD_TYPE_F z0 = SIMD_LOAD_F(in[0].z + i);

        SIMD_TYPE_F x1 = SIMD_LOAD_F(in[1].x + i);
        SIMD_TYPE_F y1 = SIMD_LOAD_F(in[1].y + i);
        SIMD_TYPE_F z1 = SIMD_LOAD_F(in[1].z + i);

        x1 = de_periodize(x1, x0, box_ext_x);
        y1 = de_periodize(y1, y0, box_ext_y);
        z1 = de_periodize(z1, z0, box_ext_z);

        const SIMD_TYPE_F x = simd::lerp(x0, x1, t);
        const SIMD_TYPE_F y = simd::lerp(y0, y1, t);
        const SIMD_TYPE_F z = simd::lerp(z0, z1, t);

        simd::store(out.x + i, x);
        simd::store(out.y + i, y);
        simd::store(out.z + i, z);
    }
}

// clang-format off
void cubic_interpolation(soa_vec3 out, const soa_vec3 in[4], i64 count, float t) {
    for (i64 i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in[0].x + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in[0].y + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in[0].z + i);
                                           
        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in[1].x + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in[1].y + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in[1].z + i);
                                    
        const SIMD_TYPE_F x2 = SIMD_LOAD_F(in[2].x + i);
        const SIMD_TYPE_F y2 = SIMD_LOAD_F(in[2].y + i);
        const SIMD_TYPE_F z2 = SIMD_LOAD_F(in[2].z + i);
                                        
        const SIMD_TYPE_F x3 = SIMD_LOAD_F(in[3].x + i);
        const SIMD_TYPE_F y3 = SIMD_LOAD_F(in[3].y + i);
        const SIMD_TYPE_F z3 = SIMD_LOAD_F(in[3].z + i);

        const SIMD_TYPE_F x = simd::cubic_spline(x0, x1, x2, x3, t);
        const SIMD_TYPE_F y = simd::cubic_spline(y0, y1, y2, y3, t);
        const SIMD_TYPE_F z = simd::cubic_spline(z0, z1, z2, z3, t);

        simd::store(out.x + i, x);
        simd::store(out.y + i, y);
        simd::store(out.z + i, z);
    }
}

void cubic_interpolation_pbc_ref(soa_vec3 out, const soa_vec3 in[4], i64 count, float t, const mat3& sim_box) {
    const vec3 box_ext = sim_box * vec3(1.0f);

    for (i64 i = 0; i < count; i++) {
        float x0 = in[0].x[i];
        float y0 = in[0].y[i];
        float z0 = in[0].z[i];

        float x1 = in[1].x[i];
        float y1 = in[1].y[i];
        float z1 = in[1].z[i];

        float x2 = in[2].x[i];
        float y2 = in[2].y[i];
        float z2 = in[2].z[i];

        float x3 = in[3].x[i];
        float y3 = in[3].y[i];
        float z3 = in[3].z[i];

        x0 = de_periodize(x0, x1, box_ext.x);
        x2 = de_periodize(x2, x1, box_ext.x);
        x3 = de_periodize(x3, x2, box_ext.x);

        y0 = de_periodize(y0, y1, box_ext.y);
        y2 = de_periodize(y2, y1, box_ext.y);
        y3 = de_periodize(y3, y2, box_ext.y);

        z0 = de_periodize(z0, z1, box_ext.z);
        z2 = de_periodize(z2, z1, box_ext.z);
        z3 = de_periodize(z3, z2, box_ext.z);

        const float x = math::cubic_spline(x0, x1, x2, x3, t);
        const float y = math::cubic_spline(y0, y1, y2, y3, t);
        const float z = math::cubic_spline(z0, z1, z2, z3, t);

        out.x[i] = x;
        out.y[i] = y;
        out.z[i] = z;
    }
}

void cubic_interpolation_pbc(soa_vec3 out, const soa_vec3 in[4], i64 count, float t, const mat3& sim_box) {
    const SIMD_TYPE_F box_ext_x = SIMD_SET_F(sim_box[0][0]);
    const SIMD_TYPE_F box_ext_y = SIMD_SET_F(sim_box[1][1]);
    const SIMD_TYPE_F box_ext_z = SIMD_SET_F(sim_box[2][2]);

    for (i64 i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in[0].x + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in[0].y + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in[0].z + i);

        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in[1].x + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in[1].y + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in[1].z + i);

        const SIMD_TYPE_F x2 = SIMD_LOAD_F(in[2].x + i);
        const SIMD_TYPE_F y2 = SIMD_LOAD_F(in[2].y + i);
        const SIMD_TYPE_F z2 = SIMD_LOAD_F(in[2].z + i);

        const SIMD_TYPE_F x3 = SIMD_LOAD_F(in[3].x + i);
        const SIMD_TYPE_F y3 = SIMD_LOAD_F(in[3].y + i);
        const SIMD_TYPE_F z3 = SIMD_LOAD_F(in[3].z + i);

        const SIMD_TYPE_F dp_x0 = de_periodize(x0, x1, box_ext_x);
        const SIMD_TYPE_F dp_x2 = de_periodize(x2, x1, box_ext_x);
        const SIMD_TYPE_F dp_x3 = de_periodize(x3, dp_x2, box_ext_x);

        const SIMD_TYPE_F dp_y0 = de_periodize(y0, y1, box_ext_y);
        const SIMD_TYPE_F dp_y2 = de_periodize(y2, y1, box_ext_y);
        const SIMD_TYPE_F dp_y3 = de_periodize(y3, dp_y2, box_ext_y);

        const SIMD_TYPE_F dp_z0 = de_periodize(z0, z1, box_ext_z);
        const SIMD_TYPE_F dp_z2 = de_periodize(z2, z1, box_ext_z);
        const SIMD_TYPE_F dp_z3 = de_periodize(z3, dp_z2, box_ext_z);

        const SIMD_TYPE_F x = simd::cubic_spline(dp_x0, x1, dp_x2, dp_x3, t);
        const SIMD_TYPE_F y = simd::cubic_spline(dp_y0, y1, dp_y2, dp_y3, t);
        const SIMD_TYPE_F z = simd::cubic_spline(dp_z0, z1, dp_z2, dp_z3, t);

        simd::store(out.x + i, x);
        simd::store(out.y + i, y);
        simd::store(out.z + i, z);
    }
}

void apply_pbc(soa_vec3 in_out, const float in_mass[], i64 count, const mat3& sim_box) {
    const vec3 box_ext = sim_box * vec3(1.0f);
    const vec3 one_over_box_ext = 1.0f / box_ext;

    const vec3 com = compute_com_periodic(in_out, in_mass, count, sim_box);
    const vec3 com_dp = math::fract(com * one_over_box_ext) * box_ext;

    for (i64 i = 0; i < count; i++) {
        in_out.x[i] = de_periodize(in_out.x[i], com_dp.x, box_ext.x);
        in_out.y[i] = de_periodize(in_out.y[i], com_dp.y, box_ext.y);
        in_out.z[i] = de_periodize(in_out.z[i], com_dp.z, box_ext.z);
    }
}

void apply_pbc(soa_vec3 in_out_pos, const float in_mass[], const AtomRange in_ranges[], i64 num_ranges, const mat3& sim_box) {
    for (i64 i = 0; i < num_ranges; i++) {
        apply_pbc(in_out_pos + in_ranges[i].beg, in_mass + in_ranges[i].beg, in_ranges[i].ext(), sim_box);
    }
}

void compute_backbone_angles(BackboneAngle out_angle[], const soa_vec3 in_pos, const BackboneSegment backbone_segments[], i64 num_segments) {
    ASSERT(out_angle);
    ASSERT(backbone_segments);
    float phi = 0, psi = 0;

    if (num_segments < 2) {
        return;
    }

    ASSERT(valid_segment(backbone_segments[0]));
    vec3 n =  { in_pos.x[backbone_segments[0].n_idx],  in_pos.y[backbone_segments[0].n_idx],  in_pos.z[backbone_segments[0].n_idx]  };
    vec3 ca = { in_pos.x[backbone_segments[0].ca_idx], in_pos.y[backbone_segments[0].ca_idx], in_pos.z[backbone_segments[0].ca_idx] };
    vec3 c =  { in_pos.x[backbone_segments[0].c_idx],  in_pos.y[backbone_segments[0].c_idx],  in_pos.z[backbone_segments[0].c_idx]  };

    vec3 c_prev = c;
    vec3 n_next = {in_pos.x[backbone_segments[1].n_idx], in_pos.y[backbone_segments[1].n_idx], in_pos.z[backbone_segments[1].n_idx]};
    phi = 0.0f;
    psi = math::dihedral_angle(n, ca, c, n_next);
    out_angle[0] = {phi, psi};

    for (i64 i = 1; i < num_segments - 1; i++) {
        ASSERT(valid_segment(backbone_segments[i]));

        c_prev = c;
        n = n_next;
        ca = {in_pos.x[backbone_segments[i].ca_idx], in_pos.y[backbone_segments[i].ca_idx], in_pos.z[backbone_segments[i].ca_idx]};
        c = {in_pos.x[backbone_segments[i].c_idx], in_pos.y[backbone_segments[i].c_idx], in_pos.z[backbone_segments[i].c_idx]};
        n_next = {in_pos.x[backbone_segments[i + 1].n_idx], in_pos.y[backbone_segments[i + 1].n_idx], in_pos.z[backbone_segments[i + 1].n_idx]};

        phi = math::dihedral_angle(c_prev, n, ca, c);
        psi = math::dihedral_angle(n, ca, c, n_next);
        out_angle[i] = {phi, psi};
    }

    auto N = num_segments - 1;
    ASSERT(valid_segment(backbone_segments[N]));

    c_prev = c;
    n = n_next;
    ca = {in_pos.x[backbone_segments[N].ca_idx], in_pos.y[backbone_segments[N].ca_idx], in_pos.z[backbone_segments[N].ca_idx]};
    c = {in_pos.x[backbone_segments[N].c_idx], in_pos.y[backbone_segments[N].c_idx], in_pos.z[backbone_segments[N].c_idx]};

    phi = math::dihedral_angle(c_prev, n, ca, c);
    psi = 0.0f;
    out_angle[N] = {phi, psi};
}

void compute_atom_radius(float out_radius[], const Element in_element[], i64 count) {
    for (i64 i = 0; i < count; i++) {
        out_radius[i] = element::vdw_radius(in_element[i]);
    }
}

void compute_atom_mass(float out_mass[], const Element in_element[], i64 count) {
    for (i64 i = 0; i < count; i++) {
        out_mass[i] = element::vdw_radius(in_element[i]);
    }
}

bool is_amino_acid(const Label& res_label) { return get_amino_acid_from_string(res_label) != AminoAcid::Unknown; }

static constexpr CStringView dna_residues[12] = {"DA", "DA3", "DA5", "DC", "DC3", "DC5", "DG", "DG3", "DG5", "DT", "DT3", "DT5"};
bool is_dna(const Label& res_label) {
    for (auto dna_res : dna_residues) {
        if (compare(res_label, dna_res)) return true;
    }
    return false;
}

DynamicArray<Label> get_unique_residue_types(const MoleculeStructure& mol) {
    DynamicArray<Label> types = {};
    Label cur_lbl = {};
    for (const auto& name : get_residue_names(mol)) {
        if (name != cur_lbl) {
            types.push_back(name);
        }
    }
    return types;
}

DynamicArray<ResIdx> get_residues_by_name(const MoleculeStructure& mol, CStringView name) {
    DynamicArray<ResIdx> residues;
    for (ResIdx i = 0; i < (ResIdx)mol.residue.count; i++) {
        if (mol.residue.name[i] == name) {
            residues.push_back(i);
        }
    }
    return residues;
}

bool atom_ranges_match(const MoleculeStructure& mol, AtomRange range_a, AtomRange range_b) {
    if (range_a.ext() == 0 || range_b.ext() == 0) return false;
    if (range_a.ext() != range_b.ext()) return false;

    const auto ele_a = get_elements(mol).subarray(range_a);
    const auto ele_b = get_elements(mol).subarray(range_b);
    const auto lbl_a = get_labels(mol).subarray(range_a);
    const auto lbl_b = get_labels(mol).subarray(range_b);
    int a_i = range_a.beg;
    int b_i = range_b.beg;
    while (a_i != range_a.end && b_i != range_b.end) {
        if (ele_a[a_i] != ele_b[b_i]) return false;
        if (lbl_a[a_i] != lbl_b[b_i]) return false;
    }
    return true;
}

DynamicArray<AtomRange> find_equivalent_structures(const MoleculeStructure& mol, AtomRange ref) {
    DynamicArray<AtomRange> matches = {};

    const auto ele = get_elements(mol);
    const auto lbl = get_labels(mol);

    const auto ele_ref = ele.subarray(ref);
    const auto lbl_ref = lbl.subarray(ref);

    int j = 0;
    for (int i = 0; i < mol.atom.count; i++) {
        if (ref.beg <= i && i < ref.end) {
            j = 0;
            i = ref.end;
        }
        if (ele[i] == ele_ref[j] && (compare(lbl[i], lbl_ref[j]))) {
            j++;
        }
        if (j == ref.ext()) {
            matches.push_back({i - j + 1, i + 1});
            j = 0;
        }
    }

    return matches;
}

// @NOTE: THIS IS STUPID
bool structure_match(const MoleculeStructure& mol, Bitfield ref_mask, int ref_offset, int test_offset, Bitfield used_mask) {

    const auto ele = get_elements(mol);
    // const auto lbl = get_labels(mol);
    // const auto res_idx = get_residue_indices(mol);

    for (int i = 0; i < ref_mask.size(); ++i) {
        if (ref_mask[i] && !used_mask[test_offset + i]) {
            const auto ref_element = ele[ref_offset + i];
            const auto element = ele[test_offset + i];
            if (element != ref_element) return false;

            // const auto& ref_label = lbl[mask_offset + i];
            // const auto& label = lbl[structure_offset + i];
            // if (compare(label, ref_label) == false) return false;

            // if (compare(mol.residues[res_idx[mask_offset + i]].name, mol.residues[res_idx[mask_offset + i]].name) == false) return false;
            // if (mol.residues[res_idx[mask_offset + i]].id != mol.residues[res_idx[mask_offset + i]].id) return false;
        }
    }

    return true;
};

DynamicArray<int> find_equivalent_structures(const MoleculeStructure& mol, Bitfield mask, int offset) {
    const auto ele = get_elements(mol);
    const auto lbl = get_labels(mol);

    Bitfield used_atoms;
    bitfield::init(&used_atoms, mol.atom.count);
    defer { bitfield::free(&used_atoms); };
    bitfield::clear_all(used_atoms);

    for (int i = 0; i < mask.size(); i++) {
        if (mask[i]) bitfield::set_bit(used_atoms, offset + i);
    }

    i64 count = bitfield::number_of_bits_set(mask);
    if (count == 0) {
        LOG_WARNING("Supplied mask was empty");
        return {};
    }

#if 0
    void* mem = TMP_MALLOC(count * sizeof(float) * 4);
    defer { TMP_FREE(mem); };
    float* x = (float*)mem + 0 * count;
    float* y = (float*)mem + 1 * count;
    float* z = (float*)mem + 2 * count;
    float* m = (float*)mem + 3 * count;

    bitfield::gather_masked(x, mol.atom.position.x, mask, offset);
    bitfield::gather_masked(y, mol.atom.position.y, mask, offset);
    bitfield::gather_masked(z, mol.atom.position.z, mask, offset);
    bitfield::gather_masked(m, mol.atom.mass, mask, offset);

    const EigenFrame ref_eigen = compute_eigen_frame(x, y, z, m, count, compute_com(x, y, z, m, count));
#endif

    DynamicArray<int> offsets = {};
    for (int i = 0; i < mol.atom.count - mask.size(); ++i) {
        if (used_atoms[i]) continue;

        if (structure_match(mol, mask, offset, i, used_atoms)) {
#if 0
            bitfield::gather_masked(x, mol.atom.position.x, mask, i);
            bitfield::gather_masked(y, mol.atom.position.y, mask, i);
            bitfield::gather_masked(z, mol.atom.position.z, mask, i);
            bitfield::gather_masked(m, mol.atom.mass, mask);

            const EigenFrame eigen = compute_eigen_frame(x, y, z, m, count, compute_com(x, y, z, m, count));

            const float ratio_x = ref_eigen.value[0] / eigen.value[0];
            const float ratio_y = ref_eigen.value[1] / eigen.value[1];
            const float ratio_z = ref_eigen.value[2] / eigen.value[2];
#endif
            for (int j = 0; j < mask.size(); j++) {
                if (mask[j]) bitfield::set_bit(used_atoms, i + j);
            }

            offsets.push_back(i);
            // i += count - 1;
        }
    }

    return offsets;
}