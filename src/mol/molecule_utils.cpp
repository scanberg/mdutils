#include "molecule_utils.h"
#include <core/common.h>
#include <core/hash.h>
#include <core/log.h>
#include <mol/trajectory_utils.h>
#include <mol/spatial_hash.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <ctype.h>

#include <core/simd.h>

//#include "compute_velocity_ispc.h"
//#include "interpolate_position_linear_ispc.h"

inline __m128 de_periodize(const __m128 a, const __m128 b, const __m128 full_ext, const __m128 half_ext) {
    const __m128 delta = simd::sub(b, a);
    const __m128 signed_mask = simd::mul(delta, simd::step(half_ext, simd::abs(delta)));
    const __m128 res = simd::sub(b, simd::mul(full_ext, signed_mask));
    return res;
}

inline __m128 apply_pbc(const __m128 x, const __m128 box_ext) {
    const __m128 add = simd::bit_and(simd::cmp_lt(x, simd::zero_float4()), box_ext);
    const __m128 sub = simd::bit_and(simd::cmp_gt(x, box_ext), box_ext);
    const __m128 res = simd::bit_and(x, simd::sub(add, sub));
    return res;
}

inline float de_periodize(float a, float b, float full_ext, float half_ext) {
    const float delta = p - ref;
    const float signed_mask = math::sign(delta) * math::step(half_ext, math::abs(delta));
    const float res = p - full_ext * signed_mask;
    return res;
}

void translate_positions(float* pos_x, float* pos_y, float* pos_z, int64 count, const vec3& translation) {
    __m128 t_x = _mm_set_ps1(translation.x);
    __m128 t_y = _mm_set_ps1(translation.y);
    __m128 t_z = _mm_set_ps1(translation.z);

    // No worries about loading out of range since we have strong guarantees on alignment
    for (int64 i = 0; i < count; i += 4) {
        __m128 p_x = _mm_load_ps(pos_x + i);
        __m128 p_y = _mm_load_ps(pos_y + i);
        __m128 p_z = _mm_load_ps(pos_z + i);

        p_x = _mm_add_ps(p_x, t_x);
        p_y = _mm_add_ps(p_y, t_y);
        p_z = _mm_add_ps(p_z, t_z);

        _mm_store_ps(pos_x + i, p_x);
        _mm_store_ps(pos_y + i, p_y);
        _mm_store_ps(pos_z + i, p_z);
    }
}

void transform_positions(float* pos_x, float* pos_y, float* pos_z, int64 count, const mat4& transformation) {
    for (int64 i = 0; i < count; i++) {
        vec4 p{pos_x[i], pos_y[i], pos_z[i], 1.0f};
        p = transformation * p;
        pos_x[i] = p.x;
        pos_y[i] = p.y;
        pos_z[i] = p.z;
    }
}

void compute_bounding_box(vec3* min_box, vec3* max_box, const float* pos_x, const float* pos_y, const float* pos_z, int64 count) {
    ASSERT(min_box);
    ASSERT(max_box);

    if (count == 0) {
        *min_box = *max_box = vec3(0);
        return;
    }

    *min_box = *max_box = vec3(pos_x[0], pos_y[0], pos_z[0]);
    for (int64 i = 1; i < count; i++) {
        const vec3 p = vec3(pos_x[i], pos_y[i], pos_z[i]);
        *min_box = math::min(*min_box, p);
        *max_box = math::max(*max_box, p);
    }
}

void compute_bounding_box(vec3* min_box, vec3* max_box, const float* pos_x, const float* pos_y, const float* pos_z, const float* radii, int64 count) {
    ASSERT(min_box);
    ASSERT(max_box);

    if (count == 0) {
        *min_box = *max_box = vec3(0);
        return;
    }

    *min_box = *max_box = vec3(pos_x[0], pos_y[0], pos_z[0]);
    for (int64 i = 1; i < count; i++) {
        const vec3 p = vec3(pos_x[i], pos_y[i], pos_z[i]);
        const float r = radii[i];
        *min_box = math::min(*min_box, p - r);
        *max_box = math::max(*max_box, p + r);
    }
}

vec3 compute_com(const float* pos_x, const float* pos_y, const float* pos_z, int64 count) {
    if (count == 0) return vec3(0);
    if (count == 1) return {pos_x[0], pos_y[0], pos_z[0]};

    vec3 sum{0};
    for (int64 i = 0; i < count; i++) {
        const vec3 pos = {pos_x[i], pos_y[i], pos_z[i]};
        sum += pos;
    }
    return sum / (float)count;
}

vec3 compute_com(const float* pos_x, const float* pos_y, const float* pos_z, const float* mass, int64 count) {
    if (count == 0) return {0, 0, 0};
    if (count == 1) return {pos_x[0], pos_y[0], pos_z[0]};

    vec3 sum{0};
    for (int32 i = 0; i < count; i++) {
        const vec3 pos = {pos_x[i], pos_y[i], pos_z[i]};
        sum += pos * mass[i];
    }

    return sum / (float)count;
}

vec3 compute_com(const float* pos_x, const float* pos_y, const float* pos_z, const Element* element, int64 count) {
    if (count == 0) return {0, 0, 0};
    if (count == 1) return {pos_x[0], pos_y[0], pos_z[0]};

    vec3 sum{0};
    for (int32 i = 0; i < count; i++) {
        const vec3 pos = {pos_x[i], pos_y[i], pos_z[i]};
        const float mass = element::atomic_mass(element[i]);
        sum += pos * mass;
    }

    return sum / (float)count;
}

/*
void recenter_trajectory(MoleculeDynamic* dynamic, ResIdx center_res_idx) {
    ASSERT(dynamic);
    if (!dynamic->operator bool()) {
        LOG_ERROR("Dynamic is not valid.");
        return;
    }

    auto& mol = dynamic->molecule;
    const auto elements = get_elements(mol);
    const auto residues = get_residues(mol);
    const auto chains = get_chains(mol);

    const int32 res_beg = residues[center_res_idx].atom_idx.beg;
    const int32 res_end = residues[center_res_idx].atom_idx.end;
    const auto res_ele = elements.subarray(res_beg, res_end - res_beg);

    for (int32 f_idx = 0; f_idx < dynamic->trajectory.num_frames; f_idx++) {
        auto frame = get_trajectory_frame(dynamic->trajectory, f_idx);
        auto positions = frame.atom_positions;
        const auto box_ext = frame.box * vec3(1);
        const auto box_center = frame.box * vec3(0.5f);

        const auto res_pos = positions.subarray(res_beg, res_end - res_beg);
        const vec3 target_center = compute_periodic_com(res_pos, res_ele, box_ext);

        const vec3 delta = box_center - target_center;

        for (int64 i = 0; i < chains.size(); i++) {
            const auto p = positions.subarray(chains[i].atom_idx.beg, chains[i].atom_idx.end - chains[i].atom_idx.beg);
            const auto e = elements.subarray(chains[i].atom_idx.beg, chains[i].atom_idx.end - chains[i].atom_idx.beg);
            const vec3 com_chain = compute_periodic_com(p, e, box_ext);
            const vec3 new_com = de_periodize(box_center, com_chain + delta, box_ext);
            const vec3 diff = new_com - com_chain;
            translate_positions(p, diff);
        }
    }
}
*/

void linear_interpolation(float* out_x, float* out_y, float* out_z, const float* p0_x, const float* p0_y, const float* p0_z, const float* p1_x, const float* p1_y, const float* p1_z, int64 count,
                          float t) {

    for (int64 i = 0; i < count; i++) {
        const float x = p0_x[i] * (1.0f - t) + p1_x[i] * t;
        const float y = p0_y[i] * (1.0f - t) + p1_y[i] * t;
        const float z = p0_z[i] * (1.0f - t) + p1_z[i] * t;
        out_x[i] = x;
        out_y[i] = y;
        out_z[i] = z;
    }
}

void linear_interpolation_pbc(float* out_x, float* out_y, float* out_z, const float* in_x0, const float* in_y0, const float* in_z0, const float* in_x1, const float* in_y1, const float* in_z1,
                              int64 count, float t, const mat3& sim_box) {

    const float full_ext_x = sim_box[0][0];
    const float full_ext_y = sim_box[1][1];
    const float full_ext_z = sim_box[2][2];

    const float half_ext_x = full_ext_x * 0.5f;
    const float half_ext_y = full_ext_y * 0.5f;
    const float half_ext_z = full_ext_z * 0.5f;

    for (int64 i = 0; i < count; i++) {
        float x0 = in_x0[i];
        float y0 = in_y0[i];
        float z0 = in_z0[i];

        float x1 = in_x1[i];
        float y1 = in_y1[i];
        float z1 = in_z1[i];

        de_periodize(x1, x0, full_ext_x, half_ext_x);
        de_periodize(y1, y0, full_ext_y, half_ext_y);
        de_periodize(z1, z0, full_ext_z, half_ext_z);

        out_x[i] = x;
        out_y[i] = y;
        out_z[i] = z;
    }
}

/*
void linear_interpolation_pbc(float* out_x, float* out_y, float* out_z, const float* in0_x, const float* in0_y, const float* in0_z, const float* in1_x, const float* in1_y, const float* in1_z,
                              int64 count, float t, const mat3& sim_box) {

    const __m128 full_box_ext_x = _mm_set_ps1(sim_box[0][0]);
    const __m128 full_box_ext_y = _mm_set_ps1(sim_box[1][1]);
    const __m128 full_box_ext_z = _mm_set_ps1(sim_box[2][2]);

    const __m128 half_box_ext_x = _mm_mul_ps(full_box_ext_x, _mm_set_ps1(0.5f));
    const __m128 half_box_ext_y = _mm_mul_ps(full_box_ext_y, _mm_set_ps1(0.5f));
    const __m128 half_box_ext_z = _mm_mul_ps(full_box_ext_z, _mm_set_ps1(0.5f));

    const __m128 t4 = _mm_set_ps1(t);
    const __m128 one_minus_t4 = _mm_set_ps1(1.0f - t);

    for (int64 i = 0; i < count; i += 4) {
        __m128 x0 = _mm_load_ps(in0_x + i);
        __m128 y0 = _mm_load_ps(in0_y + i);
        __m128 z0 = _mm_load_ps(in0_z + i);
        __m128 x1 = _mm_load_ps(in1_x + i);
        __m128 y1 = _mm_load_ps(in1_y + i);
        __m128 z1 = _mm_load_ps(in1_z + i);

        __m128 dx = _mm_sub_ps(x1, x0);
        __m128 dy = _mm_sub_ps(y1, y0);
        __m128 dz = _mm_sub_ps(z1, z0);

        __m128 abs_dx = _mm_and_ps(dx, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF)));
        __m128 abs_dy = _mm_and_ps(dy, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF)));
        __m128 abs_dz = _mm_and_ps(dz, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF)));

        __m128 pm_x = _mm_cmpgt_ps(abs_dx, half_box_ext_x);
        __m128 pm_y = _mm_cmpgt_ps(abs_dy, half_box_ext_y);
        __m128 pm_z = _mm_cmpgt_ps(abs_dz, half_box_ext_z);

        __m128 dp_x = _mm_and_ps(pm_x, full_box_ext_x);
        __m128 dp_y = _mm_and_ps(pm_y, full_box_ext_y);
        __m128 dp_z = _mm_and_ps(pm_z, full_box_ext_z);

        _mm_store_ps(out_x, x0);
        _mm_store_ps(out_y, y0);
        _mm_store_ps(out_z, z0);
    }
}
*/

void cubic_interpolation(float* out_x, float* out_y, float* out_z, const float* in_x0, const float* in_y0, const float* in_z0, const float* in_x1, const float* in_y1, const float* in_z1,
                         const float* in_x2, const float* in_y2, const float* in_z2, const float* in_x3, const float* in_y3, const float* in_z3, int64 count, float t) {
    for (int i = 0; i < count; i += 4) {
        const __m128 x0 = _mm_load_ps(in_x0);
        const __m128 y0 = _mm_load_ps(in_y0);
        const __m128 z0 = _mm_load_ps(in_z0);

        const __m128 x1 = _mm_load_ps(in_x1);
        const __m128 y1 = _mm_load_ps(in_y1);
        const __m128 z1 = _mm_load_ps(in_z1);

        const __m128 x2 = _mm_load_ps(in_x2);
        const __m128 y2 = _mm_load_ps(in_y2);
        const __m128 z2 = _mm_load_ps(in_z2);

        const __m128 x3 = _mm_load_ps(in_x3);
        const __m128 y3 = _mm_load_ps(in_y3);
        const __m128 z3 = _mm_load_ps(in_z3);

        const __m128 x = simd::cubic_spline(x0, x1, x2, x3, t);
        const __m128 y = simd::cubic_spline(y0, y1, y2, y3, t);
        const __m128 z = simd::cubic_spline(z0, z1, z2, z3, t);

        simd::store(out_x, x);
        simd::store(out_y, y);
        simd::store(out_z, z);
    }
}

void cubic_interpolation_pbc(float* out_x, float* out_y, float* out_z, const float* in_x0, const float* in_y0, const float* in_z0, const float* in_x1, const float* in_y1, const float* in_z1,
                             const float* in_x2, const float* in_y2, const float* in_z2, const float* in_x3, const float* in_y3, const float* in_z3, int64 count, float t, const mat3& sim_box) {

    const __m128 full_box_ext_x = simd::set_float4(sim_box[0][0]);
    const __m128 full_box_ext_y = simd::set_float4(sim_box[1][1]);
    const __m128 full_box_ext_z = simd::set_float4(sim_box[2][2]);

    const __m128 half_box_ext_x = simd::mul(full_box_ext_x, simd::set_float4(0.5f));
    const __m128 half_box_ext_y = simd::mul(full_box_ext_y, simd::set_float4(0.5f));
    const __m128 half_box_ext_z = simd::mul(full_box_ext_z, simd::set_float4(0.5f));

    for (int i = 0; i < count; i += 4) {
        __m128 x0 = _mm_load_ps(in_x0);
        __m128 y0 = _mm_load_ps(in_y0);
        __m128 z0 = _mm_load_ps(in_z0);

        __m128 x1 = _mm_load_ps(in_x1);
        __m128 y1 = _mm_load_ps(in_y1);
        __m128 z1 = _mm_load_ps(in_z1);

        __m128 x2 = _mm_load_ps(in_x2);
        __m128 y2 = _mm_load_ps(in_y2);
        __m128 z2 = _mm_load_ps(in_z2);

        __m128 x3 = _mm_load_ps(in_x3);
        __m128 y3 = _mm_load_ps(in_y3);
        __m128 z3 = _mm_load_ps(in_z3);

        x0 = de_periodize(x1, x0, full_box_ext_x, half_box_ext_x);
        x2 = de_periodize(x1, x2, full_box_ext_x, half_box_ext_x);
        x3 = de_periodize(x1, x3, full_box_ext_X, half_box_ext_x);

        y0 = de_periodize(y1, y0, full_box_ext_y, half_box_ext_y);
        y2 = de_periodize(y1, y2, full_box_ext_y, half_box_ext_y);
        y3 = de_periodize(y1, y3, full_box_ext_y, half_box_ext_y);

        z0 = de_periodize(z1, z0, full_box_ext_z, half_box_ext_z);
        z2 = de_periodize(z1, z2, full_box_ext_z, half_box_ext_z);
        z3 = de_periodize(z1, z3, full_box_ext_z, half_box_ext_z);

        const __m128 x = simd::cubic_spline(x0, x1, x2, x3, t);
        const __m128 y = simd::cubic_spline(y0, y1, y2, y3, t);
        const __m128 z = simd::cubic_spline(z0, z1, z2, z3, t);

        simd::store(out_x, x);
        simd::store(out_y, y);
        simd::store(out_z, z);
    }
}

void compute_velocities_pbc(Array<vec3> dst_vel, Array<const vec3> pos, Array<const vec3> old_pos, const vec3& box_ext) {
    ASSERT(dst_vel.size() == pos.size());
    ASSERT(dst_vel.size() == old_pos.size());

    const float dt = 1.f;
    // ispc::compute_velocity(dst_vel.data(), pos.data(), old_pos.data(), dt, dst_vel.size());
    /*
for (int64 i = 0; i < dst_vel.size(); i++) {
    // De-periodize previous position
    const vec3 p1 = pos[i];
    const vec3 p0 = old_pos[i];

    const vec3 delta = p1 - p0;
    const vec3 signed_mask = sign(delta) * step(box_ext * 0.5f, abs(delta));
    const vec3 dp_p0 = p0 + box_ext * signed_mask;

    dst_vel[i] = p1 - dp_p0;
}
    */
}

/*
void apply_pbc_residues(Array<vec3> positions, Array<const Residue> residues, const mat3& sim_box) {
    const glm_vec4 full_box_ext = _mm_set_ps(0.f, sim_box[2][2], sim_box[1][1], sim_box[0][0]);
    for (const auto& res : residues) {
        const auto count = res.atom_idx.end - res.atom_idx.beg;
        const glm_vec4 scl = _mm_set_ps1(1.f / (float)count);
        glm_vec4 center = _mm_setzero_ps();
        for (const auto& pos : positions.subarray(res.atom_idx.beg, count)) {
            const glm_vec4 pos_vec = *reinterpret_cast<const glm_vec4*>(&pos);
            center = glm_vec4_add(center, pos_vec);
        }
        center = glm_vec4_mul(center, scl);
        const glm_vec4 pbc_cent = apply_pbc(center, full_box_ext);
        const glm_vec4 delta = glm_vec4_sub(pbc_cent, center);

        // @TODO: if delta is zero, skip this
        // if (!all_zero(delta)) {
        for (auto& pos : positions.subarray(res.atom_idx.beg, count)) {
            glm_vec4& pos_vec = *reinterpret_cast<glm_vec4*>(&pos);
            pos_vec = glm_vec4_add(pos_vec, delta);
        }
        //}
    }
}

void apply_pbc_chains(Array<vec3> positions, Array<const Chain> chains, Array<const Residue> residues, const mat3& sim_box) {
    const glm_vec4 full_box_ext = _mm_set_ps(0.f, sim_box[2][2], sim_box[1][1], sim_box[0][0]);
    for (const auto& chain : chains) {
        const auto beg_idx = residues[chain.res_idx.beg].atom_idx.beg;
        const auto end_idx = residues[chain.res_idx.end - 1].atom_idx.end;
        const auto count = end_idx - beg_idx;
        const glm_vec4 scl = _mm_set_ps1(1.f / (float)count);
        glm_vec4 center = _mm_setzero_ps();
        for (const auto& pos : positions.subarray(beg_idx, count)) {
            const glm_vec4 pos_vec = *reinterpret_cast<const glm_vec4*>(&pos);
            center = glm_vec4_add(center, pos_vec);
        }
        center = glm_vec4_mul(center, scl);
        const glm_vec4 pbc_cent = apply_pbc(center, full_box_ext);
        const glm_vec4 delta = glm_vec4_sub(pbc_cent, center);

        // if (!all_zero(delta)) {
        for (auto& pos : positions.subarray(beg_idx, count)) {
            const glm_vec4 pos_vec = *reinterpret_cast<glm_vec4*>(&pos);
            const glm_vec4 res = glm_vec4_add(pos_vec, delta);
            pos = *reinterpret_cast<const vec3*>(&res);
        }
        //}
    }
}
*/

inline bool covelent_bond_heuristic(float x0, float y0, float z0, Element e0, float x1, float y1, float z1, Element e1) {
    const float d = element::covalent_radius(e0) + element::covalent_radius(e1);
    const float d_max = d + 0.3f;
    const float d_min = d - 0.5f;
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const float dz = z1 - z0;
    const float d2 = sqrtf(dx * dx + dy * dy + dz * dx);
    return (d_min * d_min) < d2 && d2 < (d_max * d_max);
}

bool has_covalent_bond(const Residue& res_a, const Residue& res_b) { return (res_a.bond_idx.beg < res_b.bond_idx.end && res_b.bond_idx.beg < res_a.bond_idx.end); }

bool valid_segment(const BackboneSegment& segment) { return segment.ca_idx != -1 && segment.c_idx != -1 && segment.n_idx != -1 && segment.o_idx != -1; }

// Computes covalent bonds between a set of atoms with given positions and elements.
// The approach is inspired by the technique used in NGL (https://github.com/arose/ngl)
DynamicArray<Bond> compute_covalent_bonds(Array<Residue> residues, const float* pos_x, const float* pos_y, const float* pos_z, const ResIdx* res_idx, const Element* element, int64 count) {
    if (residues.count == 0) {
        LOG_WARNING("Cannot compute covalent bonds, no residues were given.");
        return {};
    }

    constexpr float max_covelent_bond_length = 4.0f;
    spatialhash::Frame frame = spatialhash::compute_frame(pos_x, pos_y, pos_z, count, vec3(max_covelent_bond_length));
    DynamicArray<Bond> bonds;

    // @NOTE: The assumtion is that a bond is either within a single residue or between concecutive residues.

    for (int32 i = 0; i < residues.size(); i++) {
        auto& res = residues[i];

        if (i > 0) {
            // Include potential shared bonds from previous residue
            res.bond_idx.beg = residues[i - 1].bond_idx.end_internal;
            res.bond_idx.beg_internal = res.bond_idx.end_internal = res.bond_idx.end = residues[i - 1].bond_idx.end;
        } else {
            res.bond_idx.beg = res.bond_idx.end = (BondIdx)bonds.size();
            res.bond_idx.beg_internal = res.bond_idx.end_internal = res.bond_idx.end = (BondIdx)bonds.size();
        }

        // Internal bonds
        for (AtomIdx i = res.atom_idx.beg; i < res.atom_idx.end; i++) {
            const vec3& pos_xyz = {pos_x[i], pos_y[i], pos_z[i]};
            spatialhash::for_each_within(frame, pos_xyz, max_covelent_bond_length, [&bonds, &res, pos_x, pos_y, pos_z, res_idx, element, i](int j, const vec3& atom_j_pos) {
                (void)atom_j_pos;
                const bool has_bond = covelent_bond_heuristic(pos_x[i], pos_y[i], pos_z[i], element[i], pos_x[j], pos_y[j], pos_z[j], element[j]);
                if (i < j && res_idx[i] == res_idx[j] && has_bond) {
                    bonds.push_back({{i, j}});
                    res.bond_idx.end++;
                }
            });
        }
        res.bond_idx.end_internal = res.bond_idx.end;

        // Now locate external bonds to next residue
        for (AtomIdx i = res.atom_idx.beg; i < res.atom_idx.end; i++) {
            const vec3& pos_xyz = {pos_x[i], pos_y[i], pos_z[i]};
            spatialhash::for_each_within(frame, pos_xyz, max_covelent_bond_length, [&bonds, &res, pos_x, pos_y, pos_z, res_idx, element, i](int j, const vec3& atom_j_pos) {
                (void)atom_j_pos;
                const bool has_bond = covelent_bond_heuristic(pos_x[i], pos_y[i], pos_z[i], element[i], pos_x[j], pos_y[j], pos_z[j], element[j]);
                const bool consecutive_res = math::abs(res_idx[i] - res_idx[j]) == 1;
                if (i < j && consecutive_res && has_bond) {
                    bonds.push_back({{i, j}});
                    res.bond_idx.end++;
                }
            });
        }
    }

    /*
        // Old approach which does not give internal then external bonds for residues
for (int atom_i = 0; atom_i < atom_pos.count; atom_i++) {
spatialhash::for_each_within(
    frame, atom_pos[atom_i], max_covelent_bond_length,
    [&bonds, &atom_pos, &atom_elem, &atom_res_idx, atom_i](int atom_j, const vec3& atom_j_pos) {
        (void)atom_j_pos;
        if (atom_i < atom_j &&
            (math::abs(atom_res_idx[atom_i] - atom_res_idx[atom_j]) <
                2) &&  // only create bonds where i < j and res_idx is concecutive (abs(res_idx[i] - res_idx[j]) < 2)
            covelent_bond_heuristic(atom_pos[atom_i], atom_elem[atom_i], atom_pos[atom_j], atom_elem[atom_j])) {
            bonds.push_back({atom_i, atom_j});
        }
    });
}
    */

    return bonds;
}

DynamicArray<Bond> compute_covalent_bonds(const float* pos_x, const float* pos_y, const float* pos_z, const ResIdx* res_idx, const Element* element, int64 count) {
    constexpr float max_covelent_bond_length = 4.0f;
    spatialhash::Frame frame = spatialhash::compute_frame(pos_x, pos_y, pos_z, count, vec3(max_covelent_bond_length));
    DynamicArray<Bond> bonds;

    for (int i = 0; i < count; i++) {
        const vec3& pos_xyz = {pos_x[i], pos_y[i], pos_z[i]};
        spatialhash::for_each_within(frame, pos_xyz, max_covelent_bond_length, [&bonds, pos_x, pos_y, pos_z, element, i](int j, const vec3& atom_j_pos) {
            (void)atom_j_pos;
            bool has_bond = covelent_bond_heuristic(pos_x[i], pos_y[i], pos_z[i], element[i], pos_x[j], pos_y[j], pos_z[j], element[j]);
            if (i < j && has_bond) {
                bonds.push_back({{i, j}});
            }
        });
    }

    return bonds;
}

// @NOTE this method is sub-optimal and can surely be improved...
// Residues should have no more than 2 potential connections to other residues.
DynamicArray<Chain> compute_chains(Array<const Residue> residues) {

    DynamicArray<Bond> residue_bonds;
    for (ResIdx i = 0; i < (ResIdx)residues.size() - 1; i++) {
        if (has_covalent_bond(residues[i], residues[i + 1])) {
            residue_bonds.push_back({{i, i + 1}});
        }
    }

    if (residue_bonds.count == 0) {
        // No residue bonds, No chains.
        return {};
    }

    DynamicArray<int> residue_chains(residues.count, -1);

    if (residue_bonds.count > 0) {
        int curr_chain_idx = 0;
        int res_bond_idx = 0;
        for (int i = 0; i < residues.count; i++) {
            if (residue_chains[i] == -1) residue_chains[i] = curr_chain_idx++;
            for (; res_bond_idx < residue_bonds.count; res_bond_idx++) {
                const auto& res_bond = residue_bonds[res_bond_idx];
                if (i == res_bond.idx[0]) {
                    residue_chains[res_bond.idx[1]] = residue_chains[res_bond.idx[0]];
                } else if (res_bond.idx[0] > i)
                    break;
            }
        }
    }

    DynamicArray<Chain> chains;
    int curr_chain_idx = -1;
    for (int i = 0; i < residue_chains.size(); i++) {
        if (residue_chains[i] != curr_chain_idx) {
            curr_chain_idx = residue_chains[i];
            Label lbl;
            snprintf(lbl.cstr(), lbl.capacity(), "C%i", curr_chain_idx);
            chains.push_back({lbl, {(ResIdx)i, (ResIdx)i}, {}});
        }
        if (chains.size() > 0) {
            chains.back().res_idx.end++;
        }
    }

    for (auto& c : chains) {
        c.atom_idx.beg = residues[c.res_idx.beg].atom_idx.beg;
        c.atom_idx.end = residues[c.res_idx.end - 1].atom_idx.end;
    }

    return chains;
}

DynamicArray<BackboneSequence> compute_backbone_sequences(Array<const BackboneSegment> segments, Array<const Residue> residues) {
    if (segments.count == 0) return {};
    ASSERT(segments.count == residues.count);

    DynamicArray<BackboneSequence> bb_sequences;
    for (ResIdx i = 0; i < (ResIdx)residues.size(); i++) {
        if (valid_segment(segments[i])) {
            bb_sequences.push_back({i, i + 1});
            while (i < (ResIdx)residues.size() - 1 && valid_segment(segments[i + 1]) && has_covalent_bond(residues[i], residues[i + 1])) {
                bb_sequences.back().end++;
                i++;
            }
        }
    }

    return bb_sequences;
}

template <int64 N>
bool match(const Label& lbl, const char (&cstr)[N]) {
    for (int64 i = 0; i < N; i++) {
        if (tolower(lbl[i]) != tolower(cstr[i])) return false;
    }
    return true;
}

DynamicArray<BackboneSegment> compute_backbone_segments(Array<const Residue> residues, Array<const Label> atom_labels) {
    DynamicArray<BackboneSegment> segments;
    int64 invalid_segments = 0;
    constexpr int32 min_atom_count = 4;  // Must contain at least 8 atoms to be considered as an amino acid.
    for (auto& res : residues) {
        const int32 atom_count = res.atom_idx.end - res.atom_idx.beg;
        if (atom_count < min_atom_count) {
            segments.push_back({-1, -1, -1, -1});
            invalid_segments++;
            continue;
        }

        BackboneSegment seg{};
        // if (is_amino_acid(res)) {
        // find atoms
        for (int32 i = res.atom_idx.beg; i < res.atom_idx.end; i++) {
            const auto& lbl = atom_labels[i];
            if (seg.ca_idx == -1 && match(lbl, "CA")) seg.ca_idx = i;
            if (seg.n_idx == -1 && match(lbl, "N")) seg.n_idx = i;
            if (seg.c_idx == -1 && match(lbl, "C")) seg.c_idx = i;
            if (seg.o_idx == -1 && match(lbl, "O")) seg.o_idx = i;
        }

        // Could not match "O"
        if (seg.o_idx == -1) {
            // Pick first atom containing O after C atom
            for (int32 i = seg.c_idx; i < res.atom_idx.end; i++) {
                const auto& lbl = atom_labels[i];
                if (lbl[0] == 'o' || lbl[0] == 'O') seg.o_idx = i;
            }
        }

        if (!valid_segment(seg)) {
            // LOG_ERROR("Could not identify all backbone indices for residue %s.", res.name.beg());
            invalid_segments++;
        }
        //} else {
        //    invalid_segments++;
        //}
        segments.push_back(seg);
    }

    if (invalid_segments == segments.count) return {};

    return segments;
}

/*
DynamicArray<SplineSegment> compute_spline(Array<const vec3> atom_pos, Array<const uint32> colors, Array<const BackboneSegment> backbone, int32 num_subdivisions, float tension) {
    if (backbone.count < 4) return {};

    DynamicArray<vec3> p_tmp;
    DynamicArray<vec3> o_tmp;
    DynamicArray<vec3> c_tmp;
    DynamicArray<int> ca_idx;

    // Pad front with duplicated vectors
    auto d_p0 = atom_pos[backbone[1].ca_idx] - atom_pos[backbone[0].ca_idx];
    p_tmp.push_back(atom_pos[backbone[0].ca_idx] - d_p0);

    auto d_o0 = atom_pos[backbone[1].o_idx] - atom_pos[backbone[0].o_idx];
    o_tmp.push_back(atom_pos[backbone[0].o_idx] - d_o0);

    auto d_c0 = atom_pos[backbone[1].c_idx] - atom_pos[backbone[0].c_idx];
    c_tmp.push_back(atom_pos[backbone[0].c_idx] - d_c0);

    ca_idx.push_back(backbone[0].ca_idx);

    // Fetch vectors
    const int size = (int)(backbone.count);
    for (auto i = 0; i < size; i++) {
        p_tmp.push_back(atom_pos[backbone[i].ca_idx]);
        o_tmp.push_back(atom_pos[backbone[i].o_idx]);
        c_tmp.push_back(atom_pos[backbone[i].c_idx]);
        ca_idx.push_back(backbone[i].ca_idx);
    }

    // Pad back with duplicated vectors
    auto d_pn = atom_pos[backbone[size - 1].ca_idx] - atom_pos[backbone[size - 2].ca_idx];
    p_tmp.push_back(atom_pos[backbone[size - 1].ca_idx] + d_pn);
    p_tmp.push_back(p_tmp.back() + d_pn);

    auto d_on = atom_pos[backbone[size - 1].o_idx] - atom_pos[backbone[size - 2].o_idx];
    o_tmp.push_back(atom_pos[backbone[size - 1].o_idx] + d_on);
    o_tmp.push_back(o_tmp.back() + d_on);

    auto d_cn = atom_pos[backbone[size - 1].c_idx] - atom_pos[backbone[size - 2].c_idx];
    c_tmp.push_back(atom_pos[backbone[size - 1].c_idx] + d_cn);
    c_tmp.push_back(c_tmp.back() + d_cn);

    ca_idx.push_back(backbone[size - 1].ca_idx);
    ca_idx.push_back(backbone[size - 1].ca_idx);

    // @NOTE: Flip direction of support vector (O <- C) if pointing the 'wrong way' as seen from previous segment
    for (int64 i = 1; i < o_tmp.size(); i++) {
        vec3 v0 = o_tmp[i - 1] - c_tmp[i - 1];
        vec3 v1 = o_tmp[i] - c_tmp[i];

        if (glm::dot(v0, v1) < 0) {
            o_tmp[i] = c_tmp[i] - v1;
        }
    }

    DynamicArray<SplineSegment> segments;

    for (int64 i = 1; i < p_tmp.size() - 2; i++) {
        auto p0 = p_tmp[i - 1];
        auto p1 = p_tmp[i];
        auto p2 = p_tmp[i + 1];
        auto p3 = p_tmp[i + 2];

        auto o0 = o_tmp[i - 1];
        auto o1 = o_tmp[i];
        auto o2 = o_tmp[i + 1];
        auto o3 = o_tmp[i + 2];

        auto c0 = c_tmp[i - 1];
        auto c1 = c_tmp[i];
        auto c2 = c_tmp[i + 1];
        auto c3 = c_tmp[i + 2];

        uint32 idx = ca_idx[i];
        uint32 color = colors[idx];

        auto count = (i < (p_tmp.size() - 3)) ? num_subdivisions : num_subdivisions + 1;
        for (int n = 0; n < count; n++) {
            auto t = n / (float)(num_subdivisions);

            vec3 p = math::spline(p0, p1, p2, p3, t, tension);
            vec3 o = math::spline(o0, o1, o2, o3, t, tension);
            vec3 c = math::spline(c0, c1, c2, c3, t, tension);

            vec3 v_dir = math::normalize(o - c);

            vec3 tangent = math::spline_tangent(p0, p1, p2, p3, t, tension);
            vec3 normal = math::normalize(math::cross(v_dir, tangent));
            vec3 binormal = math::normalize(math::cross(tangent, normal));

            segments.push_back({p, tangent, normal, binormal, idx, color});
        }
    }

    return segments;
}
*/

DynamicArray<BackboneAngle> compute_backbone_angles(Array<const BackboneSegment> backbone, const float* pos_x, const float* pos_y, const float* pos_z, int64 count) {
    if (backbone.count == 0) return {};
    DynamicArray<BackboneAngle> angles(backbone.count);
    compute_backbone_angles(angles, backbone, pos_x, pos_y, pos_z, count);
    return angles;
}

void compute_backbone_angles(Array<BackboneAngle> dst, Array<const BackboneSegment> backbone_segments, const float* pos_x, const float* pos_y, const float* pos_z, int64 count) {
    ASSERT(dst.count >= backbone_segments.count);
    float phi, psi;

    ASSERT(valid_segment(backbone_segments[0]));
    vec3 n = {pos_x[backbone_segments[0].n_idx], pos_y[backbone_segments[0].n_idx], pos_z[backbone_segments[0].n_idx]};
    vec3 ca = {pos_x[backbone_segments[0].ca_idx], pos_y[backbone_segments[0].ca_idx], pos_z[backbone_segments[0].ca_idx]};
    vec3 c = {pos_x[backbone_segments[0].c_idx], pos_y[backbone_segments[0].c_idx], pos_z[backbone_segments[0].c_idx]};

    vec3 c_prev = c;
    vec3 n_next = {pos_x[backbone_segments[1].n_idx], pos_y[backbone_segments[1].n_idx], pos_z[backbone_segments[1].n_idx]};
    phi = 0.0f;
    psi = math::dihedral_angle(n, ca, c, n_next);
    dst[0] = {phi, psi};

    for (int64 i = 1; i < backbone_segments.count - 1; i++) {
        ASSERT(valid_segment(backbone_segments[i]));

        c_prev = c;
        n = n_next;
        ca = {pos_x[backbone_segments[i].ca_idx], pos_y[backbone_segments[i].ca_idx], pos_z[backbone_segments[i].ca_idx]};
        c = {pos_x[backbone_segments[i].c_idx], pos_y[backbone_segments[i].c_idx], pos_z[backbone_segments[i].c_idx]};
        n_next = {pos_x[backbone_segments[i + 1].n_idx], pos_y[backbone_segments[i + 1].n_idx], pos_z[backbone_segments[i + 1].n_idx]};

        phi = math::dihedral_angle(c_prev, n, ca, c);
        psi = math::dihedral_angle(n, ca, c, n_next);
        dst[i] = {phi, psi};
    }

    auto N = backbone_segments.count - 1;
    ASSERT(valid_segment(backbone_segments[N]));

    c_prev = c;
    n = n_next;
    ca = {pos_x[backbone_segments[N].ca_idx], pos_y[backbone_segments[N].ca_idx], pos_z[backbone_segments[N].ca_idx]};
    c = {pos_x[backbone_segments[N].c_idx], pos_y[backbone_segments[N].c_idx], pos_z[backbone_segments[N].c_idx]};

    phi = math::dihedral_angle(c_prev, n, ca, c);
    psi = 0.0f;
    dst[N] = {phi, psi};
}

DynamicArray<BackboneAngle> compute_backbone_angles(Array<const BackboneSegment> segments, Array<const BackboneSequence> sequences, const float* pos_x, const float* pos_y, const float* pos_z,
                                                    int64 count) {
    if (segments.size() == 0) return {};
    DynamicArray<BackboneAngle> angles(segments.count);
    compute_backbone_angles(angles, segments, sequences, pos_x, pos_y, pos_z, count);
    return angles;
}

void compute_backbone_angles(Array<BackboneAngle> dst, Array<const BackboneSegment> segments, Array<const BackboneSequence> sequences, const float* pos_x, const float* pos_y, const float* pos_z,
                             int64 count) {
    for (const auto& seq : sequences) {
        compute_backbone_angles(dst.subarray(seq.beg, seq.end - seq.beg), segments.subarray(seq.beg, seq.end - seq.beg), pos_x, pos_y, pos_z, count);
    }
}

void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(data);
    if (!dynamic.molecule || !dynamic.trajectory) return;

    if (data->angle_data) {
        FREE(data->angle_data.ptr);
    }

    int32 alloc_count = (int32)dynamic.molecule.backbone.segments.count * (int32)dynamic.trajectory.frame_buffer.count;
    data->num_segments = (int32)dynamic.molecule.backbone.segments.count;
    data->num_frames = 0;
    data->angle_data = {(vec2*)CALLOC(alloc_count, sizeof(vec2)), alloc_count};
}

void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data) {
    ASSERT(data);
    if (data->angle_data) {
        FREE(data->angle_data.ptr);
        *data = {};
    }
}

void compute_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(dynamic);
    if (dynamic.trajectory.num_frames == 0 || dynamic.molecule.backbone.segments.count == 0) return;

    //@NOTE: Trajectory may be loading while this is taking place, therefore read num_frames once and stick to that
    const int32 traj_num_frames = dynamic.trajectory.num_frames;

    // @NOTE: If we are up to date, no need to compute anything
    if (traj_num_frames == data->num_frames) {
        return;
    }

    // @TODO: parallelize?
    // @NOTE: Only compute data for indices which are new
    for (int32 f_idx = data->num_frames; f_idx < traj_num_frames; f_idx++) {
        auto pos_x = get_trajectory_position_x(dynamic.trajectory, f_idx);
        auto pos_y = get_trajectory_position_y(dynamic.trajectory, f_idx);
        auto pos_z = get_trajectory_position_z(dynamic.trajectory, f_idx);

        Array<vec2> frame_angles = get_backbone_angles(*data, f_idx);
        for (const auto& bb_seq : dynamic.molecule.backbone.sequences) {
            auto bb_segments = get_backbone(dynamic.molecule, bb_seq);
            auto bb_angles = frame_angles.subarray(bb_seq.beg, bb_seq.end - bb_seq.beg);

            if (bb_segments.size() < 2) {
                memset(bb_angles.ptr, 0, bb_angles.size_in_bytes());
            } else {
                compute_backbone_angles(bb_angles, bb_segments, pos_x.data(), pos_y.data(), pos_z.data(), dynamic.trajectory.num_atoms);
            }
        }
    }
    data->num_frames = traj_num_frames;  // update current count
}

DynamicArray<float> compute_atom_radii(Array<const Element> elements) {
    DynamicArray<float> radii(elements.size(), 0);
    compute_atom_radii(radii.data(), elements.data(), radii.size());
    return radii;
}

void compute_atom_radii(float* out_radius, const Element* element, int64 count) {
    for (int64 i = 0; i < count; i++) {
        out_radius[i] = element::vdw_radius(element[i]);
    }
}

bool is_amino_acid(const Residue& res) { return aminoacid::get_from_string(res.name) != AminoAcid::Unknown; }

bool is_dna(const Residue& res) {
    constexpr const char* dna_residues[12] = {"DA", "DA3", "DA5", "DC", "DC3", "DC5", "DG", "DG3", "DG5", "DT", "DT3", "DT5"};
    for (auto dna_res : dna_residues) {
        if (compare(res.name, dna_res)) return true;
    }
    return false;
}
