#include "structure_tracking.h"

#include <core/log.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#pragma warning( disable : 4127 ) // disable warnings about expressions which could be constexpr in Eigen
#include <Eigen/Eigen>

namespace structure_tracking {

struct Structure {
	ID id = 0;
	int32 ref_frame_idx = 0;
	int32 num_points = 0;
	int32 num_frames = 0;

	struct {
		float* x;
		float* y;
		float* z;
		float* mass;
	} reference;

	Transform* transform = nullptr;
};

struct Entry {
	ID id;
	Structure* ptr;
};

struct Context {
	DynamicArray<Entry> entries{};
};

Context* context = nullptr;

// from here https://stackoverflow.com/questions/4372224/fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
static void Diagonalize(const float(&A)[3][3], float(&Q)[3][3], float(&D)[3][3]) {
	// A must be a symmetric matrix.
	// returns Q and D such that
	// Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
	const int maxsteps = 24;  // certainly wont need that many.
	int k0, k1, k2;
	float o[3], m[3];
	float q[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float jr[4];
	float sqw, sqx, sqy, sqz;
	float tmp1, tmp2, mq;
	float AQ[3][3];
	float thet, sgn, t, c;
	for (int i = 0; i < maxsteps; ++i) {
		// quat to matrix
		sqx = q[0] * q[0];
		sqy = q[1] * q[1];
		sqz = q[2] * q[2];
		sqw = q[3] * q[3];
		Q[0][0] = (sqx - sqy - sqz + sqw);
		Q[1][1] = (-sqx + sqy - sqz + sqw);
		Q[2][2] = (-sqx - sqy + sqz + sqw);
		tmp1 = q[0] * q[1];
		tmp2 = q[2] * q[3];
		Q[1][0] = 2.0f * (tmp1 + tmp2);
		Q[0][1] = 2.0f * (tmp1 - tmp2);
		tmp1 = q[0] * q[2];
		tmp2 = q[1] * q[3];
		Q[2][0] = 2.0f * (tmp1 - tmp2);
		Q[0][2] = 2.0f * (tmp1 + tmp2);
		tmp1 = q[1] * q[2];
		tmp2 = q[0] * q[3];
		Q[2][1] = 2.0f * (tmp1 + tmp2);
		Q[1][2] = 2.0f * (tmp1 - tmp2);

		// AQ = A * Q
		AQ[0][0] = Q[0][0] * A[0][0] + Q[1][0] * A[0][1] + Q[2][0] * A[0][2];
		AQ[0][1] = Q[0][1] * A[0][0] + Q[1][1] * A[0][1] + Q[2][1] * A[0][2];
		AQ[0][2] = Q[0][2] * A[0][0] + Q[1][2] * A[0][1] + Q[2][2] * A[0][2];
		AQ[1][0] = Q[0][0] * A[0][1] + Q[1][0] * A[1][1] + Q[2][0] * A[1][2];
		AQ[1][1] = Q[0][1] * A[0][1] + Q[1][1] * A[1][1] + Q[2][1] * A[1][2];
		AQ[1][2] = Q[0][2] * A[0][1] + Q[1][2] * A[1][1] + Q[2][2] * A[1][2];
		AQ[2][0] = Q[0][0] * A[0][2] + Q[1][0] * A[1][2] + Q[2][0] * A[2][2];
		AQ[2][1] = Q[0][1] * A[0][2] + Q[1][1] * A[1][2] + Q[2][1] * A[2][2];
		AQ[2][2] = Q[0][2] * A[0][2] + Q[1][2] * A[1][2] + Q[2][2] * A[2][2];
		// D = Qt * AQ
		D[0][0] = AQ[0][0] * Q[0][0] + AQ[1][0] * Q[1][0] + AQ[2][0] * Q[2][0];
		D[0][1] = AQ[0][0] * Q[0][1] + AQ[1][0] * Q[1][1] + AQ[2][0] * Q[2][1];
		D[0][2] = AQ[0][0] * Q[0][2] + AQ[1][0] * Q[1][2] + AQ[2][0] * Q[2][2];
		D[1][0] = AQ[0][1] * Q[0][0] + AQ[1][1] * Q[1][0] + AQ[2][1] * Q[2][0];
		D[1][1] = AQ[0][1] * Q[0][1] + AQ[1][1] * Q[1][1] + AQ[2][1] * Q[2][1];
		D[1][2] = AQ[0][1] * Q[0][2] + AQ[1][1] * Q[1][2] + AQ[2][1] * Q[2][2];
		D[2][0] = AQ[0][2] * Q[0][0] + AQ[1][2] * Q[1][0] + AQ[2][2] * Q[2][0];
		D[2][1] = AQ[0][2] * Q[0][1] + AQ[1][2] * Q[1][1] + AQ[2][2] * Q[2][1];
		D[2][2] = AQ[0][2] * Q[0][2] + AQ[1][2] * Q[1][2] + AQ[2][2] * Q[2][2];
		o[0] = D[1][2];
		o[1] = D[0][2];
		o[2] = D[0][1];
		m[0] = fabs(o[0]);
		m[1] = fabs(o[1]);
		m[2] = fabs(o[2]);

		k0 = (m[0] > m[1] && m[0] > m[2]) ? 0 : (m[1] > m[2]) ? 1 : 2;  // index of largest element of offdiag
		k1 = (k0 + 1) % 3;
		k2 = (k0 + 2) % 3;
		if (o[k0] == 0.0f) {
			break;  // diagonal already
		}
		thet = (D[k2][k2] - D[k1][k1]) / (2.0f * o[k0]);
		sgn = (thet > 0.0f) ? 1.0f : -1.0f;
		thet *= sgn;                                                             // make it positive
		t = sgn / (thet + ((thet < 1.E6f) ? sqrtf(thet * thet + 1.0f) : thet));  // sign(T)/(|T|+sqrt(T^2+1))
		c = 1.0f / sqrtf(t * t + 1.0f);                                          //  c= 1/(t^2+1) , t=s/c
		if (c == 1.0f) {
			break;  // no room for improvement - reached machine precision.
		}
		jr[0] = jr[1] = jr[2] = jr[3] = 0.0f;
		jr[k0] = sgn * sqrtf((1.0f - c) / 2.0f);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
		jr[k0] *= -1.0f;                          // since our quat-to-matrix convention was for v*M instead of M*v
		jr[3] = sqrtf(1.0f - jr[k0] * jr[k0]);
		if (jr[3] == 1.0f) {
			break;  // reached limits of floating point precision
		}
		q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
		q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
		q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
		q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
		mq = sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		q[0] /= mq;
		q[1] /= mq;
		q[2] /= mq;
		q[3] /= mq;
	}
}

static void diagonalize(const mat3& M, mat3* Q, mat3* D) {
	ASSERT(Q);
	ASSERT(D);
	Diagonalize((const float(&)[3][3])M, (float(&)[3][3]) * Q, (float(&)[3][3]) * D);
}

static void decompose(const mat3& M, mat3* R, mat3* S) {
	ASSERT(R);
	ASSERT(S);
	mat3 AtA = math::transpose(M) * M;
	mat3 Q, D;
	diagonalize(AtA, &Q, &D);
	const float det = math::determinant(AtA); // For debugging
	D[0][0] = sqrtf(D[0][0]);
	D[1][1] = sqrtf(D[1][1]);
	D[2][2] = sqrtf(D[2][2]);
	*S = math::inverse(Q) * D * Q;
	*R = M * math::inverse(*S);
}

static mat3 compute_linear_transform(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
								     const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
									 const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1) {
	mat3 Apq{ 0 };
	mat3 Aqq{ 0 };

	for (int64 i = 0; i < count; i++) {
		// @TODO: Vectorize...
		const float q_x = x0[i] - com0.x;
		const float q_y = y0[i] - com0.y;
		const float q_z = z0[i] - com0.z;

		const float p_x = x1[i] - com1.x;
		const float p_y = y1[i] - com1.y;
		const float p_z = z1[i] - com1.z;

		Apq[0][0] += mass[i] * p_x * q_x;
		Apq[0][1] += mass[i] * p_y * q_x;
		Apq[0][2] += mass[i] * p_z * q_x;
		Apq[1][0] += mass[i] * p_x * q_y;
		Apq[1][1] += mass[i] * p_y * q_y;
		Apq[1][2] += mass[i] * p_z * q_y;
		Apq[2][0] += mass[i] * p_x * q_z;
		Apq[2][1] += mass[i] * p_y * q_z;
		Apq[2][2] += mass[i] * p_z * q_z;

		Aqq[0][0] += mass[i] * q_x * q_x;
		Aqq[0][1] += mass[i] * q_y * q_x;
		Aqq[0][2] += mass[i] * q_z * q_x;
		Aqq[1][0] += mass[i] * q_x * q_y;
		Aqq[1][1] += mass[i] * q_y * q_y;
		Aqq[1][2] += mass[i] * q_z * q_y;
		Aqq[2][0] += mass[i] * q_x * q_z;
		Aqq[2][1] += mass[i] * q_y * q_z;
		Aqq[2][2] += mass[i] * q_z * q_z;
	}

	return Apq / Aqq;
}

static void compute_residual_error(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								   const float* RESTRICT src_x, const float* RESTRICT src_y, const float* RESTRICT src_z,
								   const float* RESTRICT ref_x, const float* RESTRICT ref_y, const float* RESTRICT ref_z,
								   int64 count, const mat4& matrix) {
	for (int32 i = 0; i < count; i++) {
		// @TODO: Vectorize this...

		const vec4 r = { ref_x[i], ref_y[i], ref_z[i], 1.0f };
		const vec4 v = { src_x[i], src_y[i], src_z[i], 1.0f };
		const vec4 u = matrix * v;
		const vec4 d = u - r;

		out_x[i] = d.x;
		out_y[i] = d.y;
		out_z[i] = d.z;
	}
}

// RBF functions
inline float Wendland_3_1(float r) {
	const float x = 1.f - r;
	const float x2 = x * x;
	return x2 * x2 * (4.f * r + 1.f);
}

inline float Gaussian(float r) {
	const float a = 0.5f * r;
	return math::exp(-a*a);
}

#define RBF_FUNCTION Wendland_3_1

static void compute_rbf_weights(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
						 const float* RESTRICT in_pos_x, const float* RESTRICT in_pos_y, const float* RESTRICT in_pos_z,
						 const float* RESTRICT in_val_x, const float* RESTRICT in_val_y, const float* RESTRICT in_val_z,
						 int64 count, const float radial_cutoff = 10.f) {

	const float d2_max = radial_cutoff * radial_cutoff;
	const float scl = 1.f / radial_cutoff;
	const int32 N = (int32)count;

	// Create matrix A NxN
	Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N, N);
	for (int32 i = 0; i < N; i++) {
		for (int32 j = 0; j < N; j++) {
			const float dx = in_pos_x[i] - in_pos_x[j];
			const float dy = in_pos_y[i] - in_pos_y[j];
			const float dz = in_pos_z[i] - in_pos_z[j];
			const float d2 = dx * dx + dy * dy + dz * dz;
			if (d2 > d2_max) continue;

			const float d = math::sqrt(d2);
			A(i, j) = RBF_FUNCTION(d * scl);
		}
	}

	// Create Vector b
	Eigen::MatrixXf b = Eigen::MatrixXf::Zero(N, 3);
	for (int32 i = 0; i < N; i++) {
		b.row(i) = Eigen::Vector3f(in_val_x[i], in_val_y[i], in_val_z[i]);
	}

	Eigen::MatrixXf x = (A.transpose() * A).ldlt().solve(A.transpose() * b);

	for (int32 i = 0; i < N; i++) {
		out_x[i] = x(i, 0);
		out_y[i] = x(i, 1);
		out_z[i] = x(i, 2);
	};
}

static void free_structure_data(Structure* s) {
	ASSERT(s);
	if (s->reference.x) FREE(s->reference.x);
	if (s->transform) FREE(s->transform);
	/*
	if (s->data.rbf.weight.x) FREE(s->data.rbf.weight.x);
	if (s->data.rbf.pos.x) FREE(s->data.rbf.pos.x);
	*/
}

static void init_structure_data(Structure* s, ID id, int32 ref_frame_idx, int32 num_points, int32 num_frames) {
	ASSERT(s);
	free_structure_data(s);
	s->id = id;
	s->ref_frame_idx = ref_frame_idx;
	s->num_frames = num_frames;
	s->num_points = num_points;

	s->reference.x = (float*)MALLOC(sizeof(float) * num_points * 4);
	s->reference.y = s->reference.x + num_points;
	s->reference.z = s->reference.y + num_points;
	s->reference.mass = s->reference.z + num_points;
	s->transform = (Transform*)MALLOC(sizeof(Transform) * num_frames);
/*
	s->data.rbf.radial_cutoff = radial_cutoff;
	s->data.rbf.weight.x = (float*)MALLOC(sizeof(float) * num_points * num_frames * 3);
	s->data.rbf.weight.y = s->data.rbf.weight.x + num_points * num_frames;
	s->data.rbf.weight.z = s->data.rbf.weight.y + num_points * num_frames;
	s->data.rbf.pos.x = (float*)MALLOC(sizeof(float) * num_points * 3);
	s->data.rbf.pos.y = s->data.rbf.pos.x + num_points;
	s->data.rbf.pos.z = s->data.rbf.pos.y + num_points;
	*/
}

static Structure* find_structure(ID id) {
	for (const auto& e : context->entries) {
		if (e.id == id) return e.ptr;
	}
	return nullptr;
}

void initialize() {
	if (!context) {
		context = NEW(Context);
		context->entries.reserve(8);
	}
}

void shutdown() {
	if (context) {
		clear_structures();
		context->~Context();
		FREE(context);
	}
}

bool create_structure(ID id) {
	ASSERT(context);
	ASSERT(id != 0);

	if (find_structure(id) != NULL) {
		LOG_ERROR("Could not create structure: it already exists!");
		return false;
	}
	Structure* s = new (MALLOC(sizeof(Structure))) Structure();
	context->entries.push_back({id, s});
	return true;
}

bool remove_structure(ID id) {
	ASSERT(context);
	for (auto& e : context->entries) {
		if (e.id == id) {
			if (e.ptr) free_structure_data(e.ptr);
			FREE(e.ptr);
			context->entries.swap_back_and_pop(&e);
			return true;
		}
	}
	LOG_ERROR("Could not remove structure: it does not exists!");
	return false;
}

void clear_structures() {
	ASSERT(context);
	for (auto& e : context->entries) {
		if (e.ptr) free_structure_data(e.ptr);
	}
	context->entries.clear();
}

template <typename T>
void extract_data_from_indices(T* RESTRICT dst_data, const T* RESTRICT src_data, int32* RESTRICT indices, int32 count) {
	for (int32 i = 0; i < count; i++) {
		dst_data[i] = src_data[indices[i]];
	}
}

mat3 compute_rotation(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
					  const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
					  const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1)
{
	mat3 M = compute_linear_transform(x0, y0, z0, x1, y1, z1, mass, count, com0, com1);
	mat3 R, S;
	decompose(M, &R, &S);
	return R;
}

bool compute_trajectory_transform_data(ID id, Array<const bool> atom_mask, const MoleculeDynamic& dynamic, int32 target_frame_idx) {
	ASSERT(context);

	const int32 num_frames = (int32)dynamic.trajectory.num_frames;
	if (target_frame_idx >= num_frames) {
		LOG_ERROR("Supplied target frame index is out of range.");
		return false;
	}

	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Could not compute tracking data, supplied id is not valid.");
		return false;
	}

	int* indices = (int*)TMP_MALLOC(sizeof(int) * atom_mask.size());
	defer { TMP_FREE(indices); };
	int num_points = 0;
	for (int i = 0; i < (int32)atom_mask.size(); i++) {
		if (atom_mask[i]) {
			indices[num_points] = i;
			num_points++;
		}
	}

	if (num_points == 0) {
		LOG_ERROR("Supplied atom mask is empty.");
		return false;
	}

	// Allocate memory for all data
	init_structure_data(s, id, target_frame_idx, num_points, num_frames);
	s->transform[target_frame_idx] = {};

	// Scratch data
	const auto mem_size = sizeof(float) * num_points * 3;
	void* mem = TMP_MALLOC(mem_size);
	defer{ TMP_FREE(mem); };

	float* cur_x = (float*)mem;
	float* cur_y = cur_x + num_points;
	float* cur_z = cur_y + num_points;
	//float* ref_x = cur_z + num_points;
	//float* ref_y = ref_x + num_points;
	//float* ref_z = ref_y + num_points;
	//float* err_x = cur_z + num_points;
	//float* err_y = err_x + num_points;
	//float* err_z = err_y + num_points;

	float* ref_x = s->reference.x;
	float* ref_y = s->reference.y;
	float* ref_z = s->reference.z;
	float* mass = s->reference.mass;

	extract_data_from_indices(ref_x, get_trajectory_position_x(dynamic.trajectory, target_frame_idx).data(), indices, num_points);
	extract_data_from_indices(ref_y, get_trajectory_position_y(dynamic.trajectory, target_frame_idx).data(), indices, num_points);
	extract_data_from_indices(ref_z, get_trajectory_position_z(dynamic.trajectory, target_frame_idx).data(), indices, num_points);
	const vec3 ref_com = compute_com(ref_x, ref_y, ref_z, mass, num_points);

	for (int32 cur_idx = 0; cur_idx < num_frames; cur_idx++) {
		if (cur_idx == target_frame_idx) continue;
		extract_data_from_indices(cur_x, get_trajectory_position_x(dynamic.trajectory, cur_idx).data(), indices, num_points);
		extract_data_from_indices(cur_y, get_trajectory_position_y(dynamic.trajectory, cur_idx).data(), indices, num_points);
		extract_data_from_indices(cur_z, get_trajectory_position_z(dynamic.trajectory, cur_idx).data(), indices, num_points);

		const vec3 cur_com = compute_com(cur_x, cur_y, cur_z, mass, num_points);
		
		//float* rbf_x = s->data.rbf.weight.x + i * num_points;
		//float* rbf_y = s->data.rbf.weight.y + i * num_points;
		//float* rbf_z = s->data.rbf.weight.z + i * num_points;

		// @NOTE: Compute linear transformation matrix between the two sets of points.
		const mat3 cur_rot = compute_linear_transform(cur_x, cur_y, cur_z, ref_x, ref_y, ref_z, mass, num_points, cur_com, ref_com);
	
		// @NOTE: Compute residual error between the two sets of points.
		//compute_residual_error(err_x, err_y, err_z, cur_x, cur_y, cur_z, ref_x, ref_y, ref_z, num_points, *M);

		// @NOTE: Encode residual error as rbf at reference points
		//compute_rbf_weights(rbf_x, rbf_y, rbf_z, ref_x, ref_y, ref_z, err_x, err_y, err_z, num_points, rbf_radial_cutoff);

		s->transform->rotation = cur_rot;
		s->transform->translation = cur_com;
	}

	return true;
}

const Transform& get_transform_to_target_frame(ID id, int32 source_frame) {
	ASSERT(context);
	static const Transform invalid_transform{};

	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Supplied id is not valid.");
		return invalid_transform;
	}

	if (source_frame < 0 || s->num_frames <= source_frame) {
		LOG_ERROR("Supplied frame is out of range");
		return invalid_transform;
	}

	return s->transform[source_frame];
}

void apply_transform(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, int64 count, const Transform& t, TransformFlags flags) {
	const mat4 R = (flags & TransformFlag_Rotate) ? t.rotation : mat4(1);
	const mat4 T = (flags & TransformFlag_Translate) ? mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(t.translation, 1)) : mat4(1);
	const mat4 M = T * R;
	transform_positions(x, y, z, count, M);
}

} // namespace structure_tracking