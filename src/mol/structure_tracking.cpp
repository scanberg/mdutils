#include "structure_tracking.h"

#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#include <Eigen/Eigen>

namespace structure_tracking {

// from here https://stackoverflow.com/questions/4372224/fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
void Diagonalize(const float(&A)[3][3], float(&Q)[3][3], float(&D)[3][3]) {
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

void diagonalize(const mat3& M, mat3* Q, mat3* D) {
	ASSERT(Q);
	ASSERT(D);
	Diagonalize((const float(&)[3][3])M, (float(&)[3][3]) * Q, (float(&)[3][3]) * D);
}

void decompose(const mat3& M, mat3* R, mat3* S) {
	ASSERT(R);
	ASSERT(S);
	mat3 Q, D;
	mat3 AtA = math::transpose(M) * M;
	diagonalize(AtA, &Q, &D);
	float det = math::determinant(AtA);
	D[0][0] = sqrtf(D[0][0]);
	D[1][1] = sqrtf(D[1][1]);
	D[2][2] = sqrtf(D[2][2]);
	*S = math::inverse(Q) * D * Q;
	*R = M * math::inverse(*S);
}

void compute_linear_transform(mat3* dst_mat,
								const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
								const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
								const float* RESTRICT mass, int64 count) {
	
	const vec3 com_a = compute_com(x0, y0, z0, mass, count);
	const vec3 com_b = compute_com(x1, y1, z1, mass, count);

	mat3 Apq{ 0 };
	mat3 Aqq{ 0 };

	for (int32 i = 0; i < count; i++) {
		const float q_x = x0[i] - com_a.x;
		const float q_y = y0[i] - com_a.y;
		const float q_z = z0[i] - com_a.z;

		const float p_x = x1[i] - com_b.x;
		const float p_y = y1[i] - com_b.y;
		const float p_z = z1[i] - com_b.z;

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

	*dst_mat = Apq / Aqq;
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

void compute_rbf_weights(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
						 const float* in_pos_x, const float* in_pos_y, const float* in_pos_z,
						 const float* in_val_x, const float* in_val_y, const float* in_val_z,
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

struct Structure {
	ID id = 0;
	int32 ref_frame_idx = 0;
	int32 num_points = 0;
	int32 num_frames = 0;

	struct {
		mat4* matrix = nullptr;			// transform matrix from the frame to the reference

		struct {						// rbf weights encoding the residual error after linear transform (frames are stored consecutively) 0-1-2-3-4 | 0-1-2-3-4
			float* x = nullptr;
			float* y = nullptr;
			float* z = nullptr;
		} rbf;
	} data;

	StringBuffer<64> name{};
};

static struct {
	DynamicArray<Structure> structures;
} context;

static void free_structure(Structure* s) {
	if (s->data.matrix) FREE(s->data.matrix);
	if (s->data.rbf.x) FREE(s->data.rbf.x);
}

void initialize() {
	context.structures.reserve(8);
}

void shutdown() {
	context.structures.clear();
}

void compute_tracking_data(ID structure_id, Array<const bool> atom_mask, const MoleculeStructure& mol, const MoleculeTrajectory& traj, int32 ref_idx) {
	DynamicArray<int> index_list;
	for (int i = 0; i < atom_mask.size(); i++) {
		if (atom_mask[i]) index_list.push_back(i);
	}

	// Allocate memory for all data


	for (int i = 0; i < traj.num_frames; i++) {
		if (i == ref_idx) continue;
		const auto pos_x = get_trajectory_position_x(traj, i);
	}
}

void transform_coordinates_to_reference(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, ID structure_id) {

}

}