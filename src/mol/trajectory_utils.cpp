#include "trajectory_utils.h"

#include <mol/pdb_utils.h>
#include <mol/xtc_utils.h>

bool read_next_trajectory_frame(MoleculeTrajectory* traj) {
	ASSERT(traj);
	if (all_trajectory_frames_read(*traj)) return false;
	if (traj->file.tag == 0) return false;

	switch (traj->file.tag) {
	case pdb::PDB_FILE_TAG:
		return pdb::read_next_trajectory_frame(traj);
	case xtc::XTC_FILE_TAG:
		return xtc::read_next_trajectory_frame(traj);
	default:
		ASSERT(false);
	}

	return false;
}

bool close_file_handle(MoleculeTrajectory* traj) {
	ASSERT(traj);
	bool result = false;

	if (traj->file.handle) {
		switch (traj->file.tag) {
		case pdb::PDB_FILE_TAG:
			result = pdb::close_file_handle(traj);
			break;
		case xtc::XTC_FILE_TAG:
			result = xtc::close_file_handle(traj);
			break;
		default:
			ASSERT(false);
		}
		traj->file.handle = nullptr;
	}
	return result;
}