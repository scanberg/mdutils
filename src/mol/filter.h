#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/string_types.h>
#include <core/bitfield.h>
#include <mol/molecule_dynamic.h>

struct StoredSelection {
	CString name;
	Bitfield mask;
};

namespace filter {
void initialize();
void shutdown();
bool compute_filter_mask(Bitfield mask, CString filter, const MoleculeStructure& molecule, Array<const StoredSelection> stored_selectons = {});

template <typename T>
void extract_filtered_data(DynamicArray<T>* dst_data, Array<const T> src_data, Array<const bool> mask) {
    ASSERT(dst_data);
    dst_data->clear();
    for (int32 i = 0; i < mask.count; i++) {
        if (mask[i]) dst_data->push_back(src_data[i]);
    }
}

template <typename T>
DynamicArray<T> extract_filtered_data(Array<const T> data, Array<const bool> mask) {
    DynamicArray<T> result{};
    extract_filtered_data(&result, data, mask);
    return result;
}

}  // namespace filter
