#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/string_types.h>
#include <core/bitfield.h>
#include <mol/molecule_dynamic.h>

struct StoredSelection {
	CStringView name;
	Bitfield mask;
};

namespace filter {
void initialize();
void shutdown();
bool compute_filter_mask(Bitfield mask, CStringView filter, const MoleculeStructure& molecule, ArrayView<const StoredSelection> stored_selectons = {});
bool filter_uses_selection(CStringView filter, ArrayView<const StoredSelection> stored_selectons);

template <typename T>
void extract_filtered_data(DynamicArray<T>* dst_data, ArrayView<const T> src_data, ArrayView<const bool> mask) {
    ASSERT(dst_data);
    dst_data->clear();
    for (int32 i = 0; i < mask.count; i++) {
        if (mask[i]) dst_data->push_back(src_data[i]);
    }
}

template <typename T>
DynamicArray<T> extract_filtered_data(ArrayView<const T> data, ArrayView<const bool> mask) {
    DynamicArray<T> result{};
    extract_filtered_data(&result, data, mask);
    return result;
}

}  // namespace filter
