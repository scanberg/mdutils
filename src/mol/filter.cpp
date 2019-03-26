#include "filter.h"
#include <core/common.h>
#include <core/hash.h>
#include <core/log.h>
//#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>



namespace filter {

struct FilterContext {
	const MoleculeStructure& mol;
	const Array<const StoredSelection> sel;
};

typedef bool(*FilterCommandFunc)(Array<bool> mask, const FilterContext& ctx, Array<const CString> args);

struct FilterCommand {
	StringBuffer<16> keyword{};
	FilterCommandFunc func = nullptr;
};

struct Context {
	DynamicArray<FilterCommand> filter_commands{};
};

Context* context = nullptr;

static bool is_modifier(CString str) {
    if (compare_ignore_case(str, "and")) return true;
    if (compare_ignore_case(str, "or")) return true;
    return false;
}

static bool is_keyword(CString str) {
    if (compare_ignore_case(str, "and")) return true;
    if (compare_ignore_case(str, "or")) return true;
    if (compare_ignore_case(str, "not")) return true;
    return false;
}

int32 count_parentheses(CString str) {
    int beg_parentheses_count = 0;
    int end_parentheses_count = 0;

    for (int64 i = 0; i < str.count; i++) {
        if (str[i] == '(') beg_parentheses_count++;
        if (str[i] == ')') end_parentheses_count++;
    }
    return beg_parentheses_count - end_parentheses_count;
}

CString extract_parenthesis(CString str) {
    const uint8* beg = str.beg();
    while (beg != str.end() && *beg != '(') beg++;
    if (beg == str.end()) return {};

    const uint8* end = beg + 1;
    int count = 1;
    while (end++ != str.end() && count > 0) {
        if (*end == '(') count++;
        if (*end == ')') count--;
    }
    if (end == str.end()) return {};
    return CString(beg, end);
}

DynamicArray<CString> extract_chunks(CString str) {
    DynamicArray<CString> chunks;

    const uint8* beg = str.beg();
    while (beg != str.end()) {
        if (*beg == '(') {
            CString par = extract_parenthesis(CString(beg, str.end()));
            chunks.push_back({par.beg(), par.end()});  // Exclude actual parentheses
            beg = par.end();
        } else if (*beg != ' ') {
            const uint8* end = beg;
            while (end != str.end() && *end != ' ') end++;
            chunks.push_back(CString(beg, end));
            beg = end;
        } else
            beg++;
    }

    DynamicArray<CString> big_chunks;

    CString* chunk = chunks.beg();
    while (chunk != chunks.end()) {
        if (chunk->front() == '(') {
            big_chunks.push_back(*chunk);
            chunk++;
        } else if (is_keyword(*chunk)) {
            big_chunks.push_back(*chunk);
            chunk++;
        } else {
            CString* beg_chunk = chunk;
            CString* end_chunk = chunk + 1;
            while (end_chunk != chunks.end() && !is_modifier(*end_chunk) && end_chunk->front() != '(') end_chunk++;
            big_chunks.push_back({beg_chunk->beg(), (end_chunk - 1)->end()});
            chunk = end_chunk;
        }
    }

    return big_chunks;
}

FilterCommand* find_filter_command(CString command) {
	ASSERT(context);
    for (auto& f : context->filter_commands) {
        if (compare(command, f.keyword)) return &f;
    }
    return nullptr;
}

void combine_mask_and(Array<bool> dst, Array<bool> src_a, Array<bool> src_b, bool state_not) {
    if (state_not) {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] & !src_b[i];
    } else {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] & src_b[i];
    }
}

void combine_mask_or(Array<bool> dst, Array<bool> src_a, Array<bool> src_b, bool state_not) {
    if (state_not) {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] | !src_b[i];
    } else {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] | src_b[i];
    }
}

bool internal_filter_mask(Array<bool> mask, CString filter, const FilterContext& ctx) {
    DynamicArray<CString> chunks = extract_chunks(filter);
    DynamicArray<bool> chunk_mask(mask.count);

    bool state_and = true;
    bool state_or = false;
    bool state_not = false;

    for (const auto& chunk : chunks) {
        if (compare_ignore_case(chunk, "and")) {
            state_and = true;
        } else if (compare_ignore_case(chunk, "or")) {
            state_or = true;
        } else if (compare_ignore_case(chunk, "not")) {
            state_not = true;
        } else {
            if (chunk.front() == '(') {
                ASSERT(chunk.back() == ')');
                if (!internal_filter_mask(chunk_mask, CString(chunk.beg() + 1, chunk.end() - 1), ctx)) return false;
            } else {
                auto tokens = ctokenize(chunk);
                auto cmd = find_filter_command(tokens[0]);
                if (!cmd) {
                    StringBuffer<32> buf = tokens[0];
                    LOG_ERROR("Could not match command: '%s'\n", buf.beg());
                    return false;
                }

                auto args = tokens.subarray(1);

                while (args.count > 0 && compare_ignore_case(args[0], "not")) {
                    state_not = !state_not;
                    args = args.subarray(1);
                }

                if (!cmd->func(chunk_mask, ctx, args)) {
                    StringBuffer<32> buf = tokens[0];
                    LOG_ERROR("Could not parse command: '%s' with arguments: ", buf.beg());
                    for (const auto& arg : args) {
                        buf = arg;
                        printf("'%s'", buf.beg());
                    }
                    return false;
                }
            }

            if (state_and)
                combine_mask_and(mask, mask, chunk_mask, state_not);
            else if (state_or)
                combine_mask_or(mask, mask, chunk_mask, state_not);

            state_and = false;
            state_or = false;
            state_not = false;
        }

        if (state_and && state_or) {
            LOG_ERROR("Cannot use both 'and' and 'or' to combine filter options\n");
            return false;
        }
    }

    return true;
}

bool compute_filter_mask(Array<bool> mask, const CString filter, const MoleculeStructure& molecule, Array<const StoredSelection> stored_selections) {
	ASSERT(molecule);
    ASSERT(molecule.atom.count == mask.count);

    if (count_parentheses(filter) != 0) {
        LOG_ERROR("Unmatched parentheses\n");
        return false;
    }

	const FilterContext ctx{ molecule, stored_selections };

    memset(mask.data(), 1, mask.count);
    return internal_filter_mask(mask, filter, ctx);
}

void filter_colors(Array<uint32> colors, Array<bool> mask) {
    ASSERT(colors.count == mask.count);
    for (int i = 0; i < colors.count; i++) {
        if (mask[i])
            colors[i] |= 0xff000000;
        else
            colors[i] &= ~0xff000000;
    }
}

void desaturate_colors(Array<uint32> colors, Array<bool> mask, float scale) {
    ASSERT(colors.count == mask.count);
    for (int i = 0; i < colors.count; i++) {
        if (!mask[i]) continue;

        vec4 rgba = math::convert_color(colors[i]);
        vec3 hsv = math::rgb_to_hsv((vec3)rgba);
        hsv.y *= scale;
        rgba = vec4(math::hsv_to_rgb(hsv), rgba.a);
        colors[i] = math::convert_color(rgba);
    }
}

void initialize() {

	if (context) return;
	context = NEW(Context);

    /*
            all
            water
            aminoacid
            backbone?
            protein

            name
                        type (alias name)
            element
            atomicnumber
            atom
            residue
            resname
            resid
            chain
            chainid
    */

    auto filter_amino_acid = [](Array<bool> mask, const FilterContext& ctx, Array<const CString>) {
        memset(mask.data(), 0, mask.count);
        for (const auto& res : ctx.mol.residues) {
            if (is_amino_acid(res)) {
                memset(mask.data() + res.atom_range.beg, 1, res.atom_range.end - res.atom_range.beg);
            }
        }
        return true;
    };

    auto filter_atom_name = [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
        if (args.count == 0) return false;

        for (int64 i = 0; i < ctx.mol.atom.count; i++) {
            mask[i] = false;
            for (const auto& arg : args) {
                if (compare(ctx.mol.atom.label[i], arg)) {
                    mask[i] = true;
                    break;
                }
            }
        }
        return true;
    };

    context->filter_commands.push_back({"all", [](Array<bool> mask, const FilterContext&, Array<const CString>) {
                                   memset(mask.data(), 1, mask.size_in_bytes());
                                   return true;
                               }});
    context->filter_commands.push_back({"water", [](Array<bool> mask, const FilterContext& ctx, Array<const CString>) {
                                   memset(mask.data(), 0, mask.size_in_bytes());
                                   for (const auto& res : ctx.mol.residues) {
                                       const auto res_size = res.atom_range.end - res.atom_range.beg;
                                       if (res_size == 3) {
                                           int32 h_count = 0;
                                           int32 o_count = 0;
                                           for (auto e : get_elements(ctx.mol, res)) {
                                               if (e == Element::H) h_count++;
                                               if (e == Element::O) o_count++;
                                           }
                                           if (h_count == 2 && o_count == 1) {
                                               memset(mask.data() + res.atom_range.beg, 1, res_size);
                                           }
                                       }
                                   }
                                   return true;
                               }});
    context->filter_commands.push_back({"aminoacid", filter_amino_acid});
	context->filter_commands.push_back({ "backbone", [](Array<bool> mask, const FilterContext& ctx, Array<const CString>) {
									memset(mask.data(), 0, mask.size_in_bytes());
									for (int64 i = 0; i < ctx.mol.atom.count; i++) {
										if (compare_n(ctx.mol.atom.label[i], "CA", 2)) {
											mask[i] = true;
										}
									}
									return true;
								}});
    context->filter_commands.push_back({"protein", filter_amino_acid});
    context->filter_commands.push_back({"dna", [](Array<bool> mask, const FilterContext& ctx, Array<const CString>) {
                                   memset(mask.data(), 0, mask.size_in_bytes());
                                   for (const auto& res : ctx.mol.residues) {
                                       if (is_dna(res)) {
                                           memset(mask.data() + res.atom_range.beg, 1, res.atom_range.end - res.atom_range.beg);
                                       }
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"name", filter_atom_name});
    context->filter_commands.push_back({"label", filter_atom_name});
    context->filter_commands.push_back({"type", filter_atom_name});

    context->filter_commands.push_back({"element", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   Array<Element> elements = {(Element*)(TMP_MALLOC(args.count * sizeof(Element))), args.count};
                                   defer { TMP_FREE(elements.data()); };
                                   for (int64 i = 0; i < elements.count; i++) {
                                       elements[i] = element::get_from_string(args[i]);
                                       if (elements[i] == Element::Unknown) return false;
                                   }

                                   for (int64 i = 0; i < ctx.mol.atom.count; i++) {
                                       mask[i] = false;
                                       for (const auto& ele : elements) {
                                           if (ctx.mol.atom.element[i] == ele) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"atomicnumber", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   for (int64 i = 0; i < ctx.mol.atom.count; i++) {
                                       int atomnr = (int)ctx.mol.atom.element[i];
                                       mask[i] = false;
                                       for (auto range : ranges) {
                                           if (range.x == -1) range.x = 0;
                                           if (range.y == -1) range.y = element::num_elements;
                                           if (range.x <= atomnr && atomnr <= range.y) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"atom", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   memset(mask.data(), 0, mask.size_in_bytes());
                                   if (ctx.mol.atom.count == 0) return true;
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   for (auto range : ranges) {
                                       if (range.x == -1) range.x = 0;
                                       if (range.y == -1) range.y = (int32)ctx.mol.atom.count - 1;
                                       range.x = math::clamp(range.x - 1, 0, (int32)ctx.mol.atom.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)ctx.mol.atom.count - 1);
                                       if (range.x == range.y)
                                           mask[range.x] = true;
                                       else
                                           memset(mask.data() + range.x, 1, range.y - range.x);
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"residue", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   memset(mask.data(), 0, mask.size_in_bytes());
                                   if (ctx.mol.residues.count == 0) return true;
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   for (auto range : ranges) {
                                       if (range.x == -1) range.x = 0;
                                       if (range.y == -1) range.y = (int32)ctx.mol.atom.count - 1;
                                       range.x = math::clamp(range.x, 0, (int32)ctx.mol.residues.count - 1);
                                       range.y = math::clamp(range.y, 0, (int32)ctx.mol.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           const auto beg = ctx.mol.residues[i].atom_range.beg;
                                           const auto end = ctx.mol.residues[i].atom_range.end;
                                           memset(mask.data() + beg, 1, end - beg);
                                       }
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"resname", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   memset(mask.data(), 0, mask.count);
                                   for (int i = 0; i < args.count; i++) {
                                       for (const auto& res : ctx.mol.residues) {
                                           if (compare(args[i], res.name)) {
                                               const auto beg = res.atom_range.beg;
                                               const auto end = res.atom_range.end;
                                               memset(mask.data() + beg, 1, end - beg);
                                           }
                                       }
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"resid", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   memset(mask.data(), 0, mask.size_in_bytes());
                                   if (ctx.mol.residues.count == 0) return true;
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   for (auto range : ranges) {
                                       if (range.x == -1) range.x = 0;
                                       if (range.y == -1) range.y = (int32)ctx.mol.residues.count - 1;
                                       range.x = math::clamp(range.x - 1, 0, (int32)ctx.mol.residues.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)ctx.mol.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           const auto beg = ctx.mol.residues[i].atom_range.beg;
                                           const auto end = ctx.mol.residues[i].atom_range.end;
                                           memset(mask.data() + beg, 1, end - beg);
                                       }
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"chain", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   memset(mask.data(), 0, mask.size_in_bytes());
                                   if (ctx.mol.chains.count == 0) return true;
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   for (auto range : ranges) {
                                       if (range.x == -1) range.x = 0;
                                       if (range.y == -1) range.y = (int32)ctx.mol.atom.count - 1;
                                       range.x = math::clamp(range.x - 1, 0, (int32)ctx.mol.chains.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)ctx.mol.chains.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           Chain chain = get_chain(ctx.mol, (ChainIdx)i);
                                           memset(mask.data() + chain.atom_range.beg, 1, chain.atom_range.end - chain.atom_range.beg);
                                       }
                                   }
                                   return true;
                               }});

    context->filter_commands.push_back({"chainid", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
                                   memset(mask.data(), 0, mask.count);
                                   for (int i = 0; i < args.count; i++) {
                                       for (const auto& chain : ctx.mol.chains) {
                                           if (compare(args[i], chain.id)) {
                                               memset(mask.data() + chain.atom_range.beg, 1, chain.atom_range.end - chain.atom_range.beg);
											   break;
                                           }
                                       }
                                   }
                                   return true;
                               }});

	context->filter_commands.push_back({ "selection", [](Array<bool> mask, const FilterContext& ctx, Array<const CString> args) {
							   memset(mask.data(), 0, mask.count);
							   for (int i = 0; i < args.count; i++) {
								   for (const auto& s : ctx.sel) {
									   if (compare(args[i], s.name)) {
										   memcpy(mask.data(), s.mask.data(), mask.size_in_bytes());
										   break;
									   }
								   }
							   }
							   return true;
						   } });
}

void shutdown() {
	if (context) {
		context->~Context();
		FREE(context);
	}
}

}  // namespace filter
