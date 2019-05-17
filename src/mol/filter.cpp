#include "filter.h"
#include <core/common.h>
#include <core/hash.h>
#include <core/log.h>
#include <mol/molecule_utils.h>

namespace filter {

struct FilterContext {
    const MoleculeStructure& mol;
    const Array<const StoredSelection> sel;
};

typedef bool (*FilterCommandFunc)(Bitfield mask, const FilterContext& ctx, Array<const CString> args);

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
    const char* beg = str.beg();
    while (beg != str.end() && *beg != '(') beg++;
    if (beg == str.end()) return {};

    const char* end = beg + 1;
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

    const char* beg = str.beg();
    while (beg != str.end()) {
        if (*beg == '(') {
            CString par = extract_parenthesis(CString(beg, str.end()));
            chunks.push_back({par.beg(), par.end()});  // Exclude actual parentheses
            beg = par.end();
        } else if (*beg != ' ') {
            const char* end = beg;
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

bool internal_filter_mask(Bitfield mask, CString filter, const FilterContext& ctx) {
    DynamicArray<CString> chunks = extract_chunks(filter);
    Bitfield chunk_mask;
    bitfield::init(&chunk_mask, mask.count);
    defer { bitfield::free(&chunk_mask); };

    bool state_and = false;
    bool state_or = true;
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
                if (tokens.size() > 0) {
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

                    bitfield::clear_all(chunk_mask);
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
            }

            if (state_not) {
                bitfield::invert_all(chunk_mask);
            }
            if (state_and) {
                bitfield::and_field(mask, mask, chunk_mask);
            } else if (state_or) {
                bitfield::or_field(mask, mask, chunk_mask);
            }

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

bool compute_filter_mask(Bitfield mask, const CString filter, const MoleculeStructure& molecule, Array<const StoredSelection> stored_selections) {
    ASSERT(molecule);
    ASSERT(molecule.atom.count == mask.count);

    if (count_parentheses(filter) != 0) {
        LOG_ERROR("Unmatched parentheses\n");
        return false;
    }

    const FilterContext ctx{molecule, stored_selections};

    bitfield::clear_all(mask);
    return internal_filter_mask(mask, filter, ctx);
}

static inline Range<int32> fix_range(Range<int32> user_range, int32 min, int32 max) {
    if (user_range.x == -1) user_range.x = min;
    if (user_range.y == -1) user_range.y = max;
    return {math::clamp(user_range.x, min, max), math::clamp(user_range.y, min, max)};
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

    auto filter_amino_acid = [](Bitfield mask, const FilterContext& ctx, Array<const CString>) {
        for (const auto& res : ctx.mol.residues) {
            if (is_amino_acid(res)) {
                bitfield::set_range(mask, res.atom_range);
            }
        }
        return true;
    };

    auto filter_atom_name = [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
        if (args.size() == 0) return false;
        for (int64 i = 0; i < ctx.mol.atom.count; i++) {
            for (const auto& arg : args) {
                if (compare(ctx.mol.atom.label[i], arg)) {
                    bitfield::set_bit(mask, i);
                    break;
                }
            }
        }
        return true;
    };

    context->filter_commands.push_back({"all", [](Bitfield mask, const FilterContext&, Array<const CString>) {
                                            bitfield::set_all(mask);
                                            return true;
                                        }});
    context->filter_commands.push_back({"water", [](Bitfield mask, const FilterContext& ctx, Array<const CString>) {
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
                                                        bitfield::set_range(mask, res.atom_range);
                                                    }
                                                }
                                            }
                                            return true;
                                        }});
    context->filter_commands.push_back({"aminoacid", filter_amino_acid});
    context->filter_commands.push_back({"backbone", [](Bitfield mask, const FilterContext& ctx, Array<const CString>) {
                                            for (int64 i = 0; i < ctx.mol.atom.count; i++) {
                                                if (compare_n(ctx.mol.atom.label[i], "CA", 2)) {
                                                    bitfield::set_bit(mask, i);
                                                }
                                            }
                                            return true;
                                        }});
    context->filter_commands.push_back({"protein", filter_amino_acid});
    context->filter_commands.push_back({"dna", [](Bitfield mask, const FilterContext& ctx, Array<const CString>) {
                                            for (const auto& res : ctx.mol.residues) {
                                                if (is_dna(res)) {
                                                    bitfield::set_range(mask, res.atom_range);
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"name", filter_atom_name});
    context->filter_commands.push_back({"label", filter_atom_name});
    context->filter_commands.push_back({"type", filter_atom_name});

    context->filter_commands.push_back({"element", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            Array<Element> elements = {(Element*)(TMP_MALLOC(args.count * sizeof(Element))), args.count};
                                            defer { TMP_FREE(elements.data()); };
                                            for (int64 i = 0; i < elements.count; i++) {
                                                elements[i] = element::get_from_string(args[i]);
                                                if (elements[i] == Element::Unknown) return false;
                                            }

                                            for (int64 i = 0; i < ctx.mol.atom.count; i++) {
                                                for (const auto& ele : elements) {
                                                    if (ctx.mol.atom.element[i] == ele) {
                                                        bitfield::set_bit(mask, i);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"atomicnumber", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            DynamicArray<Range<int32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (int64 i = 0; i < ctx.mol.atom.count; i++) {
                                                const int atomnr = (int)ctx.mol.atom.element[i];
                                                for (auto range : ranges) {
                                                    range = fix_range(range, 1, element::num_elements - 1);
                                                    if (range.x <= atomnr && atomnr <= range.y) {
                                                        bitfield::set_bit(mask, i);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"atom", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            if (ctx.mol.atom.count == 0) return true;
                                            DynamicArray<Range<int32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto range : ranges) {
                                                range = fix_range(range, 1, (int32)ctx.mol.atom.count);
                                                bitfield::set_range(mask, Range<int32>(range.x - 1, range.y));  // @NOTE: [1, N] range to [0, N[
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"residue", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            if (ctx.mol.residues.count == 0) return true;
                                            DynamicArray<Range<int32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto range : ranges) {
                                                range = fix_range(range, 1, (int32)ctx.mol.residues.count);
                                                for (int32 i = range.x - 1; i < range.y; i++) {
                                                    const auto& res = ctx.mol.residues[i];
                                                    bitfield::set_range(mask, res.atom_range);
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"resname", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            for (int i = 0; i < args.count; i++) {
                                                for (const auto& res : ctx.mol.residues) {
                                                    if (compare(args[i], res.name)) {
                                                        bitfield::set_range(mask, res.atom_range);
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"resid", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            if (ctx.mol.residues.count == 0) return true;
                                            DynamicArray<Range<int32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto range : ranges) {
                                                for (const auto& res : ctx.mol.residues) {
                                                    if (range.x <= res.id && res.id <= range.y) {
                                                        bitfield::set_range(mask, res.atom_range);
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"chain", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            if (ctx.mol.chains.count == 0) return true;
                                            DynamicArray<Range<int32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto user_range : ranges) {
                                                auto range = fix_range(user_range, 1, (int32)ctx.mol.chains.count);
                                                for (int i = range.x - 1; i < range.y; i++) {
                                                    const Chain& chain = get_chain(ctx.mol, (ChainIdx)i);
                                                    bitfield::set_range(mask, chain.atom_range);
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"sequence", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            if (ctx.mol.sequences.count == 0) return true;
                                            DynamicArray<Range<int32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto user_range : ranges) {
                                                auto range = fix_range(user_range, 1, (int32)ctx.mol.sequences.count);
                                                for (int i = range.x - 1; i < range.y; i++) {
                                                    const Sequence& seq = get_sequence(ctx.mol, (SeqIdx)i);
                                                    bitfield::set_range(mask, seq.atom_range);
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"chainid", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            for (int i = 0; i < args.count; i++) {
                                                for (const auto& chain : ctx.mol.chains) {
                                                    if (compare(args[i], chain.id)) {
                                                        bitfield::set_range(mask, chain.atom_range);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"selection", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            for (int i = 0; i < args.count; i++) {
                                                for (const auto& s : ctx.sel) {
                                                    if (compare(args[i], s.name)) {
                                                        bitfield::copy(mask, s.mask);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"current", [](Bitfield mask, const FilterContext& ctx, Array<const CString> args) {
                                            UNUSED(args);
                                            for (const auto& s : ctx.sel) {
                                                if (compare("current", s.name)) {
                                                    bitfield::copy(mask, s.mask);
                                                    break;
                                                }
                                            }
                                            return true;
                                        }});
}

void shutdown() {
    if (context) {
        context->~Context();
        DELETE(context);
    }
}

}  // namespace filter
