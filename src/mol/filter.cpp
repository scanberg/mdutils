#include "filter.h"
#include <core/common.h>
#include <core/hash.h>
#include <core/log.h>
#include <mol/element.h>
#include <mol/element_utils.h>
#include <mol/aminoacid.h>
#include <mol/aminoacid_utils.h>
#include <mol/molecule_utils.h>

namespace filter {

struct FilterContext {
    const MoleculeStructure& mol;
    const Array<const StoredSelection> sel;
};

typedef bool (*FilterCommandFunc)(Bitfield mask, const FilterContext& ctx, Array<const CStringView> args);

struct FilterCommand {
    StringBuffer<16> keyword{};
    FilterCommandFunc func = nullptr;
};

struct Context {
    DynamicArray<FilterCommand> filter_commands{};
};

Context* context = nullptr;

static bool is_modifier(CStringView str) {
    if (compare_ignore_case(str, "and")) return true;
    if (compare_ignore_case(str, "or")) return true;
    return false;
}

static bool is_keyword(CStringView str) {
    if (compare_ignore_case(str, "and")) return true;
    if (compare_ignore_case(str, "or")) return true;
    if (compare_ignore_case(str, "not")) return true;
    return false;
}

i32 count_parentheses(CStringView str) {
    int beg_parentheses_count = 0;
    int end_parentheses_count = 0;

    for (i64 i = 0; i < str.count; i++) {
        if (str[i] == '(') beg_parentheses_count++;
        if (str[i] == ')') end_parentheses_count++;
    }
    return beg_parentheses_count - end_parentheses_count;
}

CStringView extract_parenthesis(CStringView str) {
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
    return CStringView(beg, end);
}

DynamicArray<CStringView> extract_chunks(CStringView str) {
    DynamicArray<CStringView> chunks;

    const char* beg = str.beg();
    while (beg != str.end()) {
        if (*beg == '(') {
            CStringView par = extract_parenthesis(CStringView(beg, str.end()));
            chunks.push_back({par.beg(), par.end()});  // Exclude actual parentheses
            beg = par.end();
        } else if (*beg != ' ') {
            const char* end = beg;
            while (end != str.end() && *end != ' ') end++;
            chunks.push_back(CStringView(beg, end));
            beg = end;
        } else
            beg++;
    }

    DynamicArray<CStringView> big_chunks;

    CStringView* chunk = chunks.beg();
    while (chunk != chunks.end()) {
        if (chunk->front() == '(') {
            big_chunks.push_back(*chunk);
            chunk++;
        } else if (is_keyword(*chunk)) {
            big_chunks.push_back(*chunk);
            chunk++;
        } else {
            CStringView* beg_chunk = chunk;
            CStringView* end_chunk = chunk + 1;
            while (end_chunk != chunks.end() && !is_modifier(*end_chunk) && end_chunk->front() != '(') end_chunk++;
            big_chunks.push_back({beg_chunk->beg(), (end_chunk - 1)->end()});
            chunk = end_chunk;
        }
    }

    return big_chunks;
}

FilterCommand* find_filter_command(CStringView command) {
    ASSERT(context);
    for (auto& f : context->filter_commands) {
        if (compare(command, f.keyword)) return &f;
    }
    return nullptr;
}

bool internal_filter_mask(Bitfield mask, CStringView filter, const FilterContext& ctx) {
    DynamicArray<CStringView> chunks = extract_chunks(filter);
    Bitfield chunk_mask;
    bitfield::init(&chunk_mask, ctx.mol.atom.count);
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
                if (!internal_filter_mask(chunk_mask, CStringView(chunk.beg() + 1, chunk.end() - 1), ctx)) return false;
            } else {
                auto tokens = tokenize(chunk);
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
                            LOG_NOTE("'%s'", buf.beg());
                        }
                        return false;
                    }
                }
            }

            if (state_not) {
                bitfield::invert_all(chunk_mask);
            }
            if (state_and) {
                if (mask) bitfield::and_field(mask, mask, chunk_mask);
            } else if (state_or) {
                if (mask) bitfield::or_field(mask, mask, chunk_mask);
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

bool compute_filter_mask(Bitfield mask, const CStringView filter, const MoleculeStructure& molecule, Array<const StoredSelection> stored_selections) {
    ASSERT(molecule);
    ASSERT(molecule.atom.count == mask.bit_count);

    if (count_parentheses(filter) != 0) {
        LOG_ERROR("Unmatched parentheses\n");
        return false;
    }

    const FilterContext ctx{molecule, stored_selections};

    if (mask) bitfield::clear_all(mask);
    return internal_filter_mask(mask, filter, ctx);
}

bool filter_uses_selection(CStringView filter, Array<const StoredSelection> stored_selectons) {
    if (stored_selectons.size() == 0) return false;

    if (count_parentheses(filter) != 0) {
        LOG_ERROR("Unmatched parentheses\n");
        return false;
    }

    DynamicArray<CStringView> chunks = extract_chunks(filter);
    for (const auto& chunk : chunks) {
        if (chunk.front() == '(') {
            ASSERT(chunk.back() == ')');
            if (filter_uses_selection(CStringView(chunk.beg() + 1, chunk.end() - 1), stored_selectons)) return true;
        } else {
            for (const auto& token : tokenize(chunk)) {
                for (const auto& sel : stored_selectons) {
                    if (compare(token, sel.name) == true) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

static inline Range<i32> fix_range(Range<i32> user_range, i32 min, i32 max) {
    if (user_range.beg == -1) user_range.beg = min;
    if (user_range.end == -1) user_range.end = max;
    return {math::clamp(user_range.beg, min, max), math::clamp(user_range.end, min, max)};
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

    auto filter_amino_acid = [](Bitfield mask, const FilterContext& ctx, Array<const CStringView>) {
        for (i64 i = 0; i < ctx.mol.residue.count; i++) {
            if (is_amino_acid(ctx.mol.residue.name[i])) {
                bitfield::set_range(mask, ctx.mol.residue.atom_range[i]);
            }
        }
        return true;
    };

    auto filter_atom_name = [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
        if (args.size() == 0) return false;
        for (i64 i = 0; i < ctx.mol.atom.count; i++) {
            for (const auto& arg : args) {
                if (compare(ctx.mol.atom.label[i], arg)) {
                    bitfield::set_bit(mask, i);
                    break;
                }
            }
        }
        return true;
    };

    context->filter_commands.push_back({"all", [](Bitfield mask, const FilterContext&, Array<const CStringView>) {
                                            bitfield::set_all(mask);
                                            return true;
                                        }});
    context->filter_commands.push_back({"water", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView>) {
                                            for (i64 i = 0; i < ctx.mol.residue.count; i++) {
                                                const auto res_size = ctx.mol.residue.atom_range[i].ext();
                                                if (res_size == 3) {
                                                    i32 h_count = 0;
                                                    i32 o_count = 0;
                                                    for (auto e : get_residue_elements(ctx.mol, (ResIdx)i)) {
                                                        if (e == Element::H) h_count++;
                                                        if (e == Element::O) o_count++;
                                                    }
                                                    if (h_count == 2 && o_count == 1) {
                                                        bitfield::set_range(mask, ctx.mol.residue.atom_range[i]);
                                                    }
                                                }
                                            }
                                            return true;
                                        }});
    context->filter_commands.push_back({"aminoacid", filter_amino_acid});
    context->filter_commands.push_back({"backbone", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView>) {
                                            for (i64 i = 0; i < ctx.mol.atom.count; i++) {
                                                if (compare_n(ctx.mol.atom.label[i], "CA", 2)) {
                                                    bitfield::set_bit(mask, i);
                                                }
                                            }
                                            return true;
                                        }});
    context->filter_commands.push_back({"protein", filter_amino_acid});
    context->filter_commands.push_back({"dna", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView>) {
                                            for (i64 i = 0; i < ctx.mol.residue.count; i++) {
                                                if (is_dna(ctx.mol.residue.name[i])) {
                                                    bitfield::set_range(mask, ctx.mol.residue.atom_range[i]);
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"name", filter_atom_name});
    context->filter_commands.push_back({"label", filter_atom_name});
    context->filter_commands.push_back({"type", filter_atom_name});

    context->filter_commands.push_back({"element", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            Array<Element> elements = {(Element*)(TMP_MALLOC(args.count * sizeof(Element))), args.count};
                                            defer { TMP_FREE(elements.data()); };
                                            for (i64 i = 0; i < elements.count; i++) {
                                                elements[i] = get_element_from_string(args[i]);
                                                if (elements[i] == Element::Unknown) return false;
                                            }

                                            for (i64 i = 0; i < ctx.mol.atom.count; i++) {
                                                for (const auto& ele : elements) {
                                                    if (ctx.mol.atom.element[i] == ele) {
                                                        bitfield::set_bit(mask, i);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"atomicnumber", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            DynamicArray<Range<i32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (i64 i = 0; i < ctx.mol.atom.count; i++) {
                                                const int atomnr = (int)ctx.mol.atom.element[i];
                                                for (auto range : ranges) {
                                                    range = fix_range(range, 1, element::num_elements - 1);
                                                    if (range.beg <= atomnr && atomnr <= range.end) {
                                                        bitfield::set_bit(mask, i);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"atom", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            if (ctx.mol.atom.count == 0) return true;
                                            DynamicArray<Range<i32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto range : ranges) {
                                                range = fix_range(range, 1, (i32)ctx.mol.atom.count);
                                                bitfield::set_range(mask, Range<i32>(range.beg - 1, range.end));  // @NOTE: [1, N] range to [0, N[
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"residue", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            if (ctx.mol.residue.count == 0) return true;
                                            DynamicArray<Range<i32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto range : ranges) {
                                                range = fix_range(range, 1, (i32)ctx.mol.residue.count);
                                                for (i32 i = range.beg - 1; i < range.end; i++) {
                                                    bitfield::set_range(mask, ctx.mol.residue.atom_range[i]);
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"resname", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            for (const auto& arg : args) {
                                                for (i64 i = 0; i < ctx.mol.residue.count; i++) {
                                                    if (compare(arg, ctx.mol.residue.name[i])) {
                                                        bitfield::set_range(mask, ctx.mol.residue.atom_range[i]);
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"resid", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            if (ctx.mol.residue.count == 0) return true;
                                            DynamicArray<Range<i32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto range : ranges) {
                                                for (i64 i = 0; i < ctx.mol.residue.count; i++) {
                                                    if (range.beg <= ctx.mol.residue.id[i] && ctx.mol.residue.id[i] <= range.end) {
                                                        bitfield::set_range(mask, ctx.mol.residue.atom_range[i]);
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"chain", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            if (ctx.mol.chain.count == 0) return true;
                                            DynamicArray<Range<i32>> ranges;
                                            if (!extract_ranges(&ranges, args)) return false;
                                            for (auto user_range : ranges) {
                                                auto range = fix_range(user_range, 1, (i32)ctx.mol.chain.count);
                                                for (int i = range.beg - 1; i < range.end; i++) {
                                                    bitfield::set_range(mask, ctx.mol.chain.atom_range[i]);
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"chainid", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            for (const auto& arg : args) {
                                                for (i64 i = 0; i < ctx.mol.chain.count; i++) {
                                                    if (compare(arg, ctx.mol.chain.id[i])) {
                                                        bitfield::set_range(mask, ctx.mol.chain.atom_range[i]);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"selection", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
                                            for (const auto& arg : args) {
                                                for (const auto& s : ctx.sel) {
                                                    if (compare(arg, s.name)) {
                                                        bitfield::copy(mask, s.mask);
                                                        break;
                                                    }
                                                }
                                            }
                                            return true;
                                        }});

    context->filter_commands.push_back({"current", [](Bitfield mask, const FilterContext& ctx, Array<const CStringView> args) {
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
