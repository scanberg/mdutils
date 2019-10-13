#include "bitfield.h"

namespace bitfield {

#define B2C(b, i) (b & (1U << i) ? '1' : '0')

void print(const Bitfield field) {
	for (int64 i = 0; i < field.size_in_bytes(); i++) {
		const uint8 b = ((uint8*)field.block_ptr)[i];
		if (i % 4 == 0) printf("\n");
		printf("%c%c%c%c%c%c%c%c ", B2C(b, 0), B2C(b, 1), B2C(b, 2), B2C(b, 3), B2C(b, 4), B2C(b, 5), B2C(b, 6), B2C(b, 7));
	}
}

}