#include "bit.h"

extern inline bool bmm_bit_test(unsigned char, unsigned char);

extern inline unsigned char bmm_bit_set(unsigned char, unsigned char);

extern inline unsigned char bmm_bit_clear(unsigned char, unsigned char);

extern inline unsigned char bmm_bit_compl(unsigned char, unsigned char);

extern inline void bmm_bit_pset(unsigned char*, unsigned char);

extern inline void bmm_bit_pclear(unsigned char*, unsigned char);

extern inline void bmm_bit_pcompl(unsigned char*, unsigned char);
