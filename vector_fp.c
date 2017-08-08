extern inline void $(zero, A, D)(A *);
extern inline void $(add, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void $(neg, A, D)(A *restrict, A const *restrict);
extern inline void $(sub, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void $(one, A, D)(A *);
extern inline void $(mul, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void $(recip, A, D)(A *restrict, A const *restrict);
extern inline void $(div, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void $(quott, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void $(remt, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void $(quote, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void $(reme, A, D)(A *restrict,
    A const *restrict, A const *restrict);

extern inline void $(add_mut, A, D)(A *restrict, A const *restrict);
extern inline void $(neg_mut, A, D)(A *restrict);
extern inline void $(sub_mut, A, D)(A *restrict, A const *restrict);
extern inline void $(mul_mut, A, D)(A *restrict, A const *restrict);
extern inline void $(recip_mut, A, D)(A *restrict);
extern inline void $(div_mut, A, D)(A *restrict, A const *restrict);
extern inline void $(quott_mut, A, D)(A *restrict, A const *restrict);
extern inline void $(remt_mut, A, D)(A *restrict, A const *restrict);
extern inline void $(quote_mut, A, D)(A *restrict, A const *restrict);
extern inline void $(reme_mut, A, D)(A *restrict, A const *restrict);
