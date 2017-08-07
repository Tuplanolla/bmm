extern inline void type(zero, A, D)(A *);
extern inline void type(add, A, D)(A *restrict,
    A const *restrict, A const *restrict);
inline void type(neg, A, D)(A *restrict, A const *restrict);
extern inline void type(sub, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void type(one, A, D)(A *);
extern inline void type(mul, A, D)(A *restrict,
    A const *restrict, A const *restrict);
inline void type(recip, A, D)(A *restrict, A const *restrict);
extern inline void type(div, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void type(quot, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void type(rem, A, D)(A *restrict,
    A const *restrict, A const *restrict);

extern inline void type(add_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(neg_mut, A, D)(A *restrict);
extern inline void type(sub_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(mul_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(recip_mut, A, D)(A *restrict);
extern inline void type(div_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(quot_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(rem_mut, A, D)(A *restrict, A const *restrict);
