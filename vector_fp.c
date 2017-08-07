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
extern inline void type(quott, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void type(remt, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void type(quote, A, D)(A *restrict,
    A const *restrict, A const *restrict);
extern inline void type(reme, A, D)(A *restrict,
    A const *restrict, A const *restrict);

extern inline void type(add_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(neg_mut, A, D)(A *restrict);
extern inline void type(sub_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(mul_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(recip_mut, A, D)(A *restrict);
extern inline void type(div_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(quott_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(remt_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(quote_mut, A, D)(A *restrict, A const *restrict);
extern inline void type(reme_mut, A, D)(A *restrict, A const *restrict);
