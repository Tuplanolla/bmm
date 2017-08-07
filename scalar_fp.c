extern inline A type(add, A)(A, A);
extern inline A type(neg, A)(A);
extern inline A type(sub, A)(A, A);
extern inline A type(mul, A)(A, A);
extern inline A type(recip, A)(A);
extern inline A type(div, A)(A, A);
extern inline A type(quote, A)(A, A);
extern inline A type(reme, A)(A, A);

extern inline void type(add_mut, A)(A *, A);
extern inline void type(neg_mut, A)(A *);
extern inline void type(sub_mut, A)(A *, A);
extern inline void type(mul_mut, A)(A *, A);
extern inline void type(recip_mut, A)(A *);
extern inline void type(div_mut, A)(A *, A);
extern inline void type(quott_mut, A)(A *, A);
extern inline void type(remt_mut, A)(A *, A);
extern inline void type(quote_mut, A)(A *, A);
extern inline void type(reme_mut, A)(A *, A);
