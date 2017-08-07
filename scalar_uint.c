extern inline bool type(add_ovf, A)(A, A);
extern inline bool type(sub_ovf, A)(A, A);
extern inline bool type(mul_ovf, A)(A, A);
extern inline bool type(quott_ovf, A)(A, A);
extern inline bool type(remt_ovf, A)(A, A);
extern inline bool type(quote_ovf, A)(A, A);
extern inline bool type(reme_ovf, A)(A, A);

extern inline A type(add, A)(A, A);
extern inline A type(sub, A)(A, A);
extern inline A type(mul, A)(A, A);
extern inline A type(quott, A)(A, A);
extern inline A type(remt, A)(A, A);
extern inline A type(quote, A)(A, A);
extern inline A type(reme, A)(A, A);

extern inline void type(add_mut, A)(A *, A);
extern inline void type(sub_mut, A)(A *, A);
extern inline void type(mul_mut, A)(A *, A);
extern inline void type(quott_mut, A)(A *, A);
extern inline void type(remt_mut, A)(A *, A);
extern inline void type(quote_mut, A)(A *, A);
extern inline void type(reme_mut, A)(A *, A);
