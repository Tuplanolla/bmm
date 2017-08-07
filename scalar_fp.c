extern inline A $(add, A)(A, A);
extern inline A $(neg, A)(A);
extern inline A $(sub, A)(A, A);
extern inline A $(mul, A)(A, A);
extern inline A $(recip, A)(A);
extern inline A $(div, A)(A, A);
extern inline A $(quote, A)(A, A);
extern inline A $(reme, A)(A, A);

extern inline void $(add_mut, A)(A *, A);
extern inline void $(neg_mut, A)(A *);
extern inline void $(sub_mut, A)(A *, A);
extern inline void $(mul_mut, A)(A *, A);
extern inline void $(recip_mut, A)(A *);
extern inline void $(div_mut, A)(A *, A);
extern inline void $(quott_mut, A)(A *, A);
extern inline void $(remt_mut, A)(A *, A);
extern inline void $(quote_mut, A)(A *, A);
extern inline void $(reme_mut, A)(A *, A);
