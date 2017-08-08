extern inline bool $(add_ovf, A)(A, A);
extern inline bool $(sub_ovf, A)(A, A);
extern inline bool $(mul_ovf, A)(A, A);
extern inline bool $(quott_ovf, A)(A, A);
extern inline bool $(remt_ovf, A)(A, A);

extern inline A $(add, A)(A, A);
extern inline A $(sub, A)(A, A);
extern inline A $(mul, A)(A, A);
extern inline A $(quott, A)(A, A);
extern inline A $(remt, A)(A, A);
extern inline void $(quotremt, A)(A *, A, A);
extern inline A $(quote, A)(A, A);
extern inline A $(reme, A)(A, A);
extern inline void $(quotreme, A)(A *, A, A);

extern inline void $(add_mut, A)(A *, A);
extern inline void $(sub_mut, A)(A *, A);
extern inline void $(mul_mut, A)(A *, A);
extern inline void $(quott_mut, A)(A *, A);
extern inline void $(remt_mut, A)(A *, A);
extern inline void $(quotremt_mut, A)(A *);
extern inline void $(quote_mut, A)(A *, A);
extern inline void $(reme_mut, A)(A *, A);
extern inline void $(quotreme_mut, A)(A *);
