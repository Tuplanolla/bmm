extern inline signed_char type(zero, signed_char)(void);
extern inline unsigned_char type(zero, unsigned_char)(void);
extern inline int type(zero, int)(void);
extern inline double type(zero, double)(void);
extern inline size_t type(zero, size_t)(void);

extern inline signed_char type(one, signed_char)(void);
extern inline unsigned_char type(one, unsigned_char)(void);
extern inline int type(one, int)(void);
extern inline double type(one, double)(void);
extern inline size_t type(one, size_t)(void);

extern inline signed_char type(minval, signed_char)(void);
extern inline unsigned_char type(minval, unsigned_char)(void);
extern inline int type(minval, int)(void);
extern inline double type(minval, double)(void);
extern inline size_t type(minval, size_t)(void);

extern inline signed_char type(maxval, signed_char)(void);
extern inline unsigned_char type(maxval, unsigned_char)(void);
extern inline int type(maxval, int)(void);
extern inline double type(maxval, double)(void);
extern inline size_t type(maxval, size_t)(void);

extern inline signed_char type(trunc, signed_char)(signed_char, signed_char);
extern inline unsigned_char type(trunc, unsigned_char)(unsigned_char,
    unsigned_char);
extern inline int type(trunc, int)(int, int);
extern inline double type(trunc, double)(double, double);
extern inline size_t type(trunc, size_t)(size_t, size_t);

extern inline signed_char type(mod, signed_char)(signed_char, signed_char);
extern inline unsigned_char type(mod, unsigned_char)(unsigned_char,
    unsigned_char);
extern inline int type(mod, int)(int, int);
extern inline double type(mod, double)(double, double);
extern inline size_t type(mod, size_t)(size_t, size_t);
