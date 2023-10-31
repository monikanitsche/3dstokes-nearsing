#define THIRD
#define FOURTH

#define DBLE
!#define QUAD

#ifdef DBLE
#define MODE real*8
#endif

#ifdef QUAD
#define MODE real*16
#endif
