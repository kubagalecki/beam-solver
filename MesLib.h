#ifndef METNUM_MESLIB_H
#define METNUM_MESLIB_H

#include <stdexcept>
#include <cstddef>
#include <vector>

extern std::uint_fast16_t mx;
extern std::uint_fast16_t my;
extern bool* fix;
extern double* F;

using std::size_t;

size_t P(const size_t x, const size_t y, const size_t z) noexcept;

size_t Q(const size_t x, const size_t y) noexcept;

size_t DOF(const size_t elidx, const size_t elidy, const size_t locdofid);

constexpr double skala = 1.0;

//------------------------------- MACIERZ SZTYWNOSCI ELEMENTU --------------------------
constexpr double E = 100.;
constexpr double nu = 0.3;
constexpr double elems[] = {
	0, 1. / 2. - nu / 6., 1. / 8. + nu / 8., -1. / 4. - nu / 12., -1. / 8. + 3. * nu / 8.,
	-1. / 4. + nu / 12., -1. / 8. - nu / 8., nu / 6., 1. / 8. - 3. * nu / 8.
};

constexpr double Md = 1.;

constexpr double K[8][8] = { { elems[1], elems[2], elems[3], elems[4], elems[5], elems[6], elems[7], elems[8] },
{ elems[2], elems[1], elems[8], elems[7], elems[6], elems[5], elems[4], elems[3] },
{ elems[3], elems[8], elems[1], elems[6], elems[7], elems[4], elems[5], elems[2] },
{ elems[4], elems[7], elems[6], elems[1], elems[8], elems[3], elems[2], elems[5] },
{ elems[5], elems[6], elems[7], elems[8], elems[1], elems[2], elems[3], elems[4] },
{ elems[6], elems[5], elems[4], elems[3], elems[2], elems[1], elems[8], elems[7] },
{ elems[7], elems[4], elems[5], elems[2], elems[3], elems[8], elems[1], elems[6] },
{ elems[8], elems[3], elems[2], elems[5], elems[4], elems[7], elems[6], elems[1] } };

#endif // !METNUM_MESLIB_H
