#include "MesLib.h"

size_t P(const size_t x, const size_t y, const size_t z) noexcept
{
    return (y * (mx + 1) + x) * 2 + z;
}

size_t Q(const size_t x, const size_t y) noexcept
{
    return y * mx + x;
}

size_t DOF(const size_t elidx, const size_t elidy, const size_t locdofid)
{
    if (locdofid < 4)
    {
        return 2 * (elidy * (mx + 1) + elidx) + locdofid;
    }
    if (locdofid == 4 || locdofid == 5)
    {
        return 2 * ((elidy + 1) * (mx + 1) + elidx - 1) + locdofid;
    }
    if (locdofid == 6 || locdofid == 7)
    {
        return 2 * ((elidy + 1) * (mx + 1) + elidx - 3) + locdofid;
    }
    throw std::logic_error{ "Lokalny stopien swobody musi nalezec do zakresu 0-7" };
}