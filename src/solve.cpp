#include "solve.h"

extern std::uint_fast16_t mx;
extern std::uint_fast16_t my;
extern bool*              fix;
extern double*            F;

using std::size_t;

// Algo for accumulation of consecutively repeating entries
template < class Iterator, class Equality, class Accumulate >
Iterator accumulate_consecutive(const Iterator&   first_,
                                const Iterator&   last_,
                                const Equality&   eq_op,
                                const Accumulate& ac_op)
{
    // Return if empty
    if (first_ == last_)
        return last_;

    auto first = first_;

    bool moved = false; // Indicate if first move has been made (avoid redundant in-place move)
    auto end   = first; // Iterator to the end of processed block

    while (first != last_)
    {
        if (std::next(first) != last_ && eq_op(*first, *std::next(first)))
        {
            auto current = ac_op(*first, *std::next(first));
            ++first;
            while (std::next(first) != last_ && eq_op(*first, *std::next(first)))
            {
                ++first;
                current = ac_op(current, *first);
            }

            *end++ = std::move(current);
            ++first;
            moved = true;
        }
        else
        {
            if (moved)
            {
                *end = std::move(*first);
            }
            ++first;
            ++end;
        }
    }

    return end;
}

// Block type for sparse matrix
template < class T, std::int_fast8_t SIZE >
struct Block2D
{
    // Block element access
    T& operator()(std::int_fast8_t i, std::int_fast8_t j) { return m_data[i][j]; }

    const T& operator()(std::int_fast8_t i, std::int_fast8_t j) const { return m_data[i][j]; }

    // Add block
    Block2D& operator+=(const Block2D& b)
    {
        for (std::int_fast8_t i = 0; i < SIZE; ++i)
        {
            for (std::int_fast8_t j = 0; j < SIZE; ++j)
            {
                m_data[i][j] += b.m_data[i][j];
            }
        }
        return *this;
    }

    // Block data
    std::array< std::array< T, SIZE >, SIZE > m_data;

    static constexpr std::int_fast8_t size_v = SIZE;
};

template < class T, std::int_fast8_t SIZE >
std::ostream& operator<<(std::ostream& io, const Block2D< T, SIZE >& b)
{
    std::for_each(b.m_data.cbegin(), b.m_data.cend(), [](const std::array< T, SIZE > b) {
        std::for_each(b.cbegin(), b.cend(), [](const T& e) { std::cout << e << "\t"; });
        std::cout << '\n';
    });

    return std::cout;
}

// Block type for vector (easier to keep dimension compatibility than with raw vector)
template < class T, std::int_fast8_t SIZE >
struct Block1D
{
    // Element access
    T&       operator()(std::int_fast8_t i) { return m_data[i]; }
    const T& operator()(std::int_fast8_t i) const { return m_data[i]; }

    // Add block
    Block1D& operator+=(const Block1D& b)
    {
        for (std::int_fast8_t i = 0; i < SIZE; ++i)
        {
            m_data[i] += b.m_data[i];
        }
        return *this;
    }

    // Vector data
    std::array< T, SIZE > m_data;
};

// Block matrix-vector product
template < class T, std::int_fast8_t SIZE >
Block1D< T, SIZE > operator*(const Block2D< T, SIZE >& mat, const Block1D< T, SIZE >& vec)
{
    Block1D< T, SIZE > prod_vec{}; // The default constructor fills arrays with zeros

    std::transform(mat.m_data.cbegin(),
                   mat.m_data.cend(),
                   prod_vec.m_data.begin(),
                   [&vec](const std::array< T, SIZE >& mat_row) {
                       return std::inner_product(
                           mat_row.cbegin(), mat_row.cend(), vec.m_data.cbegin(), 0.);
                   });

    return prod_vec;
}

// Block vector-vector product
template < class T, std::int_fast8_t SIZE >
T operator*(const Block1D< T, SIZE >& a, const Block1D< T, SIZE >& b)
{
    return std::transform_reduce(a.m_data.cbegin(), a.m_data.cend(), b.m_data.cbegin(), 0.);
}

// Block vector-scalar product
template < class T, std::int_fast8_t SIZE >
Block1D< T, SIZE > operator*(Block1D< T, SIZE > a, const double& s)
{
    for (auto& e : a.m_data)
        e *= s;

    return a;
}

// Block vector-vector addition
template < class T, std::int_fast8_t SIZE >
Block1D< T, SIZE > operator+(const Block1D< T, SIZE >& a, const Block1D< T, SIZE >& b)
{
    Block1D< T, SIZE > sum{};
    std::transform(a.m_data.cbegin(),
                   a.m_data.cend(),
                   b.m_data.cbegin(),
                   sum.m_data.begin(),
                   [](const T& v1, const T& v2) { return v1 + v2; });
    return sum;
}

// Block vector-vector subtraction
template < class T, std::int_fast8_t SIZE >
Block1D< T, SIZE > operator-(const Block1D< T, SIZE >& a, const Block1D< T, SIZE >& b)
{
    Block1D< T, SIZE > sum{};
    std::transform(a.m_data.cbegin(),
                   a.m_data.cend(),
                   b.m_data.cbegin(),
                   sum.m_data.begin(),
                   [](const T& v1, const T& v2) { return v1 - v2; });
    return sum;
}

// Block sparse matrix type. T denotes the block type. For simplicity, matlab storage format is
// used, but with row-major order.
template < class T >
class BlockSparseMatrix
{
public:
    using index_t     = std::int_fast64_t;
    using entry_t     = std::pair< std::pair< index_t, index_t >, T >;
    using data_t      = std::vector< entry_t >;
    using data_iter_t = typename data_t::const_iterator;
    using part_t      = std::vector< data_iter_t >;

    static BlockSparseMatrix
                             assembleBSMatrix(data_t&&, const index_t&, const index_t&, const size_t&);
    static BlockSparseMatrix assembleIdentBSMatrix(const size_t);

    BlockSparseMatrix()                             = delete;
    BlockSparseMatrix(const BlockSparseMatrix&)     = default;
    BlockSparseMatrix(BlockSparseMatrix&&) noexcept = default;
    BlockSparseMatrix& operator=(const BlockSparseMatrix&) = default;
    BlockSparseMatrix& operator=(BlockSparseMatrix&&) noexcept = default;

    void print() const;
    void print_sparsity_pattern() const;
    void print_to_file(const char*) const;

    [[nodiscard]] typename data_t::const_iterator begin() const { return m_data.cbegin(); }
    [[nodiscard]] typename data_t::const_iterator end() const { return m_data.cend(); }

    [[nodiscard]] index_t n_rows() const { return m_rows; }
    [[nodiscard]] index_t n_cols() const { return m_cols; }

    [[nodiscard]] const part_t& get_partition() const { return m_partition; }

private:
    BlockSparseMatrix(data_t&& d, const index_t& r, const index_t& c, const size_t& n_parts);

    data_t  m_data;
    index_t m_rows, m_cols;
    part_t  m_partition;
};

template < class T >
BlockSparseMatrix< T > BlockSparseMatrix< T >::assembleBSMatrix(data_t&&       data,
                                                                const index_t& r,
                                                                const index_t& c,
                                                                const size_t&  n_parts)
{
    return BlockSparseMatrix(std::forward< data_t >(data), r, c, n_parts);
}

template < class T >
void BlockSparseMatrix< T >::print_sparsity_pattern() const
{
    auto next     = m_data.cbegin();
    auto any_left = true;
    auto m        = 0u;
    auto n        = 0u;

    std::cout << "| ";
    for (size_t i = 0; i < m_cols; ++i)
        std::cout << "--";
    std::cout << " |\n| ";

    while (true)
    {
        if (any_left && (*next).first.first == m && (*next).first.second == n)
        {
            std::cout << "X ";
            if (++next == m_data.cend())
                any_left = false;
        }
        else
        {
            std::cout << "  ";
        }

        if (++n == m_cols)
        {
            std::cout << " |\n| ";
            n = 0;
            ++m;
            if (m == m_rows)
                break;
        }
    }

    for (size_t i = 0; i < m_cols; ++i)
        std::cout << "--";
    std::cout << " |\n";
}

template < class T >
void BlockSparseMatrix< T >::print_to_file(const char* file_name) const
{
    std::fstream fs(file_name, std::ios_base::out);

    auto next     = m_data.cbegin();
    auto any_left = true;
    auto m        = 0u;
    auto n        = 0u;
    auto b_ind    = 0u;

    while (true)
    {
        if (any_left && (*next).first.first == m && (*next).first.second == n)
        {
            fs << (*next).second(b_ind, 0) << ' ' << (*next).second(b_ind, 1) << ' ';
            if (++next == m_data.cend())
                any_left = false;
        }
        else
        {
            fs << "0. 0. ";
        }

        if (++n == m_cols)
        {
            fs << '\n';
            n = 0;
            if (b_ind == 0)
            {
                ++b_ind;
                while (next != m_data.cbegin() && (*std::prev(next)).first.first == m)
                {
                    --next;
                    any_left = true;
                }
            }
            else
            {
                ++m;

                if (m == m_rows)
                    break;

                b_ind = 0;
            }
        }
    }
}

template < class T >
void BlockSparseMatrix< T >::print() const
{
    std::for_each(m_data.cbegin(), m_data.cend(), [](const entry_t& entry) {
        std::cout << "( " << entry.first.first << ", " << entry.first.second << " ):\n"
                  << entry.second << '\n';
    });
}

template < class T >
BlockSparseMatrix< T > BlockSparseMatrix< T >::assembleIdentBSMatrix(const size_t size)
{
    data_t entries;
    entries.reserve(size);
    const auto eye_block = []() {
        T block{};
        for (size_t i = 0; i < T::size_v; ++i)
        {
            block(i, i) = 1.;
        }
        return block;
    }();

    for (size_t i = 0; i < size; ++i)
        entries.emplace_back(std::make_pair(i, i), eye_block);

    return BlockSparseMatrix(std::move(entries), size, size);
}

template < class T >
BlockSparseMatrix< T >::BlockSparseMatrix(BlockSparseMatrix::data_t&&       d,
                                          const BlockSparseMatrix::index_t& r,
                                          const BlockSparseMatrix::index_t& c,
                                          const size_t&                     n_parts)
    : m_data(std::forward< BlockSparseMatrix::data_t >(d)), m_rows(r), m_cols(c)
{
    // Partition
    m_partition.reserve(n_parts + 1);
    m_partition.push_back(m_data.cbegin());

    const index_t row_block_size = r / n_parts;
    const index_t mod            = r % n_parts;

    for (size_t i = 1; i < n_parts; ++i)
    {
        const index_t s_val = i * row_block_size - 1 - static_cast< index_t >(i <= mod);
        m_partition.push_back(std::lower_bound(
            m_data.cbegin(), m_data.cend(), s_val, [](const entry_t& b, const index_t& ind) {
                return b.first.first < ind;
            }));
    }

    m_partition.push_back(m_data.cend());
}

// Sparse matrix-vector multiplication
template < class T2D, class T1D >
std::vector< T1D > operator*(const BlockSparseMatrix< T2D >& mat, const std::vector< T1D >& vec)
{
    auto mat_iter = mat.begin();

    // Assume 0-initialization for blocks
    auto ret_vec = std::vector< T1D >(vec.size());

    while (mat_iter != mat.end())
    {
        ret_vec[(*mat_iter).first.first] += (((*mat_iter).second) * vec[(*mat_iter).first.second]);
        ++mat_iter;
    }

    return ret_vec;
}

template < class T2D, class T1D >
void mat_vec_prod(const BlockSparseMatrix< T2D >&      mat,
                  const std::vector< T1D >&            vec,
                  std::vector< T1D >&                  ret_vec,
                  const std::shared_ptr< ThreadPool >& thread_pool)
{
    // Assume 0-initialization for blocks
    std::fill(ret_vec.begin(), ret_vec.end(), T1D{});

    const auto mult = [&](auto iter_begin, const auto& iter_end) {
        while (iter_begin != iter_end)
        {
            ret_vec[iter_begin->first.first] +=
                ((iter_begin->second) * vec[iter_begin->first.second]);
            ++iter_begin;
        }
    };

    const size_t tp_size = thread_pool->size();

    for (size_t i = 0; i < tp_size; ++i)
    {
        thread_pool->push_job(mult, mat.get_partition()[i], mat.get_partition()[i + 1]);
    }

    mult(mat.get_partition()[tp_size], mat.get_partition()[tp_size + 1]);

    thread_pool->complete();
}

// Rewrite local K as 4x4 block array
constexpr std::array< std::array< Block2D< double, 2 >, 4 >, 4 > makeBlockK()
{
    auto a = std::array< std::array< Block2D< double, 2 >, 4 >, 4 >{};

    for (size_t i_b = 0; i_b < 4; ++i_b)
    {
        for (size_t j_b = 0; j_b < 4; ++j_b)
        {
            for (size_t i = 0; i < 2; i++)
            {
                for (size_t j = 0; j < 2; j++)
                {
                    a[i_b][j_b].m_data[i][j] = K[2 * i_b + i][2 * j_b + j];
                }
            }
        }
    }

    return a;
}

BlockSparseMatrix< Block2D< double, 2 > > assembleStiffnessMatrix(const size_t n_part = 1)
{
    // Rewrite K to blocks
    constexpr auto K_block = makeBlockK();

    typename BlockSparseMatrix< Block2D< double, 2 > >::data_t entries;
    entries.reserve(mx * my * 16);

    // Lambda for calculating block DOF
    constexpr auto block_DOF = [](const size_t& row_ind, const size_t& col_ind, const size_t& i) {
        if (i < 2)
            return col_ind * (mx + 1) + row_ind + i;
        else if (i == 2)
            return (col_ind + 1) * (mx + 1) + row_ind + 1;
        else if (i == 3)
            return (col_ind + 1) * (mx + 1) + row_ind;
        else
            throw std::logic_error("bad dim");
    };

    for (size_t row_ind = 0; row_ind < mx; row_ind++)
    {
        for (size_t col_ind = 0; col_ind < my; col_ind++)
        {
            for (size_t i = 0; i < 4; i++)
            {
                const auto dof_i = block_DOF(row_ind, col_ind, i);

                for (size_t j = 0; j < 4; j++)
                {
                    const auto dof_j = block_DOF(row_ind, col_ind, j);

                    if (dof_i != dof_j && ((fix[2 * dof_i] && fix[2 * dof_i + 1]) ||
                                           (fix[2 * dof_j] && fix[2 * dof_j + 1])))
                        continue;

                    entries.emplace_back(std::make_pair(dof_i, dof_j), K_block[i][j]);

                    if (dof_i != dof_j)
                    {
                        if (fix[2 * dof_i])
                            entries.back().second.m_data[0].fill(0.);

                        if (fix[2 * dof_i + 1])
                            entries.back().second.m_data[1].fill(0.);

                        if (fix[2 * dof_j])
                        {
                            entries.back().second.m_data[0][0] = 0.;
                            entries.back().second.m_data[1][0] = 0.;
                        }

                        if (fix[2 * dof_j + 1])
                        {
                            entries.back().second.m_data[0][1] = 0.;
                            entries.back().second.m_data[1][1] = 0.;
                        }
                    }
                }
            }
        }
    }

    std::sort(std::execution::par,
              entries.begin(),
              entries.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    entries.erase(accumulate_consecutive(
                      entries.begin(),
                      entries.end(),
                      [](const auto& a, const auto& b) { return a.first == b.first; },
                      [](const auto& a, const auto& b) {
                          auto ret = a;
                          ret.second += b.second;
                          return ret;
                      }),
                  entries.end());

    entries.shrink_to_fit();

    std::for_each(entries.begin(), entries.end(), [](auto& entry) {
        if (entry.first.first == entry.first.second)
        {
            if (fix[2 * entry.first.first])
            {
                entry.second.m_data[0][0] = 1.;
                entry.second.m_data[0][1] = 0.;
                entry.second.m_data[1][0] = 0.;
            }

            if (fix[2 * entry.first.first + 1])
            {
                entry.second.m_data[0][1] = 0.;
                entry.second.m_data[1][1] = 1.;
                entry.second.m_data[1][0] = 0.;
            }
        }
    });

    return BlockSparseMatrix< Block2D< double, 2 > >::assembleBSMatrix(
        std::move(entries), (mx + 1) * (my + 1), (mx + 1) * (my + 1), n_part);
}

// conjugate gradient
template < class T2D, class T1D >
std::vector< T1D > cg(const BlockSparseMatrix< T2D >&      A,
                      const std::vector< T1D >&            b,
                      const std::vector< T1D >&            x0,
                      const std::shared_ptr< ThreadPool >& thread_pool,
                      const size_t                         maxit = 10000,
                      const double                         tol   = 1e-6)
{
    const auto minus = [](const std::vector< T1D >& a, const std::vector< T1D >& b) {
        std::vector< T1D > sum(a.size());
        std::transform(
            a.cbegin(), a.cend(), b.cbegin(), sum.begin(), [](const T1D& b1, const T1D& b2) {
                return b1 - b2;
            });
        return sum;
    };

    const auto times = [](const std::vector< T1D >& a, const std::vector< T1D >& b) {
        return std::transform_reduce(a.cbegin(), a.cend(), b.cbegin(), 0.);
    };

    const auto increment_by_scaled_v =
        [](std::vector< T1D >& a, const std::vector< T1D >& b, const double& s) {
            std::transform(
                a.cbegin(), a.cend(), b.cbegin(), a.begin(), [&s](const T1D& b_a, const T1D& b_b) {
                    return b_b * s + b_a;
                });
        };

    const auto residual = [&minus](const BlockSparseMatrix< T2D >& A,
                                   const std::vector< T1D >&       b,
                                   const std::vector< T1D >&       x) {
        return minus(b, A * x);
    };

    // Assume 0-initialized blocks
    auto x  = x0;
    auto r  = residual(A, b, x0);
    auto p  = r;
    auto Ap = std::vector< T1D >(x0.size());

    double rxr   = times(r, r);
    double rxr_p = rxr;

    size_t iter;
    for (iter = 0; iter < maxit; ++iter)
    {
        // Ap = A * p;
        mat_vec_prod(A, p, Ap, thread_pool);

        const double alpha = rxr / times(p, Ap);
        increment_by_scaled_v(x, p, alpha);
        increment_by_scaled_v(r, Ap, -alpha);

        rxr = times(r, r);

        if (sqrt(rxr) < tol)
            break;

        std::transform(
            p.cbegin(), p.cend(), r.cbegin(), p.begin(), [s = rxr / rxr_p](auto p, auto r) {
                return p * s + r;
            });
        rxr_p = rxr;
    }

    // std::cout << "Performed " << (iter < maxit ? iter + 1 : iter) << " iterations\n";
    // std::cout << "CG returned residual: " << sqrt(rxr) << '\n';
    return x;
}

void solve(double* wynik, const size_t nt)
{
    const auto                          K_global = assembleStiffnessMatrix(nt);
    std::vector< Block1D< double, 2 > > F_block{static_cast< size_t >(K_global.n_rows())};
    for (size_t i = 0; i < (mx + 1) * (my + 1); ++i)
    {
        F_block[i](0) = F[i * 2];
        F_block[i](1) = F[i * 2 + 1];
    }

    const auto d = cg(K_global,
                      F_block,
                      decltype(F_block){F_block.size()},
                      std::make_shared< ThreadPool >(nt - 1));

    for (size_t i = 0; const auto& d_b : d)
    {
        wynik[i * 2]     = d_b(0);
        wynik[i * 2 + 1] = d_b(1);
        ++i;
    }
    int dummy = 0;
}
