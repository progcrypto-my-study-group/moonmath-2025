# Elliptic Curve Performance Benchmark Results

This document summarizes the performance benchmarks of group law operations for three different elliptic curve implementations:
- **Weierstrass curves** (Short Weierstrass form)
- **Montgomery curves** 
- **Twisted Edwards curves**

## Benchmark Configuration

- **Prime field**: F₁₃ (13)
- **Curve parameters**:
  - Weierstrass: y² = x³ + x + 1 (mod 13)
  - Montgomery: By² = x³ + Ax² + x (mod 13) with A=1, B=1
  - Twisted Edwards: ax² + y² = 1 + dx²y² (mod 13) with a=1, d=2
- **Benchmark tool**: Criterion.rs with optimized release builds

## Results Summary

### Point Addition Performance

| Curve Type | Average Time | Performance Rank |
|------------|-------------|------------------|
| **Weierstrass** | **2.13 ns** | **🥇 1st (FASTEST)** |
| Montgomery | 8.08 ns | 🥈 2nd |
| Twisted Edwards | 8.66 ns | 🥉 3rd |

**Winner: Weierstrass** - Nearly **4x faster** than Montgomery and Twisted Edwards!

### Point Doubling Performance

| Curve Type | Average Time | Performance Rank |
|------------|-------------|------------------|
| **Montgomery** | **1.73 ns** | **🥇 1st (FASTEST)** |
| Twisted Edwards | 6.82 ns | 🥈 2nd |
| Weierstrass | 8.00 ns | 🥉 3rd |

**Winner: Montgomery** - **4x faster** than Weierstrass and Twisted Edwards!

### Scalar Multiplication Performance

#### Scalar = 3
| Curve Type | Average Time | Performance Rank |
|------------|-------------|------------------|
| **Montgomery** | **7.60 ns** | **🥇 1st (FASTEST)** |
| Twisted Edwards | 12.04 ns | 🥈 2nd |
| Weierstrass | 45.68 ns | 🥉 3rd |

#### Scalar = 5
| Curve Type | Average Time | Performance Rank |
|------------|-------------|------------------|
| **Montgomery** | **9.76 ns** | **🥇 1st (FASTEST)** |
| Twisted Edwards | 14.23 ns | 🥈 2nd |
| Weierstrass | 80.42 ns | 🥉 3rd |

#### Scalar = 7
| Curve Type | Average Time | Performance Rank |
|------------|-------------|------------------|
| **Montgomery** | **9.70 ns** | **🥇 1st (FASTEST)** |
| Twisted Edwards | 14.57 ns | 🥈 2nd |
| Weierstrass | 80.57 ns | 🥉 3rd |

#### Scalar = 11
| Curve Type | Average Time | Performance Rank |
|------------|-------------|------------------|
| **Montgomery** | **11.88 ns** | **🥇 1st (FASTEST)** |
| Twisted Edwards | 16.53 ns | 🥈 2nd |
| Weierstrass | 109.05 ns | 🥉 3rd |

**Winner: Montgomery** - Consistently **6-9x faster** than Weierstrass for scalar multiplication!

### Comprehensive Mixed Operations

| Curve Type | Average Time | Performance Rank |
|------------|-------------|------------------|
| **Twisted Edwards** | **30.41 ns** | **🥇 1st (FASTEST)** |
| Montgomery | 34.13 ns | 🥈 2nd |
| Weierstrass | 89.93 ns | 🥉 3rd |

**Winner: Twisted Edwards** - Best overall performance for mixed operations!

## Analysis and Conclusions

### Overall Performance Ranking

1. **🏆 Montgomery Curves** - **BEST OVERALL PERFORMANCE**
   - ✅ **Fastest** point doubling (1.73 ns)
   - ✅ **Fastest** scalar multiplication (consistently 6-9x faster than Weierstrass)
   - ⚠️ Moderate point addition (8.08 ns)

2. **🥈 Twisted Edwards Curves** - **BALANCED PERFORMANCE**
   - ✅ **Fastest** comprehensive mixed operations (30.41 ns)
   - ✅ Consistent performance across operations
   - ⚠️ Moderate point addition and doubling

3. **🥉 Weierstrass Curves** - **SPECIALIZED FOR ADDITION**
   - ✅ **Fastest** point addition (2.13 ns)
   - ❌ **Slowest** scalar multiplication (6-9x slower than Montgomery)
   - ❌ Poor overall performance for most operations

### Key Insights

1. **Montgomery curves excel at scalar multiplication** - This makes them ideal for cryptographic operations that rely heavily on scalar multiplication (like ECDH key exchange).

2. **Weierstrass curves are optimized for single point additions** - They shine in scenarios requiring many individual point additions but struggle with repeated operations.

3. **Twisted Edwards curves provide balanced performance** - They offer good all-around performance and are excellent for applications requiring mixed operations.

4. **Performance differences are significant** - The fastest implementation can be 4-9x faster than the slowest for the same operation.

### Recommendations

- **For ECDH/ECDSA and scalar-heavy operations**: Use **Montgomery curves**
- **For mixed cryptographic operations**: Use **Twisted Edwards curves**  
- **For specialized point addition workloads**: Use **Weierstrass curves**
- **For general-purpose cryptography**: **Montgomery curves** offer the best overall performance

### Technical Notes

- All benchmarks use the same finite field F₁₃ for fair comparison
- Results may vary with different field sizes and curve parameters
- Montgomery curves' superiority in scalar multiplication is due to efficient ladder algorithms
- The benchmark uses optimized release builds with Rust's criterion framework

### How to Run Benchmarks

```bash
# Run all benchmarks
cargo bench

# View detailed HTML reports
open target/criterion/report/index.html
```

The HTML reports provide detailed performance graphs, statistical analysis, and comparison charts for all curve operations.