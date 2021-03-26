//! `morton_encoding` is a crate that helps interleave the bits of numbers,
//! (called “co-ordinates” in the source code)
//! thereby creating a so-called “Morton encoding” (also known as a
//! “Z-order encoding”) in whose final result (“Key”)
//! all bits of the same significance are together.
//! This is helpful for linearising the points 
//! of data structures whilst preserving
//! some measure of locality. Z-order encoding isn't as good at preserving
//! locality as Hilbert encoding, but it is much easier to compute.
//! (For Hilbert encoding please refer to the 
//! ([`lindel`](https://crates.io/crates/lindel)
//! crate instead.)
//!
//! # Usage
//! The user is cordially advised to look for solutions in the following order:
//! 1. By far the most ergonomic, performant, and panic-free usage
//! is to use the [`morton_encode`](fn@morton_encode) and
//! [`morton_decode`](fn@morton_decode) functions with arrays of
//! primitive unsigned integers, as follows:
//! ```
//! # use morton_encoding::*;
//! let input = [5u8, 4, 12, 129];
//! let output = morton_encode(input);
//! assert_eq!(output, 268447241);
//! let reassembled_input = morton_decode(output);
//! assert_eq!(input, reassembled_input);
//! ```
//! 
//! 2. If 128 bits are not enough, the user should use
//! [`lindel`](https://crates.io/crates/lindel),
//! which contains the `create_lineariseable_data_type` macro. This
//! will unfortunately require manual implementation of any trait that
//! the macro doesn't implement automatically, such as subtraction.
//! 
//! 3. If the co-ordinates come at the end of a long iterator chain of length
//! `<=8`, the user should use the `_generic` functions, taking care to not
//! compare Keys resulting from unequal-length iterators.
//! 
//! 4. Finally, if performance is of no concern at all, the user can make use
//! of the `arbitrary_bit_size` module or of the `_generic` functions for
//! arbitrary iterator length.
//!
//! # Contents
//!
//! The functions contained herein are broadly split into four categories:
//! * “Bloating functions”, that interleave a number's bits with zeroes;
//! * “Shrinking functions”, that perform the opposite operation;
//! * “Morton encoding functions”, that interleave some numbers'
//! bits with each other;
//! * “Morton decoding functions”, that perform the opposite operation.
//!
//! For more detailed information,
//! please refer to each individual function's documentation.
//!
//! # `std` vs `no_std`
//!
//! `std` is entirely optional for this crate, and is only useful if
//! the user needs to use `BigUint`s. Currently, the crate
//! (or its documentation, in any event) has currently been compiled with `std`
#![cfg_attr(
    feature = "std",
    doc = " **on**, meaning that it can interleave arbitrary amounts of arbitrarily large integers, as long as the user only runs it within an operating system."
)]
#![cfg_attr(
    not(feature = "std"),
    doc = " **off**, meaning that it can be used even in embedded environments, but will demand home-made `uint` types if 128 bits are not enough for the user."
)]
//!
//!
//! # Similar crates
//! A quick sift through crates.rs comes up with the following 3 similar crates:
//! * [`morton`](https://crates.io/crates/morton): Just as fast as ours, but only works for `[u16; 2] <-> u32`.
//! * [`bitwise`](https://crates.io/crates/bitwise):
//!     1. On my machine it doesn't even compile.
//!     2. It hasn't been implemented for more than 3 dimensions.
//!     3. Its portable implementation, when copy-pasted, is 3.5x
//!     slower than ours.
//!     4. Its non-portable implementation uses the
//!     `core::arch::x86_64::_pdep_u64` primitive, so we wrote some extra
//!     implementations that used it and benchmarked them. To our surprise,
//!     it was consistently slightly slower than our own code!
//! * [`zdex`](https://crates.io/crates/zdex):
//!     Despite not built with performance in mind, its speed
//!     is within a few orders of magnitude from ours, its results are almost
//!     correct, and the API ought to be usable with enough dedication and
//!     experimentation. It does, however, need `std` whereas ours doesn't.
//!
//! # Operation
//! This code works by splitting each coordinate's bits in half, and pushing
//! them apart; then each half is split into quarters, and pushed apart again;
//! and so on until each bit is `n-1` bits apart from its neighbours. From
//! there, the complete key is computed by a mere shift-and-bitwise-or.
//! A simplified version of this algorithm (and, in fact, the inspiration
//! for this crate) can be found [here](https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN).
//!
//! # Performance
//! This crate's performace depends *crucially* on three things:
//! 1. The CPU's ability to shift numbers by more than one bit at a time
//! (which some microcontrollers can't do)
//! 2. The compiler's ability to utilise compile-time information in order to
//! [fold all the constants](https://en.wikipedia.org/wiki/Constant_folding)
//! 3. The result of the operations fitting within a single register.
//!
//! Assuming, then, that the compiler *can* optimise to the best of its ability,
//! the code is reasonably fast, at least for an operation whose hardware
//! implementation only needs wires. For the transformation between
//! `[u{N}; M] <-> u{N*M}`, the time needed is `O(M*log(N))`.
//! The program *generally* needs about 20 machine instructions per coordinate,
//! though this can certainly vary. <details>
//! <summary>A full list of the machine instructions
//! needed for each transformation can be found below; each column lists
//! the operation first, then the total amount of machine instructions,
//! then the amount of machine instructions per co-ordinate rounded up.</summary>
//!
//! ```should_not_compile
//!    u16   -> [u8; 2] :  23     12
//! [u8; 2]  ->   u16   :  27     14
//!    u32   -> [u16; 2]:  28     14
//!    u64   -> [u32; 2]:  32     16
//! [u16; 2] ->   u32   :  34     17
//! [u8; 3]  ->   u32   :  41     14
//!    u32   -> [u8;  3]:  48     16
//! [u32; 2] ->   u64   :  48     24
//! [u8; 4]  ->   u32   :  55     14
//!    u32   -> [u8;  4]:  55     14
//! [u16; 3] ->   u64   :  56     19
//!    u64   -> [u16; 3]:  62     21
//! [u8; 5]  ->   u64   :  64     13
//!    u64   -> [u16; 4]:  68     17
//! [u16; 4] ->   u64   :  71     18
//!    u64   -> [u8;  5]:  72     15
//! [u8; 6]  ->   u64   :  79     14
//!    u64   -> [u8;  6]:  87     15
//! [u64; 2] ->   u128  :  93     47
//! [u8; 7]  ->   u64   :  94     14
//!    u64   -> [u8;  7]: 102     15
//! [u8; 8]  ->   u64   : 102     13
//!   u128   -> [u64; 2]: 110     55
//!   u128   -> [u16; 5]: 131     27
//!    u64   -> [u8;  8]: 132     17
//! [u32; 3] ->   u128  : 142     48
//!   u128   -> [u32; 3]: 142     48
//! [u8; 9]  ->   u128  : 149     17
//! [u32; 4] ->   u128  : 152     38
//!   u128   -> [u32; 4]: 153     39
//! [u16; 5] ->   u128  : 159     32
//!   u128   -> [u16; 6]: 185     31
//! [u8; 10] ->   u128  : 190     19
//!   u128   -> [u8; 9]:  194     22
//! [u16; 6] ->   u128  : 213     36
//!   u128   -> [u8; 10]: 216     22
//!   u128   -> [u16; 7]: 217     31
//! [u16; 7] ->   u128  : 230     33
//! [u8; 11] ->   u128  : 244     23
//!   u128   -> [u8; 11]: 247     23
//! [u16; 8] ->   u128  : 253     32
//!   u128   -> [u16; 8]: 255     32
//!   u128   -> [u8; 12]: 266     23
//! [u8; 12] ->   u128  : 268     23
//! [u8; 13] ->   u128  : 285     22
//!   u128   -> [u8; 13]: 287     23
//!   u128   -> [u8; 14]: 308     22
//! [u8; 14] ->   u128  : 318     23
//!   u128   -> [u8; 15]: 332     23
//! [u8; 15] ->   u128  : 361     25
//!   u128   -> [u8; 16]: 382     24
//! [u8; 16] ->   u128  : 419     27
//! ```
//! As we can see, anything that doesn't involve `u128`s generally takes fewer
//! than 20 machine code instructions per co-ordinate; the actual spread
//! ranges between 12 and 24 instructions. For calculations that include
//! `u128`s and above, the instructions needed increase linearly with the
//! bit length.
//! </details>

#![cfg_attr(not(feature = "std"), no_std)]

type SizeRatio = core::num::NonZeroUsize;

use core::convert::From;
use core::ops::BitAndAssign;
use core::ops::BitOrAssign;
use core::ops::ShlAssign;
use num::traits::int::PrimInt;
use num_traits::ToPrimitive;

macro_rules! sizeof {
    ($t: ty) => {
        core::mem::size_of::<$t>()
    };
}

/// A convenience function.
/// 
/// Persuades the compiler that a non-zero value is, in fact, non-zero.
///
/// Will (unsurprisingly) panic if called with a zero.
pub fn nz(x: usize) -> SizeRatio {
    SizeRatio::new(x).unwrap()
}

fn fits_n_times<A, B>(n: SizeRatio) -> bool {
    sizeof!(A) >= n.get() * sizeof!(B)
}

/// An 8-bit co-ordinate requires 3 steps, regardless of the size ratio.
/// For every doubling of its bits, another step is required.
fn ceil_of_log2_of_coor_bits<Coor>() -> u32 {
    (sizeof!(Coor) * 8).next_power_of_two().trailing_zeros()
}

/// Let `x = 1<<step_number`, and `y = x * (siz_rat - 1)`.
/// Essentially, this function outputs a number that has `y` zeroes,
/// then `x` ones, then `y` zeroes again, then `x` ones again,
/// repeating that pattern for the whole number.
/// In case the masking isn't really needed,
/// all ones are returned instead so the compiler understands to
/// not actually perform the masking.
///
/// The “bit twiddling hacks” web-site calls these “magic numbers”,
/// but they're essentially just what I mentioned.
/// https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
///
/// Side-note: You really, really want this function to be
/// optimised away.
fn get_mask<Key>(step_number: usize, siz_rat: SizeRatio) -> Key
where
    Key: PrimInt + BitOrAssign + ShlAssign<usize>,
{
    let siz_rat = siz_rat.get();
    let tentative_mask = {
        let key_bits = sizeof!(Key) * 8;
        let one_1 = Key::one();
        let amount_of_1s_per_pattern = 1 << step_number;
        let pattern_length = siz_rat * amount_of_1s_per_pattern;
        let pattern = (one_1 << amount_of_1s_per_pattern) - one_1;
        let mut insert_patterns_here: Key = pattern;
        let amt_of_patterns_that_fit_in_key = key_bits / pattern_length;
        let amt_of_patterns_we_still_need_to_insert = amt_of_patterns_that_fit_in_key - 1;
        for _ in 0..amt_of_patterns_we_still_need_to_insert {
            insert_patterns_here <<= pattern_length;
            insert_patterns_here |= pattern;
        }
        insert_patterns_here
    };

    let masking_is_necessary = {
        let usize_bits = 8 * (sizeof!(usize) as u32);
        let floor_of_log2_of_siz_rat = usize_bits - siz_rat.leading_zeros() - 1;
        ((step_number as u32) % floor_of_log2_of_siz_rat) == 0
    };
    if masking_is_necessary {
        tentative_mask
    } else {
        Key::max_value()
    }
}

/// Guarantees that a given data structure is suitable for a Morton Key.
///
/// Essentially, it guarantees that a) it behaves as if it's an integer
/// and b) that it's larger than the relevant coordinate.
///
/// Ideally, we would also like a `const N: usize` here that would
/// assert that a `Key` value is at least `N` times larger than a `Coor`
/// value, but as of current (1.50) Rust this is not possible.
pub trait ValidKey<Coor>:
    PrimInt + From<Coor> + BitOrAssign + BitAndAssign + ShlAssign<usize>
{
}
/// All primitive int types are suitable for Morton Keys.
impl<Coor, Key: PrimInt + From<Coor> + BitOrAssign + BitAndAssign + ShlAssign<usize>> ValidKey<Coor>
    for Key
{
}

/// Given a data type and a number, it yields the smallest unsigned
/// integer that's at least `N` times larger.
/// 
/// Implemented by brute force.
pub trait IdealKey<const N: usize> {
    type Key;
}

// UGLY WORK-AROUND, HOOOOOOOO!
impl IdealKey<1> for u8 {
    type Key = u8;
}
impl IdealKey<2> for u8 {
    type Key = u16;
}
impl IdealKey<3> for u8 {
    type Key = u32;
}
impl IdealKey<4> for u8 {
    type Key = u32;
}
impl IdealKey<5> for u8 {
    type Key = u64;
}
impl IdealKey<6> for u8 {
    type Key = u64;
}
impl IdealKey<7> for u8 {
    type Key = u64;
}
impl IdealKey<8> for u8 {
    type Key = u64;
}
impl IdealKey<9> for u8 {
    type Key = u128;
}
impl IdealKey<10> for u8 {
    type Key = u128;
}
impl IdealKey<11> for u8 {
    type Key = u128;
}
impl IdealKey<12> for u8 {
    type Key = u128;
}
impl IdealKey<13> for u8 {
    type Key = u128;
}
impl IdealKey<14> for u8 {
    type Key = u128;
}
impl IdealKey<15> for u8 {
    type Key = u128;
}
impl IdealKey<16> for u8 {
    type Key = u128;
}

impl IdealKey<1> for u16 {
    type Key = u16;
}
impl IdealKey<2> for u16 {
    type Key = u32;
}
impl IdealKey<3> for u16 {
    type Key = u64;
}
impl IdealKey<4> for u16 {
    type Key = u64;
}
impl IdealKey<5> for u16 {
    type Key = u128;
}
impl IdealKey<6> for u16 {
    type Key = u128;
}
impl IdealKey<7> for u16 {
    type Key = u128;
}
impl IdealKey<8> for u16 {
    type Key = u128;
}

impl IdealKey<1> for u32 {
    type Key = u32;
}
impl IdealKey<2> for u32 {
    type Key = u64;
}
impl IdealKey<3> for u32 {
    type Key = u128;
}
impl IdealKey<4> for u32 {
    type Key = u128;
}

impl IdealKey<1> for u64 {
    type Key = u64;
}
impl IdealKey<2> for u64 {
    type Key = u128;
}

impl IdealKey<1> for u128 {
    type Key = u128;
}

/// “Bloats” a given number by interleaving its bits with zeroes.
///
/// Each input bit is interleaved with `N - 1` zeroes.
///
/// # Examples
/// ```
/// # use morton_encoding::bloat_custom;
/// assert_eq!(bloat_custom::<_, u32, 3>(0x55u8), 0x041041);
/// ```
/// # Panics
///
/// With the advent of `min_const_generics` in this crate, it was
/// hoped that this function could recognise when a `Key` value
/// is not large enough for `N` `Coor` values, and fail to compile
/// in that case. Sadly, the compiler isn't quite clever enough for that yet.
/// Thus, we have to make do with run-time panics instead, as seen in the
/// following example:
/// ```should_panic
/// # use morton_encoding::bloat_custom;
/// bloat_custom::<u8, u16, 3>(0x55u8); // Compiles, but panics.
/// ```
/// Therefore, the user is solely responsible for maintaining the invariants
/// needed.
pub fn bloat_custom<Coor, Key, const N: usize>(x: Coor) -> Key
where
    Coor: ToPrimitive,
    Key: ValidKey<Coor>,
{
    //static_assertions::const_assert!(fitsntimes::<Key, Coor, {N}>());
    assert!(N > 0, "N must not equal zero!");
    assert!(
        sizeof!(Key) >= N * sizeof!(Coor),
        "Key value is not big enough for {} Coordinate values!",
        N
    );
    bloat_custom_checked(x, nz(N)).unwrap()
}

/// Fallibly “bloats” a given number, interleaving its bits with zeroes.
/// 
/// Returns an `Option`.
///
/// Each input bit is interleaved with `siz_rat - 1` zeroes.
/// Alternatively, if the provided `Key` type is too small to
/// fit `siz_rat` `Coor` numbers inside it, a `None` is returned.
///
/// This function is usable when the size ratio is computed at run-time,
/// but bear in mind that it might run much more slowly than its
/// compile-time counterparts.
///  
/// # Examples
/// ```
/// # use morton_encoding::bloat_custom_checked;
/// # use morton_encoding::nz;
/// assert_eq!(bloat_custom_checked(0xFFFFu16, nz(2)), Some(0x55555555u32));
/// assert_eq!(bloat_custom_checked::<u16, u32>(0xFFFF, nz(3)), None);
/// ```
pub fn bloat_custom_checked<Coor, Key>(x: Coor, siz_rat: SizeRatio) -> Option<Key>
where
    Coor: ToPrimitive,
    Key: ValidKey<Coor>,
{
    if !fits_n_times::<Key, Coor>(siz_rat) {
        return None;
    }

    let mut result: Key = From::from(x);
    if siz_rat.get() == 1 {
        return Some(result);
    }

    let shift_bitor_mask = |x: &mut Key, y| {
        let shift_amt = (siz_rat.get() - 1) << y;
        *x |= *x << shift_amt;
        (*x) &= get_mask(y, siz_rat)
    };

    let step_amount = ceil_of_log2_of_coor_bits::<Coor>();

    (0..step_amount)
        .rev()
        .for_each(|q| shift_bitor_mask(&mut result, q as usize));

    Some(result)
}

/// “Shrinks” a given number, only keeping a user-provided amount of its bits.
///
/// The 0th bit is always kept, and so is every `siz_rat`th bit thereafter. All
/// other bits are lost.
///
/// This function sanitises its input, so the bits that are thrown away do not
/// need to be cleared by the user.
///
/// # Examples
/// ```
/// # use morton_encoding::shrink_custom;
/// assert_eq!(shrink_custom::<u8, u32, 3>(0x145145), 0x55);
/// ```
///
/// # Panics
///
/// With the advent of min_const_generics in this crate, it was
/// hoped that this function could recognise when a `Key` value
/// is not large enough for `N` `Coor` values, and fail to compile
/// in that case. Sadly, the compiler isn't quite clever enough for that yet.
/// Thus, we have to make do with run-time panics instead, as seen in the
/// following example:
/// ```should_panic
/// # use morton_encoding::shrink_custom;
/// shrink_custom::<u8, u16, 3>(0x5555u16); // Compiles, but panics.
/// ```
/// Therefore, the user is solely responsible for maintaining the invariants
/// needed.
///
pub fn shrink_custom<Coor, Key, const N: usize>(x: Key) -> Coor
where
    Coor: ToPrimitive + PrimInt,
    Key: ValidKey<Coor>,
{
    if N == 0 || sizeof!(Key) < N * sizeof!(Coor) {
        panic!();
    }
    shrink_custom_checked(x, nz(N)).unwrap()
}

/// Fallibly “shrinks” a given number, only keeping a certain amount of its bits.
///
/// The 0th bit is always kept, and so is every `siz_rat`th bit thereafter. All
/// other bits are lost. Alternatively, if the provided `Key` type is too small
/// to fit `siz_rat` `Coor` numbers inside it, a `None` is returned .
///
/// This function is usable when the size ratio is computed at run-time,
/// but bear in mind that it might run much more slowly than its
/// compile-time counterparts.
///
/// This function sanitises its input, so the bits that are thrown away do not
/// need to be cleared by the user.
///
/// # Examples
/// ```
/// # use morton_encoding::shrink_custom_checked;
/// # use morton_encoding::nz;
/// assert_eq!(shrink_custom_checked::<u16, u32>(0xFFFF_FFFF, nz(2)), Some(0xFFFF));
/// assert_eq!(shrink_custom_checked::<u16, u32>(0x185185, nz(3)), None);
/// ```
pub fn shrink_custom_checked<Coor, Key>(x: Key, siz_rat: SizeRatio) -> Option<Coor>
where
    Coor: ToPrimitive + PrimInt,
    Key: ValidKey<Coor>,
{
    if !fits_n_times::<Key, Coor>(siz_rat) {
        return None;
    }

    let mut result: Key = x;
    if siz_rat.get() == 1 {
        return Coor::from(result);
    }

    let shift_bitor_mask = |x: &mut Key, y| {
        (*x) &= get_mask(y, siz_rat);
        let shift_amt = (siz_rat.get() - 1) << y;
        *x |= *x >> shift_amt;
    };

    let step_amount = ceil_of_log2_of_coor_bits::<Coor>();

    (0..step_amount).for_each(|q| shift_bitor_mask(&mut result, q as usize));

    Coor::from(result & (From::from(Coor::max_value())))
}

/// Receives an iterator of `Coor` values and encodes them all in a `Key` value.
///
/// Meant to be used with statically-known iterator lengths. If that is not
/// the case, the checked version of this function (`morton_encode_checked`)
/// should be used instead.
///
/// # Panics
///
/// This function will panic if the `Key` value provided does not have enough
/// length for all the values in the iterator.
///
/// In the future, we would like to constrain it for iterators that have exactly
/// `const N` elements, but that's still quite far away.
///
/// # Examples
/// ```
/// # use morton_encoding::morton_encode_generic;
/// let input = vec!(0u8, 255);
/// let result: u16 = morton_encode_generic(input);
/// assert_eq!(result, 0x5555);
/// ```
/// In the event that one `Key` is strictly larger than all `Coor`s put together, all significant bits are gathered at the end, and extraneous bits at the beginning are cleared. The following example illustrates the difference:
/// ```
/// # use morton_encoding::morton_encode_generic;
/// let input = vec!(255u8; 3);
/// let result: u32 = morton_encode_generic(input);
/// assert_eq!(result, 0b_0000_0000_111_111_111_111_111_111_111_111u32);
/// let input = vec!(0u8, 255, 255, 255);
/// let result: u32 = morton_encode_generic(input);
/// assert_eq!(result, 0b_0111_0111_0111_0111_0111_0111_0111_0111u32);
/// ```
pub fn morton_encode_generic<Coor, Key, Coors>(coors: Coors) -> Key
where
    Coor: ToPrimitive,
    Key: ValidKey<Coor>,
    Coors: IntoIterator<Item = Coor>,
    <Coors as core::iter::IntoIterator>::IntoIter: ExactSizeIterator,
{
    morton_encode_generic_checked(coors).unwrap()
}

/// Receives an iterator of `Coor` values and encodes them all in a `Option<Key>` value.
///
/// Returns a `None` if the iterator is empty, or if it is too large for all the `Coor` values contained within to fit inside a `Key` type.
///
/// # Examples
/// ```
/// # use morton_encoding::morton_encode_generic_checked;
/// assert_eq!(morton_encode_generic_checked(vec!(1u8, 1)), Some(3u16));
/// assert_eq!(morton_encode_generic_checked::<u16, u32, Vec<u16>>(vec!()), None);
/// assert_eq!(morton_encode_generic_checked::<_, u16, _>(vec!(0u8, 255, 200)), None);
/// ```
pub fn morton_encode_generic_checked<Coor, Key, Coors>(coors: Coors) -> Option<Key>
where
    Coor: ToPrimitive,
    Key: ValidKey<Coor>,
    Coors: IntoIterator<Item = Coor>,
    <Coors as core::iter::IntoIterator>::IntoIter: ExactSizeIterator,
{
    let iterator = coors.into_iter();
    let siz_rat = SizeRatio::new(iterator.len())?;
    if !fits_n_times::<Key, Coor>(siz_rat) {
        None
    } else {
        Some({
            let bloat_fn = |x| match siz_rat.get() {
                1 => bloat_custom::<Coor, Key, 1>(x),
                2 => bloat_custom::<Coor, Key, 2>(x),
                3 => bloat_custom::<Coor, Key, 3>(x),
                4 => bloat_custom::<Coor, Key, 4>(x),
                5 => bloat_custom::<Coor, Key, 5>(x),
                6 => bloat_custom::<Coor, Key, 6>(x),
                7 => bloat_custom::<Coor, Key, 7>(x),
                8 => bloat_custom::<Coor, Key, 8>(x),
                _ => bloat_custom_checked::<Coor, Key>(x, siz_rat).unwrap(),
            };
            iterator
                .map(bloat_fn)
                .fold(Key::zero(), |acc, x| (acc << 1) | x)
        })
    }
}

/// Receives a `Key` value and returns an iterator of `Coor` values that were decoded from it.
///
/// Returns an empty iterator if the `Key` value is too small for so many `Coor` values.
///
/// # Examples
/// ```
/// # use morton_encoding::morton_decode_generic;
/// # use morton_encoding::morton_encode_generic;
/// # use morton_encoding::nz;
/// assert_eq!(morton_decode_generic::<_, _, Vec<u8>>(3u16, nz(2)), vec!(1, 1));
/// assert_eq!(morton_decode_generic::<u16, u32, Vec<u16>>(0, nz(4)), vec!());
/// let input = vec!(1u8, 2);
/// let encoded_input: u16 = morton_encode_generic(input.clone());
/// let reassembled_input: Vec<u8> = morton_decode_generic(encoded_input, nz(2));
/// assert_eq!(input, reassembled_input);
/// ```
pub fn morton_decode_generic<Coor, Key, Coors>(key: Key, siz_rat: SizeRatio) -> Coors
where
    Coor: ToPrimitive + PrimInt,
    Key: ValidKey<Coor>,
    Coors: core::iter::FromIterator<Coor>,
{
    if !fits_n_times::<Key, Coor>(siz_rat) {
        core::iter::empty().collect::<Coors>()
    } else {
        let shrink_fn = |x| match siz_rat.get() {
            1 => shrink_custom::<Coor, Key, 1>(x),
            2 => shrink_custom::<Coor, Key, 2>(x),
            3 => shrink_custom::<Coor, Key, 3>(x),
            4 => shrink_custom::<Coor, Key, 4>(x),
            5 => shrink_custom::<Coor, Key, 5>(x),
            6 => shrink_custom::<Coor, Key, 6>(x),
            7 => shrink_custom::<Coor, Key, 7>(x),
            8 => shrink_custom::<Coor, Key, 8>(x),
            _ => shrink_custom_checked::<Coor, Key>(x, siz_rat).unwrap(),
        };
        let siz_rat = siz_rat.get();
        (0..siz_rat)
            .map(|x| key >> (siz_rat - (x + 1)))
            .map(shrink_fn)
            .collect::<Coors>()
    }
}

/// Encodes an array of `N` `Coordinates` into a `Key`.
/// 
/// Receives an array of `N` values of `Coordinate` type (`[Coordinate; N]`),
/// and encodes them all in a `Key` value. In the event that one `Key` is
/// strictly larger than `N` `Coordinate`s, all significant bits are gathered
/// at the end, and extraneous bits at the beginning are cleared.
///
/// # Panics
///
/// This function will panic if the `Key` value provided does not have enough
/// length for all the values in the iterator.
///
/// In the future, we would like for it to fail to compile for invalid
/// parametres, but that's still quite far away.
/// ```should_panic
/// # use morton_encoding::morton_encode_array;
/// let _: u32 = morton_encode_array([15u16, 256, 10000]); // Compiles, but panics.
/// ```
/// # Examples
/// ```
/// # use morton_encoding::morton_encode_array;
/// let input = [255u8, 0];
/// let result: u16 = morton_encode_array(input);
/// assert_eq!(result, 0xAAAA);
/// ```
/// In the event that one `Key` is strictly larger than all `Coor`s put
/// together, all significant bits are gathered at the end, and extraneous
/// bits at the beginning are cleared. The following example illustrates
/// the difference:
/// ```
/// # use morton_encoding::morton_encode_array;
/// let input = [255u8; 3];
/// let result: u32 = morton_encode_array(input);
/// assert_eq!(result, 0b_0000_0000_111_111_111_111_111_111_111_111u32);
/// let input = [0u8, 255, 255, 255];
/// let result: u32 = morton_encode_array(input);
/// assert_eq!(result, 0b_0111_0111_0111_0111_0111_0111_0111_0111u32);
/// ```
pub fn morton_encode_array<Coordinate, Key, const N: usize>(input: [Coordinate; N]) -> Key
where
    Coordinate: ToPrimitive + PrimInt,
    Key: ValidKey<Coordinate>,
{
    input
        .iter()
        .map(|&m| bloat_custom::<Coordinate, Key, N>(m))
        .fold(Key::zero(), |acc, x| (acc << 1) | x)
}

/// The most ergonomic way to perform the Morton encoding operation.
/// 
/// Works for all primitive unsigned integers, and never panics for them.
/// Will work for others if the user implements the `IdealKey`
/// trait for them, but will panic if the trait is misimplemented.
///
/// # Examples
/// ```
/// # use morton_encoding::morton_encode;
/// # use morton_encoding::nz;
/// assert_eq!(morton_encode([0u8; 2]), 0u16);
/// assert_eq!(morton_encode([0u8; 3]), 0u32);
/// assert_eq!(morton_encode([0u32; 3]), 0u128);
/// ```
pub fn morton_encode<Coordinate, const N: usize>(
    input: [Coordinate; N],
) -> <Coordinate as IdealKey<N>>::Key
where
    Coordinate: IdealKey<N> + ToPrimitive + PrimInt,
    <Coordinate as IdealKey<N>>::Key: ValidKey<Coordinate>,
{
    morton_encode_array(input)
}

/// Receives a `Key` value and unscrambles it into an array.
/// 
/// Returns an array of `Coor` values that were decoded from the input.
///
/// Panics if the `Key` value is too small for so many `Coor` values.
/// ```should_panic
/// # use morton_encoding::morton_decode_array;
/// let _: [u16; 5] = morton_decode_array(6000000u64);
/// ```
///
/// # Examples
/// ```
/// # use morton_encoding::morton_decode_array;
/// # use morton_encoding::morton_encode_array;
/// # use morton_encoding::nz;
/// assert_eq!(morton_decode_array::<u8, u16, 2>(3u16), [1u8, 1]);
/// let input = [1u32, 2];
/// let encoded_input: u64 = morton_encode_array(input);
/// let reassembled_input: [u32; 2] = morton_decode_array(encoded_input);
/// assert_eq!(input, reassembled_input);
/// ```
pub fn morton_decode_array<Coordinate, Key, const N: usize>(input: Key) -> [Coordinate; N]
where
    Coordinate: ToPrimitive + PrimInt,
    Key: ValidKey<Coordinate>,
{
    let mut result = [Coordinate::zero(); N];
    for (i, n) in result.iter_mut().rev().enumerate() {
        *n = shrink_custom::<Coordinate, Key, N>(input >> i);
    }
    result
}

/// The most ergonomic way to perform the Morton decoding operation.
/// 
/// Works for all primitive unsigned integers, and never panics for them.
/// Will work for others if the user implements the `IdealKey`
/// trait for them, but will panic if the trait is misimplemented.
///
/// # Examples
/// ```
/// # use morton_encoding::morton_encode;
/// # use morton_encoding::morton_decode;
/// let input = 992;
/// let output: [u8; 5] = morton_decode(input);
/// assert_eq!(output, [2u8; 5]);
/// let input = [543u32, 23765];
/// let encoded_input: u64 = morton_encode(input);
/// let reassembled_input: [u32; 2] = morton_decode(encoded_input);
/// assert_eq!(input, reassembled_input);
/// ```
pub fn morton_decode<Coordinate, const N: usize>(
    input: <Coordinate as IdealKey<N>>::Key,
) -> [Coordinate; N]
where
    Coordinate: IdealKey<N> + ToPrimitive + PrimInt,
    <Coordinate as IdealKey<N>>::Key: ValidKey<Coordinate>,
{
    morton_decode_array(input)
}

#[cfg(feature = "std")]
#[cfg(test)]
mod tests {
    use super::*;

    fn fmt_bin<T>(x: T) -> std::string::String
    where
        T: std::fmt::Binary,
    {
        let size = sizeof!(T);
        if size == 16 {
            format!("{:#0130b}", x)
        } else if size == 8 {
            format!("{:#066b}", x)
        } else if size == 4 {
            format!("{:#034b}", x)
        } else if size == 2 {
            format!("{:#018b}", x)
        } else if size == 1 {
            format!("{:#010b}", x)
        } else {
            format!("{:#b}", x)
        }
    }

    fn bloat_slow_custom<Coor, Key>(x: Coor, siz_rat: SizeRatio) -> Key
    where
        Coor: ToPrimitive + PrimInt + std::fmt::Binary + core::ops::BitAnd,
        Key: ValidKey<Coor>,
        <Key as num_traits::Num>::FromStrRadixErr: std::fmt::Debug,
    {
        let coor_siz = sizeof!(Coor);
        let coor_bits = coor_siz * 8;
        let key_siz = sizeof!(Key);
        let siz_rat = siz_rat.get();
        if siz_rat == 1 {
            return core::convert::From::from(x);
        }
        assert!(siz_rat >= 1 && coor_siz * siz_rat <= key_siz);
        let key_zero = Key::zero();
        let coor_zero = Coor::zero();
        let mut result = key_zero;
        let get_mask_key = |b| (Key::one()) << (b * siz_rat);
        let get_mask_coor = |b| (Coor::one()) << b;
        let tst_bit = |a: Coor, b| (a & get_mask_coor(b)) != coor_zero;
        let set_bit = |a: &mut Key, b| (*a) |= get_mask_key(b);
        for bit_examined in 0..coor_bits {
            if tst_bit(x, bit_examined) {
                set_bit(&mut result, bit_examined)
            }
        }
        result
    }

    const MAX_BITS: u8 = 20;

    fn size_ratio<A, B>() -> SizeRatio
    where
        A: From<B>,
    {
        nz(sizeof!(A) / sizeof!(B))
        // The "From" trait ensures that the result is never zero; consequently, this function can never panic.
    }

    fn test_all_possible_values<Coor, Key>(siz_rat: Option<usize>)
    where
        Coor: num_traits::ToPrimitive
            + std::fmt::Binary
            + core::ops::BitAnd
            + num::traits::PrimInt
            + std::fmt::Display
            + std::fmt::Debug,
        rand::distributions::Standard: rand::distributions::Distribution<Coor>,
        Key: num::traits::int::PrimInt
            + core::convert::From<Coor>
            + core::ops::BitOrAssign
            + core::ops::BitAndAssign
            + std::fmt::Binary
            + core::ops::ShlAssign<usize>,
        <Key as num_traits::Num>::FromStrRadixErr: std::fmt::Debug,
        u128: core::convert::From<Coor>,
    {
        let siz_rat = siz_rat
            .and_then(|x| SizeRatio::new(x))
            .unwrap_or(size_ratio::<Key, Coor>());
        let testing_function = |x| bloat_slow_custom::<Coor, Key>(x, siz_rat);
        let function_under_test = |x| bloat_custom_checked::<Coor, Key>(x, siz_rat).unwrap();
        let limit = core::cmp::min(u128::from(Coor::max_value()), 1u128 << MAX_BITS);
        let limit = limit as usize;
        for x in 0..=limit {
            let x_ = Coor::from(x).expect("Coor and usize incompatible.");
            if x % (1 << 26) == 0 && x != 0 {
                dbg!(x >> 26);
            }
            let fn2 = function_under_test(x_);
            let fn1 = testing_function(x_);
            if fn1 != fn2 {
                panic!(
                    "x = {} \ntesting_function = {}, \nfunction_under_test = {}",
                    x,
                    fmt_bin(fn1),
                    fmt_bin(fn2)
                )
            }
            if x_ != shrink_custom_checked(fn2, siz_rat).unwrap() {
                panic!(
                    "x = {} \nBloated = {}, \nShrunk = {}",
                    x,
                    fmt_bin(fn2),
                    shrink_custom_checked(fn2, siz_rat).unwrap()
                )
            }
            if x % 16 == 0 {
                let input = (0..siz_rat.get())
                    .map(|_| rand::random::<Coor>())
                    .collect::<Vec<Coor>>();
                let encoded_input: Key = morton_encode_generic(input.clone());
                let reinstated_input: Vec<Coor> = morton_decode_generic(encoded_input, siz_rat);
                assert_eq!(input, reinstated_input);
            }
        }
    }

    #[test]
    fn test_u8_u16() {
        test_all_possible_values::<u8, u16>(None);
    }
    #[test]
    fn test_u8_u24() {
        test_all_possible_values::<u8, u32>(Some(3));
    }
    #[test]
    fn test_u8_u32() {
        test_all_possible_values::<u8, u32>(None);
    }
    #[test]
    fn test_u8_u40() {
        test_all_possible_values::<u8, u64>(Some(5));
    }
    #[test]
    fn test_u8_u48() {
        test_all_possible_values::<u8, u64>(Some(6));
    }
    #[test]
    fn test_u8_u56() {
        test_all_possible_values::<u8, u64>(Some(7));
    }
    #[test]
    fn test_u8_u64() {
        test_all_possible_values::<u8, u64>(None);
    }
    #[test]
    fn test_u8_u72() {
        test_all_possible_values::<u8, u128>(Some(9));
    }
    #[test]
    fn test_u8_u80() {
        test_all_possible_values::<u8, u128>(Some(10));
    }
    #[test]
    fn test_u8_u88() {
        test_all_possible_values::<u8, u128>(Some(11));
    }
    #[test]
    fn test_u8_u96() {
        test_all_possible_values::<u8, u128>(Some(12));
    }
    #[test]
    fn test_u8_u104() {
        test_all_possible_values::<u8, u128>(Some(13));
    }
    #[test]
    fn test_u8_u112() {
        test_all_possible_values::<u8, u128>(Some(14));
    }
    #[test]
    fn test_u8_u120() {
        test_all_possible_values::<u8, u128>(Some(15));
    }
    #[test]
    fn test_u8_u128() {
        test_all_possible_values::<u8, u128>(None);
    }
    #[test]
    fn test_u16_u32() {
        test_all_possible_values::<u16, u32>(None);
    }
    #[test]
    fn test_u16_u48() {
        test_all_possible_values::<u16, u64>(Some(3));
    }
    #[test]
    fn test_u16_u64() {
        test_all_possible_values::<u16, u64>(None);
    }
    #[test]
    fn test_u16_u80() {
        test_all_possible_values::<u16, u128>(Some(5));
    }
    #[test]
    fn test_u16_u96() {
        test_all_possible_values::<u16, u128>(Some(6));
    }
    #[test]
    fn test_u16_u112() {
        test_all_possible_values::<u16, u128>(Some(7));
    }
    #[test]
    fn test_u16_u128() {
        test_all_possible_values::<u16, u128>(None);
    }
    #[test]
    fn test_u32_u64() {
        test_all_possible_values::<u32, u64>(None);
    }
    #[test]
    fn test_u32_u96() {
        test_all_possible_values::<u32, u128>(Some(3));
    }
    #[test]
    fn test_u32_u128() {
        test_all_possible_values::<u32, u128>(None);
    }
    #[test]
    fn test_u64_u128() {
        test_all_possible_values::<u64, u128>(None);
    }

    macro_rules! test_functions {
        ($encod: ident, $decod: ident, $t1: ty, $t2: ty, $siz_rat: expr) => {
            let max_bits_for_this = (MAX_BITS - 3) as usize;
            let tmp_bit_limit = $siz_rat * sizeof!($t1) * 8;
            if tmp_bit_limit > max_bits_for_this {
                for _ in 0..(1 << max_bits_for_this) {
                    let mut input: [$t1; $siz_rat] = [0; $siz_rat];
                    for m in 0..$siz_rat {
                        input[m] = rand::random::<$t1>();
                    }
                    let tmp = morton_encode(input);
                    assert_eq!(input, morton_decode(tmp));
                }
            } else {
                let limit = ((1 as u128) << tmp_bit_limit).wrapping_sub(1);
                let limit = limit as $t2;
                for input in 0..=limit {
                    let tmp: [$t1; $siz_rat] = morton_decode(input);
                    assert_eq!(input, morton_encode(tmp));
                }
            }
        };
    }
    #[test]
    fn test_all_spec_funcs() {
        test_functions!(morton_encode_u8_2d, morton_decode_u8_2d, u8, u16, 2);
        test_functions!(morton_encode_u8_3d, morton_decode_u8_3d, u8, u32, 3);
        test_functions!(morton_encode_u8_4d, morton_decode_u8_4d, u8, u32, 4);
        test_functions!(morton_encode_u8_5d, morton_decode_u8_5d, u8, u64, 5);
        test_functions!(morton_encode_u8_6d, morton_decode_u8_6d, u8, u64, 6);
        test_functions!(morton_encode_u8_7d, morton_decode_u8_7d, u8, u64, 7);
        test_functions!(morton_encode_u8_8d, morton_decode_u8_8d, u8, u64, 8);
        test_functions!(morton_encode_u8_9d, morton_decode_u8_9d, u8, u128, 9);
        test_functions!(morton_encode_u8_10d, morton_decode_u8_10d, u8, u128, 10);
        test_functions!(morton_encode_u8_11d, morton_decode_u8_11d, u8, u128, 11);
        test_functions!(morton_encode_u8_12d, morton_decode_u8_12d, u8, u128, 12);
        test_functions!(morton_encode_u8_13d, morton_decode_u8_13d, u8, u128, 13);
        test_functions!(morton_encode_u8_14d, morton_decode_u8_14d, u8, u128, 14);
        test_functions!(morton_encode_u8_15d, morton_decode_u8_15d, u8, u128, 15);
        test_functions!(morton_encode_u8_16d, morton_decode_u8_16d, u8, u128, 16);
        test_functions!(morton_encode_u16_2d, morton_decode_u16_2d, u16, u32, 2);
        test_functions!(morton_encode_u16_3d, morton_decode_u16_3d, u16, u64, 3);
        test_functions!(morton_encode_u16_4d, morton_decode_u16_4d, u16, u64, 4);
        test_functions!(morton_encode_u16_5d, morton_decode_u16_5d, u16, u128, 5);
        test_functions!(morton_encode_u16_6d, morton_decode_u16_6d, u16, u128, 6);
        test_functions!(morton_encode_u16_7d, morton_decode_u16_7d, u16, u128, 7);
        test_functions!(morton_encode_u16_8d, morton_decode_u16_8d, u16, u128, 8);
        test_functions!(morton_encode_u32_2d, morton_decode_u32_2d, u32, u64, 2);
        test_functions!(morton_encode_u32_3d, morton_decode_u32_3d, u32, u128, 3);
        test_functions!(morton_encode_u32_4d, morton_decode_u32_4d, u32, u128, 4);
        test_functions!(morton_encode_u64_2d, morton_decode_u64_2d, u64, u128, 2);
    }
}

#[cfg(feature = "std")]
/// The most general, but least performant, ways to perform
/// Morton encoding and decoding. The functions contained herein work
/// correctly, but the author makes absolutely no guarantees about
/// how performant they are.
pub mod arbitrary_bit_size {
    use super::*;
    use num::BigUint;
    /// A convenience function for converting primitive values to BigUints.
    /// Will fail to compile if called with a literal that's not explicitly
    /// declared as unsigned.
    ///
    /// Unavailable if compiled with `no_std`.
    pub fn tobuint<Coor>(x: Coor) -> BigUint
    where
        Coor: num::Unsigned,
        num::BigUint: core::convert::From<Coor>,
    {
        BigUint::from(x)
    }

    fn get_mask_biguint(
        step_number: usize,
        siz_rat: SizeRatio,
        key_bits: usize,
    ) -> Option<num::BigUint> {
        use num_traits::identities::One;

        let siz_rat = siz_rat.get();
        let tentative_mask = {
            let one_1 = BigUint::one();
            let amount_of_1s_per_pattern = 1 << step_number;
            let pattern = (one_1 << amount_of_1s_per_pattern) - BigUint::one();
            let mut insert_patterns_here = pattern.clone();
            let pattern_length = siz_rat * amount_of_1s_per_pattern;
            let amt_of_patterns_that_fit_in_key = key_bits / pattern_length;
            let amt_of_patterns_we_still_need_to_insert = amt_of_patterns_that_fit_in_key - 1;
            for _ in 0..amt_of_patterns_we_still_need_to_insert {
                insert_patterns_here <<= pattern_length;
                insert_patterns_here |= pattern.clone();
            }
            insert_patterns_here
        };

        let masking_is_necessary = {
            let floor_of_log2_of_siz_rat =
                8 * (sizeof!(usize) as u32) - siz_rat.leading_zeros() - 1;
            (step_number as u32) % floor_of_log2_of_siz_rat == 0
        };
        if masking_is_necessary {
            Some(tentative_mask)
        } else {
            None
        }
    }

    /// “Bloats” a given number to an arbitrarily large BigUint.
    ///
    /// Each bit of the input is interleaved with `siz_rat - 1` zeroes.
    ///
    /// This function assumes that the user will not need to use numbers larger than
    ///`usize::max_value()` bytes in size.
    ///
    /// Unavailable if compiled with `no_std`.
    ///
    /// # Examples
    /// ```
    /// # use morton_encoding::arbitrary_bit_size::bloat_custom_biguint;
    /// # use morton_encoding::nz;
    /// # use morton_encoding::arbitrary_bit_size::tobuint;
    /// # use num::BigUint;
    /// assert_eq!(bloat_custom_biguint(1u32, nz(2)), tobuint(1u8));
    /// assert_eq!(bloat_custom_biguint(u128::max_value(), nz(32)), BigUint::new(vec!(1u32; 128)));
    /// ```
    pub fn bloat_custom_biguint<Coor>(x: Coor, siz_rat: SizeRatio) -> num::BigUint
    where
        Coor: num_traits::ToPrimitive,
        num::BigUint: core::convert::From<Coor>,
    {
        let coor_siz = sizeof!(Coor);
        let key_siz = coor_siz * siz_rat.get();

        let mut result = num::BigUint::from(x);
        if siz_rat.get() == 1 {
            return result;
        }

        let shift_bitor_mask = |x: &mut num::BigUint, y| {
            let shift_amt = (siz_rat.get() - 1) << y;
            *x |= ((*x).clone()) << shift_amt;
            if let Some(mask) = get_mask_biguint(y, siz_rat, key_siz * 8) {
                (*x) &= mask;
            }
        };

        let op_amt = coor_siz.next_power_of_two().trailing_zeros() + 3;
        (0..op_amt)
            .rev()
            .for_each(|q| shift_bitor_mask(&mut result, q as usize));
        result
    }

    /// “Shrinks” a given number from an arbitrarily large BigUint.
    ///
    /// The 0th bit is always kept, as is every `siz_rat`th bit thereafter.
    ///
    /// This function sanitises its input, such that the bits that are thrown away do not need to be cleared by the user. It is also assumed that the user will not need to use numbers larger than `usize::max_value()` bytes in size.
    ///
    /// Unavailable if compiled with `no_std`.
    ///
    /// # Examples
    /// ```
    /// # use morton_encoding::arbitrary_bit_size::shrink_custom_biguint;
    /// # use morton_encoding::nz;
    /// # use morton_encoding::arbitrary_bit_size::tobuint;
    /// assert_eq!(shrink_custom_biguint(tobuint(1u8), nz(2)), tobuint(1u8));
    /// assert_eq!(shrink_custom_biguint(num::BigUint::new(vec!(3u32; 128)), nz(32)), num::BigUint::new(vec!(u32::max_value(); 4)));
    /// ```
    pub fn shrink_custom_biguint(x: BigUint, siz_rat: SizeRatio) -> num::BigUint {
        use num_traits::identities::One;

        let key_bits = x.bits() + siz_rat.get() - 1;
        let coor_bits = key_bits / siz_rat.get();
        let key_bits = coor_bits * siz_rat.get();

        let mut result = x;
        if siz_rat.get() == 1 {
            return result;
        }

        let shift_bitor_mask = |x: &mut BigUint, y| {
            if let Some(mask) = get_mask_biguint(y, siz_rat, key_bits) {
                (*x) &= mask;
            }
            //dbg!(x);
            let shift_amt = (siz_rat.get() - 1) << y;
            *x |= ((*x).clone()) >> shift_amt;
        };

        let step_amount = coor_bits.next_power_of_two().trailing_zeros();
        (0..step_amount).for_each(|q| shift_bitor_mask(&mut result, q as usize));

        let final_mask = {
            let mut tmpmask = BigUint::one();
            tmpmask <<= coor_bits;
            tmpmask - BigUint::one()
        };

        result & final_mask
    }

    /// Receives an iterator of `Coor` values and encodes them all in a `BigUint`.
    ///
    /// Unavailable if compiled with `no_std`.
    ///
    /// # Examples
    /// ```
    /// # use morton_encoding::arbitrary_bit_size::morton_encode_biguint;
    /// # use morton_encoding::nz;
    /// # use morton_encoding::arbitrary_bit_size::tobuint;
    /// let input = vec!(0u8, 255u8);
    /// let result = morton_encode_biguint(input);
    /// assert_eq!(result, tobuint(0x5555u16));
    /// ```
    pub fn morton_encode_biguint<Coor, Coors>(coors: Coors) -> BigUint
    where
        Coor: ToPrimitive,
        num::BigUint: core::convert::From<Coor>,
        Coors: IntoIterator<Item = Coor>,
        <Coors as core::iter::IntoIterator>::IntoIter: ExactSizeIterator,
    {
        let zero = BigUint::new(vec![0u32]);
        let iterator = coors.into_iter();
        if let Some(siz_rat) = SizeRatio::new(iterator.len()) {
            let bloat_fn = |x: Coor| bloat_custom_biguint::<Coor>(x, siz_rat);
            iterator.map(bloat_fn).fold(zero, |acc, x| (acc << 1) | x)
        } else {
            zero
        }
    }

    /// Receives a `BigUint` value and returns an iterator of `BigUint` values that were decoded from it.
    ///
    /// Unavailable if compiled with `no_std`.
    ///
    /// # Examples
    /// ```
    /// # use morton_encoding::arbitrary_bit_size::morton_decode_biguint;
    /// # use morton_encoding::arbitrary_bit_size::morton_encode_biguint;
    /// # use morton_encoding::nz;
    /// # use morton_encoding::arbitrary_bit_size::tobuint;
    /// assert_eq!(morton_decode_biguint::<Vec<_>>(tobuint(3u8), nz(2)), vec!(tobuint(1u8); 2));
    /// let input = vec!(tobuint(1u8), tobuint(2u8));
    /// let encoded_input = morton_encode_biguint(input.clone());
    /// let reassembled_input: Vec<_> = morton_decode_biguint(encoded_input, nz(2));
    /// assert_eq!(input, reassembled_input);
    /// ```
    pub fn morton_decode_biguint<Coors>(key: BigUint, siz_rat: SizeRatio) -> Coors
    where
        Coors: core::iter::FromIterator<BigUint> + IntoIterator<Item = BigUint>,
        <Coors as core::iter::IntoIterator>::IntoIter: ExactSizeIterator,
    {
        let shrink_fn = |x: BigUint| shrink_custom_biguint(x, siz_rat);
        let siz_rat = siz_rat.get();
        (0..siz_rat)
            .map(|x| key.clone() >> (siz_rat - (x + 1)))
            .map(shrink_fn)
            .collect::<Coors>()
    }
}
