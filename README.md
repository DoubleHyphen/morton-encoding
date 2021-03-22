# Morton encoding (Z-order encoding)

## Introduction
The `morton-encoding` crate offers convenience functions for transforming arrays of primitive unsigned integers to Morton keys and back, via the eponymous encoding process, which basically groups all same-order bits together. This helps linearise data points while preserving some measure of locality.

This crate was originally built with fractal analysis in mind. Nevertheless, Morton encoding can also be used for other use-cases, such as the efficient searching of data-bases.

## Usage
The `morton_encode` and `morton_decode` functions ought to be sufficient for most use-cases. They are used as follows:

/// ```rust
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

For more detailed information, as well as information on other similar crates, please look at the [documentation](https://docs.rs/morton-encoding/2.0.0/morton_encoding/).

## Advantages
The Morton encoding can be computed very efficiently and branchlessly. For use-cases where one needs to keep splitting the available space recursively into halves, it's unmatched.

## Disadvantages
Morton encoding in general isn't very good at preserving locality, at least when compared to eg Hilbert encoding. Furthermore, the present crate only implements it for cases where each coordinate has the same amount of significant bits; for cases where that's not true, the user is urged to use the `lindel` crate, which will hopefully come out soon:tm:.
