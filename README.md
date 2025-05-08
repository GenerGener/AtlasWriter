# AtlasWriter
Write transcripts with locus segmentation, alternate start sites, and polyA tails if you feel like it!

TODO: Fix bug in locus segments preceding locus segment 9. Segment 8 was missing.

Corrected segment 7 and segment 8:

```
>7
GTAGGATCTCTACAGTACTTGGCACTAGCAGCATTAATAAAACCAAAACAGATAAAGCCACCTTTGCCTAGTGTTAGGAAACTGACAGAGGACAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCATACAATGAATGGACACTAGAGCTTTTAGAGGAACTTAAGAGTGAAGCTGTTAGACATTTTCCTAGGATATGGCTCCATAACTTAGGACAACATATCTATGAAACTTACGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTCAG
>8
AATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATGACAAAAGCCTTAGGCATCTCCTATGGCAG
>9
GAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCA
```

# Documentation

AtlasWriter is a Python tool for extracting and manipulating HIV locus segments from FASTA files. It allows you to select specific genomic segments by index and supports sequence truncation at specified positions.

## Locus Segment Clarification

The segmentation of the HIV genome in this tool is based on standardized segment inclusion from deep long-read cDNA RNA-seq datasets. The script identifies locus segments based on the order they appear in the FASTA file. It assigns an index (1, 2, 3, etc.) to each FASTA record in sequential order. When using the tool, segment indices directly correspond to the sequence order in the input file, not to any predefined genomic feature names.

## Installation

No installation is needed beyond having Python 3.6 or newer. Simply download the `AtlasWriter.py` script to your local machine.

## Basic Usage

```bash
python AtlasWriter.py --input <fasta_file> --include_locus_segment <segments> [options]
```

## Command Line Arguments

| Argument | Description |
|----------|-------------|
| `--input` | Path to the input FASTA file containing locus segments |
| `--include_locus_segment` | Comma-separated list of segment indices or ranges to include |
| `--truncate` | Truncate sequence at specified positions |
| `--polyA` | Add a polyA tail after a specified position |
| `--alt_start` | Specify an alternative start site |
| `--output` | Output file path (if not specified, output is printed to stdout) |
| `--list` | List available segments and exit |

## Segment Selection

You can select segments to include using comma-separated indices, ranges, or a combination of both:

- Individual segments: `1,3,5,7`
- Range of segments: `1-5`
- Mixed format: `1-3,5,7-9`

## Alternative Start Site

The `--alt_start` parameter allows you to specify an alternative start site for the sequence, which is useful for modeling alternative transcription start sites or promoters:

### 1. Global Position Format (Recommended)

```
--alt_start /position
```

Where:
- The slash (`/`) comes before the position number to indicate "start from this position"
- `position` is the position in the combined sequence that will become the new first base (inclusive)

This format is more intuitive as it works with the absolute position in the final concatenated sequence, regardless of which segments are included.

### 2. Segment-Specific Format

```
--alt_start segment:/position
```

Where:
- `segment` is the segment index
- The slash (`/`) comes before the position number to indicate "start from this position"
- `position` is the position within that specific segment that will become the new first base (inclusive)

**Important Note**: When using the segment-specific format, the position must be within the valid range for that segment. For example, if segment 2 is 4169 bases long, the position must be between 1 and 4169.

Alternative start sites are useful for modeling different transcript isoforms, alternative promoters, or creating constructs with specific 5' ends.

## Truncation Options

The `--truncate` parameter allows you to cut the sequence at specific positions:

1. **Cutting after a specific base**: 
   ```
   --truncate segment:position/
   ```
   This keeps all bases up to and including the specified position and removes downstream bases. The slash (`/`) visually indicates the cutting point.

2. **Keeping a range between two positions**:
   ```
   --truncate upstream_segment:upstream_base-downstream_segment:downstream_base
   ```
   This keeps only the sequence between the two specified positions (inclusive).

## PolyA Tail Addition

The `--polyA` parameter allows you to add a polyA tail at a specific position in the sequence using two different formats:

### 1. Global Position Format (Recommended)

```
--polyA position/count
```

Where:
- `position` is the position in the combined sequence after which to add the polyA tail
- `count` is the number of A's to add

This format is more intuitive as it works with the absolute position in the final concatenated sequence, regardless of which segments are included.

### 2. Segment-Specific Format

```
--polyA segment:position/count
```

Where:
- `segment` is the segment index
- `position` is the position within that specific segment after which to add the polyA tail
- `count` is the number of A's to add

**Important Note**: When using the segment-specific format, the position must be within the valid range for that segment. For example, if segment 10 is 1258 bases long, the position must be between 1 and 1258.

In both formats, the slash notation (`/`) visually indicates the point after which the polyA tail will be added. The polyA tail replaces any downstream sequence that would have been included.

## Combining Multiple Operations

AtlasWriter allows you to combine multiple operations to create complex sequence modifications. The polyA tail can be useful for various modeling purposes, including mimicking natural mRNA processing or creating constructs for RNA stability studies.

### Order of Operations

When multiple operations are combined, they are applied in the following order:
1. Segment selection (as specified by `--include_locus_segment`)
2. Alternative start site (as specified by `--alt_start`)
3. Truncation (as specified by `--truncate`)
4. PolyA tail addition (as specified by `--polyA`)

### Combining PolyA Addition with Alternative Start Site

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-11 --alt_start /397 --polyA 3000/50
```

This will:
1. Combine segments 1 through 11
2. Start the sequence at position 397
3. Add a polyA tail of 50 A's after position 3000 of the resulting sequence
4. Discard any sequence after position 3000

### Combining Truncation with Alternative Start Site

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-11 --alt_start /397 --truncate 2:500-5:200
```

This will:
1. Combine segments 1 through 11
2. Start the sequence at position 397
3. Keep only the sequence between position 500 in segment 2 and position 200 in segment 5 (after applying the alternative start)

## Examples

Here are several examples that demonstrate different ways to use AtlasWriter:

### Listing Available Segments

To see what segments are available in your FASTA file:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --list
```

Output:
```
Loaded 11 locus segments from MZ242719_NL4-3_locus_segments.fasta

Available Locus Segments:
--------------------------------------------------------------------------------
Index |  Accession  |      Range      | Length  | Description
--------------------------------------------------------------------------------
  1   | MZ242719.1  |    1   -  290   |   290   | Mutant HIV-1 clone MSTRG.2.4 from USA genomic sequence
  2   | MZ242719.1  |   291  - 4459   |  4169   | Mutant HIV-1 clone MSTRG.2.4 from USA genomic sequence
  3   | MZ242719.1  |  4460  - 4509   |   50    | Mutant HIV-1 clone MSTRG.2.4 from USA genomic sequence
  ...
```

### Basic Segment Selection

Extract segments 1 through 3 and print to stdout:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-3
```

Extract specific segments and save to a file:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1,3,5 --output output.fasta
```

### Truncation Examples

To extract segments 1-5 but keep only the sequence between position 10 in segment 2 and position 20 in segment 4:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-5 --truncate 2:10-4:20
```

To cut the sequence after position 150 in segment 3:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-5 --truncate 3:150/
```

### Specify an Alternative Start Site

#### Using Global Position (Recommended)

To specify an alternative start site at position 397 in the combined sequence:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-11 --alt_start /397
```

This will:
1. Combine segments 1 through 11
2. Remove bases 1-396
3. Start the sequence at what was originally position 397
4. The first base in the output will be the base that was at position 397 in the original sequence

#### Using Segment-Specific Position

To specify an alternative start site at position 107 in segment 2:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-11 --alt_start 2:/107
```

This will:
1. Combine segments 1 through 11
2. Remove all bases before position 107 in segment 2
3. Start the sequence at what was position 107 in segment 2
4. The first base in the output will be the base that was at position 107 in segment 2

### Add a PolyA Tail

#### Using Global Position (Recommended)

To add a polyA tail after position 5697 in the combined sequence:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-10 --polyA 5697/50
```

This will:
1. Combine segments 1 through 10
2. Add 50 adenine bases (A's) after position 5697 in the combined sequence
3. Discard any sequence that would have followed position 5697

#### Using Segment-Specific Position

To add a polyA tail after position 200 in segment 10:

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-10 --polyA 10:200/50
```

This will:
1. Combine segments 1 through 10
2. Add 50 adenine bases (A's) after position 200 in segment 10
3. Discard any sequence that would have followed that position

### Working with HIV Genome Regions

Extract the gag gene region (approximately segments 2-3):

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 2-3
```

Extract the env gene region (approximately segments 7-8):

```bash
python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 7-8
```

## Output Format

The output is a properly formatted FASTA file with a header containing:
- The accession number
- Start and end positions
- Description including which segments were combined and if truncation was applied

Example output header:
```
>MZ242719.1:1-4509 Combined locus segments 1,2,3
```

## Handling Errors

AtlasWriter performs validation checks and provides clear error messages:

- If segments are out of range
- If truncation positions are invalid
- If file operations fail

## Notes on Coordinates

- Coordinates in the input FASTA file are 1-based
- Positions for truncation are also 1-based
- When truncating with `segment:position/`, the base at that position is included in the output

## Caveats

- The script assumes that the segments in the input FASTA file are in the correct order
- Segment indices in the tool (1,2,3,...) correspond to the order of sequences in the FASTA file

## Technical Details

AtlasWriter is written in Python and uses the following standard libraries:
- `argparse` for command-line argument parsing
- `re` for regular expression matching of FASTA headers
- `sys` for system operations

No external dependencies are required.
