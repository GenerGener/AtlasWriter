#!/usr/bin/env python3
"""
AtlasWriter - A tool for extracting and manipulating HIV locus segments from FASTA files.

This script allows selection of specific locus segments by index and supports truncation
of sequences at specified positions.
"""

import argparse
import re
import sys
from typing import List, Tuple, Dict, Optional


class AtlasWriter:
    def __init__(self, fasta_file: str):
        """
        Initialize the AtlasWriter with a FASTA file containing locus segments.
        
        Args:
            fasta_file: Path to the FASTA file containing locus segments
        """
        self.fasta_file = fasta_file
        self.segments = []
        self.segment_data = {}
        self.load_segments()
    
    def load_segments(self) -> None:
        """Load and parse the locus segments from the FASTA file."""
        try:
            with open(self.fasta_file, 'r') as file:
                content = file.read()
                
            # Extract header and sequence pairs from the FASTA file
            entries = re.findall(r'(>.+?\n)([^>]+)', content, re.DOTALL)
            
            for i, (header, sequence) in enumerate(entries, 1):
                # Clean up header and sequence
                header = header.strip()
                sequence = sequence.replace('\n', '').strip()
                
                # Extract position information from header (e.g., MZ242719.1:1-290)
                match = re.search(r'([^:]+):(\d+)-(\d+)', header)
                if match:
                    accession = match.group(1)
                    start = int(match.group(2))
                    end = int(match.group(3))
                    
                    # Store segment information
                    self.segments.append((i, accession, start, end, header, sequence))
                    self.segment_data[i] = {
                        'accession': accession,
                        'start': start,
                        'end': end,
                        'header': header,
                        'sequence': sequence
                    }
                else:
                    print(f"Warning: Could not parse position information from header: {header}")
            
            # Sort segments by start position
            self.segments.sort(key=lambda x: x[2])
            
            print(f"Loaded {len(self.segments)} locus segments from {self.fasta_file}")
        
        except FileNotFoundError:
            sys.exit(f"Error: FASTA file '{self.fasta_file}' not found")
        except Exception as e:
            sys.exit(f"Error loading FASTA file: {str(e)}")
    
    def parse_segment_selection(self, selection: str) -> List[int]:
        """
        Parse a segment selection string into a list of segment indices.
        
        Args:
            selection: A string like "1,2,3,4,5" or "1-5,7,9-11"
            
        Returns:
            A list of segment indices
        """
        selected_indices = []
        
        # Split the selection string by commas
        parts = selection.split(',')
        
        for part in parts:
            if '-' in part:
                # Handle ranges (e.g., "1-5")
                try:
                    start, end = map(int, part.split('-'))
                    if start < 1 or end > len(self.segments) or start > end:
                        raise ValueError(f"Invalid range: {part}")
                    selected_indices.extend(range(start, end + 1))
                except ValueError:
                    sys.exit(f"Error: Invalid range format in '{part}'")
            else:
                # Handle single indices
                try:
                    index = int(part)
                    if index < 1 or index > len(self.segments):
                        raise ValueError(f"Invalid segment index: {index}")
                    selected_indices.append(index)
                except ValueError:
                    sys.exit(f"Error: Invalid segment index: {part}")
        
        # Remove duplicates and sort
        return sorted(set(selected_indices))
    
    def parse_truncation(self, truncation: str, selected_indices: List[int]) -> Tuple[Optional[int], Optional[int], Optional[int]]:
        """
        Parse a truncation string to determine start and end positions.
        
        Args:
            truncation: A string specifying where to truncate the sequence
                - For a range: "upstream_segment:upstream_base-downstream_segment:downstream_base"
                - For a single cut: "segment:position/downstream_base"
            selected_indices: List of selected segment indices
            
        Returns:
            A tuple of (global_start, global_end, poly_a_count)
        """
        if not truncation:
            return None, None, None
        
        # Check for polyA specification
        poly_a_count = None
        
        if "-" in truncation:
            # Handle range truncation (keep sequence between two positions)
            try:
                start_spec, end_spec = truncation.split('-')
                
                start_segment, start_pos = map(int, start_spec.split(':'))
                end_segment, end_pos = map(int, end_spec.split(':'))
                
                if start_segment not in selected_indices or end_segment not in selected_indices:
                    sys.exit(f"Error: Truncation segments must be included in the selected segments")
                
                # Calculate global positions
                global_start = sum(len(self.segment_data[i]['sequence']) for i in selected_indices 
                                  if i < start_segment)
                global_start += start_pos - 1  # Convert to 0-based index (keep the base at start_pos)
                
                global_end = sum(len(self.segment_data[i]['sequence']) for i in selected_indices 
                                if i < end_segment)
                global_end += end_pos  # Include the base at end_pos
                
                return global_start, global_end, poly_a_count
            
            except ValueError:
                sys.exit(f"Error: Invalid truncation format. Use 'upstream_segment:upstream_base-downstream_segment:downstream_base'")
        elif "/" in truncation:
            # Handle single position truncation with explicit upstream/downstream specification
            try:
                segment_pos_part = truncation.split('/')[0]
                poly_a_part = truncation.split('/')[1] if len(truncation.split('/')) > 1 else None
                
                segment, position = map(int, segment_pos_part.split(':'))
                
                if poly_a_part and poly_a_part.isdigit():
                    poly_a_count = int(poly_a_part)
                
                if segment not in selected_indices:
                    sys.exit(f"Error: Truncation segment must be included in the selected segments")
                
                # Calculate global position
                global_pos = sum(len(self.segment_data[i]['sequence']) for i in selected_indices 
                                if i < segment)
                global_pos += position  # This is the position AFTER which we'll cut
                
                return None, global_pos, poly_a_count  # Truncate up to this position
            
            except ValueError:
                sys.exit(f"Error: Invalid truncation format. Use 'segment:position'")
        else:
            # Backward compatibility for segment:position format
            try:
                segment, position = map(int, truncation.split(':'))
                
                if segment not in selected_indices:
                    sys.exit(f"Error: Truncation segment must be included in the selected segments")
                
                # Calculate global position
                global_pos = sum(len(self.segment_data[i]['sequence']) for i in selected_indices 
                                if i < segment)
                global_pos += position
                
                print(f"Warning: Deprecated truncation format. Consider using 'segment:position/' for clarity.")
                return None, global_pos, poly_a_count  # Truncate up to this position
            
            except ValueError:
                sys.exit(f"Error: Invalid truncation format. Use 'segment:position/' or 'upstream_segment:upstream_base-downstream_segment:downstream_base'")
    
    def parse_poly_a(self, poly_a: str, selected_indices: List[int]) -> Tuple[Optional[int], Optional[int]]:
        """
        Parse a polyA string to determine position and count.
        
        Args:
            poly_a: A string like "segment:position/count" or "position/count"
            selected_indices: List of selected segment indices
            
        Returns:
            A tuple of (global_position, count)
        """
        if not poly_a:
            return None, None
        
        try:
            if "/" in poly_a:
                position_part, count_str = poly_a.split("/")
                count = int(count_str)
                
                # Check if format is segment:position or just a global position
                if ":" in position_part:
                    segment, position = map(int, position_part.split(':'))
                    
                    if segment not in selected_indices:
                        sys.exit(f"Error: PolyA segment must be included in the selected segments")
                    
                    # Check if position is valid for the specified segment
                    if position > len(self.segment_data[segment]['sequence']):
                        sys.exit(f"Error: Position {position} is beyond the length of segment {segment} (length: {len(self.segment_data[segment]['sequence'])})")
                    
                    # Calculate the correct global position by summing only selected segments
                    global_pos = 0
                    for idx in selected_indices:
                        if idx < segment:
                            global_pos += len(self.segment_data[idx]['sequence'])
                        elif idx == segment:
                            global_pos += position
                            break
                else:
                    # Treat as a global position directly
                    global_pos = int(position_part)
                
                return global_pos, count
            else:
                sys.exit(f"Error: Invalid polyA format. Use 'segment:position/count' or 'position/count'")
        except ValueError:
            sys.exit(f"Error: Invalid polyA format. Use 'segment:position/count' or 'position/count'")
    
    def combine_segments(self, indices: List[int], 
                         global_start: Optional[int] = None, 
                         global_end: Optional[int] = None,
                         poly_a_pos: Optional[int] = None,
                         poly_a_count: Optional[int] = None) -> Tuple[str, str]:
        """
        Combine selected segments into a single sequence, with optional truncation and polyA addition.
        
        Args:
            indices: List of segment indices to combine
            global_start: Optional starting position for truncation (inclusive)
            global_end: Optional ending position for truncation (exclusive)
            poly_a_pos: Optional position to add polyA tail
            poly_a_count: Optional number of A's to add
            
        Returns:
            A tuple of (header, sequence)
        """
        if not indices:
            sys.exit("Error: No segments selected")
        
        # Get the first and last segments for header information
        first_segment = self.segment_data[indices[0]]
        last_segment = self.segment_data[indices[-1]]
        
        # Create a new header
        accession = first_segment['accession'].split('.')[0]
        start_pos = first_segment['start']
        end_pos = last_segment['end']
        
        # Combine sequences
        combined_sequence = ''
        for idx in indices:
            combined_sequence += self.segment_data[idx]['sequence']
        
        # Apply truncation if specified
        truncated_sequence = combined_sequence
        if global_start is not None or global_end is not None:
            start_idx = global_start if global_start is not None else 0
            end_idx = global_end if global_end is not None else len(combined_sequence)
            
            if start_idx < 0 or end_idx > len(combined_sequence) or start_idx >= end_idx:
                sys.exit(f"Error: Invalid truncation range: {start_idx}-{end_idx}")
            
            truncated_sequence = combined_sequence[start_idx:end_idx]
            
            # Update the positions in the header
            if global_start is not None:
                segment_idx = 0
                cumulative_length = 0
                
                # Find which segment contains the start position
                for i, idx in enumerate(indices):
                    seg_length = len(self.segment_data[idx]['sequence'])
                    if cumulative_length + seg_length > global_start:
                        segment_idx = i
                        break
                    cumulative_length += seg_length
                
                # Calculate the new start position
                start_seg = self.segment_data[indices[segment_idx]]
                relative_pos = global_start - cumulative_length
                start_pos = start_seg['start'] + relative_pos
            
            if global_end is not None:
                segment_idx = 0
                cumulative_length = 0
                
                # Find which segment contains the end position
                for i, idx in enumerate(indices):
                    seg_length = len(self.segment_data[idx]['sequence'])
                    if cumulative_length + seg_length >= global_end:
                        segment_idx = i
                        break
                    cumulative_length += seg_length
                
                # Calculate the new end position
                end_seg = self.segment_data[indices[segment_idx]]
                relative_pos = global_end - cumulative_length
                end_pos = end_seg['start'] + relative_pos - 1
        
        # Construct the final header
        description = f"Combined locus segments {','.join(map(str, indices))}"
        if global_start is not None or global_end is not None:
            description += f" (truncated)"
        
        # Apply polyA tail if specified
        final_sequence = truncated_sequence
        if poly_a_pos is not None and poly_a_count is not None:
            if poly_a_pos > len(truncated_sequence):
                sys.exit(f"Error: PolyA position {poly_a_pos} is beyond the length of the combined sequence {len(truncated_sequence)}")
            
            final_sequence = truncated_sequence[:poly_a_pos] + 'A' * poly_a_count
            description += f" (with {poly_a_count} polyA tail)"
        else:
            final_sequence = truncated_sequence
        
        header = f">{accession}.1:{start_pos}-{end_pos} {description}"
        
        return header, final_sequence
    
    def format_fasta(self, header: str, sequence: str, width: int = 70) -> str:
        """
        Format a sequence as FASTA with specified line width.
        
        Args:
            header: The FASTA header
            sequence: The sequence to format
            width: Line width for the sequence
            
        Returns:
            Formatted FASTA string
        """
        # Add header
        fasta = header + '\n'
        
        # Add sequence with specified line width
        for i in range(0, len(sequence), width):
            fasta += sequence[i:i+width] + '\n'
        
        return fasta
    
    def write_output(self, header: str, sequence: str, output_file: Optional[str] = None) -> None:
        """
        Write the formatted FASTA output to a file or stdout.
        
        Args:
            header: The FASTA header
            sequence: The sequence to write
            output_file: Optional file path for output
        """
        formatted_fasta = self.format_fasta(header, sequence)
        
        if output_file:
            try:
                with open(output_file, 'w') as file:
                    file.write(formatted_fasta)
                print(f"Output written to {output_file}")
            except Exception as e:
                sys.exit(f"Error writing to file: {str(e)}")
        else:
            print(formatted_fasta)
    
    def print_segment_info(self) -> None:
        """Print information about available segments."""
        print("\nAvailable Locus Segments:")
        print("-" * 80)
        print(f"{'Index':^5} | {'Accession':^10} | {'Range':^15} | {'Length':^8} | {'Description'}")
        print("-" * 80)
        
        for idx, segment in sorted(self.segment_data.items()):
            accession = segment['accession']
            start = segment['start']
            end = segment['end']
            length = len(segment['sequence'])
            
            # Extract description from header if available
            header = segment['header']
            description = header.split(' ', 1)[1] if ' ' in header else "No description"
            
            print(f"{idx:^5} | {accession:^10} | {start:^7}-{end:^7} | {length:^8} | {description}")


def main():
    parser = argparse.ArgumentParser(
        description="AtlasWriter - A tool for extracting and manipulating HIV locus segments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Select segments 1-3 and output to stdout
  python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-3
  
  # Select segments 1,3,5 and save to output.fasta
  python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1,3,5 --output output.fasta
  
  # Select segments 1-5 and keep sequence between position 10 of segment 2 and position 20 of segment 4
  python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-5 --truncate 2:10-4:20
  
  # Cut after position 150 in segment 3 (keeps bases upstream of the cut, removes downstream bases)
  python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-5 --truncate 3:150/
  
  # Cut after position 150 in segment 3 and add 20 A's
  python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-5 --polyA 3:150/20
  
  # Add 50 A's after global position 5000 in the combined sequence
  python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --include_locus_segment 1-10 --polyA 5000/50
  
  # List available segments
  python AtlasWriter.py --input MZ242719_NL4-3_locus_segments.fasta --list
"""
    )
    
    parser.add_argument("--input", type=str, required=True,
                        help="Input FASTA file containing locus segments")
    
    parser.add_argument("--include_locus_segment", type=str, 
                        help="Comma-separated list of segment indices or ranges (e.g., '1,2,3' or '1-3,5,7-9')")
    
    parser.add_argument("--truncate", type=str,
                        help="Truncate sequence at specified positions. Use 'segment:position/' to cut after position (keep upstream bases) or "
                             "'upstream_segment:upstream_base-downstream_segment:downstream_base' to keep sequence between the two bases")
    
    parser.add_argument("--polyA", type=str, 
                        help="Add a polyA tail after a specified position. Use 'segment:position/N' where N is the number of A's to add, "
                             "or simply 'position/N' to use a global position in the combined sequence")
    
    parser.add_argument("--output", type=str,
                        help="Output file path (if not specified, output is printed to stdout)")
    
    parser.add_argument("--list", action="store_true",
                        help="List available segments and exit")
    
    args = parser.parse_args()
    
    # Initialize AtlasWriter
    atlas = AtlasWriter(args.input)
    
    # If --list flag is provided, print segment info and exit
    if args.list:
        atlas.print_segment_info()
        return
    
    # Check if segment selection is provided
    if not args.include_locus_segment:
        sys.exit("Error: No segment selection provided. Use --include_locus_segment or --list")
    
    # Parse segment selection
    selected_indices = atlas.parse_segment_selection(args.include_locus_segment)
    
    # Parse truncation if provided
    global_start, global_end, poly_a_count_from_truncate = None, None, None
    if args.truncate:
        global_start, global_end, poly_a_count_from_truncate = atlas.parse_truncation(args.truncate, selected_indices)
    
    # Parse polyA if provided
    poly_a_pos, poly_a_count = None, None
    if args.polyA:
        poly_a_pos, poly_a_count = atlas.parse_poly_a(args.polyA, selected_indices)
    elif poly_a_count_from_truncate is not None:
        # If polyA was specified in truncate parameter
        poly_a_pos = global_end
        poly_a_count = poly_a_count_from_truncate
    
    # Combine selected segments
    header, sequence = atlas.combine_segments(selected_indices, global_start, global_end, poly_a_pos, poly_a_count)
    
    # Write output
    atlas.write_output(header, sequence, args.output)


if __name__ == "__main__":
    main()
