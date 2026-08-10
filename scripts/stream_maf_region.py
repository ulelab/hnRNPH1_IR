#!/usr/bin/env python3
"""
Stream MAF file from UCSC multiz470way server and extract a specific region
Usage: python stream_maf_region.py <chr> <start> <end> [output.maf]
"""

import sys
import urllib.request
import gzip

def parse_maf_src(src):
    """Parse MAF source line to extract chromosome
    Format: <db>.<chrom> (e.g., hg38.chr5)
    """
    if '.' in src:
        return src.split('.', 1)[1]
    return None

def region_overlaps(start, end, target_start, target_end):
    """Check if a MAF block overlaps with target region"""
    return not (end <= target_start or start >= target_end)

def stream_maf_region(maf_url, chrom, target_start, target_end, output_file=None):
    """
    Stream MAF file from URL and extract matching region
    """
    # Open URL stream
    try:
        print(f"Connecting to: {maf_url}", file=sys.stderr)
        if maf_url.endswith('.gz'):
            response = urllib.request.urlopen(maf_url)
            maf_stream = gzip.open(response, 'rt')
        else:
            response = urllib.request.urlopen(maf_url)
            maf_stream = response
            # Decode bytes to string if needed
            if hasattr(maf_stream, 'read'):
                import io
                maf_stream = io.TextIOWrapper(maf_stream, encoding='utf-8')
    except Exception as e:
        print(f"ERROR: Failed to open URL {maf_url}: {e}", file=sys.stderr)
        return False
    
    # Set up output
    if output_file:
        out = open(output_file, 'w')
    else:
        out = sys.stdout
    
    # Process MAF file
    maf_header = None
    current_block = []
    in_block = False
    block_matches = False
    blocks_found = 0
    
    try:
        for line in maf_stream:
            # Handle bytes vs string
            if isinstance(line, bytes):
                line_str = line.decode('utf-8')
            else:
                line_str = line
            
            # Capture header
            if line_str.startswith('##maf'):
                maf_header = line_str
                out.write(line_str)
                continue
            
            # Start of alignment block
            if line_str.startswith('a '):
                # Process previous block if it matched
                if in_block and block_matches:
                    for block_line in current_block:
                        if isinstance(block_line, bytes):
                            block_line = block_line.decode('utf-8')
                        out.write(block_line)
                    out.write('\n')
                    blocks_found += 1
                
                # Reset for new block
                current_block = [line_str]
                in_block = False
                block_matches = False
                continue
            
            # Sequence line in alignment block
            if line_str.startswith('s '):
                current_block.append(line_str)
                fields = line_str.strip().split()
                if len(fields) >= 6:
                    src = fields[1]  # e.g., "hg38.chr5"
                    parsed_chrom = parse_maf_src(src)
                    
                    # Check if this is the reference species (hg38) and matches our chromosome
                    if src.startswith('hg38.') and parsed_chrom == chrom:
                        start = int(fields[2])
                        size = int(fields[3])
                        end = start + size
                        
                        # Check if this block overlaps with our target region
                        if region_overlaps(start, end, target_start, target_end):
                            in_block = True
                            block_matches = True
                continue
            
            # Other lines in block (i, e, q, etc.)
            if current_block:
                current_block.append(line_str)
            
            # Empty line ends block
            if line_str.strip() == '':
                if in_block and block_matches:
                    for block_line in current_block:
                        if isinstance(block_line, bytes):
                            block_line = block_line.decode('utf-8')
                        out.write(block_line)
                    out.write('\n')
                    blocks_found += 1
                current_block = []
                in_block = False
                block_matches = False
        
        # Handle last block if needed
        if in_block and block_matches and current_block:
            for block_line in current_block:
                if isinstance(block_line, bytes):
                    block_line = block_line.decode('utf-8')
                out.write(block_line)
            out.write('\n')
            blocks_found += 1
    
    except Exception as e:
        print(f"ERROR: Failed to process MAF stream: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False
    finally:
        maf_stream.close()
        if output_file:
            out.close()
    
    print(f"Extracted {blocks_found} MAF blocks for region {chrom}:{target_start}-{target_end}", file=sys.stderr)
    return True

def get_maf_url(chrom, base_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf"):
    """Construct MAF file URL for a chromosome"""
    # Handle special chromosomes
    if chrom.startswith('chr'):
        maf_name = f"{chrom}.maf"
    else:
        maf_name = f"chr{chrom}.maf"
    
    return f"{base_url}/{maf_name}"

def main():
    if len(sys.argv) < 4:
        print("Usage: python stream_maf_region.py <chr> <start> <end> [output.maf]", file=sys.stderr)
        print("Example: python stream_maf_region.py chr5 179619358 179620941 output.maf", file=sys.stderr)
        sys.exit(1)
    
    chrom = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    output_file = sys.argv[4] if len(sys.argv) > 4 else None
    
    # Construct MAF URL
    maf_url = get_maf_url(chrom)
    
    print(f"Extracting region: {chrom}:{start}-{end}", file=sys.stderr)
    print(f"MAF URL: {maf_url}", file=sys.stderr)
    
    # Stream and extract
    success = stream_maf_region(maf_url, chrom, start, end, output_file)
    
    if success:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
