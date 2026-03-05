process EXTRACT_KASP_MARKERS {
    tag "$meta.id"
    publishDir "${params.outdir}/kasp_markers", mode: 'copy'
    
    input:
    tuple val(meta), path(filtered_vcf), path(original_vcf)
    path reference_genome
    val flanking_region

    output:
    tuple val(meta), path("${meta.id}.kasp_markers.txt"), emit: kasp_markers

    script:
    """
    #!/usr/bin/env python3
    
    import pysam
    import re
    import sys
    
    def extract_kasp_markers(filtered_vcf_file, original_vcf_file, reference_file, flanking_bp, output_file):
        # Extract KASP markers from filtered VCF file
        
        try:
            # Load reference genome
            ref_fasta = pysam.FastaFile(reference_file)
            print(f"Loaded reference genome: {reference_file}")
        except Exception as e:
            print(f"Error loading reference genome: {e}")
            sys.exit(1)
        
        try:
            # Parse filtered VCF for KASP markers
            vcf = pysam.VariantFile(filtered_vcf_file)
            print(f"Loaded filtered VCF file: {filtered_vcf_file}")
        except Exception as e:
            print(f"Error loading filtered VCF file: {e}")
            sys.exit(1)
        
        kasp_markers = []
        skipped_variants = 0
        
        for record in vcf:
            chrom = record.chrom
            pos = record.pos - 1  # Convert to 0-based coordinates
            ref_allele = record.ref
            alt_allele = str(record.alts[0]) if record.alts else ""
            
            # Skip if not a simple SNP
            if len(ref_allele) != 1 or len(alt_allele) != 1:
                skipped_variants += 1
                print(f"Skipping complex variant at {chrom}:{pos+1} ({ref_allele}->{alt_allele})")
                continue
                
            # Get chromosome length
            try:
                chrom_length = ref_fasta.get_reference_length(chrom)
            except KeyError:
                print(f"Warning: Chromosome {chrom} not found in reference")
                skipped_variants += 1
                continue
            
            # Calculate flanking region boundaries
            start_pos = max(0, pos - flanking_bp)
            end_pos = min(chrom_length, pos + flanking_bp + 1)
            
            # Extract sequence
            try:
                sequence = ref_fasta.fetch(chrom, start_pos, end_pos).upper()
            except Exception as e:
                print(f"Warning: Could not extract sequence for {chrom}:{pos+1} - {e}")
                skipped_variants += 1
                continue
            
            # Find the SNP position within the extracted sequence
            snp_pos_in_seq = pos - start_pos
            
            # Verify we have the correct reference base and sufficient flanking sequence
            if (snp_pos_in_seq < len(sequence) and 
                snp_pos_in_seq >= 0 and 
                sequence[snp_pos_in_seq] == ref_allele):
                
                # Create KASP marker format
                left_flank = sequence[:snp_pos_in_seq]
                right_flank = sequence[snp_pos_in_seq + 1:]
                
                kasp_sequence = f"{left_flank}[{ref_allele}/{alt_allele}]{right_flank}"
                
                kasp_markers.append({
                    'id': f"{chrom}_{pos+1}_{ref_allele}_{alt_allele}",
                    'sequence': kasp_sequence,
                    'chrom': chrom,
                    'pos': pos + 1,
                    'ref': ref_allele,
                    'alt': alt_allele,
                    'genomic_start': start_pos,
                    'genomic_end': end_pos,
                    'snp_pos_in_seq': snp_pos_in_seq
                })
            else:
                print(f"Warning: Reference mismatch at {chrom}:{pos+1} - expected {ref_allele}, got {sequence[snp_pos_in_seq] if snp_pos_in_seq < len(sequence) else 'N/A'}")
                skipped_variants += 1
        
        print(f"Processed {len(kasp_markers)} valid gene-located SNPs, skipped {skipped_variants} variants")
        
        # Now mask other SNPs as 'N' in the flanking regions using original VCF
        if kasp_markers:
            masked_markers = mask_other_variants(kasp_markers, original_vcf_file, flanking_bp)
        else:
            masked_markers = []
        
        # Write output
        with open(output_file, 'w') as f:
            f.write("ID\\tSequence\\tChromosome\\tPosition\\tRef\\tAlt\\tLength\\tMaskedPositions\\tGeneLocation\\n")
            for marker in masked_markers:
                seq_length = len(marker['sequence'].replace('[', '').replace(']', '').replace('/', ''))
                masked_count = marker.get('masked_positions', 0)
                f.write(f"{marker['id']}\\t{marker['sequence']}\\t{marker['chrom']}\\t{marker['pos']}\\t{marker['ref']}\\t{marker['alt']}\\t{seq_length}\\t{masked_count}\\tGene\\n")
        
        ref_fasta.close()
        return len(masked_markers)
    
    def mask_other_variants(kasp_markers, original_vcf_file, flanking_bp):
        # Mask other variants as 'N' within the flanking regions
        
        print("Masking other variants in flanking regions...")
        print(f"Reading all variants from original VCF: {original_vcf_file}")
        
        # First pass: collect ALL variants from the original VCF (not just filtered ones)
        vcf = pysam.VariantFile(original_vcf_file)
        all_variants = {}
        variant_count = 0
        
        for record in vcf:
            chrom = record.chrom
            pos = record.pos - 1  # Convert to 0-based
            ref_len = len(record.ref)
            alt_len = len(str(record.alts[0])) if record.alts else 0
            
            if chrom not in all_variants:
                all_variants[chrom] = {}
            
            # For SNPs, mark just the position
            if ref_len == 1 and alt_len == 1:
                all_variants[chrom][pos] = {
                    'type': 'snp',
                    'ref': record.ref,
                    'alt': str(record.alts[0]) if record.alts else '',
                    'start_pos': pos
                }
                variant_count += 1
            else:
                # For MNPs and indels, mark all affected positions
                # For deletions, mark all deleted positions
                # For insertions, mark the insertion point
                # For MNPs, mark all substituted positions
                max_len = max(ref_len, alt_len)
                for i in range(max_len):
                    variant_pos = pos + i
                    if variant_pos not in all_variants[chrom]:  # Don't overwrite existing
                        all_variants[chrom][variant_pos] = {
                            'type': 'indel' if ref_len != alt_len else 'mnp',
                            'ref': record.ref,
                            'alt': str(record.alts[0]) if record.alts else '',
                            'start_pos': pos
                        }
                        variant_count += 1
        
        print(f"Found {variant_count} variant positions across {len(all_variants)} chromosomes")
        
        # Process each KASP marker
        masked_markers = []
        
        for marker in kasp_markers:
            chrom = marker['chrom']
            target_pos = marker['pos'] - 1  # Convert to 0-based
            sequence = marker['sequence']
            genomic_start = marker['genomic_start']
            snp_pos_in_seq = marker['snp_pos_in_seq']
            
            # Find the SNP position in the sequence (look for [REF/ALT] pattern)
            snp_match = re.search(r'\\[([ATCG])/([ATCG])\\]', sequence)
            if not snp_match:
                print(f"Warning: Could not find SNP pattern in sequence for {marker['id']}")
                masked_markers.append(marker)
                continue
            
            snp_marker_start = snp_match.start()
            snp_marker_end = snp_match.end()
            
            masked_positions = 0
            debug_info = []
            
            # Process the entire sequence and mask variants
            if chrom in all_variants:
                sequence_chars = list(sequence)
                
                # Check each position in the sequence (excluding the SNP marker itself)
                for seq_pos in range(len(sequence)):
                    # Skip the SNP marker part [REF/ALT]
                    if seq_pos >= snp_marker_start and seq_pos < snp_marker_end:
                        continue
                    
                    # Calculate genomic position for this sequence position
                    if seq_pos < snp_marker_start:
                        # Left flank
                        genomic_pos = genomic_start + seq_pos
                    else:
                        # Right flank - account for the SNP marker being replaced by single base
                        offset = seq_pos - snp_marker_end
                        genomic_pos = target_pos + 1 + offset
                    
                    # Check if this genomic position has a variant (excluding our target SNP)
                    if (genomic_pos in all_variants[chrom] and 
                        genomic_pos != target_pos):
                        
                        variant_info = all_variants[chrom][genomic_pos]
                        original_base = sequence_chars[seq_pos]
                        sequence_chars[seq_pos] = 'N'
                        masked_positions += 1
                        debug_info.append(f"pos {genomic_pos+1}({original_base}→N,{variant_info['type']})")
                
                # Reconstruct the sequence
                masked_sequence = ''.join(sequence_chars)
            else:
                masked_sequence = sequence
            
            if debug_info:
                print(f"Gene marker {marker['id']}: masked {debug_info}")
            elif masked_positions == 0:
                print(f"Gene marker {marker['id']}: no positions masked")
            
            masked_marker = marker.copy()
            masked_marker['sequence'] = masked_sequence
            masked_marker['masked_positions'] = masked_positions
            masked_markers.append(masked_marker)
        
        total_masked = sum(m.get('masked_positions', 0) for m in masked_markers)
        print(f"Total masked positions across all gene markers: {total_masked}")
        
        return masked_markers
    
    # Main execution
    try:
        num_markers = extract_kasp_markers("${filtered_vcf}", "${original_vcf}", "${reference_genome}", ${flanking_region}, "${meta.id}.kasp_markers.txt")
        print(f"Successfully extracted {num_markers} gene-located KASP markers for sample ${meta.id}")
        
        if num_markers == 0:
            print("Warning: No gene-located KASP markers were extracted. Check your VCF file, BED file, and reference genome.")
            # Create empty output file
            with open("${meta.id}.kasp_markers.txt", 'w') as f:
                f.write("ID\\tSequence\\tChromosome\\tPosition\\tRef\\tAlt\\tLength\\tMaskedPositions\\tGeneLocation\\n")
                
    except Exception as e:
        print(f"Error in gene-located KASP marker extraction: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    """
}