
# This method outputs information for each SNP, e.g. synonymous and non-synonymous etc.
# Its is called in lib/snp-search.rb

def information()

  strains = Strain.all

  output = File.open("#{out}", "w")

  output.puts "start_cds_in_ref\tend_cds_in_ref\tpos_of_SNP_in_ref\tSNP_base\tsynonymous or non-synonymous\tpossible_pseudogene?\tamino_acid_original\tamino_acid_change\tchange_in_hydrophobicity_of_AA?\tchange_in_polarisation_of_AA?\tchange_in_size_of_AA?\t#{strains.map{|strain| strain.name}.join("\t")}"

  snp_info = 0

  snps = Snp.find_by_sql("SELECT snps.* from snps inner join alleles on snps.id = alleles.snp_id")

  snps.each do |snp|
    snp.alleles.each do |allele|

      features = Feature.joins(:snps).where("snps.id = ?", snp.id)

      features.each do |feature|
        if feature.feature == "CDS" || "RNA"
          annotation = Annotation.where("annotations.qualifier = 'product' and annotations.feature_id = ?", feature.id).first
        else
          annotation = Annotation.where("annotations.qualifier = 'note' and annotations.feature_id = ?", feature.id).first
        end
 
  
    all_seqs_mutated = genome_sequence.seq
    mutated_seq_translated = []
    original_seq_translated = []

    # Mutate the sequence with the SNP
    all_seqs_mutated[snp.ref_pos.to_i-1] = allele.base

    # The reason why we use all the sequence from the genome and not just the feature sequence is because the ref_pos is based on the genome sequence and when we come to mutate it we have to have the full sequence to find the position.

    mutated_seq = Bio::Sequence.auto(all_seqs_mutated[feature.start-1..feature.end-1])

    original_seq =  Bio::Sequence.auto(all_seqs_original[feature.start-1..feature.end-1])

    #If the strand is negative then reverse complement

    if feature.strand == -1
      mutated_seq_translated << mutated_seq.reverse_complement.translate
      original_seq_translated << original_seq.reverse_complement.translate

    else
      mutated_seq_translated << mutated_seq.translate
      original_seq_translated << original_seq.translate

    end

    # Remove the star at the end of each translated sequence.

    mutated_seq_translated.zip(original_seq_translated).each do |mut, org|
      mutated_seq_translated_clean = mut.gsub(/\*$/,"")
      original_seq_translated_clean = org.gsub(/\*$/,"")

      # Amino acid properties

      hydrophobic = ["I", "L", "V", "C", "A", "G", "M", "F", "Y", "W", "H", "T"]
      non_hydrophobic = ["K", "E", "Q", "D", "N", "S", "P", "B"]

      polar = ["R", "N", "D", "E", "Q", "H", "K", "S", "T", "Y"]
      non_polar = ["A", "C", "G", "I", "L", "M", "F", "P", "W", "V"]

      small = ["V","C","A","G","D","N","S","T","P"]
      non_small = ["I","L","M","F","Y","W","H","K","R","E","Q"]

      alleles_array = []
      strains.each do |strain|
        allele = Allele.joins(:genotypes => :strain).where("strains.id = ? AND alleles.snp_id = ?", strain.id, snp.id).first
        alleles_array << allele.base
      end

      snp_info +=1

      if original_seq_translated_clean == mutated_seq_translated_clean
        if mutated_seq_translated_clean =~ /\*/
          output.puts "#{variant.start}\t#{variant.end}\t#{snp.ref_pos}\t#{Allele.find(snp.reference_allele_id).base}\tsynonymous\tYes\t\t\t\t\t\t#{alleles_array.join("\t")}"
        else
          output.puts "#{variant.start}\t#{variant.end}\t#{snp.ref_pos}\t#{Allele.find(snp.reference_allele_id).base}\tsynonymous\t\t\t\t\t\t\t#{alleles_array.join("\t")}"
        end
      else

        diffs = Diff::LCS.diff(original_seq_translated_clean, mutated_seq_translated_clean)

        if mutated_seq_translated_clean =~ /\*/
          output.puts "#{variant.start}\t#{variant.end}\t#{snp.ref_pos}\t#{Allele.find(snp.reference_allele_id).base}\tnon-synonymous\tYes\t#{diffs[0][0].element}\t#{diffs[0][1].element}\t#{'Yes' if (hydrophobic.include? diffs[0][0].element) == (non_hydrophobic.include? diffs[0][1].element)}\t#{'Yes' if (polar.include? diffs[0][0].element) == (non_polar.include? diffs[0][1].element)}\t#{'Yes' if (small.include? diffs[0][0].element) == (non_small.include? diffs[0][1].element)}\t#{alleles_array.join("\t")}"
          puts "#{variant.start}\t#{variant.end}\t#{snp.ref_pos}\t#{Allele.find(snp.reference_allele_id).base}\tnon-synonymous\tYes\t#{diffs[0][0].element}\t#{diffs[0][1].element}\t#{'Yes' if (hydrophobic.include? diffs[0][0].element) == (non_hydrophobic.include? diffs[0][1].element)}\t#{'Yes' if (polar.include? diffs[0][0].element) == (non_polar.include? diffs[0][1].element)}\t#{'Yes' if (small.include? diffs[0][0].element) == (non_small.include? diffs[0][1].element)}\t#{alleles_array.join("\t")}"
        else
          output.puts "#{variant.start}\t#{variant.end}\t#{snp.ref_pos}\t#{Allele.find(snp.reference_allele_id).base}\tnon-synonymous\t\t#{diffs[0][0].element}\t#{diffs[0][1].element}\t#{'Yes' if (hydrophobic.include? diffs[0][0].element) == (non_hydrophobic.include? diffs[0][1].element)}\t#{'Yes' if (polar.include? diffs[0][0].element) == (non_polar.include? diffs[0][1].element)}\t#{'Yes' if (small.include? diffs[0][0].element) == (non_small.include? diffs[0][1].element)}\t#{alleles_array.join("\t")}"
          puts "#{variant.start}\t#{variant.end}\t#{snp.ref_pos}\t#{Allele.find(snp.reference_allele_id).base}\tnon-synonymous\t\t#{diffs[0][0].element}\t#{diffs[0][1].element}\t#{'Yes' if (hydrophobic.include? diffs[0][0].element) == (non_hydrophobic.include? diffs[0][1].element)}\t#{'Yes' if (polar.include? diffs[0][0].element) == (non_polar.include? diffs[0][1].element)}\t#{'Yes' if (small.include? diffs[0][0].element) == (non_small.include? diffs[0][1].element)}\t#{alleles_array.join("\t")}"
        end
      end
      puts "Total SNPs outputted so far: #{snp_info}" if snp_info % 50 == 0 
    end
  end
end

#Take all SNP positions in ref genome
# snp_positions = Feature.find_by_sql("select snps.ref_pos from features inner join snps on features.id = snps.feature_id inner join alleles on snps.id = alleles.snp_id where alleles.id <> snps.reference_allele_id and features.name = 'CDS'").map{|snp| snp.ref_pos}

# # Take all SNP nucleotide
# snps = Feature.find_by_sql("select alleles.base from features inner join snps on features.id = snps.feature_id inner join alleles on snps.id = alleles.snp_id where alleles.id <> snps.reference_allele_id and features.name = 'CDS'").map{|allele| allele.base}

# # Mutate (substitute) the original sequence with the SNPs

# # Here all_seqs_original are all the nucelotide sequences but with the snps subsituted in them

# #Get start position of CDS with SNP
# coordinates_start = Feature.find_by_sql("select start from features inner join snps on features.id = snps.feature_id inner join alleles on snps.id = alleles.snp_id where features.name = 'CDS' and alleles.id <> snps.reference_allele_id").map{|feature| feature.start}

# #Get end position of CDS with SNP
# coordinates_end = Feature.find_by_sql("select end from features inner join snps on features.id = snps.feature_id inner join alleles on snps.id = alleles.snp_id where features.name = 'CDS' and alleles.id <> snps.reference_allele_id").map{|feature| feature.end} 
