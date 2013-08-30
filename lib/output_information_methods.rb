
#This method finds info_snps from a list of strains given by user.
# Its is called in lib/snp-search.rb

def output_information_methods(snps, outfile, cuttoff_genotype, cuttoff_snp, info = true)

  strains = Strain.all

  outfile.puts "pos_of_SNP_in_ref\tref_base\tSNP_base\tsynonymous or non-synonymous\tGene_annotation\tpossible_pseudogene?\tamino_acid_original\tamino_acid_change\tchange_in_hydrophobicity_of_AA?\tchange_in_polarisation_of_AA?\tchange_in_size_of_AA?\t#{strains.map{|strain| strain.name}.join("\t") if info}"

  snps_counter = 0
  cds_snps_counter = 0
  total_number_of_syn_snps = 0
  total_number_of_non_syn_snps = 0
  total_number_of_pseudo = 0
  snps.each do |snp|

    ActiveRecord::Base.transaction do
      snp.alleles.each do |allele|
        next if snp.alleles.any?{|allele| allele.base.length > 1} # indel
        if allele.id != snp.reference_allele_id
          
          # get annotation (if there is any) for each SNP
          features = Feature.joins(:snps).where("snps.id = ?", snp.id)
          
          # get snp quality for each snp
          snp_qual = Snp.find_by_sql("select qual from snps where snps.id = #{snp.id}")
          next if snp_qual.any?{|snps_quality| snps_quality.qual < cuttoff_snp.to_i}
          # ignore snp if the snp qual is less than cuttoff.
          # next if snp.snp_qual < cuttoff_snp.to_i

          # get all genotype qualities for each snp.
          gqs = Genotype.find_by_sql("select geno_qual from genotypes inner join alleles on alleles.id = genotypes.allele_id inner join snps on snps.id = alleles.snp_id where snps.id = #{snp.id}")
          # ignore snp if any of its genotype qualities is lower than the cuttoff.          
          next if gqs.any?{|genotype_quality| genotype_quality.geno_qual < cuttoff_genotype.to_i}
          
          ref_base = Bio::Sequence.auto(Allele.find(snp.reference_allele_id).base)
          snp_base = Bio::Sequence.auto(allele.base)
          # count snps now: after you have selected the snps with gqs and snp_qual greater than the threshold.
          snps_counter += 1
          # If the feature is empty then just output basic information about the snp.

          if features.empty?
            outfile.puts "#{snp.ref_pos}\t#{features.map{|feature| feature.strand == 1} ? "#{ref_base.upcase}" : "#{ref_base.reverse_complement.upcase}"}\t#{features.map{|feature| feature.strand == 1} ? "#{snp_base.upcase}" : "#{snp_base.reverse_complement.upcase}"}"
          
          else  
            features.each do |feature|
              if feature.name == "CDS"

                cds_snps_counter +=1

                annotation = Annotation.where("annotations.qualifier = 'product' and annotations.feature_id = ?", feature.id).first
                #if annotation is nil, or empty
                if annotation.nil?
                  outfile.puts "#{snp.ref_pos}\t#{feature.strand == 1 ? "#{ref_base.upcase}" : "#{ref_base.reverse_complement.upcase}"}\t#{feature.strand == 1 ? "#{snp_base.upcase}" : "#{snp_base.reverse_complement.upcase}"}"
                else

                  feature_sequence = feature.sequence

                  feature_sequence_bio = Bio::Sequence::NA.new(feature_sequence)

                  #Mutate sequence with SNP
                  feature_sequence_mutated = feature.sequence
                  feature_sequence_snp_pos = (snp.ref_pos-1) - (feature.start-1)
                  feature_sequence_mutated[feature_sequence_snp_pos] = allele.base
                  feature_sequence_mutated_bio = Bio::Sequence::NA.new(feature_sequence_mutated)

                  # Translate the sequences
                  if feature.strand == -1
                    mutated_seq_translated = feature_sequence_mutated_bio.reverse_complement.translate
                    original_seq_translated = feature_sequence_bio.reverse_complement.translate

                  else
                    mutated_seq_translated = feature_sequence_mutated_bio.translate
                    original_seq_translated = feature_sequence_bio.translate

                  end

                  # Remove the star at the end of each translated sequence.
                  mutated_seq_translated_clean = mutated_seq_translated.gsub(/\*$/,"")
                  original_seq_translated_clean = original_seq_translated.gsub(/\*$/,"")

                  # Amino acid properties
                  hydrophobic = ["I", "L", "V", "C", "A", "G", "M", "F", "Y", "W", "H", "T"]
                  non_hydrophobic = ["K", "E", "Q", "D", "N", "S", "P", "B"]

                  polar = ["Y", "W", "H", "K", "R", "E", "Q", "D", "N", "S", "P", "B"]
                  non_polar = ["I", "L", "V", "C", "A", "G", "M", "F", "T"]

                  small = ["V","C","A","G","D","N","S","T","P"]
                  non_small = ["I","L","M","F","Y","W","H","K","R","E","Q"]

                  # Get alleles for each strain
                  alleles_array = []
                    strains.each do |strain|
                      allele_for_strains = Allele.joins(:genotypes => :strain).where("strains.id = ? AND alleles.snp_id = ?", strain.id, snp.id).first
                      alleles_array << allele_for_strains.base
                    end

                  # If no difference between the amino acids then its synonymous SNP, if different then its non-synonymous.
                  if original_seq_translated_clean == mutated_seq_translated_clean
                    total_number_of_syn_snps +=1
                    if mutated_seq_translated_clean =~ /\*/
                      total_number_of_pseudo +=1
                      outfile.puts "#{snp.ref_pos}\t#{features.map{|feature| feature.strand == 1} ? "#{ref_base.upcase}" : "#{ref_base.reverse_complement.upcase}"}\t#{features.map{|feature| feature.strand == 1} ? "#{snp_base.upcase}" : "#{snp_base.reverse_complement.upcase}"}\tsynonymous\t#{annotation.value}\tYes\tN/A\tN/A\tN/A\tN/A\tN/A\t#{alleles_array.join("\t") if info}"
                    else
                      outfile.puts "#{snp.ref_pos}\t#{features.map{|feature| feature.strand == 1} ? "#{ref_base.upcase}" : "#{ref_base.reverse_complement.upcase}"}\t#{features.map{|feature| feature.strand == 1} ? "#{snp_base.upcase}" : "#{snp_base.reverse_complement.upcase}"}\tsynonymous\t#{annotation.value}\tNo\tN/A\tN/A\tN/A\tN/A\tN/A\t#{alleles_array.join("\t") if info}"
                    end
                  else
                    total_number_of_non_syn_snps +=1
                    diffs = Diff::LCS.diff(original_seq_translated_clean, mutated_seq_translated_clean)

                    if mutated_seq_translated_clean =~ /\*/
                      total_number_of_pseudo +=1
                      outfile.puts "#{snp.ref_pos}\t#{features.map{|feature| feature.strand == 1} ? "#{ref_base.upcase}" : "#{ref_base.reverse_complement.upcase}"}\t#{features.map{|feature| feature.strand == 1} ? "#{snp_base.upcase}" : "#{snp_base.reverse_complement.upcase}"}\tnon-synonymous\t#{annotation.value}\tYes\t#{diffs[0][0].element}\t#{diffs[0][1].element}\t#{'Yes' if (hydrophobic.include? diffs[0][0].element) == (non_hydrophobic.include? diffs[0][1].element)}#{'No' if (hydrophobic.include? diffs[0][0].element) != (non_hydrophobic.include? diffs[0][1].element)}\t#{'Yes' if (polar.include? diffs[0][0].element) == (non_polar.include? diffs[0][1].element)}#{'No' if (polar.include? diffs[0][0].element) != (non_polar.include? diffs[0][1].element)}\t#{'Yes' if (small.include? diffs[0][0].element) == (non_small.include? diffs[0][1].element)}#{'No' if (small.include? diffs[0][0].element) != (non_small.include? diffs[0][1].element)}\t#{alleles_array.join("\t") if info}"
                    else
                      outfile.puts "#{snp.ref_pos-1}\t#{features.map{|feature| feature.strand == 1} ? "#{ref_base.upcase}" : "#{ref_base.reverse_complement.upcase}"}\t#{features.map{|feature| feature.strand == 1} ? "#{snp_base.upcase}" : "#{snp_base.reverse_complement.upcase}"}\tnon-synonymous\t#{annotation.value}\tNo\t#{diffs[0][0].element}\t#{diffs[0][1].element}\t#{'Yes' if (hydrophobic.include? diffs[0][0].element) == (non_hydrophobic.include? diffs[0][1].element)}#{'No' if (hydrophobic.include? diffs[0][0].element) != (non_hydrophobic.include? diffs[0][1].element)}\t#{'Yes' if (polar.include? diffs[0][0].element) == (non_polar.include? diffs[0][1].element)}#{'No' if (polar.include? diffs[0][0].element) != (non_polar.include? diffs[0][1].element)}\t#{'Yes' if (small.include? diffs[0][0].element) == (non_small.include? diffs[0][1].element)}#{'No' if (small.include? diffs[0][0].element) != (non_small.include? diffs[0][1].element)}\t#{alleles_array.join("\t") if info}"
                    end
                  end
                end
              end
            end
          end
        end
      end
      puts "Total SNPs added so far: #{snps_counter}" if snps_counter % 100 == 0 
    end
  end
  puts "Total number of snps: #{snps_counter}"
  puts "Total number of snps in CDS region: #{cds_snps_counter}"
  puts "Total number of synonymous SNPs: #{total_number_of_syn_snps}"
  puts "Total number of non-synonymous SNPs: #{total_number_of_non_syn_snps}"
  puts "Total number of pseudogenes: #{total_number_of_pseudo}"
  outfile.puts "Total number of snps: #{snps_counter}"
  outfile.puts "Total number of snps in CDS region: #{cds_snps_counter}"
  outfile.puts "Total number of synonymous SNPs: #{total_number_of_syn_snps}"
  outfile.puts "Total number of non-synonymous SNPs: #{total_number_of_non_syn_snps}"
  outfile.puts "Total number of possible pseudogenes: #{total_number_of_pseudo}"
end
